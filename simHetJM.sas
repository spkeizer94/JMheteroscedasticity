%macro simHetJM(nsim=, n=, t=, ve=, u0=, u1=, u2=, u3=, 
             alpha0=,alpha1=, gamma=, s0=, s1=, shape=, scale=,
               v00=, v11=, v22=, v33=, v44=, cov01=, cov02=, cov03=, cov04=, cov12=, cov13=, cov14=, cov23=, cov24=, cov34=,
				startsim =,endsim=);
%let seed1 = 12345;
%let seed2 = 54321;
%let seed3 = 24680; 
%let seed4 = 13579;
%let seed5 = 112233
/***************************/;
/*Simulate latent variables*/;
/***************************/;
proc iml;
    call randseed(12345);
    sim     =   colvec(repeat(T(1:&nsim),1, &n));
    id      =   colvec(repeat(T(1:&n), &nsim));
    mean    =   {0          0          0       0    0};
    cov     =   {&v00      &cov01    &cov02    &cov03  &cov04,
                 &cov01    &v11      &cov12    &cov13  &cov14,
                 &cov02    &cov12    &v22      &cov23  &cov24,
                 &cov03    &cov13    &cov23    &v33    &cov34,
                 &cov04    &cov14    &cov24    &cov34  &v44};
    parms   =   RandNormal(&nsim*&n, mean, cov);
        u   =   j(&n, 1);
    call randgen(u, "Uniform");
    uu      =   colvec(repeat(T(u), &nsim));
    z       =   sim || id || parms || uu;
    create parameters from z[c={"sim" "id" "b0" "b1" "b2" "b3" "ci" "u"}];
    append from z;
    close parameters;
quit;
/********************************/;
/*Simulate longitudinal profiles*/;
/********************************/;

data longitudinal;
	set parameters;
    retain _zaad1 &seed1 _zaad2 &seed2 _zaad3 &seed3 _zaad4 &seed4 _zaad5 &seed5;
    ti0 = 0;
/*  ci = sqrt(&v44) * rannor(_zaad4);*/
    do period = 1 to &t;
        retain yij;
        retain zij;
    	y_prev = yij;
        z_prev = zij;
        tij = period + ranuni(_zaad2);
    	start = lag1(tij);
        if (period = 1) then do;
      		start = 0;
            zij = 0;
      		z_prev = 0;
      		eij = sqrt(&ve * exp(ci)) * rannor(_zaad1);
      		yij = &u0 + b0 + (&u1 + b1) * tij + eij;
            ti0 = start;
        end;
        else do;
			
			if (z_prev = 0) then do;
        		pij = logistic(&alpha0 + &alpha1 * yij);
        		if (pij = 1) then zij = 1;
        		else if (pij = 0) then zij = 0;
        		else zij = ranbin(_zaad4, 1, pij);
      		end;
      		else zij = 1;
        	if (zij - z_prev = 1) then ti0 =tij; *;
			*eij = sqrt(&ve * exp(&gamma * z_prev + ci)) * rannor(_zaad1);
        	*yij = &u0 + b0 + (&u1 + b1) * start + z_prev * (&u2 + b2 + (&u3 + b3) * (start - ti0)) + eij;
			eij = sqrt(&ve * exp(&gamma * z_prev + ci)) * rannor(_zaad1);
        	yij = &u0 + b0 + (&u1 + b1) * tij + z_prev * (&u2 + b2 + (&u3 + b3) * (tij - ti0)) + eij;
        end; 
    output; 
    end;
run;
/************************/;
/*Simulate survival data*/;
/************************/;
data longitudinal;
set longitudinal;
where sim>=&startsim and sim<=&endsim;
m_log_u = -log(u);
run;
proc iml;
	/*functions used to calculate cumulative hazard*/
	start integrand(t) global(k, a, b);
		res = a*t**(k-1)*exp(b*t);
		return(res);
	finish;
	start Tevent(t) global(t_start, t_end, k,a,b , cum_start , mlogu);  *cum_start + interval_cum_haz = -logU to find Tevent;
		call quad(int, "integrand", t_start||t);
		return( int+cum_start-mlogu);
	finish;

	use longitudinal;
		read all var {b0 b1 b2 b3 ci start tij ti0 z_prev id sim period m_log_u};
	close longitudinal;

	mu_c = &u0 +b0 +  (&u2+b2)#z_prev -ti0#z_prev#(&u3+b3); /*constant part of mu*/
	mu_t = &u1 +b1 + (&u3+b3)#z_prev; /*time dependent part of mu;  mu(t) = mu_c + mu_t*t*/
	logsig = log(&ve) + &gamma#z_prev+ci; /*log sigma*/

	exp_t = mu_t*&s0; /*time dependent part inside the exponential*/
	constant = &shape * (&scale**(-&shape)); /*constant of the baseline hazard*/
	exp_c = exp(mu_c*&s0 +logsig*&s1)*constant; /*constant of the hazard ( h = exp_c * t^(k-1)*exp(exp_t*t)*/


	x = j(nrow(id),1);
	event = j(nrow(id),1);
	individual = 0;

	do i=1 to nrow(id);
		k = &shape;
		a = exp_c[i];
		b = exp_t[i];
		
		t_start =start[i];
		t_end =tij[i];
		call quad(temp,"integrand",t_start||t_end); /*calculate the increase in cumulative hazard on the interval t(i,j-1) to t(i,j)*/
		if(individual=id[i] & i>1)then do;/*if individual is not ID[i,1] then this is the first measurement of this individual*/
			x[i] = x[i-1]+temp;
			event[i] = (x[i]>m_log_u[i] & x[i-1]<m_log_u[i]);/*is this the first time the cumsum is above mlogu*/
			if(event[i-1]=1) then event[i]=1;
			end;
		else do;/*if first measurement, cumhaz is just the integral, and id is updated*/
			x[i]=temp;
			event[i] = (x[i]>m_log_u[i]);
			
			individual =id[i];
			end;
	end;
	RES = sim||id||period||m_log_u||x||event||start||tij||exp_c||exp_t;*||tnew;
	create survival2 from RES[c={"sim" "id" "period" "m_log_u" "cumsum" "event" "start" "tij" "exp_c" "exp_t"}];* "cumsum" "event" "st"}];
	append from RES;
	close survival2;
quit;

proc iml;
start integrand(t) global(k, a, b);
		res = a*t**(k-1)*exp(b*t);
		return(res);
	finish;
	start Tevent(t) global(t_start, t_end, k,a,b , cum_start , mlogu);  *cum_start + interval_cum_haz = -logU to find Tevent;
		call quad(int, "integrand", t_start||t);
		return( int+cum_start-mlogu);
	finish;

use survival2;
read all var{sim id period m_log_u cumsum event start tij exp_c exp_t};
close survival2;

tnew = j(nrow(id),1);
do i=1 to nrow(id);
		if(event[i]=1)then do;
			k = &shape;
			a = exp_c[i];
			b = exp_t[i];
		
			t_start =start[i];
			t_end =tij[i];
		
			if(period[i]>1) then cum_start=cumsum[i-1];
			else cum_start=0;
			mlogu = m_log_u[i];
			tnew[i] = froot("Tevent",t_start||t_end);/*if the event happens when does it actually happen?*/
		end;
		else tnew[i] = 0;
end;
	RES = sim||id||period||m_log_u||cumsum||event||tnew;
	create survival2 from RES[c={"sim" "id" "period" "mlogu" "cumsum" "event" "st"}];
	append from RES;
	close survival2;
quit;

data survival;
merge longitudinal survival2;
by sim id period;
cum_pre = lag1(cumsum);
if(period=1) then cum_pre=0;
run;

data simdata;
  set survival;
  time2 = z_prev * (tij - ti0);																					
  if (mlogu < cum_pre) then delete; *delete all measurement taken after event;
  if (st > 0) then tij=st; *st !=0 only at event time;
  if (st > 0) then yij = &u0 + b0 + (&u1 + b1) * tij + z_prev * (&u2 + b2 + (&u3 + b3) * (tij - ti0)) + eij; *update the yij at event time;
run;

data analysis;
  set simdata;
  keep sim--id period--start event time2 ti0;
run;
%mend simHetJM;
