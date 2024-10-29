options mautolocdisplay mautosource 
        sasautos = ("C:\Users\s130398\Documents\JM\code\SAS\macros\" sasautos);
Libname SAS 'C:\Users\s130398\JM\MixedJump\data';
/*nsim = 500, n in 100,500, t=10,50*/
*********************************************************************************************************************;
%simHetJM(nsim = 100, n = 100, t =30, ve = 2, u0 = 12, u1 = 0, u2 =-0.3, u3 = 0,
		   alpha0 = -3, alpha1 = 0.05, gamma = 0.5, s0 = 0.05, s1 = 0.2, shape =1.5, scale =150,
		   v00 = 4, v11 = 0.005, v22 = 1, v33 = 0.015, v44 = 0.3, cov01 = -0.02, cov02 = -0.2, cov03 = -0.05, cov04=0.6,
		   cov12 = -0.015, cov13 = 0.005, cov14=0.01, cov23 = -0.02, cov24=-0.1, cov34=-0.025);

proc export data=work.Simdata
outfile='...\simulated_data.csv'
dbms=csv replace;
putnames=yes;
run;
*********************************************************************************************************************;

/******************************************************/;
/*Calculate percentage of patients that have treatment*/;
/******************************************************/;

data analysis_last;
	set analysis;
	by sim id;
	if LAST.id;
	output;
run;


proc means data=analysis_last N MEAN;
	var event;
run;
proc means data=analysis_last N MEAN;
	var zij;
run;
