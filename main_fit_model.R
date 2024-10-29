library("rstan") # observe startup messages

rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPP = "-mtune = corei7")
options(mc.cores = parallel::detectCores())

initf1 <- function(){list(shape  = 1.5,
                          alpha  = c(0, 0.05, 0.1),
                          beta   = c(12, 0.1, -0.3, -0.05),
                          gamma  = 0.5,
                          sigma0 = 1.4)}
#compile stan model
stm <- stan_model(file = '.../JM.stan')

# read in from sas data
bigdata <- read.csv(file = "../simulated_data.csv", header = TRUE)
big_split <- split(bigdata, bigdata$sim)
ndata <- length(big_split)

runtime<-c()
for (i in 1:ndata) {  
  # prepare dataset
  sim_dat <- big_split[[i]]
  N <- dim(sim_dat)[1]
  M <- length(unique(sim_dat$id))
  # s <- as.vector(table(sim_dat$id))
  yij <- sim_dat$yij
  dij <- sim_dat$event
  z_prev <- sim_dat$z_pre
  tij <- sim_dat$tij
  time2 <- sim_dat$time2
  start <- sim_dat$start
  subj <- sim_dat$id

  df <- list(N=N, M=M, yij=yij, dij=dij, z_prev=z_prev, tij=tij,time2=time2,start=start,subj=subj)
  # fit with stan
  fit <- sampling(stm, init = initf1, iter = 1000, chains = 2, pars = c("shape", "alpha", "beta", "gamma", "sigma0", "tau", "L_Omega"),
                  include = TRUE, seed = 19881110, data=df)

  save(fit,file = paste("../stanfit_jm", i, ".RData", sep = ""))
}
