### Cumulative Incidence Function Estimation Based on Population-Based Biobank Data
### Malka Gorfine, David M. Zucker, and Shoval Shoham

### Code by David Zucker
### Version of 26 March 2024

### R code for carrying out a simulation study of a new method
### for estimating the cumulative incidence function for disease
### in an illness-death model

main.path = 'C:/Users/owner/Dropbox/WORK/SHARED FOLDERS/NewMarginalSurvival-forEMR-Data/'
code.path = paste0(main.path, 'Codes2/')
sim.path = paste0(main.path, 'Simulation Results/Results15JAN24/')
setwd(code.path)

library(etm)
library(tidyverse)
library(truncdist)
library(survival)
library(beepr)
library(Rcpp)
library(ggplot2)
library(reshape2)
library(ggpubr)

#OPTIONS
options(max.print = 999999)
options(width=250)
set.seed(7903718)

### LEGEND FOR DATA GENERATING SETTINGS ########################################
# T1 = age at disease diagnosis, T1_setting:
# 1: Truncated Weibull (shape=4, scale=115, a=40), a = minimal recruitment age, scale set for age-80 CIF of 7.5%
# 2: Truncated Weibull (shape=4, k=130, a=40), a = minimal recruitment age, scale set for age-80 CIF of 15.0% 
# 3: Weibull (shape=3.5, scale=200) #parameters set for an age-40 CIF of 5% and an age-80 CIF of about 25%
#
# T2 = age at death, T2_setting:
# 0 - The disease does not affect life expectancy
# Types 1 and 2 are settings with death distribution after disease different from
#   the general death distribution  
# Both types are based on setting the expected survival after disease
# For T1_setting = 1 or 2
#  type 1: expected survival = 2.5 years
#  type 2: expected survival = 7.5 years
# For T1_setting = 3
#  type 1: expected survival = 5 years
#  type 2: expected survival = 10 years
#
# R = recruitment age, R_setting:
# 0: no left truncation, set all the observations to be R=L_R (lower bound)
# 1: uniform distribution, between L_R to U_R
# 2: from the UKB recruitment distribution
#
# C = censoring age, C_setting:
# 1: Recruitment age + uniform(11,15) - similar to current UKB data
# 2: Recruitment age + uniform(11,25) - longer follow-up times
#
# L_R=40, U_R=69 are the default lower and upper recruitment ages in UKB
################################################################################

### BEGIN SETTINGS ##

#BASIC SETTINGS
n = 7500                #sample size
nrep = 1000             #number of simulation replications
nresam = 250            #number of replications for bootstrap
covpr = 0.95            #target coverage probability for confidence intervals and band
bwts = 1                #bootstrap weights: 1=normal, 2=exponential, 3=Poisson
trans.opt = 3           #1=no transformation, 2=-log(1-u), 3=0.5*pi-arcsin(sqrt(1-u))
pt.ci = 1               #for pointwise confidence intervals: 1=normal-theory, 2=bootstrap
auxflg = FALSE          #flag for computing auxiliary terms
plot.to.file = TRUE     #flag for printing plots in a file

#DATA GENERATING SETTINGS
scenario = c(3,2,2,2)           #simulation scenario
popcurv = F                     #flag for generating true CIF
npop = 1e6                      #initial size of popn for generating true CIF
rnddat = TRUE                   #flag for rounding
rnd.unit = 1/365.24             #unit for rounding
L_R = 40                        #minimum recruitment age
U_R = 69                        #maximum recruitment age
pr.enter = 0.87                 #probability that a random person has T2 > R
T1_setting = scenario[1]        #setting for age at disease   3 choices: 1, 2, 3
T2_setting = scenario[2]        #setting for RR for death     2 choices: 1, 2
R_setting = scenario[3]         #setting for recruitment age  2 choices: 1, 2
C_setting = scenario[4]         #setting for censoring age    2 choices: 1, 2
ipre = (T1_setting == 3)        #flag for events before L_R

#TIME GRIDS #

#grid over which CIF estimates will be computed
tgrd = seq(40, 80, 1)
if (T1_setting == 3) {tgrd = seq(30, 80, 1)}
ntgrd = length(tgrd)

#grid over which pointwise CI's will be displayed
ixcb1 = 6:ntgrd
#for Scenarios 1 and 2 this is 45 to 80
#for Scenario 3 this is 35 to 80
tgrd1 = tgrd[ixcb1]         
ntgrd1 = length(tgrd1)

#grid over which confidence band will be computed
ixcb2 = 11:41       #50 to 80
if ((T2_setting == 2) & (C_setting == 1)) {
  ixcb2 = 11:36     #50 to 75
}
if (T1_setting == 3) {ixcb2 = 6:41}     #35 to 80
tgrd2 = tgrd[ixcb2]         
ntgrd2 = length(tgrd2)

#WEIGHTS FOR COMBO ESTIMATOR
adpwt = FALSE                         #flag for adaptive weights
#cmb.wgt = 0.5 - 0.2*(tgrd-40)/40     #weight on new estimator in combo estimate (NA for adaptive)
cmb.wgt = rep(0.5,ntgrd)              #weight on new estimator in combo estimate (NA for adaptive)

### END SETTINGS ###

#LABELS AND FILENAMES
scenario.pr = sum(scenario*(10^(3:0)))
label0 = paste0('Scenario ',scenario.pr)
label0a = paste0('Scenario_',scenario.pr)
label = paste0(label0, ', n=', n)
ciftru.file.name = paste0('ciftru-',label0a,'.csv')     #filename for true CIF 
plot.pdf.name = paste0(sim.path,label0a, '-n=',         #filename for plots
  n, '-',  Sys.Date(), '-1.pdf')                        #filename for plots continued

#PRINT OUT MAIN SETTINGS
print(cbind(n,T1_setting,T2_setting,R_setting,C_setting))
cat('\n')
settings = c('n','T1_setting','T2_setting', 'R_setting', 'C_setting')

#flag for computing the combination estimator
#this estimator is relevant only for the case where there no disease events
#prior to the lowest recruitment age L_R
cmb.flg = !ipre

ibreak = FALSE

### MAIN ROUTINE FOR SIMULATING DATA ###########################################

#RETURNS A DATAFRAME WITH VARIABLES: T1, T2, R, C

load("death_group1.Rdata")  #UK Office of National Statistics age at death distribution
load("R_group1.Rdata")      #Age at recruitment distribution in UK Biobank

sim_data <- function(nn,T1_setting,T2_setting,R_setting,C_setting,L_R=40,U_R=69) {

  # Generate T1 = age at disease
  if (T1_setting==1) {T1<- rtrunc(nn, spec="weibull", shape = 4, scale=140, a = 40)} 
  if (T1_setting==2) {T1<- rtrunc(nn, spec="weibull", shape = 4, scale=116, a = 40)} 
  #if (T1_setting==3) {T1<- rweibull(nn, shape=3.5, scale=200)}
  if (T1_setting==3) {T1<- 20 + rweibull(nn, shape=0.85, scale=170)}

  # Generate T2 = age at death based on the UK Office of National Statistics
  dth_ch = death_group$ch #this the cuml hazard for death
  dth_ch_fun = approxfun(0:105, dth_ch, rule=2)
  dth_inv_ch = approxfun(dth_ch, 0:105, rule=2)
  eee = rexp(nn)
  T2 = dth_inv_ch(eee)

  data <- data.frame(cbind(T1,T2,eee))
  data <- data %>% 
    mutate(delta1 = T1 <= T2) %>% 
    arrange(delta1)
  n_disease <- sum(data$delta1) 
  
  #for those with the disease, a new death time is sampled based on a death distribution after diagnosis
  #we have two choices for expected survival after disease
  if (T2_setting > 0) {
    temp_data <- filter(data, delta1==TRUE) 
    dur.shape = 4
    dur.es.vec = c(2.5,7.5)
    dur.es = dur.es.vec[T2_setting]
    dur.scale = dur.es / gamma(1+1/dur.shape)
    ddur = rweibull(n_disease, shape=dur.shape, scale=dur.scale)
    temp_data$T2 = temp_data$T1 + ddur
    dshift = ifelse(temp_data$T1<40,5,0)
    temp_data$T2 = temp_data$T2 + dshift
    data <- rbind(data[1:(nn-n_disease),],temp_data)
  }  
  
  #Generate R = recruitment age
  if (R_setting==0) data$R = L_R
  if (R_setting==1) data$R = runif(nn,L_R,U_R)  #uniform, independent of disease status
  if (R_setting==2) {
    #independent of disease status - based on the marginal recruitment age distribution of UKB
    rct_ch = R_group$ch
    rct_inv_ch = approxfun(rct_ch, 0:29, rule=2)
    rrr = rexp(nn)
    data$R = 40 + rct_inv_ch(rrr)
  }
 
  # Generate C = age at censoring
  if (C_setting == 0) data$C <- rep(Inf,nn) 
  if (C_setting == 1) data$C <- data$R + runif(nn,11,15)  #as in UKB   
  if (C_setting == 2) data$C <- data$R + runif(nn,11,25)  #longer follow-up
  
  data <- data %>% select(-delta1)
  
  data = data[,-3]
  
  return(data)
  
}
### END OF MAIN DATA SIMULATION ROUTINE ########################################

### ADDITIONAL DATA GENERATION CODE ############################################

#FUNCTION FOR ROUNDING
my.round = function(a,unit) {
  return(unit*round(a/unit))
}  

#GENERATE TRUE CIF CURVE IF REQUESTED
if (popcurv) {
  
  #GENERATE A LARGE POPULATION  
  popdata0 <- sim_data(nn=npop, T1_setting, T2_setting, R_setting = 0,
    C_setting = 0, L_R, U_R)

  #REMOVE SUBJECTS WHO DIED BEFORE L_R
  popdata1 = filter(popdata0, T2 > L_R)
  npd1 = nrow(popdata1)
  
  #TRUE CIF COMPUTATION - BEFORE ROUNDING
  #THE CIF IS CONDITIONAL ON T2>L_R
  ciftru = rep(0,100)
  for (ig in 1:100) {
    ciftru[ig] = mean((popdata1$T1 <= ig) & (popdata1$T1 <= popdata1$T2))
  }
  
  #TRUE CIF COMPUTATION FOR AJ ESTIMATOR
  #THE CIF IS CONDITIONAL ON T1>L_R AND T2>T_R
  ciftru.aj = ciftru
  if (ipre) {
    aj.select = (popdata1$T1 > L_R) 
    popdata1a = filter(popdata1, aj.select)
    for (ig in 1:100) {
      ciftru.aj[ig] = mean((popdata1a$T1 <= ig) & (popdata1a$T1 <= popdata1a$T2))
    }
    rm(popdata1a)
  }
   
  pr.enter <<- mean(popdata0$T2>popdata0$R)
  print(noquote(paste0('pr.enter = ', pr.enter)))
  
  ciftru.fin = cbind(ciftru,ciftru.aj)
  write.csv(ciftru.fin,ciftru.file.name)

  rm(popdata0)
  rm(popdata1)
  stop('Population curve has been generated')

}
#END OF CODE FOR GENERATING THE TRUE CIF CURVE

#IF USING EXISTING DATA FILES FOR THE TRUE CIF, READ THEM IN
if (!popcurv) {
  ciftru.dat = read.csv(ciftru.file.name)
  ciftru = ciftru.dat[,2]
  ciftru.aj = ciftru.dat[,3]
}

#FUNCTION FOR GENERATING A NEW SAMPLE
samgen = function(nsam, pr.enter, T1_setting, T2_setting, R_setting,
  C_setting, L_R, U_R, round.flg) {

  #GENERATE THE SAMPLE
  ngen = round(1.30 * nsam / pr.enter)  
  sam.ini = sim_data(ngen, T1_setting, T2_setting, R_setting, 
    C_setting, L_R, U_R)
  sam.fltr = filter(sam.ini, T2 > R)
  select1 = sample.int(nrow(sam.fltr), nsam)
  sam.fltr = sam.fltr[select1,]
  
  #ROUND DATA IF DESIRED
  sam.fin = sam.fltr
  if (round.flg) {
   adjust1 = 0.005*rnd.unit
   adjust2 = 0.01*rnd.unit
   sam.fin = my.round(sam.fltr,rnd.unit)
   for (i in 1:nsam) {
     if (sam.fin$R[i] == sam.fin$T2[i]) {
       if (sam.fltr$R[i] < sam.fltr$T2[i]) {
         sam.fin$R[i] = sam.fin$R[i] - adjust2
       }
       else {
         sam.fin$R[i] = sam.fin$R[i] + adjust1
       }
     }
     if (sam.fin$C[i] == sam.fin$T1[i]) {
       if (sam.fltr$C[i] > sam.fltr$T1[i]) {
         sam.fin$C[i] = sam.fin$C[i] + adjust2
       }
       else {
          sam.fin$C[i] = sam.fin$C[i] - adjust2
        }
      }
      if (sam.fin$C[i] == sam.fin$T2[i]) {
        if (sam.fltr$C[i] > sam.fltr$T2[i]) {
          sam.fin$C[i] = sam.fin$C[i] + adjust2
        }
        else {
          sam.fin$C[i] = sam.fin$C[i] - adjust2
        }
      }
      if (sam.fin$T1[i] == sam.fin$T2[i]) {
        if (sam.fltr$T2[i] > sam.fltr$T1[i]) {
          sam.fin$T2[i] = sam.fin$T2[i] + adjust1
        }
        else {
          sam.fin$T2[i] = sam.fin$T2[i] - adjust1
        }
      }
   }
  }
  
  #MAIN "AGE AT" VARIABLES
  sam.fin$age_recr = sam.fin$R
  sam.fin$age_diag = ifelse(sam.fin$T1 <= sam.fin$T2, sam.fin$T1, Inf)
  sam.fin$age_diag = ifelse(sam.fin$T1 <= sam.fin$C, sam.fin$age_diag, Inf)
  sam.fin$age_death = ifelse(sam.fin$T2 <= sam.fin$C, sam.fin$T2, Inf)
  sam.fin$age_end_fu = ifelse(sam.fin$T2 <= sam.fin$C, sam.fin$T2, sam.fin$C)

  #FINAL STATUS
  ix.1 = which((sam.fin$T2 < sam.fin$T1) & (sam.fin$T2 <= sam.fin$C))
  ix.2 = which((sam.fin$age_diag < Inf) & (sam.fin$C < sam.fin$T2))
  ix.3 = which((sam.fin$age_diag < Inf) & (sam.fin$T2 <= sam.fin$C))
  ix.0 = (1:nsam)[-c(ix.1,ix.2,ix.3)]
  status_end = rep(NA,nsam)
  status_end[ix.0] = 0
  status_end[ix.1] = 1 
  status_end[ix.2] = 2 
  status_end[ix.3] = 3 
  stat.end = rep(NA,nsam)
  stat.end[ix.0] = '0_censored'
  stat.end[ix.1] = '1_died_without_disease'
  stat.end[ix.2] = '2_alive_with_disease'
  stat.end[ix.3] = '3_died_with_disease'
  sam.fin$status_end = status_end
  sam.fin$stat.end = stat.end
  return(sam.fin)
  
} #END OF SAMPLE GENERATION

### END OF DATA GENERATION CODE ################################################

### FUNCTIONS ##################################################################

trans = function(u,t.opt) {
  if (t.opt == 1) ans = u
  if (t.opt == 2) ans = -log(1-u)
  if (t.opt == 3) ans = 0.5*pi - asin(sqrt(1-u))
  nnn = sum(is.nan(ans))
  if (nnn>0) {ibreak <<- TRUE}
  return(ans)
}

trans.deriv = function(u,t.opt) {
  didl = 1e-8
  if (t.opt == 1) ans = rep(1,length(u))
  if (t.opt == 2) ans = 1/(1-u)
  if (t.opt == 3) ans = 0.5 / (sqrt(u+didl)*sqrt(1-u+didl))
  return(ans)
}  

trans.inv = function(v,t.opt) {
  if (t.opt == 1) ans = v
  if (t.opt == 2) ans = 1-exp(-v)
  if (t.opt == 3) ans = 1 - (sin(0.5*pi-v)^2)
  return(ans)
}  

#KAPLAN-MEIER FUNCTION
KM = function(trecr,tfu,del,wt) {
  
  #IDENTIFY UNIQUE EVENT TIMES
  evtim = tfu[which(del==1)]
  time = unique(evtim)
  time = sort(time)
  ndist = length(time)

  #CALCULATIONS
  atrskvec = rep(0,ndist)
  km_surv = rep(0,ndist)
  skmcur = 1
  for (m in 1:ndist) {
    cur.tim = time[m] 
    d = sum(wt*(trecr < cur.tim)*(tfu==cur.tim)*del)
    atrsk = sum(wt*(trecr < cur.tim)*(tfu >= cur.tim))
    p = d/atrsk
    skmcur = skmcur * (1-p)
    atrskvec[m] = atrsk
    km_surv[m] = skmcur
  }  
  ans = list(time=time, surv=km_surv, atrsk=atrskvec)
  
  return(ans)

}
#END KAPLAN-MEIER FUNCTION

#C BACKEND FOR AJ ESTIMATOR
cppFunction('NumericVector x_aux_aj(
  int nsam, int n_aj, int nfst, NumericVector noncens, NumericVector prev, 
  NumericVector case_vec, NumericVector V1, NumericVector km_first_times,
  NumericVector ycal_circ, NumericVector ycal_circ_recip,
  NumericVector ycal_star, NumericVector ycal_star_recip,
  NumericVector sfac, NumericVector age_recr, int ntgrd, NumericVector tgrd,
  NumericVector age_first, NumericVector afun, NumericMatrix at_rsk_i,
  IntegerVector ixf_vec, int ncase, IntegerVector ix_case) {
	
NumericMatrix xa_mat(ntgrd,nsam);

double n_aj_d = n_aj;
double nsam_d = nsam;
double ppi = n_aj_d/nsam_d;

for (int it = 0; it < ntgrd; it++) {
  for (int ii = 0; ii < nsam; ii++) {
    xa_mat(it,ii) = 0;
  }
}  

for (int ic = 0; ic < ncase; ic++) {
  int ij = ix_case(ic) - 1;
  int ixfv_ij = ixf_vec(ij) - 1;
  double atrsk_cur = ycal_circ(ixfv_ij); 
	double atrsk_recip_cur = ycal_circ_recip(ixfv_ij);
  for (int ii = 0; ii < nsam; ii++) {
    if (prev(ii) == 0) {
      double trm1 = 0;
      double trm2 = 0;
      double trm3 = 0;
      double V3i = age_first(ii);
      if (noncens(ii) == 1) {
        int ixfv_ii = ixf_vec(ii) - 1;
        if (V3i < V1(ij)) {
           trm1 = sfac(ixfv_ii) * ycal_star_recip(ixfv_ii);
        }
      }
      for (int ixf = 0; ixf < nfst; ixf++) {
        double ftime = km_first_times[ixf];
        if ((ftime >= age_recr(ii)) && (ftime <= V3i) && (ftime < V1(ij))) {
          trm2 += sfac(ixf) * pow(ycal_star_recip(ixf),2);
        }
      }  
      trm2 = trm2 / n_aj;
      trm3 = atrsk_recip_cur*(at_rsk_i(ii,ixfv_ij) - atrsk_cur);
    	for (int it = 0; it < ntgrd; it++) {
      	if (V1(ij) <= tgrd(it)) {
	        xa_mat(it,ii) += -afun(ixfv_ij)*((trm1 - trm2)/ppi + trm3) / nsam;
      	}
    	}	
    }
  }    
}

return(xa_mat);

}')
#END C BACKEND FOR AJ ESTIMATOR
        
#FUNCTION FOR CIF COMPUTATION FOR AALEN-JOHANSEN ESTIMATOR
cifcmp.aj = function(age_recr, age_diag, age_death, age_end_fu, status_end, tgrd, wts) {
  
  #status_end
  #0 = alive without disease (i.e. censored)  
  #1 = died without disease
  #2 = alive with disease
  #3 = died with disease 
  
  #tgrd = grid of timepoints over which the CIF will be computed
  ntgrd = length(tgrd)
  
  #sample size  
  nsam = length(age_recr)

  #SETUPS
  noncens = (status_end != 0)
  ix.cen = which(status_end == 0)
  diseased = (status_end >= 2)
  prev = ((age_diag <= age_recr) & diseased)
  ix.prev = which(prev)
  case = (age_diag <= age_death) & (!prev) & (noncens)
  ix.case = which(case)
  ncase = length(ix.case)
  age_diag[ix.cen] = Inf
  age_death[ix.cen] = Inf
  cifest = rep(NA,ntgrd)
  cifest.sd1 = rep(NA,ntgrd)
  cifest.sd2 = rep(NA,ntgrd)
  age_first = pmin(age_diag, age_death, age_end_fu)
  
  #REMOVE PREVALENT CASES FOR COMPUTATION OF KM FOR FIRST TRANSITION
  nprev = length(ix.prev)
  n.aj = nsam - length(ix.prev)
  age_recr.aj = age_recr
  age_diag.aj = age_diag
  age_death.aj = age_death
  age_end_fu.aj = age_end_fu
  status_end.aj = status_end
  age_first.aj = age_first
  wts.aj = wts
  noncens.aj = noncens
  if (nprev > 0) {
    age_recr.aj = age_recr[-ix.prev]
    age_diag.aj = age_diag[-ix.prev]
    age_death.aj = age_death[-ix.prev]
    age_end_fu.aj = age_end_fu[-ix.prev]
    status_end.aj = status_end[-ix.prev]
    age_first.aj = age_first[-ix.prev]
    wts.aj = wts[-ix.prev]
    noncens.aj = noncens[-ix.prev]
  }  

  #KAPLAN-MEIER ESTIMATE FOR TIME TO FIRST TRANSITION
  #OMITTING PREVALENT CASES
  km_fst = KM(age_recr.aj, age_first.aj, noncens.aj, wts.aj)
  km_fst_times = km_fst$time
  km_fst_surv = km_fst$surv
  km_fst_surv1 = c(1,km_fst_surv)
  km_fst_at_rsk = km_fst$atrsk
  ycal.circ = km_fst_at_rsk/nsam
  ycal.star = km_fst_at_rsk/n.aj
  ycal.circ.recip = 1/ycal.circ
  ycal.star.recip = 1/ycal.star
  nfst = length(km_fst_times)
  sfac = km_fst_surv1[1:nfst]/km_fst_surv
  sfac[which(is.infinite(sfac))] = 0
  sfac[which(is.nan(sfac))] = 0
  afun = km_fst_surv1[1:nfst]/ycal.circ
  
  #PRELIMINARIES
  at_rsk_i = matrix(0,nsam,nfst)
  ixf_vec = rep(0,nsam)
  for (ix.f in 1:nfst) {
    icur = which((age_first==km_fst_times[ix.f]) & (!prev) & (noncens))
    ixf_vec[icur] = ix.f
    icur1 = which((age_recr < km_fst_times[ix.f]) &
      (age_first >= km_fst_times[ix.f]))
    at_rsk_i[icur1,ix.f] = 1
  }

  #CIF COMPUTATION
  x.main = matrix(0,ntgrd,nsam) 
  for (ix.g in 1:ntgrd) {
    x.main[ix.g,ix.case] = afun[ixf_vec[ix.case]] * (age_diag[ix.case] <= tgrd[ix.g])
    cifest[ix.g] = sum(x.main[ix.g,ix.case])/nsam 
    cifest.sd1[ix.g] = sqrt((sum(x.main[ix.g,ix.case]^2)/nsam - cifest[ix.g]^2)/nsam)  
  }
  
  #COMPUTATION OF TERMS GOING INTO IID REPRESENTATION FOR VARIANCE CALCULATION
  if (auxflg) {
    x.aux = x_aux_aj(nsam, n.aj, nfst, noncens, prev, case, age_diag, km_fst_times,
      ycal.circ, ycal.circ.recip, ycal.star, ycal.star.recip, sfac, age_recr, 
      ntgrd, tgrd, age_first, afun, at_rsk_i, ixf_vec, ncase, ix.case)
  }
  
  #WRAP UP IID REPRESENTATION AND VARIANCE CALCULATION
  eps = matrix(0,ntgrd,nsam)
  if (auxflg) {
    for (ix.g in 1:ntgrd) {
      eps[ix.g,] = x.main[ix.g,] + x.aux[ix.g,] - cifest[ix.g]
      cifest.sd2[ix.g] = sqrt(mean(eps[ix.g,]^2)/nsam)
    }  
  } 
  else {
    for (ix.g in 1:ntgrd) {
      eps[ix.g,] = x.main[ix.g,] - cifest[ix.g]
    }
    cifest.sd2 = cifest.sd1
  } 
 
  #COMPUTE ESTIMATOR USING etmCIF AND EXTRACT VARIANCE
  fstatus = rep(0,n.aj)
  fstatus[which(status_end.aj==1)] = 2
  fstatus[which(status_end.aj==2)] = 1
  fstatus[which(status_end.aj==3)] = 1
  mydat = data.frame(start=age_recr.aj,finish=age_first.aj,fstatus=fstatus)
  ajr = try(etmCIF(Surv(start,finish,fstatus != 0) ~ 1, data=mydat, etype=fstatus, failcode=1))
  ajr1 = summary(ajr, ci.fun='linear')
  ajr1 = ajr1[[1]]$'CIF 1'
  if (is.null(ajr1)) {
    aj_est_fin = cifest
    cifest.sd3=cifest.sd2
  }
  else {
  aj.time = ajr1$time
  aj.est = ajr1$P
  aj.sd = sqrt(ajr1$var)

  #ADAPT TO TIME GRID
  aj_est_fin = NULL
  cifest.sd3 = NULL
  ajcur = 0
  sdcur = 0
  for (ig in 1:ntgrd) {
    ixt = which(aj.time <= tgrd[ig])
    if (length(ixt)>0) {
      ixt1 = max(ixt)
      ajcur = aj.est[ixt1]
      sdcur = aj.sd[ixt1]
    }
    aj_est_fin = c(aj_est_fin, ajcur)
    cifest.sd3 = c(cifest.sd3, sdcur)
  }
  }
  
  #RETURN RESULT
  ans = list(cifest=cifest, cifest.sd1=cifest.sd1, cifest.sd2=cifest.sd2,
    cifest.sd3=cifest.sd3, eps=eps)
  return(ans)
  
}#END OF FUNCTION TO COMPUTE AJ ESTIMATOR

#C BACKEND FOR NEW ESTIMATOR
cppFunction('NumericVector x_aux_new(
  int nsam, int ndth, NumericVector died, NumericVector case_vec, NumericVector V1,
  NumericVector V2, NumericVector km_at_rsk_nrm, NumericVector km_at_rsk_nrm_recip,
  NumericVector sfac, NumericVector age_recr, int ntgrd, NumericVector tgrd, 
  NumericVector dth_vec, NumericVector afun, NumericMatrix at_rsk_i,
  IntegerVector ixd_vec, int ncase, IntegerVector ix_case) {
  
NumericMatrix xa_mat(ntgrd,nsam);

for (int ic = 0; ic < ncase; ic++) {
  int ij = ix_case(ic) - 1;
  int ixdv_ij = ixd_vec(ij) - 1;
  for (int ii = 0; ii < nsam; ii++) {
    double trm1 = 0;
    if ((died(ii) == 1) && (V2(ii) < V2(ij))) {
      int ixdv_ii = ixd_vec(ii) - 1;
      trm1 = sfac(ixdv_ii)*km_at_rsk_nrm_recip(ixdv_ii);
    }  
    double trm2 = 0;
    for (int ixd = 0; ixd < ndth; ixd++) {
      double dtime = dth_vec(ixd);
      if ((dtime >= age_recr(ii)) && (dtime <= V2(ii)) && (dtime < V2(ij))) {
        trm2 += sfac(ixd) * pow(km_at_rsk_nrm_recip(ixd),2);
      }       
    }
    trm2 = trm2 / nsam;
    double atrsk_cur = km_at_rsk_nrm(ixdv_ij); 
    double atrsk_recip_cur = km_at_rsk_nrm_recip(ixdv_ij);
    double trm3 = atrsk_recip_cur*(at_rsk_i(ii,ixdv_ij) - atrsk_cur);
    for (int it = 0; it < ntgrd; it++) {
      if (V1(ij) <= tgrd(it)) {
        xa_mat(it,ii) += -afun(ixdv_ij)*(trm1 - trm2 + trm3) / nsam;
      } 
    }         
  }    
}

return(xa_mat);
}')
#END C BACKEND FOR NEW ESTIMATOR 

#FUNCTION FOR CIF COMPUTATION FOR PROPOSED ESTIMATOR
cifcmp.new = function(age_recr, age_diag, age_death, age_end_fu, 
  status_end, tgrd, wts) {
  
  #status_end
  #0 = alive without disease (i.e. censored)  
  #1 = died without disease
  #2 = alive with disease
  #3 = died with disease 
  
  #sample size  
  nsam = length(age_recr)
  
  #tgrd = grid of timepoints over which the CIF will be computed
  ntgrd = length(tgrd)
  
  #SETUPS
  ix.cen = which(status_end == 0)
  died = (status_end == 1) | (status_end == 3)
  ix.died = which(died)
  case = (status_end == 3)
  ix.case = which(case)
  ncase = length(ix.case)
  cifest = rep(NA,ntgrd)
  cifest.sd1 = rep(NA,ntgrd)
  cifest.sd2 = rep(NA,ntgrd)
  age_diag[ix.cen] = Inf
  age_death[ix.cen] = Inf
  
  #KAPLAN-MEIER ESTIMATE OF DEATH TIME DISTN
  km_dth = KM(age_recr, age_end_fu, died, wts)
  km_dth_times = km_dth$time
  km_dth_surv = km_dth$surv
  km_dth_surv1 = c(1, km_dth_surv)
  km_at_rsk = km_dth$atrsk
  km_at_rsk_nrm = km_at_rsk/nsam   #estimate of script Y2 at death times
  km_at_rsk_nrm_recip = 1/km_at_rsk_nrm 
  ndth = length(km_dth_times)
  sfac = km_dth_surv1[1:ndth]/km_dth_surv
  sfac[which(is.infinite(sfac))] = 0
  sfac[which(is.nan(sfac))] = 0
  intwts.numer = -diff(km_dth_surv1)
  
  #PRELIMINARIES
  ixd_vec = rep(0,nsam)
  at_rsk_i = matrix(0,nsam,ndth)
  afun = rep(0,ndth)
  for (ix.d in 1:ndth) {
    icur = which(age_death == km_dth_times[ix.d])
    ixd_vec[icur] = ix.d
    ndth.cur.nrm = length(icur)/nsam
    afun[ix.d] = intwts.numer[ix.d] / ndth.cur.nrm 
    icur1 = which((age_recr < km_dth_times[ix.d]) &
      (age_death >= km_dth_times[ix.d]))
    at_rsk_i[icur1,ix.d] = 1
  } 
 
  #CIF COMPUTATION
  x.main = matrix(0,ntgrd,nsam) 
  for (ix.g in 1:ntgrd) {
    x.main[ix.g,ix.case] = afun[ixd_vec[ix.case]] * (age_diag[ix.case] <= tgrd[ix.g])
    cifest[ix.g] = sum(x.main[ix.g,ix.case])/nsam 
    cifest.sd1[ix.g] = sqrt((sum(x.main[ix.g,ix.case]^2)/nsam - cifest[ix.g]^2)/nsam)  
  }
  
  #COMPUTATION OF ADDITIONAL TERM GOING INTO IID REPRESENTATION FOR VARIANCE CALCULATION
  if (auxflg) {
    x.aux = x_aux_new(nsam, ndth, died, case, age_diag, age_death, km_at_rsk_nrm,
      km_at_rsk_nrm_recip, sfac, age_recr, ntgrd, tgrd, km_dth_times, afun, 
      at_rsk_i, ixd_vec, ncase, ix.case)
  }  
  
  #WRAP UP IID REPRESENTATION AND VARIANCE CALCULATION
  eps = matrix(0,ntgrd,nsam)
  if (auxflg) {
    for (ix.g in 1:ntgrd) {
      eps[ix.g,] = x.main[ix.g,] + x.aux[ix.g,] - cifest[ix.g]
      cifest.sd2[ix.g] = sqrt(mean(eps[ix.g,]^2)/nsam)
    }  
  } 
  else {
    for (ix.g in 1:ntgrd) {
      eps[ix.g,] = x.main[ix.g,] - cifest[ix.g]
    }
    cifest.sd2 = cifest.sd1
  }  

  #RETURN RESULT
  ans = list(cifest=cifest, cifest.sd1=cifest.sd1, cifest.sd2=cifest.sd2, eps=eps)
  return(ans)
  
}
# END OF FUNCTION TO COMPUTE NEW ESTIMATOR

#FUNCTION COMPUTE ALL ESTIMATORS WITH CONFIDENCE INTERVALS AND BANDS
cifcmp.full = function(age_recr, age_diag, age_death, age_end_fu, status_end, 
  tgrd, tgrd1, tgrd2, covpr, nresam, adpwt, cmb.wgt, t.opt, pt.ci) {

#PRELIMINARIES
zcrit = qnorm(1-(1-covpr)/2)
didl = 1e-8 #a small number
n = length(age_recr)
sq.n = sqrt(n)
ntgrd1 = length(tgrd1)
bwts = 1 #bootstrap weights: 1=normal, 2=exponential, 3=Poisson
ixcb2 = which(tgrd %in% tgrd2)

#AALEN-JOHANSEN ESTMATOR
print(noquote('Computing AJ estimator ...'))
cifest.aj = cifcmp.aj(age_recr, age_diag, age_death, age_end_fu,
  status_end, tgrd, rep(1,n))
aj.est = cifest.aj$cifest
aj.sd1 = cifest.aj$cifest.sd1
aj.sd2 = cifest.aj$cifest.sd2
aj.sd3 = cifest.aj$cifest.sd3
aj.est.tr = trans(aj.est,t.opt)
ajderv = trans.deriv(aj.est,t.opt)
aj.sd1.tr = aj.sd1*ajderv
aj.sd2.tr = aj.sd2*ajderv
aj.sd3.tr = aj.sd3*ajderv
eps.aj = cifest.aj$eps

#NEW ESTIMATOR
print(noquote('Computing new estimator ...'))
cifest.new = cifcmp.new(age_recr, age_diag, age_death, age_end_fu,
  status_end, tgrd, rep(1,n))
new.est = cifest.new$cifest
new.sd1 = cifest.new$cifest.sd1
new.sd2 = cifest.new$cifest.sd2
new.est.tr = trans(new.est,t.opt)
newderv = trans.deriv(new.est,t.opt)
new.sd1.tr = new.sd1*newderv
new.sd2.tr = new.sd2*newderv
eps.new = cifest.new$eps

#COMBINED ESTIMATOR
cmb.est = NULL
cmb.var = NULL
cmb.sd = NULL
cmb.ptwise.ci.width = NULL
cmb.ptwise.ci.lo = NULL
cmb.ptwise.ci.hi = NULL
cmb.wts = NULL
cmb.band.lo = NULL
cmb.band.hi = NULL
if (cmb.flg) {
  var.aj = aj.sd2^2
  var.new = new.sd2^2
  covar = rep(NA,ntgrd)
  for (ix.g in 1:ntgrd) {
    covar[ix.g] = cov(eps.aj[ix.g,],eps.new[ix.g,])/n
  }
  if (adpwt) {
    #adaptive weights
    wnumer = var.aj - covar
    wdenom = var.aj + var.new - 2*covar
    cmb.wts = ifelse(wdenom > 0, wnumer/wdenom, 1)
    cmb.wts = ifelse(cmb.wts < 1, cmb.wts, 1)
  }
  else {
    #fixed weights
    cmb.wts = cmb.wgt
  }
  wts1 = 1 - cmb.wts
  wts2 = cmb.wts
  cmb.est = wts1*aj.est + wts2*new.est
  cmb.var = (wts1^2)*var.aj + (wts2^2)*var.new + 2*wts1*wts2*covar
  cmb.sd = sqrt(pmax(cmb.var,1e-8))
  cmb.est.tr = trans(cmb.est,t.opt)
  cmbderv = trans.deriv(cmb.est,t.opt)
  cmb.sd.tr = cmb.sd*cmbderv
  eps.cmb = wts1*eps.aj + wts2*eps.new
}
  
#BOOTSTRAP RESAMPLING SCHEME
print(noquote('Bootstrap Reps for Confidence Bands ...'))
eps.mean.aj = apply(eps.aj,1,mean)
eps.mean.new = apply(eps.new,1,mean)
if (cmb.flg) {eps.mean.cmb = wts1*eps.mean.aj + wts2*eps.mean.new}  
if (pt.ci == 2) {
  aj.zstat.boot = matrix(NA,nresam,ntgrd)
  new.zstat.boot = matrix(NA,nresam,ntgrd)
  if (cmb.flg) {cmb.zstat.boot = matrix(NA,nresam,ntgrd)}
}
supvec.aj = rep(0,nresam)
supvec.new = rep(0,nresam)
if (cmb.flg) {supvec.cmb = rep(0,nresam)}
if (bwts == 1) {qqmat = matrix(rnorm(n*nresam),nresam,n)}
if (bwts == 2) {qqmat = matrix(rexp(n*nresam,1),nresam,n)}
if (bwts == 3) {qqmat = matrix(rpois(n*nresam,1),nresam,n)}
for (b in 1:nresam) {
  if (b %% 25 == 0) print(b)
  qq = qqmat[b,]
  qq = matrix(qq,1,n) %x% matrix(1,ntgrd,1)
  #AJ
  qe.aj = qq*eps.aj
  qe.mean.aj = apply(qe.aj,1,mean)
  dif.raw.aj = qe.mean.aj - eps.mean.aj
  aj.est.cur = aj.est + dif.raw.aj
  aj.est.cur = ifelse(aj.est.cur>0, aj.est.cur, 0)
  aj.est.cur = ifelse(aj.est.cur<1, aj.est.cur, 1)
  dif.aj = (trans(aj.est.cur,t.opt) - aj.est.tr) / 
    (aj.sd3*trans.deriv(aj.est.cur,t.opt)+didl)
  if (pt.ci == 2) {aj.zstat.boot[b,] = dif.aj}
  supvec.aj[b] = max(abs(dif.aj[ixcb2]))
  #NEW
  qe.new = qq*eps.new
  qe.mean.new = apply(qe.new,1,mean)
  dif.raw.new = qe.mean.new - eps.mean.new
  new.est.cur = new.est + dif.raw.new
  new.est.cur = ifelse(new.est.cur>0, new.est.cur, 0)
  new.est.cur = ifelse(new.est.cur<1, new.est.cur, 1)
  dif.new = (trans(new.est.cur,t.opt) - new.est.tr) / 
    (new.sd2*trans.deriv(new.est.cur,t.opt)+didl)
  if (pt.ci == 2) {new.zstat.boot[b,] = dif.new}
  supvec.new[b] = max(abs(dif.new[ixcb2]))
  #CMB
  if (cmb.flg) {
    qe.cmb = qq*eps.cmb
    qe.mean.cmb = apply(qe.cmb,1,mean)
    dif.raw.cmb = qe.mean.cmb - eps.mean.cmb
    cmb.est.cur = cmb.est + dif.raw.cmb
    cmb.est.cur = ifelse(cmb.est.cur>0, cmb.est.cur, 0)
    cmb.est.cur = ifelse(cmb.est.cur<1, cmb.est.cur, 1)
    dif.cmb = (trans(cmb.est.cur,t.opt) - cmb.est.tr) / 
      (cmb.sd*trans.deriv(cmb.est.cur,t.opt)+didl)
    if (pt.ci == 2) {cmb.zstat.boot[b,] = dif.cmb}
    supvec.cmb[b] = max(abs(dif.cmb[ixcb2]))
  }
#if (ibreak) {
#  beep()
#  browser()
#}    
}

#POINTWISE CONFIDENCE INTERVALS
if (pt.ci == 1) {
  aj.pt.crit = zcrit
  new.pt.crit = zcrit
  if (cmb.flg) {cmb.pt.crit = zcrit}
}
else {
  aj.pt.crit = apply(abs(aj.zstat.boot), 2, quantile, probs=covpr)
  new.pt.crit = apply(abs(aj.zstat.boot), 2, quantile, probs=covpr)
  if (cmb.flg) {cmb.pt.crit = apply(abs(cmb.zstat.boot), 2, quantile, probs=covpr)}
}
#AJ
aj.ptwise.ci.hwd.tr.1 = aj.pt.crit*aj.sd1.tr
aj.ptwise.ci.lo.1 = trans.inv(aj.est.tr - aj.ptwise.ci.hwd.tr.1, t.opt)
aj.ptwise.ci.hi.1 = trans.inv(aj.est.tr + aj.ptwise.ci.hwd.tr.1, t.opt)
aj.ptwise.ci.width.1 = aj.ptwise.ci.hi.1 - aj.ptwise.ci.lo.1
aj.ptwise.ci.hwd.tr.2 = aj.pt.crit*aj.sd2.tr
aj.ptwise.ci.lo.2 = trans.inv(aj.est.tr - aj.ptwise.ci.hwd.tr.1, t.opt)
aj.ptwise.ci.hi.2 = trans.inv(aj.est.tr + aj.ptwise.ci.hwd.tr.1, t.opt)
aj.ptwise.ci.width.2 = aj.ptwise.ci.hi.2 - aj.ptwise.ci.lo.2
aj.ptwise.ci.hwd.tr.3 = aj.pt.crit*aj.sd3.tr
aj.ptwise.ci.lo.3 = trans.inv(aj.est.tr - aj.ptwise.ci.hwd.tr.3, t.opt)
aj.ptwise.ci.hi.3 = trans.inv(aj.est.tr + aj.ptwise.ci.hwd.tr.3, t.opt)
aj.ptwise.ci.width.3 = aj.ptwise.ci.hi.3 - aj.ptwise.ci.lo.3
#NEW
new.ptwise.ci.hwd.tr = new.pt.crit*new.sd2.tr
new.ptwise.ci.lo = trans.inv(new.est.tr - new.ptwise.ci.hwd.tr, t.opt)
new.ptwise.ci.hi = trans.inv(new.est.tr + new.ptwise.ci.hwd.tr, t.opt)
new.ptwise.ci.width = new.ptwise.ci.hi - new.ptwise.ci.lo
#COMBO
if (cmb.flg) {
  cmb.ptwise.ci.hwd.tr = cmb.pt.crit*cmb.sd.tr
  cmb.ptwise.ci.lo = trans.inv(cmb.est.tr - cmb.ptwise.ci.hwd.tr, t.opt)
  cmb.ptwise.ci.hi = trans.inv(cmb.est.tr + cmb.ptwise.ci.hwd.tr, t.opt)
  cmb.ptwise.ci.width = cmb.ptwise.ci.hi - cmb.ptwise.ci.lo
}

#CONFIDENCE BANDS
#AJ
aj.band.crit = quantile(supvec.aj,probs=covpr,na.rm=T)
aj.band.hwd.tr = aj.band.crit*aj.sd3[ixcb2]*ajderv[ixcb2]
aj.band.lo = trans.inv(aj.est.tr[ixcb2] - aj.band.hwd.tr, t.opt)
aj.band.hi = trans.inv(aj.est.tr[ixcb2] + aj.band.hwd.tr, t.opt)
aj.band.width = aj.band.hi - aj.band.lo
#NEW
new.band.crit = quantile(supvec.new,probs=covpr,na.rm=T)
new.band.hwd.tr = new.band.crit*new.sd2[ixcb2]*newderv[ixcb2]
new.band.lo = trans.inv(new.est.tr[ixcb2] - new.band.hwd.tr, t.opt)
new.band.hi = trans.inv(new.est.tr[ixcb2] + new.band.hwd.tr, t.opt)
new.band.width = new.band.hi - new.band.lo
#CMB
if (cmb.flg) {
  cmb.band.crit = quantile(supvec.cmb,probs=covpr,na.rm=T)
  cmb.band.hwd.tr = cmb.band.crit*cmb.sd[ixcb2]*cmbderv[ixcb2]
  cmb.band.lo = trans.inv(cmb.est.tr[ixcb2] - cmb.band.hwd.tr, t.opt)
  cmb.band.hi = trans.inv(cmb.est.tr[ixcb2] + cmb.band.hwd.tr, t.opt)
  cmb.band.width = cmb.band.hi - cmb.band.lo
}  
else {
  cmb.band.width = 0
  cmb.band.lo = 0
  cmb.band.hi = 0
}

ans = list(
  aj.est = aj.est,
  aj.sd1 = aj.sd1,
  aj.sd2 = aj.sd2,
  aj.sd3 = aj.sd3,
  aj.ptwise.ci.width.1 = aj.ptwise.ci.width.1,
  aj.ptwise.ci.lo.1 = aj.ptwise.ci.lo.1,
  aj.ptwise.ci.hi.1 = aj.ptwise.ci.hi.1, 
  aj.ptwise.ci.width.2 = aj.ptwise.ci.width.2,
  aj.ptwise.ci.lo.2 = aj.ptwise.ci.lo.2,
  aj.ptwise.ci.hi.2 = aj.ptwise.ci.hi.2, 
  aj.ptwise.ci.width.3 = aj.ptwise.ci.width.3,
  aj.ptwise.ci.lo.3 = aj.ptwise.ci.lo.3,
  aj.ptwise.ci.hi.3 = aj.ptwise.ci.hi.3, 
  new.est = new.est,
  new.sd1 = new.sd1,
  new.sd2 = new.sd2,
  new.ptwise.ci.width = new.ptwise.ci.width,
  new.ptwise.ci.lo = new.ptwise.ci.lo,
  new.ptwise.ci.hi = new.ptwise.ci.hi,
  cmb.est = cmb.est,
  cmb.sd1 = cmb.sd,
  cmb.sd2 = cmb.sd,
  cmb.wts = cmb.wts,
  cmb.ptwise.ci.width = cmb.ptwise.ci.width,
  cmb.ptwise.ci.lo = cmb.ptwise.ci.lo,
  cmb.ptwise.ci.hi = cmb.ptwise.ci.hi,
  aj.band.width = aj.band.width,
  aj.band.lo = aj.band.lo,
  aj.band.hi = aj.band.hi,
  new.band.width = new.band.width,
  new.band.lo = new.band.lo,
  new.band.hi = new.band.hi,
  cmb.band.width = cmb.band.width,
  cmb.band.lo = cmb.band.lo,
  cmb.band.hi = cmb.band.hi)
#ALTOGETHER 35 OUTPUT ITEMS

return(ans)

}
#END OF FUNCTION COMPUTE ALL ESTIMATORS WITH CONFIDENCE INTERVALS AND BANDS

### MAIN PROGRAM ###############################################################

#SET UP MATRICES TO STORE SIMULATION RESULTS
cif.res.mat.aj.est = matrix(NA,nrep,ntgrd)
cif.res.mat.aj.sd1 = matrix(NA,nrep,ntgrd)
cif.res.mat.aj.sd2 = matrix(NA,nrep,ntgrd)
cif.res.mat.aj.sd3 = matrix(NA,nrep,ntgrd)
cif.res.mat.aj.width = matrix(NA,nrep,ntgrd)
cif.res.mat.aj.band.width = matrix(NA,nrep,ntgrd2)
cif.res.mat.aj.cvr.0 = matrix(NA,nrep,ntgrd)
cif.res.mat.aj.cvr.1 = matrix(NA,nrep,ntgrd)
cif.res.mat.aj.cvr.2 = matrix(NA,nrep,ntgrd)
cif.res.mat.aj.cvr.3 = matrix(NA,nrep,ntgrd)
cif.res.mat.new.est = matrix(NA,nrep,ntgrd)
cif.res.mat.new.sd1 = matrix(NA,nrep,ntgrd)
cif.res.mat.new.sd2 = matrix(NA,nrep,ntgrd)
cif.res.mat.new.width = matrix(NA,nrep,ntgrd)
cif.res.mat.new.band.width = matrix(NA,nrep,ntgrd2)
cif.res.mat.new.cvr = matrix(NA,nrep,ntgrd)
if (cmb.flg) {
  cif.res.mat.cmb.est = matrix(NA,nrep,ntgrd)
  cif.res.mat.cmb.sd1 = matrix(NA,nrep,ntgrd)
  cif.res.mat.cmb.sd2 = matrix(NA,nrep,ntgrd)
  cif.res.mat.cmb.width = matrix(NA,nrep,ntgrd)
  cif.res.mat.cmb.band.width = matrix(NA,nrep,ntgrd2)
  cif.res.mat.cmb.cvr = matrix(NA,nrep,ntgrd)
  cif.res.mat.cmb.wts = matrix(NA,nrep,ntgrd)
}
cif.res.aj.band.cvr.0 = rep(NA,ntgrd2)
cif.res.aj.band.cvr = rep(NA,ntgrd2)
cif.res.new.band.cvr = rep(NA,ntgrd2)
if (cmb.flg) {cif.res.cmb.band.cvr = rep(NA,ntgrd2)}
aj.band.cvr.0 = rep(NA,nrep)
aj.band.cvr = rep(NA,nrep)
new.band.cvr = rep(NA,nrep)
if (cmb.flg) {cmb.band.cvr = rep(NA,nrep)}

#SIMULATION LOOP
tstart = proc.time()
for (irep in 1:nrep) {
  
  print(irep)
  #if ((irep %% 10) == 0) print(irep)

  #GENERATE DATA
  datcur = samgen(n, pr.enter, T1_setting, T2_setting, R_setting,
    C_setting, L_R, U_R, rnddat) 
  age_recr = datcur$age_recr
  age_diag = datcur$age_diag
  age_death = datcur$age_death
  age_end_fu = datcur$age_end_fu
  status_end = datcur$status_end

  #COMPUTE ESTIMATORS
  cifest.rslt = cifcmp.full(age_recr, age_diag, age_death, age_end_fu, status_end,
    tgrd, tgrd1, tgrd2, covpr, nresam, adpwt, cmb.wgt, trans.opt, pt.ci)
  
  #STORE RESULTS FOR THIS REPLICATION 
  
  #FIRST BLOCK OF RESULTS
  cif.res.mat.aj.est[irep,] = cifest.rslt$aj.est  
  cif.res.mat.aj.sd1[irep,] = cifest.rslt$aj.sd1 
  cif.res.mat.aj.sd2[irep,] = cifest.rslt$aj.sd2
  cif.res.mat.aj.sd3[irep,] = cifest.rslt$aj.sd3
  cif.res.mat.aj.width[irep,] = cifest.rslt$aj.ptwise.ci.width.3
  cif.res.mat.aj.band.width[irep,] = cifest.rslt$aj.band.width
  cif.res.mat.new.est[irep,] = cifest.rslt$new.est  
  cif.res.mat.new.sd1[irep,] = cifest.rslt$new.sd1 
  cif.res.mat.new.sd2[irep,] = cifest.rslt$new.sd2
  cif.res.mat.new.width[irep,] = cifest.rslt$new.ptwise.ci.width
  cif.res.mat.new.band.width[irep,] = cifest.rslt$new.band.width
  if (cmb.flg) {
    cif.res.mat.cmb.est[irep,] = cifest.rslt$cmb.est  
    cif.res.mat.cmb.sd1[irep,] = cifest.rslt$cmb.sd1 
    cif.res.mat.cmb.sd2[irep,] = cifest.rslt$cmb.sd2
    cif.res.mat.cmb.width[irep,] = cifest.rslt$cmb.ptwise.ci.width 
    cif.res.mat.cmb.band.width[irep,] = cifest.rslt$cmb.band.width
    cif.res.mat.cmb.wts[irep,] = cifest.rslt$cmb.wts
  }  

  #POINTWISE CI COVERAGE
  for (ix.g in 1:ntgrd) {
    if (ipre) {
      cif.res.mat.aj.cvr.0[irep,ix.g] =
        (ciftru[tgrd[ix.g]] >= cifest.rslt$aj.ptwise.ci.lo.1[ix.g]) & 
        (ciftru[tgrd[ix.g]] <= cifest.rslt$aj.ptwise.ci.hi.1[ix.g])
    }
    cif.res.mat.aj.cvr.1[irep,ix.g] =
      (ciftru.aj[tgrd[ix.g]] >= cifest.rslt$aj.ptwise.ci.lo.1[ix.g]) & 
      (ciftru.aj[tgrd[ix.g]] <= cifest.rslt$aj.ptwise.ci.hi.1[ix.g])
    cif.res.mat.aj.cvr.2[irep,ix.g] =
      (ciftru.aj[tgrd[ix.g]] >= cifest.rslt$aj.ptwise.ci.lo.2[ix.g]) & 
      (ciftru.aj[tgrd[ix.g]] <= cifest.rslt$aj.ptwise.ci.hi.2[ix.g])
    cif.res.mat.aj.cvr.3[irep,ix.g] =
      (ciftru.aj[tgrd[ix.g]] >= cifest.rslt$aj.ptwise.ci.lo.3[ix.g]) & 
      (ciftru.aj[tgrd[ix.g]] <= cifest.rslt$aj.ptwise.ci.hi.3[ix.g])
    cif.res.mat.new.cvr[irep,ix.g] =
      (ciftru[tgrd[ix.g]] >= cifest.rslt$new.ptwise.ci.lo[ix.g]) & 
      (ciftru[tgrd[ix.g]] <= cifest.rslt$new.ptwise.ci.hi[ix.g])
    if (cmb.flg) {
      cif.res.mat.cmb.cvr[irep,ix.g] =
        (ciftru[tgrd[ix.g]] >= cifest.rslt$cmb.ptwise.ci.lo[ix.g]) & 
        (ciftru[tgrd[ix.g]] <= cifest.rslt$cmb.ptwise.ci.hi[ix.g])
    }
  }    
  
  #CONFIDENCE BAND COVERAGE
  for (ix.g in 1:ntgrd2) {
    if (ipre) {
      cif.res.aj.band.cvr.0[ix.g] =
        (ciftru[tgrd2[ix.g]] >= cifest.rslt$aj.band.lo[ix.g]) & 
        (ciftru[tgrd2[ix.g]] <= cifest.rslt$aj.band.hi[ix.g])
    } 
    cif.res.aj.band.cvr[ix.g] =
      (ciftru.aj[tgrd2[ix.g]] >= cifest.rslt$aj.band.lo[ix.g]) & 
      (ciftru.aj[tgrd2[ix.g]] <= cifest.rslt$aj.band.hi[ix.g])
    cif.res.new.band.cvr[ix.g] =
      (ciftru[tgrd2[ix.g]] >= cifest.rslt$new.band.lo[ix.g]) & 
      (ciftru[tgrd2[ix.g]] <= cifest.rslt$new.band.hi[ix.g])
    if (cmb.flg) {
      cif.res.cmb.band.cvr[ix.g] =
        (ciftru[tgrd2[ix.g]] >= cifest.rslt$cmb.band.lo[ix.g]) & 
        (ciftru[tgrd2[ix.g]] <= cifest.rslt$cmb.band.hi[ix.g])
    }
  }
  aj.band.cvr.0[irep] = prod(cif.res.aj.band.cvr.0)
  aj.band.cvr[irep] = prod(cif.res.aj.band.cvr)
  new.band.cvr[irep] = prod(cif.res.new.band.cvr)
  if (cmb.flg) {cmb.band.cvr[irep] = prod(cif.res.cmb.band.cvr)}
  #browser()
  
}
#END OF SIMULATION LOOP
tstop = proc.time()
print(tstop-tstart)

#MEANS OVER THE SIMULATIONS
cif.mean.aj.est =  colMeans(cif.res.mat.aj.est,na.rm=T)
cif.mean.aj.sd1 =  colMeans(cif.res.mat.aj.sd1,na.rm=T)
cif.mean.aj.sd2 =  colMeans(cif.res.mat.aj.sd2,na.rm=T)
cif.mean.aj.sd3 =  colMeans(cif.res.mat.aj.sd3,na.rm=T)
cif.mean.new.est = colMeans(cif.res.mat.new.est,na.rm=T)
cif.mean.new.sd1 = colMeans(cif.res.mat.new.sd1,na.rm=T)
cif.mean.new.sd2 = colMeans(cif.res.mat.new.sd2,na.rm=T)
if (cmb.flg) {
  cif.mean.cmb.est = colMeans(cif.res.mat.cmb.est,na.rm=T)
  cif.mean.cmb.sd1 = colMeans(cif.res.mat.cmb.sd1,na.rm=T)
  cif.mean.cmb.sd2 = colMeans(cif.res.mat.cmb.sd2,na.rm=T)
  cif.mean.cmb.wts = colMeans(cif.res.mat.cmb.wts,na.rm=T)
}  

#MEDIANS OVER THE SIMULATIONS
cif.med.aj.est =  apply(cif.res.mat.aj.est,2,median,na.rm=T)
cif.med.new.est = apply(cif.res.mat.new.est,2,median,na.rm=T)
if (cmb.flg) {cif.med.cmb.est = apply(cif.res.mat.cmb.est,2,median,na.rm=T)}

#EMPIRICAL SD'S OVER THE SIMULATIONS
cif.sd.aj =  apply(cif.res.mat.aj.est,2,sd,na.rm=T)
cif.sd.new = apply(cif.res.mat.new.est,2,sd,na.rm=T)
if (cmb.flg) {cif.sd.cmb = apply(cif.res.mat.cmb.est,2,sd,na.rm=T)}
cif.sd.aj.1 =  round(cif.sd.aj,digits=4)  
cif.sd.new.1 = round(cif.sd.new,digits=4)
if (cmb.flg) {cif.sd.cmb.1 = round(cif.sd.cmb,digits=4)}

#IQR'S OVER THE SIMULATIONS
cif.iqr.aj =  apply(cif.res.mat.aj.est,2,IQR,na.rm=T)
cif.iqr.new = apply(cif.res.mat.new.est,2,IQR,na.rm=T)
if (cmb.flg) {cif.iqr.cmb = apply(cif.res.mat.cmb.est,2,IQR,na.rm=T)}

#POINTWISE CONFIDENCE INTERVAL COVERAGE RATES
cvr.rte.ptwise.aj.0 =  colMeans(cif.res.mat.aj.cvr.0,na.rm=T)
cvr.rte.ptwise.aj.1 =  colMeans(cif.res.mat.aj.cvr.1,na.rm=T)
cvr.rte.ptwise.aj.2 =  colMeans(cif.res.mat.aj.cvr.2,na.rm=T)
cvr.rte.ptwise.aj.3 =  colMeans(cif.res.mat.aj.cvr.3,na.rm=T)
cvr.rte.ptwise.new = colMeans(cif.res.mat.new.cvr,na.rm=T)
if (cmb.flg) {cvr.rte.ptwise.cmb = colMeans(cif.res.mat.cmb.cvr,na.rm=T)}

#POINTWISE CONFIDENCE INTERVAL WIDTHS
cif.mean.width.aj =  colMeans(cif.res.mat.aj.width,na.rm=T)
cif.mean.width.new = colMeans(cif.res.mat.new.width,na.rm=T)
if (cmb.flg) {cif.mean.width.cmb = colMeans(cif.res.mat.cmb.width,na.rm=T)}

#CONFIDENCE BAND COVERAGE RATES
cvr.rte.band.aj.0 =  mean(aj.band.cvr.0,na.rm=T)
cvr.rte.band.aj =  mean(aj.band.cvr,na.rm=T)
cvr.rte.band.new = mean(new.band.cvr,na.rm=T)
if (cmb.flg) {cvr.rte.band.cmb = mean(cmb.band.cvr,na.rm=T)}

#CONFIDENCE BAND WIDTHS
cif.mean.band.width.aj =  colMeans(cif.res.mat.aj.band.width,na.rm=T)
cif.mean.band.width.new = colMeans(cif.res.mat.new.band.width,na.rm=T)
if (cmb.flg) {cif.mean.band.width.cmb = colMeans(cif.res.mat.cmb.band.width,na.rm=T)}

### PLOTS ###

#LABELS FOR PLOTS
tru.lab = rep("True",ntgrd)
tru0.lab = rep("True0",ntgrd)
tru.aj.lab = rep("AJ-Target",ntgrd)
aj.lab = rep("AJ",ntgrd)
new.lab = rep("New",ntgrd)
comb.lab = rep("Comb",ntgrd)
ajmod.lab = rep("AJ-Modified",ntgrd)
emp.lab = rep("Empirical",ntgrd)
est.lab = rep("Estimated",ntgrd)
est1.lab = rep("Estimated-1",ntgrd)
est2.lab = rep("Estimated-2",ntgrd)
est3.lab = rep("Estimated-etm",ntgrd)
cif.lab = rep("CIF",ntgrd)

if (plot.to.file) {pdf(plot.pdf.name)}

#LINE TYPES AND COLORS FOR PLOTS
#  “solid”, “dashed”, “dotted”, “dotdash”, “longdash”, “twodash”
#   black    blue      red       green       purple

#LINE WIDTHS
lwdot = 1.5
lwrest = 0.8

#COVER PAGE
trns.lbl = c('none', '-log(1-u)', '0.5*pi - asin(sqrt(1-u))' )
pt.ci.lbl = c('normal-theory', 'bootstrap')
bwt.lbl = c('normal','exponential','Poisson')
txt0 = noquote(paste0('Scenario: ', scenario.pr))
txt1 = noquote(paste0('sample size = ', n))
txt2 = noquote(paste0('number of simulation replications = ', nrep))
txt3 = noquote(paste0('number of bootstrap replications = ', nresam))
txt4 = noquote(paste0('transformation: ', trns.lbl[trans.opt]))
txt5 = noquote(paste0("pointwise CI's done by: ", pt.ci.lbl[pt.ci]))
txt6 = noquote(paste0('auxflg = ', auxflg))
txt7 = noquote(paste0('bootstrap weights: ', bwt.lbl[bwts]))
txt8 = noquote(paste0('Date/Time: ', Sys.time()))
plot.new()
text(x=.1, y=1, 'SETTINGS', adj=c(0,1))
text(x=.1, y=.9, txt0, adj=c(0,1))
text(x=.1, y=.8, txt1, adj=c(0,1))
text(x=.1, y=.7, txt2, adj=c(0,1))
text(x=.1, y=.6, txt3, adj=c(0,1))
text(x=.1, y=.5, txt4, adj=c(0,1))
text(x=.1, y=.4, txt5, adj=c(0,1))
text(x=.1, y=.3, txt6, adj=c(0,1))
text(x=.1, y=.2, txt7, adj=c(0,1))
text(x=.1, y=.1, txt8, adj=c(0,1))

#PLOT MEANS
if (!ipre) {
  temp = cbind(c(tru.lab,aj.lab,new.lab,comb.lab), c(tgrd,tgrd,tgrd,tgrd), 
    c(ciftru[tgrd],cif.mean.aj.est,cif.mean.new.est,cif.mean.cmb.est))
  lev.vec = c('True', 'AJ', 'New', 'Comb')
  col.vec = c('black', 'blue', 'red', 'purple')
  typ.vec = c('solid', 'dashed', 'dotted', 'longdash')
  lw.vec = c(lwrest, lwrest, lwdot, lwrest)
}
if (ipre) {
  temp = cbind(c(tru.lab,tru.aj.lab,aj.lab,new.lab), c(tgrd,tgrd,tgrd,tgrd), 
    c(ciftru[tgrd],ciftru.aj[tgrd],cif.mean.aj.est,cif.mean.new.est))
  lev.vec = c('True', 'AJ-Target', 'AJ', 'New')
  col.vec = c('black', 'purple', 'blue', 'red')
  typ.vec = c('solid', 'longdash', 'dashed', 'dotted')
  lw.vec = c(lwrest, lwrest, lwrest, lwdot)
}
df1 = as.data.frame(temp)
colnames(df1) = c("Method","Age","CIF")
df1$Age = as.numeric(df1$Age)
df1$CIF = as.numeric(df1$CIF)
df1$Method = factor(df1$Method, levels=lev.vec)
plt.means = ggplot(df1, aes(x = Age, y = CIF)) + 
  theme_bw() +
  geom_line(aes(color = Method, linetype = Method, linewidth = Method)) + 
  scale_color_manual(values = col.vec) +
  scale_linetype_manual(values = typ.vec) +
  scale_linewidth_manual(values = lw.vec) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  guides(linewidth='none') +
  ggtitle(paste0(label, ", Means"))
plot(plt.means)

#PLOT MEDIANS
if (!ipre) {
  temp = cbind(c(tru.lab,aj.lab,new.lab,comb.lab), c(tgrd,tgrd,tgrd,tgrd), 
    c(ciftru[tgrd],cif.med.aj.est,cif.med.new.est,cif.med.cmb.est))
  lev.vec = c('True', 'AJ', 'New', 'Comb')
  col.vec = c('black', 'blue', 'red', 'purple')
  typ.vec = c('solid', 'dashed', 'dotted', 'longdash')
  lw.vec = c(lwrest, lwrest, lwdot, lwrest)
}
if (ipre) {
  temp = cbind(c(tru.lab,tru.aj.lab,aj.lab,new.lab), c(tgrd,tgrd,tgrd,tgrd), 
    c(ciftru[tgrd],ciftru.aj[tgrd],cif.med.aj.est,cif.med.new.est))
  lev.vec = c('True', 'AJ-Target', 'AJ', 'New')
  col.vec = c('black', 'purple', 'blue', 'red')
  typ.vec = c('solid', 'longdash', 'dashed', 'dotted')
  lw.vec = c(lwrest, lwrest, lwrest, lwdot)
  
}
df1 = as.data.frame(temp)
colnames(df1) = c("Method","Age","CIF")
df1$Age = as.numeric(df1$Age)
df1$CIF = as.numeric(df1$CIF)
df1$Method = factor(df1$Method, levels=lev.vec)
plt.medians = ggplot(df1, aes(x = Age, y = CIF)) + 
  theme_bw() +
  geom_line(aes(color = Method, linetype = Method, linewidth = Method)) + 
  scale_color_manual(values = col.vec) +
  scale_linetype_manual(values = typ.vec) +
  scale_linewidth_manual(values = lw.vec) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  guides(linewidth='none') +
  ggtitle(paste0(label, ", Medians"))
plot(plt.medians)

#PLOT SD'S
if (!ipre) {
  temp = cbind(c(aj.lab,new.lab,comb.lab), c(tgrd,tgrd,tgrd), 
    c(cif.sd.aj,cif.sd.new,cif.sd.cmb))
  lev.vec = c('AJ', 'New', 'Comb')
  col.vec = c('blue', 'red', 'purple')
  typ.vec = c('dashed', 'dotted', 'longdash')
  lw.vec = c(lwrest, lwdot, lwrest)
}
if (ipre) {
  temp = cbind(c(aj.lab,new.lab), c(tgrd,tgrd), 
    c(cif.sd.aj,cif.sd.new))
  lev.vec = c('AJ', 'New')
  col.vec = c('blue', 'red')
  typ.vec = c('dashed', 'dotted')
  lw.vec = c(lwrest, lwdot)
}  
df1 = as.data.frame(temp)
colnames(df1) = c("Method","Age","SD")
df1$Age = as.numeric(df1$Age)
df1$SD = as.numeric(df1$SD)
df1$Method = factor(df1$Method, levels=lev.vec)
plt.sds = ggplot(df1, aes(x = Age, y = SD)) + 
  theme_bw() +
  geom_line(aes(color = Method, linetype = Method, linewidth = Method)) + 
  scale_color_manual(values = col.vec) +
  scale_linetype_manual(values = typ.vec) +
  scale_linewidth_manual(values = lw.vec) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  guides(linewidth='none') +
  ggtitle(paste0(label, ", SD'S"))
plot(plt.sds)

#PLOT IQR'S
if (!ipre) {
  temp = cbind(c(aj.lab,new.lab,comb.lab), c(tgrd,tgrd,tgrd), 
    c(cif.iqr.aj,cif.iqr.new,cif.iqr.cmb))
  lev.vec = c('AJ', 'New', 'Comb')
  col.vec = c('blue', 'red', 'purple')
  typ.vec = c('dashed', 'dotted', 'longdash')
  lw.vec = c(lwrest, lwdot, lwrest)
}
if (ipre) {
  temp = cbind(c(aj.lab,new.lab), c(tgrd,tgrd), 
    c(cif.iqr.aj,cif.iqr.new))
  lev.vec = c('AJ', 'New')
  col.vec = c('blue', 'red')
  typ.vec = c('dashed', 'dotted')
  lw.vec = c(lwrest, lwdot)
}  
df1 = as.data.frame(temp)
colnames(df1) = c("Method","Age","IQR")
df1$Age = as.numeric(df1$Age)
df1$IQR = as.numeric(df1$IQR)
df1$Method = factor(df1$Method, levels=lev.vec)
plt.iqr = ggplot(df1, aes(x = Age, y = IQR)) + 
  theme_bw() +
  geom_line(aes(color = Method, linetype = Method, linewidth = Method)) + 
  scale_color_manual(values = col.vec) +
  scale_linetype_manual(values = typ.vec) +
  scale_linewidth_manual(values = lw.vec) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  guides(linewidth='none') +
  ggtitle(paste0(label, ", IQR'S"))
plot(plt.iqr)

#PLOT MEAN ESTIMATED SD AND EMPIRICAL SD FOR AJ ESTIMATOR
if (!auxflg) {
  temp = cbind(c(emp.lab,est.lab,est3.lab), c(tgrd,tgrd,tgrd),
    c(cif.sd.aj,cif.mean.aj.sd1,cif.mean.aj.sd3))
  lev.vec = c('Empirical', 'Estimated', 'Estimated-etm')
  col.vec = c('black', 'blue', 'red')
  typ.vec = c('solid', 'dashed', 'dotted')
  lw.vec = c(lwrest, lwrest, lwdot)
}
if (auxflg) {
  temp = cbind(c(emp.lab,est1.lab,est2.lab,est3.lab), c(tgrd,tgrd,tgrd,tgrd),
    c(cif.sd.aj,cif.mean.aj.sd1,cif.mean.aj.sd2,cif.mean.aj.sd3))
  lev.vec = c('Empirical', 'Estimated-1', 'Estimated-2', 'Estimated-etm')
  col.vec = c('black', 'blue', 'green', 'red')
  typ.vec = c('solid', 'dashed', 'dotdash', 'dotted')
  lw.vec = c(lwrest, lwrest, lwrest, lwdot)
}
df1 = as.data.frame(temp)
colnames(df1) = c("Method","Age","SD")
df1$Age = as.numeric(df1$Age)
df1$SD = as.numeric(df1$SD)
df1$Method = factor(df1$Method, levels=lev.vec)
plt.sdsd = ggplot(df1, aes(x = Age, y = SD)) + 
  theme_bw() +
  geom_line(aes(color = Method, linetype = Method, linewidth = Method)) + 
  scale_color_manual(values = col.vec) +
  scale_linetype_manual(values = typ.vec) +
  scale_linewidth_manual(values = lw.vec) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  guides(linewidth='none') +
  ggtitle(paste0(label, ", AJ Estimator, Empirical vs. Estimated SD's"))
plot(plt.sdsd)  

#PLOT MEAN ESTIMATED SD AND EMPIRICAL SD FOR NEW ESTIMATOR
if (!auxflg) {
  temp = cbind(c(emp.lab,est.lab), c(tgrd,tgrd),
    c(cif.sd.new,cif.mean.new.sd1))
  lev.vec = c('Empirical', 'Estimated')
  col.vec = c('black', 'blue')
  typ.vec = c('solid', 'dashed')
  lw.vec = c(lwrest, lwrest)
}
if (auxflg) {
  temp = cbind(c(emp.lab,est1.lab,est2.lab), c(tgrd,tgrd,tgrd),
    c(cif.sd.new,cif.mean.new.sd1,cif.mean.new.sd2))
  lev.vec = c('Empirical', 'Estimated-1', 'Estimated-2')
  col.vec = c('black', 'blue', 'green')
  typ.vec = c('solid', 'dashed', 'dotdash')
  lw.vec = w.vec = c(lwrest, lwrest, lwrest)
}
df1 = as.data.frame(temp)
colnames(df1) = c("Method","Age","SD")
df1$Age = as.numeric(df1$Age)
df1$SD = as.numeric(df1$SD)
df1$Method = factor(df1$Method, levels=lev.vec)
plt.sdsd = ggplot(df1, aes(x = Age, y = SD)) + 
  theme_bw() +
  geom_line(aes(color = Method, linetype = Method, linewidth = Method)) + 
  scale_color_manual(values = col.vec) +
  scale_linetype_manual(values = typ.vec) +
  scale_linewidth_manual(values = lw.vec) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  guides(linewidth='none') +
  ggtitle(paste0(label, ", New Estimator, Empirical vs. Estimated SD's"))
plot(plt.sdsd)  

#PLOT MEAN ESTIMATED SD AND EMPIRICAL SD FOR COMBINED ESTIMATOR
if (cmb.flg) {
if (!auxflg) {
  temp = cbind(c(emp.lab,est.lab), c(tgrd,tgrd),
    c(cif.sd.cmb,cif.mean.cmb.sd1))
  lev.vec = c('Empirical', 'Estimated')
  col.vec = c('black', 'blue')
  typ.vec = c('solid', 'dashed')
  lw.vec = c(lwrest, lwrest)
}
if (auxflg) {
  temp = cbind(c(emp.lab,est1.lab,est2.lab), c(tgrd,tgrd,tgrd),
    c(cif.sd.cmb,cif.mean.cmb.sd1,cif.mean.cmb.sd2))
  lev.vec = c('Empirical', 'Estimated-1', 'Estimated-2')
  col.vec = c('black', 'blue', 'red')
  typ.vec = c('solid', 'dashed', 'dotted')
  lw.vec = c(lwrest, lwrest, lwdot)
}
df1 = as.data.frame(temp)
colnames(df1) = c("Method","Age","SD")
df1$Age = as.numeric(df1$Age)
df1$SD = as.numeric(df1$SD)
df1$Method = factor(df1$Method, levels=lev.vec)
plt.sdsd = ggplot(df1, aes(x = Age, y = SD)) + 
  theme_bw() +
  geom_line(aes(color = Method, linetype = Method, linewidth = Method)) + 
  scale_color_manual(values = col.vec) +
  scale_linetype_manual(values = typ.vec) +
  scale_linewidth_manual(values = lw.vec) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  guides(linewidth='none') +
  ggtitle(paste0(label, ", Combined Estimator, Empirical vs. Estimated SD's"))
plot(plt.sdsd)  
}

#PLOT POINTWISE CI COVERAGE RATES
if (!ipre) {
  temp = cbind(c(aj.lab[ixcb1],new.lab[ixcb1],comb.lab[ixcb1]), c(tgrd1,tgrd1,tgrd1), 
    c(cvr.rte.ptwise.aj.3[ixcb1],cvr.rte.ptwise.new[ixcb1],cvr.rte.ptwise.cmb[ixcb1]))
  lev.vec = c('AJ', 'New', 'Comb')
  col.vec = c('blue', 'red', 'purple')
  typ.vec = c('dashed', 'dotted', 'longdash')
  lw.vec = c(lwrest, lwdot, lwrest)
  df1 = as.data.frame(temp)
  colnames(df1) = c("Method","Age","CI_Coverage_Rate")
  df1$Age = as.numeric(df1$Age)
  df1$CI_Coverage_Rate = as.numeric(df1$CI_Coverage_Rate)
  df1$Method = factor(df1$Method, levels=lev.vec)
  plt.cvr = ggplot(df1, aes(x = Age, y = CI_Coverage_Rate)) + 
    ylim(0.70, 1.00) +
    theme_bw() +
    geom_line(aes(color = Method, linetype = Method, linewidth = Method)) + 
    scale_color_manual(values = col.vec) +
    scale_linetype_manual(values = typ.vec) +
    scale_linewidth_manual(values = lw.vec) +
    guides(color = guide_legend(override.aes = list(size=3))) +
    guides(linewidth='none') +
    ggtitle(paste0(label, ", CICR'S")) + 
    geom_hline(yintercept=0.95, color='green')
  plot(plt.cvr)
}
if (ipre) {
  temp = cbind(tgrd1,cvr.rte.ptwise.new[ixcb1])
  df1 = as.data.frame(temp)
  colnames(df1) = c("Age","CI_Coverage_Rate")
  df1$Age = as.numeric(df1$Age)
  df1$CI_Coverage_Rate = as.numeric(df1$CI_Coverage_Rate)
  plt.cvr = ggplot(df1, aes(x = Age, y = CI_Coverage_Rate)) + 
    geom_line() + 
    ylim(0.70, 1.00) +
    theme_bw() +
    ggtitle(paste0(label, ", CI Coverage Rate for New Method")) + 
    geom_hline(yintercept=0.95, color='green')
  plot(plt.cvr)
}

#PLOT POINTWISE CI WIDTHS
if (!ipre) {
  temp = cbind(c(aj.lab[ixcb1],new.lab[ixcb1],comb.lab[ixcb1]), c(tgrd1,tgrd1,tgrd1), 
    c(cif.mean.width.aj[ixcb1],cif.mean.width.new[ixcb1],cif.mean.width.cmb[ixcb1]))
  lev.vec = c('AJ', 'New', 'Comb')
  col.vec = c('blue', 'red', 'purple')
  typ.vec = c('dashed', 'dotted', 'longdash')
  lw.vec = c(lwrest, lwdot, lwrest)
df1 = as.data.frame(temp)
colnames(df1) = c("Method","Age","CI_Width")
df1$Age = as.numeric(df1$Age)
df1$CI_Width = as.numeric(df1$CI_Width)
df1$Method = factor(df1$Method, levels=lev.vec)
plt.width = ggplot(df1, aes(x = Age, y = CI_Width)) + 
  theme_bw() +
  geom_line(aes(color = Method, linetype = Method, linewidth = Method)) + 
  scale_color_manual(values = col.vec) +
  scale_linetype_manual(values = typ.vec) +
  scale_linewidth_manual(values = lw.vec) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  guides(linewidth='none') +
  ggtitle(paste0(label, ", CI Width"))
plot(plt.width)
}
if (ipre) {
  temp = cbind(tgrd1,cif.mean.width.new[ixcb1])
  df1 = as.data.frame(temp)
  colnames(df1) = c("Age","CI_Width")
  df1$Age = as.numeric(df1$Age)
  df1$CI_Width = as.numeric(df1$CI_Width)
  plt.width = ggplot(df1, aes(x = Age, y = CI_Width)) + 
    theme_bw() +
    geom_line() + 
    ggtitle(paste0(label, ", CI Width for New Estimator"))
  plot(plt.width)
}

#CONFIDENCE BAND COVERAGE RATES
txt0 = noquote(paste0('Scenario: ', scenario.pr))
if (T1_setting != 3) {
  txt1 = noquote(paste0('AJ: ', round(cvr.rte.band.aj, digits=3)))
  txt2 = noquote(paste0('new: ', round(cvr.rte.band.new, digits=3)))
  txt3 = noquote(paste0('Combo: ', round(cvr.rte.band.cmb, digits=3)))
}
if (T1_setting == 3) {
  txt1 = noquote(paste0('AJ0: ', round(cvr.rte.band.aj.0, digits=3)))
  txt2 = noquote(paste0('AJ: ', round(cvr.rte.band.aj, digits=3)))
  txt3 = noquote(paste0('New: ', round(cvr.rte.band.new, digits=3)))
}
plot.new()
text(x=.1, y=1, 'CONFIDENCE BAND COVERAGE RATES', adj=c(0,1))
text(x=.1, y=.9, txt0, adj=c(0,1))
text(x=.1, y=.8, txt1, adj=c(0,1))
text(x=.1, y=.7, txt2, adj=c(0,1))
text(x=.1, y=.6, txt3, adj=c(0,1))

#PLOT CONFIDENCE BAND WIDTHS
if (!ipre) {
  temp = cbind(c(aj.lab[ixcb2],new.lab[ixcb2],comb.lab[ixcb2]), c(tgrd2,tgrd2,tgrd2), 
    c(cif.mean.band.width.aj,cif.mean.band.width.new,cif.mean.band.width.cmb))
  lev.vec = c('AJ', 'New', 'Comb')
  col.vec = c('blue', 'red', 'purple')
  typ.vec = c('dashed', 'dotted', 'longdash')
  lw.vec = c(lwrest, lwdot, lwrest)
  df1 = as.data.frame(temp)
  colnames(df1) = c("Method","Age","Conf_Band_Width")
  df1$Age = as.numeric(df1$Age)
  df1$Conf_Band_Width = as.numeric(df1$Conf_Band_Width)
  df1$Method = factor(df1$Method, levels=lev.vec)
  plt.band.width = ggplot(df1, aes(x = Age, y = Conf_Band_Width)) + 
    theme_bw() +
    geom_line(aes(color = Method, linetype = Method, linewidth = Method)) +
    scale_color_manual(values = col.vec) +
    scale_linetype_manual(values = typ.vec) +
    scale_linewidth_manual(values = lw.vec) +
    guides(color = guide_legend(override.aes = list(size=3))) +
    guides(linewidth='none') +
    ggtitle(paste0(label, ", Confidence Band Width"))
  plot(plt.band.width)
}
if (ipre) {
  temp = cbind(tgrd2,cif.mean.band.width.new)
  df1 = as.data.frame(temp)
  colnames(df1) = c("Age","Conf_Band_Width")
  df1$Age = as.numeric(df1$Age)
  df1$Conf_Band_Width = as.numeric(df1$Conf_Band_Width)
  plt.band.width = ggplot(df1, aes(x = Age, y = Conf_Band_Width)) + 
    theme_bw() +
    geom_line() +
  ggtitle(paste0(label, ", Confidence Band Width for New Method"))
  plot(plt.band.width)
}

if (plot.to.file) {dev.off()}

cat('\n')
print(table(aj.band.cvr,useNA='ifany'))
cat('\n')
print(table(new.band.cvr,useNA='ifany'))

cat('\n')
if (cmb.flg) {
  cvr.rte.band.all = cbind(cvr.rte.band.aj, cvr.rte.band.new, cvr.rte.band.cmb)
}
if (!cmb.flg) { 
  cvr.rte.band.all = cbind(cvr.rte.band.aj.0, cvr.rte.band.aj, cvr.rte.band.new)
}
cvr.rte.band.all = round(cvr.rte.band.all,digits=3)
cvr.rte.band.all = cbind(scenario.pr,cvr.rte.band.all)
print(cvr.rte.band.all)

#CREATE PLOT OBJECT
lgnd = TRUE
if (lgnd) {
pprplt = ggarrange(plt.means, plt.sds, plt.cvr, nrow=1, ncol=3,
  common.legend=TRUE, legend='bottom')
}
if (!lgnd) {
  pprplt = ggarrange(plt.means, plt.sds, plt.cvr, nrow=1, ncol=3, legend='none')
}

#GIVE NAME TO PLOT OBJECT
pname0 = paste0('pprplt_new_',label)
assign(pname0,pprplt)

#SAVE PLOT OBJECT
pl.name = paste0(pname0,'.RData')
pname = paste0(sim.path,pl.name)
save(list=pname0, file=pname)

#SAVE WORKSPACE
wsn = paste0('/',label,'.RData')
imnam = paste0(sim.path,wsn)
save.image(imnam)

beep()
