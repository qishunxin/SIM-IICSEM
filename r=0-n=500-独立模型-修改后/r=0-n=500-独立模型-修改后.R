####IC-r = 0
rm(list=ls());
#install.packages("MCMCpack")
library(MASS)
library(HI)
library(coda)
library(MCMCpack)
set.seed(11)

   t0 <- Sys.time()
   loop  = 100 
 
   p=6  
   q=2   
   m=4  
   kh=10
   klam=10
   n = 500;

sig.mh = 3

### for the measurment equation

   Id.Delta <- matrix(rep(1, q*(1 + p)), nrow  = q)
   Id.Delta <- (Id.Delta > 0)

##  for the structual equation
   Id.se <- cbind(Id.Delta)
   n.se.row <- rowSums(Id.se)  
   n.se <- sum(Id.se)

##########prior setting

alpha.psd <- 9      # 
beta.psd <- 4       # 
Delta0 <- matrix(rep(0, 1 + p), nrow = q, ncol = (1 + p), byrow = TRUE)
sig.delta <- rep(0.001,n.se)         
L.se0 <- cbind(Delta0)


   loop_L.se    = array(rep(0,loop*n.se),dim=c(loop,n.se))
   loop_psi    = array(rep(0,loop*m),dim=c(loop,m))
   loop_At    = array(rep(0,loop*p),dim=c(loop,p))
#   loop_Ag    = array(rep(0,loop*p),dim=c(loop,p))
   loop_betat = array(rep(0,loop*q),dim=c(loop,q))
#   loop_betag = array(rep(0,loop*q),dim=c(loop,q))
   loop_gam   = array(rep(0,loop*kh),dim=c(loop,kh))
   loop_alpha = array(rep(0,loop*klam),dim=c(loop,klam))
   loop_eta   = array(rep(0,loop),dim=c(loop,1))
   loop_pho   = array(rep(0,loop),dim=c(loop,1))
   loop_Lamda = array(rep(0,loop*(m-2)),dim=c(loop,(m-2)))
   loop_phy   = array(rep(0,loop*q),dim=c(loop,q))
   loop_omega = array(0,dim=c(n,q,loop))


   parsee_L.se    = array(rep(0,loop*n.se),dim=c(loop,n.se))
   parsee_psi    = array(rep(0,loop*m),dim=c(loop,m))
   parsee_At    = array(rep(0,loop*p),dim=c(loop,p))
#   parsee_Ag    = array(rep(0,loop*p),dim=c(loop,p))
   parsee_betat = array(rep(0,loop*q),dim=c(loop,q))
#   parsee_betag = array(rep(0,loop*q),dim=c(loop,q))
   parsee_Lamda = array(rep(0,loop*(m-2)),dim=c(loop,(m-2)))
   parsee_phy   = array(rep(0,loop*q),dim=c(loop,q))
   accept.rate  = array(rep(0,loop*n),dim=c(loop,n))
   pd=NULL

   loop_surv.low_1 = array(rep(0,loop*100),dim=c(loop,100))
   loop_surv.upp_1 = array(rep(0,loop*100),dim=c(loop,100))
   loop_surv.med_1 = array(rep(0,loop*100),dim=c(loop,100))

   loop_surv.low_2 = array(rep(0,loop*100),dim=c(loop,100))
   loop_surv.upp_2 = array(rep(0,loop*100),dim=c(loop,100))
   loop_surv.med_2 = array(rep(0,loop*100),dim=c(loop,100))

   effect_names <- c("Total_Effect", "Direct_Effect", "IE_via_M1", "IE_via_M2",
                     "Check_TE_minus_sum")
   contrast_names <- c("Z1_vs_ref", "Z2_vs_ref")
   loop_effect_t_hat  <- array(NA_real_, dim = c(loop, 2, length(effect_names)),
                               dimnames = list(NULL, contrast_names, effect_names))
   loop_effect_t_true <- array(NA_real_, dim = c(loop, 2, length(effect_names)),
                               dimnames = list(NULL, contrast_names, effect_names))

make_H0_y_interval <- function(HL, HR, status, eps = 1e-8) {
   status <- as.numeric(status)
   H0 <- ifelse(status == 0, HR,
                ifelse(status == 1, HR - HL, HL))
   pmax(as.vector(H0), eps)
}

calc_cf_mean_q2 <- function(z_out, z_med, Zbase, H0, gamma, Delta_vec,
                            betaM, phy_vec, std_norm, eps = 1e-8) {
   Delta_mat <- matrix(Delta_vec, nrow = q, ncol = 1 + p)
   mc <- dim(std_norm)[3]
   out <- numeric(mc)

   for (bb in 1:mc) {
      Mu <- matrix(0, nrow = nrow(Zbase), ncol = q)
      for (jj in 1:q) {
         Mu[, jj] <- Delta_mat[jj, 1] +
            Zbase %*% Delta_mat[jj, 2:5] +
            z_med[[jj]][1] * Delta_mat[jj, 6] +
            z_med[[jj]][2] * Delta_mat[jj, 7]
      }

      Omega_cf <- Mu + sweep(std_norm[, , bb], 2, sqrt(phy_vec), "*")
      eta <- as.vector(Zbase %*% gamma[1:4]) +
         z_out[1] * gamma[5] + z_out[2] * gamma[6] +
         as.vector(Omega_cf %*% betaM)

      if (abs(r) < eps) {
         out[bb] <- mean(exp(-H0 * exp(eta)))
      } else {
         out[bb] <- mean((1 + r * H0 * exp(eta))^(-1 / r))
      }
   }

   mean(out)
}

calc_effect_q2 <- function(z_trt, Zbase, H0, gamma, Delta_vec, betaM,
                           phy_vec, std_norm) {
   z_ref <- c(0, 0)

   e00 <- calc_cf_mean_q2(z_ref, list(z_ref, z_ref), Zbase, H0, gamma,
                          Delta_vec, betaM, phy_vec, std_norm)
   e10 <- calc_cf_mean_q2(z_trt, list(z_ref, z_ref), Zbase, H0, gamma,
                          Delta_vec, betaM, phy_vec, std_norm)
   e11 <- calc_cf_mean_q2(z_trt, list(z_trt, z_ref), Zbase, H0, gamma,
                          Delta_vec, betaM, phy_vec, std_norm)
   e12 <- calc_cf_mean_q2(z_trt, list(z_trt, z_trt), Zbase, H0, gamma,
                          Delta_vec, betaM, phy_vec, std_norm)

   c(Total_Effect = e12 - e00,
     Direct_Effect = e10 - e00,
     IE_via_M1 = e11 - e10,
     IE_via_M2 = e12 - e11,
     Check_TE_minus_sum = (e12 - e00) - ((e10 - e00) + (e11 - e10) + (e12 - e11)))
}

summarize_effects <- function(est_array, true_array, family_name) {
   rows <- list()
   rr <- 1
   for (cc in seq_along(contrast_names)) {
      for (ee in seq_along(effect_names)) {
         est <- est_array[, cc, ee]
         tru <- true_array[, cc, ee]
         rows[[rr]] <- data.frame(
            Family = family_name,
            Contrast = contrast_names[cc],
            Effect = effect_names[ee],
            True_Value = mean(tru, na.rm = TRUE),
            Estimate = mean(est, na.rm = TRUE),
            Bias = mean(est - tru, na.rm = TRUE),
            RMSE = sqrt(mean((est - tru)^2, na.rm = TRUE)),
            SD = sd(est, na.rm = TRUE)
         )
         rr <- rr + 1
      }
   }
   do.call(rbind, rows)
}

for(loopp in 1:loop){

## True coeffs
gam.z.true.x <- c(0.1, 0.1, 0.5, 0.5)
gam.s.true.x <- c(0.5, 0.1)
#gam.z.true.g <- c(0.5, 0.5, 0.1, 0.1)
#gam.s.true.g <- c(0.1, 0.5)
AT <- c(gam.z.true.x, gam.s.true.x)
#AG <- c(gam.z.true.g, gam.s.true.g)
 
  betaT = c(0.5,0.1);       #gam.m.true.t
#  betaG = c(1,1);       #gam.m.true.g
  r  = 0

Delta.true <- matrix(c(0, 0.2, 0, 0.2, 0, 0, 0.2, 
                     0.2, 0, 0.2, 0,0,0.2,0.2), 
              nrow = q, ncol = (1 + p), byrow = TRUE)

L.se.true <- as.vector(Delta.true)
psi.true <- rep(0.2, m)  #GAM = psi.true
 
## generate x

  z <- 4  
  S <- t(rmultinom(n, 1, prob = c(0.5, 0.2, 0.3)))
  S <- S[, 2:3]
  Z <- mvrnorm(n, rep(0,z), diag(1,z))
  X <- cbind(Z, S);

## generate omega

 Dat <- cbind(rep(1,n), Z, S)
 Coe <- cbind(Delta.true)

phy.true    =  rep(1, q)        ###閻╃缍嬫禍宸攕dmatrix(c(1,0,0,0,1,0,0,0,1),q,q)
OmegaT = Dat %*% t(Coe) + mvrnorm(n,rep(0,q),diag(phy.true))

#########################################Generate B: model 3
#鐎规矮绠焞amda閻晠妯€閿涙碍缍旈崣姗€鍣簑閻ㄥ嫮閮撮 ?
   lamda1 = c(1,0.9,0,0)     #lamda閻╃缍嬫禍 ?
   lamda2 = c(0,0,1,0.8)
   LamdaT = matrix(c(lamda1,lamda2),ncol=m,byrow=TRUE)  #L.me.true <- LamdaT
   B <- OmegaT %*% LamdaT + mvrnorm(n, rep(0, m), diag(psi.true))  

########################################################
                             #model two閻╅晲绶?
########################################################
## Baseline Survival
#S0ofg = function(t) 2*((1+t)^2-1);

## The Survival function:
#Siofg = function(t,x,omegaT,v=0) exp((-S0ofg(t))*exp(sum(x*AG)+sum(omegaT*betaG)+v)) ;
#Fiofg = function(t,x,omegaT,v=0) 1-Siofg(t,x,omegaT,v);

## The inverse for Fioft
#Finvg = function(u, x, omegaT, v=0) uniroot(function (t) Fiofg(t,x,omegaT,v)-u,lower=1e-100, upper=1e100,extendInt ="yes")$root

## generate survival times
#u  = runif(n);
#tG = rep(0, n);
#for (i in 1:n){
#  tG[i] = Finvg(u[i], X[i,], OmegaT[i,]);
#}

########################################################
                             #model one
########################################################
## Baseline Survival
S0oft = function(t) log(1+t)+exp(t)-1;

## The Survival function:
Sioft = function(t,x,omegaT,v=0) exp( (-S0oft(t))*exp(sum(x*AT)+sum(omegaT*betaT)+v) );
Fioft = function(t,x,omegaT,v=0) 1-Sioft(t,x,omegaT,v);

## The inverse for Fioft
Finv  = function(u, x, omegaT, v=0) uniroot(function (t) Fioft(t,x,omegaT,v)-u,lower=1e-10, upper=1e100,extendInt ="yes")$root

##-------------Generate data-------------------##
## generate survival times  model one
u  = runif(n);
tT = rep(0, n);
for (i in 1:n){
  tT[i] = Finv(u[i], X[i,], OmegaT[i,]);
}

### ----------- interval-censored -------------###
# interval-censored part 
   t1_int=rep(NA, n);t2_int=rep(NA, n); delta_int=rep(NA, n);
   npois = rpois(n, 2)+1;
   for(i in 1:n){
   #    tGi    = as.matrix(seq(0.5,npois[i],1))     
   #   gaptime = (tGi%*%(tG))[,i];
      gaptime = cumsum(rexp(npois[i], 1));
      pp  = Fioft(gaptime, X[i,], OmegaT[i,]);
      ind = sum(u[i]>pp);
      if (ind==0){
         delta_int[i] = 0;
         t2_int[i] = gaptime[1];
      }else if (ind==npois[i]){
            delta_int[i] = 2;
            t1_int[i] = gaptime[ind];
         }else{
             delta_int[i] = 1;
             t1_int[i] = gaptime[ind];
             t2_int[i] = gaptime[ind+1];
          }
    }
    t1=t1_int; t2=t2_int; delta = delta_int;

## make a data frame
   #d = data.frame(t1=t1, t2=t2, Z=Z, S=S, delta=delta, tT=tT, tG=tG);
   d = data.frame(t1=t1, t2=t2, Z=Z, S=S, delta=delta, tT=tT);
   table(d$delta)/n;
   ppd = d
   pd  = rbind(pd,ppd)

#################################################################

Ispline<-function(x,order,knots){

k=order+1
m=length(knots)
n=m-2+k # number of parameters
t=c(rep(1,k)*knots[1], knots[2:(m-1)], rep(1,k)*knots[m]) # newknots

yy1=array(rep(0,(n+k-1)*length(x)),dim=c(n+k-1, length(x)))
for (l in k:n){
    yy1[l,]=(x>=t[l] & x<t[l+1])/(t[l+1]-t[l])
}

yytem1=yy1
for (ii in 1:order){
   yytem2=array(rep(0,(n+k-1-ii)*length(x)),dim=c(n+k-1-ii, length(x)))
   for (i in (k-ii):n){
      yytem2[i,]=(ii+1)*((x-t[i])*yytem1[i,]+(t[i+ii+1]-x)*yytem1[i+1,])/(t[i+ii+1]-t[i])/ii
   }
   yytem1=yytem2
}

index=rep(0,length(x))
for (i in 1:length(x)){
    index[i]=sum(t<=x[i])
}

yy=array(rep(0,(n-1)*length(x)),dim=c(n-1,length(x)))

if (order==1){
   for (i in 2:n){
      yy[i-1,]=(i<index-order+1)+(i==index)*(t[i+order+1]-t[i])*yytem2[i,]/(order+1)
   }
}else{
   for (j in 1:length(x)){
      for (i in 2:n){
         if (i<(index[j]-order+1)){
            yy[i-1,j]=1
         }else if ((i<=index[j]) && (i>=(index[j]-order+1))){
            yy[i-1,j]=(t[(i+order+1):(index[j]+order+1)]-t[i:index[j]])%*%yytem2[i:index[j],j]/(order+1)
         }else{
            yy[i-1,j]=0
         }
      }
   }
}
return(yy)
}


### get Mspline bases ###
Mspline<-function(x,order,knots){

k1=order
m=length(knots)
n1=m-2+k1 # number of parameters
t1=c(rep(1,k1)*knots[1], knots[2:(m-1)], rep(1,k1)*knots[m]) # new knots

tem1=array(rep(0,(n1+k1-1)*length(x)),dim=c(n1+k1-1, length(x)))
for (l in k1:n1){
    tem1[l,]=(x>=t1[l] & x<t1[l+1])/(t1[l+1]-t1[l])
}

if (order==1){
   mbases=tem1
}else{
   mbases=tem1
   for (ii in 1:(order-1)){
      tem=array(rep(0,(n1+k1-1-ii)*length(x)),dim=c(n1+k1-1-ii, length(x)))
      for (i in (k1-ii):n1){
         tem[i,]=(ii+1)*((x-t1[i])*mbases[i,]+(t1[i+ii+1]-x)*mbases[i+1,])/(t1[i+ii+1]-t1[i])/ii
      }
      mbases=tem
   }
}
return(mbases)
}

####################################################################
positivepoissonrnd<-function(lambda){ 
   samp = rpois(1, lambda)
   while (samp==0) {
      samp = rpois(1, lambda)
   }
   return(samp)
}

At_fun <- function(x,j,At,betat,Omega,xx,te1,te2,sig0) 
{
  At[j]<-x
  tt<-sum(xx[,j]*x*te1)-sum(exp(xx%*%At+Omega%*%betat)*te2)-x^2/sig0^2/2  
  return(tt)     
}
ind_fun_At <- function(x,j,At,betat,Omega,xx,te1,te2,sig0) (x>-coef_range)*(x<coef_range)

betat_fun <- function(x,j,At,betat,Omega,xx,te1,te2,sig0) 
{
  betat[j]<-x
  tt<-sum(Omega[,j]*x*te1)-sum(exp(xx%*%At+Omega%*%betat)*te2)-x^2/sig0^2/2  
  return(tt)     
}
ind_fun_betat <- function(x,j,At,betat,Omega,xx,te1,te2,sig0) (x>-coef_range)*(x<coef_range)

#Ag_fun <- function(x,j,Ag,betag,Omega,xx,teG,sig0) 
#{
#  Ag[j]<-x
#  tt1<-sum(xx[,j]*x-teG*exp(xx%*%Ag+Omega%*%betag))-x^2/sig0^2/2
#  return(tt1)     
#}
#ind_fun_Ag <- function(x,j,Ag,betag,Omega,xx,teG,sig0) (x>-coef_range)*(x<coef_range)


#betag_fun <- function(x,j,Ag,betag,Omega,xx,teG,sig0) 
#{
#  betag[j]<-x
#  tt1<-sum(Omega[,j]*x-teG*exp(xx%*%Ag+Omega%*%betag))-x^2/sig0^2/2
#  return(tt1)     
#}
#ind_fun_betag <- function(x,j,Ag,betag,Omega,xx,teG,sig0) (x>-coef_range)*(x<coef_range)

  
#omega_fun <- function(x,At,betat,Ag,betag,Lamda,phy,psi,B,xx,Delta,mcen,te1,te2,teG) 
#{
#  tt2 <- x%*%betat*te1-exp(xx%*%At+x%*%betat)*te2
#  tt3 <- x%*%betag-teG*exp(xx%*%Ag+x%*%betag)
#  tt4 <- sum(-0.5*(1/psi)*(B-x%*%Lamda)^2)-0.5*(x-mcen%*%invPI0)%*%ISIG%*%t(x-mcen%*%invPI0)
#  tt  <- tt2+tt3+tt4
#  return(tt)     
#}
omega_fun <- function(x,At,betat,Lamda,phy,psi,B,xx,Delta,mcen,te1,te2) 
{
  tt2 <- x%*%betat*te1-exp(xx%*%At+x%*%betat)*te2
  #tt3 <- x%*%betag-teG*exp(xx%*%Ag+x%*%betag)
  tt4 <- sum(-0.5*(1/psi)*(B-x%*%Lamda)^2)-0.5*(x-mcen%*%invPI0)%*%ISIG%*%t(x-mcen%*%invPI0)
  tt  <- tt2+tt4
  return(tt)     
}
#x=Omega[i,]
#Mcen=Mcen[i,]
### main routine ###
   L  = matrix(d$t1,ncol=1,n)   ###data
   R  = matrix(d$t2,ncol=1,n)   ###data

   L  = ifelse(is.na(L),0,L)
   R2 = ifelse(is.na(R),0,R)

   status = matrix(d$delta,ncol=1)   ###data 

   order = 2       ###data
   ## generate basis functions
   knots   <- seq(min(c(L,R2)),max(c(L,R2))+0.001,length=kh)
   #knots_g <- seq(0,max(tG)+0.001,length=klam)
   kh      <- length(knots)-2+order
   #klam    <- length(knots_g)-2+order

   grids    <- seq(min(c(L,R2)),max(c(L,R2)),length=100)
   #grids_g  <- seq(0,max(tG),length=100)
   kgrids   <- length(grids)
   #kgrids_g <- length(grids_g)

   bisL  = Ispline(L,order,knots)
   bisR  = Ispline(R2,order,knots)
   bisTI = Ispline(L,order,knots)    ##knot
   bisT  = Mspline(L,order,knots)    ##knot   
   #bisG  = Ispline(tG,order,knots_g)
   #bisGM = Mspline(tG,order,knots_g)

   gridss_1 = seq(0,4,length=100)
   #gridss_2 = seq(0,4,length=100)
      bgs_1 = Ispline(gridss_1,order,knots)
      #bgs_2 = Ispline(gridss_2,order,knots_g)

###data
   burn_in    = 10000
   niter      = 20000
   coef_range = 5
   sig0       = 10
   a_eta      = 1
   b_eta      = 1
   a_pho      = 1
   b_pho      = 1

   eta = rgamma(1,a_eta,rate=b_eta)
   PHO = rgamma(1,a_pho,rate=b_pho)

   ## initial value

Delta <- matrix(rep(0, q*(1 + p)), nrow = q, ncol = (1 + p), byrow = TRUE)
L.se <- cbind(Delta)
psi <- rep(0.5, m) 
iv.psi <- 1/psi
iv.sqrt.psi <- sqrt(iv.psi)
phy =  rep(1, q) 
iv.phy <- 1/phy
iv.sqrt.phy <- sqrt(iv.phy)

if (q > 0) {
  const <- rep(1,n)
  Dat <- cbind(const, Z, S)
  Coe <- cbind(Delta)
  Omega <- Dat %*% t(Coe) + mvrnorm(n, rep(0, q), diag(phy, q))
  Omega <- Omega %*% t(solve((diag(1, q))))   
}
 
   gamcoef     = matrix(rgamma(kh, 1, 1),ncol=kh)
   alpha       = matrix(rgamma(klam, 1, 1),ncol=klam)
   At          = matrix(rep(1,p),ncol=1)
   #Ag          = matrix(rep(1,p),ncol=1)
   betat       = matrix(rep(1,q),ncol=1)
   #betag       = matrix(rep(1,q),ncol=1)
   Lamda       = matrix(c(
                 1,1,0,0,
                 0,0,1,1
                 ),ncol=m,byr=T)

    invPI0 <- solve(diag(1, q))   #dim(invPI0) 3   3  
    const <- rep(1,n)
    Dat1 <- cbind(const, Z, S)    #dim(Dat1)  200  7
    Coe1 <- cbind(Delta)          #dim(Delta) 3    7
    Mcen <- Dat1 %*% t(Coe1)      #dim(Mcen)  200  3
    SIG <- invPI0 %*% diag(phy) %*% t(invPI0)
    ISIG <- solve(SIG)
    mcen=Mcen%*%invPI0

   ##Lmabda ini
   IDY<-matrix(c(
    0,1,0,0,
    0,0,0,1
   ),ncol=m,byr=T)
   PLY<-matrix(c(
    1,0.5,0,0,
    0,0,1,0.5
   ),ncol=m,byr=T)

IDY <- (IDY > 0)
## free loadings in each row of the loading matrix B
n.b.row <- colSums(IDY)  
n.B <- sum(IDY)

###for the measurment equation
Id.me <- IDY
n.me.row <- colSums(Id.me)  
n.me <- sum(Id.me)

sig.b <- 0.001  
alpha.psi <- 9      #  ?_ ?
beta.psi <- 4       #  ?_ ? 

   LambdatL = t(gamcoef%*%bisL) # n x 1
   LambdatR = t(gamcoef%*%bisR) # n x 1
   LambdatT = t(gamcoef%*%bisTI) # n x 1
   #LambdatG = t(alpha%*%bisG) # n x 1
   #teG      = LambdatG
  
   parL.se  = array(0, dim = c(niter, n.se))
   parpsi   = array(0, dim = c(niter, m))
   parAt    = array(rep(0,niter*p),dim=c(niter,p))
   #parAg    = array(rep(0,niter*p),dim=c(niter,p))
   parbetat = array(rep(0,niter*q),dim=c(niter,q))
   #parbetag = array(rep(0,niter*q),dim=c(niter,q))
   pargam   = array(rep(0,niter*kh),dim=c(niter,kh))
   paralpha = array(rep(0,niter*klam),dim=c(niter,klam))
   pareta   = array(rep(0,niter),dim=c(niter,1))
   parpho   = array(rep(0,niter),dim=c(niter,1))
   parLamda = array(rep(0,niter*(m-2)),dim=c(niter,(m-2))) 
   parphy   = array(rep(0,niter*(q)),dim=c(niter,(q)))  
   paromega = array(0,dim=c(n,q,niter))
   parsurv_1= array(rep(0,niter*kgrids),dim=c(niter,kgrids))
   #parsurv_2= array(rep(0,niter*kgrids_g),dim=c(niter,kgrids_g))
   acc.prob = rep(0,n)
   
   ## iteration
   iter=1
   while (iter<niter+1)
   { 
    cat("loopp", loopp, "Iteration", iter, fill=TRUE);
      # sample z, zz, w and ww
      z=array(rep(0,n),dim=c(n,1)); w=z
      zz=array(rep(0,n*kh),dim=c(n,kh)); ww=zz

      for (i in 1:n){
         if (status[i]==0){
            templam1=LambdatR[i]*exp(X[i,]%*%At+Omega[i,]%*%betat)  ###缁楊兛绔村顨€i閻ㄥ嫬寮弫?
            z[i]=positivepoissonrnd(templam1)
            zz[i,]=rmultinom(1,z[i],gamcoef*t(bisR[,i]))
         }else if (status[i]==1){
            templam1=(LambdatR[i]-LambdatL[i])*exp(X[i,]%*%At+Omega[i,]%*%betat)
            w[i]=positivepoissonrnd(templam1)
            ww[i,]=rmultinom(1,w[i],gamcoef*t(bisR[,i]-bisL[,i]))
         }
      }

      # sample pai
      #pai=array(rep(0,n*klam),dim=c(n,klam)); 
      #for (i in 1:n){
      #      pai[i,]=rmultinom(1,1, alpha*t(bisGM[,i]))
      #}

      # sample At
      te1=z*as.numeric(status==0)+w*as.numeric(status==1)
      te2=LambdatR*as.numeric(status==0)+LambdatR*as.numeric(status==1)+LambdatL*as.numeric(status==2)

      for(j in 1:p){
      At[j]<-arms(At[j],At_fun,ind_fun_At,1,j=j,At=At,betat=betat,Omega=Omega,xx=X,te1=te1,te2=te2,sig0=sig0)
      }

      # sample betat
      for(j in 1:q){
      betat[j]<-arms(betat[j],betat_fun,ind_fun_betat,1,j=j,At=At,betat=betat,Omega=Omega,xx=X,te1=te1,te2=te2,sig0=sig0)
      }

      # sample Ag
      #for(j in 1:p){
      #Ag[j]<-arms(Ag[j],Ag_fun,ind_fun_Ag,1,j=j,Ag=Ag,betag=betag,Omega=Omega,xx=X,teG=teG,sig0=sig0)
      #}

      # sample betag
      #for(j in 1:q){
      #betag[j]<-arms(betag[j],betag_fun,ind_fun_betag,1,j=j,Ag=Ag,betag=betag,Omega=Omega,xx=X,teG=teG,sig0=sig0)
      #}
    
      # sample gamcoef
      for (l in 1:kh){
         tempa=1+sum(zz[,l]*as.numeric(status==0)*(bisR[l,]>0)+ww[,l]*as.numeric(status==1)*((bisR[l,]-bisL[l,])>0))
         tempb=eta+sum((bisR[l,]*as.numeric(status==0)+bisR[l,]*as.numeric(status==1)+bisL[l,]*as.numeric(status==2))*exp(X%*%At+Omega%*%betat)) 
         gamcoef[l]=rgamma(1,tempa,rate=tempb)
      }

      #for (l in 1:klam){
      #   alpha_a=1+sum(pai[,l])                                 ##1+sum(pai[,l]*(bisGM[l,]>0))
      #   alpha_b=PHO+sum(bisG[l,]*exp(X%*%Ag+Omega%*%betag)) 
      #   alpha[l]=rgamma(1,alpha_a,rate=alpha_b)
      #}

      LambdatL=t(gamcoef%*%bisL) 
      LambdatR=t(gamcoef%*%bisR)
      LambdatT=t(gamcoef%*%bisTI)
      #LambdatG=t(alpha%*%bisG)
      #teG     =LambdatG

      #sample eta
      eta=rgamma(1,a_eta+kh, rate=b_eta+sum(gamcoef))

      #sample PHO
      PHO=rgamma(1,a_pho+klam, rate=b_pho+sum(alpha))


  count.B <- 1
  L.me <- Lamda
  omg.me <- Omega
  for (k in 1:m) {
    free <- IDY[, k]                    # which column on row K is free/fixed
    len <- n.me.row[k]
    Vcen <- t(t(L.me[!free, k, drop = F]) %*% t(omg.me)[!free, , drop = F])
    Psiginv <- rep(sig.b, len)          # H_0yk
 
    Vk.star <- B[ , k] - Vcen
    alpha.psi.star <- alpha.psi + 0.5*n
    beta.psi.star <- beta.psi + 0.5*sum(Vk.star^2)
    if (len > 0) {
	if(len==1){Mk=matrix(omg.me[,free],nrow=1)}
      if(len>1){Mk<-omg.me[,free]}
      A_vk <- chol2inv(chol(diag(Psiginv, len) + tcrossprod(Mk)))
      temp <- Psiginv*PLY[free, k] + Mk %*% Vk.star
      a_vk <- A_vk %*% temp
      beta.psi.star <- beta.psi.star + 0.5*(sum(PLY[free, k]*Psiginv*PLY[free, k]) - sum(temp*a_vk))
    }

    iv.psi[k] <- rgamma(1, shape = alpha.psi.star, rate = beta.psi.star)
    psi[k] <- 1/iv.psi[k]                  
    iv.sqrt.psi[k] <- sqrt(iv.psi[k])

    if (len > 0) {
      L.me[free, k] <- mvrnorm(1, a_vk, psi[k]*A_vk)
      if (n.b.row[k] > 0) {
        Lamda[, k] <- L.me[, k]
      }
 
        if (n.b.row[k] > 0) {
          parLamda[iter, count.B:(count.B + n.b.row[k] - 1)] <- Lamda[IDY[, k], k]
        }
 
      count.B <- count.B + n.b.row[k]
    }   
  }
    Lamda_j = c(Lamda[1,2],Lamda[2,4])
 
  # step1: update latent mediators M
  if (m > 0) {
    PI0 <- (diag(1, m))
    ISG <- crossprod(iv.sqrt.psi * Lamda) + crossprod(iv.sqrt.phy * PI0)     # 閳叡锠區-1
    SIG <- chol2inv(chol(ISG))                     
    SIG <- sig.mh*SIG
    cSIG <- chol(SIG)
  }

      #sample omega
      for(i in 1:n){
        k_X=Omega[i,]
        pi_X=omega_fun(k_X,At=At,betat=betat,Lamda=Lamda,
                       phy=phy,psi=psi,B=B[i,],xx=X[i,],Delta=Delta,
                        mcen=mcen[i,],te1=te1[i],te2=te2[i]);
        k_Y=mvrnorm(1,k_X,0.2*diag(q));   # proposal distribution
        pi_Y=omega_fun(k_Y,At=At,betat=betat,Lamda=Lamda,
                       phy=phy,psi=psi,B=B[i,],xx=X[i,],Delta=Delta,
                        mcen=mcen[i,],te1=te1[i],te2=te2[i]);
        a_X_Y=((pi_Y)-(pi_X))
        if( log(runif(1,0,1)) < a_X_Y)
        {     Omega[i,] = k_Y
              acc.prob[i] <- acc.prob[i]+1
        }      
      }

####

  ## Note!! Need revise if the number of latent mediator changed!

  count.se <- 1
  L.se <- cbind(Delta)
  omg.se <- cbind(rep(1,n), Z, S)
  for (k in 1:q) {
    free <- Id.se[k, ]  # which column on row K is free/fixed
    len <- n.se.row[k]
    Mcen <- t(L.se[k, !free, drop = F] %*% t(omg.se)[!free, , drop = F])
    Mk.star <- Omega[ , k] - as.vector(Mcen)
    alpha.psd.star <- alpha.psd + 0.5*n
    beta.psd.star <- beta.psd + 0.5*sum(Mk.star^2)

    if (len > 0) {
      Yk <- omg.se[, free, drop = F]  
      iH0dk <- diag(sig.delta, len)
      L.se0k <- L.se0[k, free]
      A_dk <- chol2inv(chol(iH0dk + crossprod(Yk)))
      temp <- iH0dk %*% L.se0k + t(Yk) %*% Mk.star
      a_dk <- A_dk %*% temp
      beta.psd.star <- beta.psd.star + 0.5*(crossprod(crossprod(iH0dk, L.se0k), L.se0k)
                                    - crossprod(a_dk, temp))
    }

    iv.phy[k] <- rgamma(1, shape = alpha.psd.star, rate = beta.psd.star)
    phy[k] <- 1/iv.phy[k]
    iv.sqrt.phy[k] <- sqrt(iv.phy[k])
 
    if (len > 0) {
      L.se[k, free] <- mvrnorm(1, a_dk, phy[k]*A_dk)

       parL.se[iter, count.se:(count.se + len - 1)] <- L.se[k, free]
       count.se <- count.se + len
      }
      } 
   
L.se_e <- as.vector(L.se)
 
      
      parL.se[iter,]  = L.se_e	  #dim( parL.se)
      parpsi[iter,]   = psi																																
      pargam[iter,]   = gamcoef
      parAt[iter,]    = At
      #parAg[iter,]    = Ag
      parbetat[iter,] = betat
      #parbetag[iter,] = betag
      pareta[iter]    = eta
      parpho[iter]    = PHO
      parLamda[iter,] = Lamda_j       #dim(parLamda)
      parphy[iter,]   = phy
      paromega[,,iter]= Omega
      parsurv_1[iter,]= exp(-gamcoef%*%bgs_1)
      #parsurv_2[iter,]= exp(-alpha%*%bgs_2)
 
      iter=iter+1
   } 

   accept.rate[loopp,]<-(acc.prob/niter)

   loop_L.se[loopp,]   = colMeans(parL.se[burn_in:niter,])
   loop_psi[loopp,]    = colMeans(parpsi[burn_in:niter,])
   loop_At[loopp,]     = colMeans(parAt[burn_in:niter,])
   #loop_Ag[loopp,]     = colMeans(parAg[burn_in:niter,])
   loop_betat[loopp,]  = colMeans(parbetat[burn_in:niter,])
   #loop_betag[loopp,]  = colMeans(parbetag[burn_in:niter,])
   loop_gam[loopp,]    = colMeans(pargam[burn_in:niter,])
   loop_eta[loopp]     = mean(pareta[burn_in:niter])
   loop_pho[loopp]     = mean(parpho[burn_in:niter])
   loop_Lamda[loopp,]  = colMeans(parLamda[burn_in:niter,])
   loop_phy[loopp,]    = colMeans(parphy[burn_in:niter,])
  loop_omega[,,loopp] = apply(paromega[,,burn_in:niter],c(1,2),mean)

   loop_surv.med_1[loopp,] = colMeans(parsurv_1[burn_in:niter,])
   #loop_surv.med_2[loopp,] = colMeans(parsurv_2[burn_in:niter,])

   parsee_L.se[loopp,]    = apply(parL.se[burn_in:niter,], 2, sd)
   parsee_psi[loopp,]    = apply(parpsi[burn_in:niter,], 2, sd)
   parsee_At[loopp,]    = apply(parAt[burn_in:niter,], 2, sd)
   #parsee_Ag[loopp,]    = apply(parAg[burn_in:niter,], 2, sd)
   parsee_betat[loopp,] = apply(parbetat[burn_in:niter,], 2, sd)
   #parsee_betag[loopp,] = apply(parbetag[burn_in:niter,], 2, sd)
   parsee_Lamda[loopp,] = apply(parLamda[burn_in:niter,], 2, sd)
   parsee_phy[loopp,]   = apply(parphy[burn_in:niter,], 2, sd)

   for(i in 1:kgrids){
      loop_surv.low_1[loopp,i] <- quantile(parsurv_1[burn_in:niter,i],  0.025)
      loop_surv.upp_1[loopp,i] <- quantile(parsurv_1[burn_in:niter,i],  0.975)
   }

   #for(i in 1:kgrids_g){
   #   loop_surv.low_2[loopp,i] <- quantile(parsurv_2[burn_in:niter,i],  0.025)
   #   loop_surv.upp_2[loopp,i] <- quantile(parsurv_2[burn_in:niter,i],  0.975)
   #}

   mc_effect <- 50
   std_norm_effect <- array(rnorm(n * q * mc_effect), dim = c(n, q, mc_effect))

   gam_hat <- loop_gam[loopp,]
   At_hat <- loop_At[loopp,]
   betat_hat <- loop_betat[loopp,]
   Delta_hat <- loop_L.se[loopp,]
   phy_hat <- loop_phy[loopp,]

   H0_t_L_hat <- as.vector(t(gam_hat %*% bisL))
   H0_t_R_hat <- as.vector(t(gam_hat %*% bisR))
   H0_t_hat <- make_H0_y_interval(H0_t_L_hat, H0_t_R_hat, status)
   H0_t_true <- make_H0_y_interval(S0oft(L), S0oft(R2), status)

   z_contrasts <- list(c(1, 0), c(0, 1))
   for (cc in seq_along(z_contrasts)) {
      z_trt <- z_contrasts[[cc]]

      loop_effect_t_hat[loopp, cc, ] <- calc_effect_q2(
         z_trt = z_trt, Zbase = Z, H0 = H0_t_hat, gamma = At_hat,
         Delta_vec = Delta_hat, betaM = betat_hat, phy_vec = phy_hat,
         std_norm = std_norm_effect
      )
      loop_effect_t_true[loopp, cc, ] <- calc_effect_q2(
         z_trt = z_trt, Zbase = Z, H0 = H0_t_true, gamma = AT,
         Delta_vec = L.se.true, betaM = betaT, phy_vec = phy.true,
         std_norm = std_norm_effect
      )
   }

#plot(grids,exp(-S0oft(grids)),lty=1,col="yellow") 
#lines(grids,exp(-colMeans(pargam)%*%bgs),lty=3,col="red")

#plot(grids_g,exp(-S0ofg(grids_g)),lty=1,col="yellow") 
#lines(grids_g,exp(-colMeans(alpha)%*%bgs_g),lty=3,col="red")


}
 
effect_summary_t <- summarize_effects(loop_effect_t_hat, loop_effect_t_true, "Y")
effect_summary <- effect_summary_t

Lamda.true  = c(0.9,0.8)

bias.L.se     = as.matrix(colMeans(loop_L.se)-L.se.true)
bias.psi     = as.matrix(colMeans(loop_psi)-psi.true)
bias.AT     = as.matrix(colMeans(loop_At)-AT)
#bias.AG     = as.matrix(colMeans(loop_Ag)-AG)
bias.betaT  = as.matrix(colMeans(loop_betat)-betaT)
#bias.betaG  = as.matrix(colMeans(loop_betag)-betaG)
bias.Lamda  = as.matrix(colMeans(loop_Lamda)-Lamda.true)
bias.phy    = as.matrix(colMeans(loop_phy)-phy.true)


RMSE.L.se=matrix(NA,n.se,1)
for(i in 1:n.se){
RMSE.L.se[i,1] = sqrt(mean((L.se.true[i]-loop_L.se[,i])**2))
}
SEE.L.se = t(t(apply((parsee_L.se),2,mean)))


RMSE.psi=matrix(NA,m,1)
for(i in 1:m){
RMSE.psi[i,1] = sqrt(mean((psi.true[i]-loop_psi[,i])**2))
}
SEE.psi = t(t(apply((parsee_psi),2,mean)))


RMSE.AT=matrix(NA,p,1)
#RMSE.AG=matrix(NA,p,1)
for(i in 1:p){
RMSE.AT[i,1] = sqrt(mean((AT[i]-loop_At[,i])**2))
#RMSE.AG[i,1] = sqrt(mean((AG[i]-loop_Ag[,i])**2))
}
SEE.AT = t(t(apply((parsee_At),2,mean)))
#SEE.AG = t(t(apply((parsee_Ag),2,mean)))

RMSE.betaT=matrix(NA,q,1)
#RMSE.betaG=matrix(NA,q,1)
for(i in 1:q){
RMSE.betaT[i] = sqrt(mean((betaT[i]-loop_betat[,i])**2))
#RMSE.betaG[i] = sqrt(mean((betaG[i]-loop_betag[,i])**2))
}
SEE.betaT = t(t(apply((parsee_betat),2,mean)))
#SEE.betaG = t(t(apply((parsee_betag),2,mean)))

RMSE.phy =matrix(NA,q,1)
for(i in 1:(q)){
RMSE.phy[i] = sqrt(mean((phy.true[i]-loop_phy[,i])**2))
}
SEE.phy = t(t(apply((parsee_phy),2,mean)))

RMSE.Lamda =matrix(NA,m-2,1)
for(i in 1:(m-2)){
RMSE.Lamda[i] = sqrt(mean((Lamda.true[i]-loop_Lamda[,i])**2))
}
SEE.Lamda = t(t(apply((parsee_Lamda),2,mean)))
 
######################################################
final.bias = rbind(bias.AT,bias.betaT,bias.Lamda,bias.psi)
final.RMSE = rbind(RMSE.AT,RMSE.betaT,RMSE.Lamda,RMSE.psi)
final.SEE  = rbind(SEE.AT, SEE.betaT,SEE.Lamda,SEE.psi)
final      = cbind(final.bias,final.RMSE,final.SEE)
dim2       = c("BIAS","RMSE","SEE")

vector1 <- c()
for (i in 1:m) {
  vector1 <- c(vector1, paste("psi", i, sep = ""))
}

dim1       = c("AT1","AT2","AT3","AT4","AT5","AT6",
               "betaT1","betaT2","Lamda1","Lamda2",vector1)
final      = matrix(final,nrow=14,ncol=3,byrow=FALSE,dimnames=list(dim1,dim2))

######################################################
final2.bias = rbind(bias.L.se,bias.phy)
final2.RMSE = rbind(RMSE.L.se,RMSE.phy)
final2.SEE  = rbind(SEE.L.se,SEE.phy)
final2      = cbind(final2.bias,final2.RMSE,final2.SEE)
dim4      = c("BIAS","RMSE","SEE")

vector <- c()
for (i in 1:14) {
  vector <- c(vector, paste("L.se", i, sep = ""))
}


dim3       = c(vector,"phy1","phy2")
final2      = matrix(final2,nrow=16,ncol=3,byrow=FALSE,dimnames=list(dim3,dim4))

if (any(final[, "RMSE"] + 1e-12 < abs(final[, "BIAS"])) ||
    any(final2[, "RMSE"] + 1e-12 < abs(final2[, "BIAS"]))) {
  stop("Invalid summary: RMSE is smaller than absolute bias.")
}
if (any(effect_summary$RMSE + 1e-12 < abs(effect_summary$Bias))) {
  stop("Invalid effect summary: RMSE is smaller than absolute bias.")
}

scenario.tag <- paste0("r0-n", n, "-independent")
write.csv(final2,paste0("sem-structural-",scenario.tag,".csv"))
write.csv(final,paste0("sem-other-parameters-",scenario.tag,".csv"))
write.csv(effect_summary_t,paste0("effect-Y-",scenario.tag,".csv"), row.names = FALSE)
write.csv(effect_summary,paste0("effect-all-",scenario.tag,".csv"), row.names = FALSE)
save.image(paste0("simulation-",scenario.tag,"-loop",loop,".RData"))


#var.omega=apply(loop_omega,c(1,2),mean)
#var.omega
#plot(parL.se[,1],type="l")
#plot(parpsi[,1],type="l")
#plot(parphy[,1],type="l")
#plot(parLamda[,2],type="l")
#plot(parAt[,2],type="l")
#plot(parbetat[,2],type="l")
#plot(parAg[,2],type="l")
#plot(parbetag[,2],type="l")

Sys.time()-t0;
#dev.off()

#write.csv(pd,"C:\\Users\\Administrator\\Desktop\\data\\r0.5-p2-ex40-n400-6.csv")
##################plot
#plot(gridss,exp(-S0oft(gridss)),xlim=c(0.1,1),ylim=c(-1.5,1.5),type="l",col="black") 
#lines(gridss,colMeans(loop_surv.med),lty=2,col="red")
#lines(gridss,apply(loop_surv.low,2,mean),lty=2,col="black")
#lines(gridss,apply(loop_surv.upp,2,mean),lty=2,col="black")
#write.csv(pd,"C:\\Users\\Administrator\\Desktop\\data\\r0.5-p2-ex40-n400-6.csv")
