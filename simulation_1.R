
get.mat = function(x){
  return(do.call(rbind, lapply(x,as.vector)) )
}
get.mat <- compiler::cmpfun(get.mat)
P_Omega = function(a,entri){
  a[entri] = 0
  return(a)
}
P_Omega <- compiler::cmpfun(P_Omega)
library(softImpute)
m = 100
p = 100
J = 5
mp = m*p
missfrac = 0.8

###uses regular matrix method for matrices with NAs
sig = 1
sigma = 1

gibb = mala = lmc = list()
for(ss in 1:50){
  x = matrix(rnorm(m*J),nr = m)%*%matrix(rnorm(J*p),nc = p) #+ 0.1*matrix(rnorm(m*50),nr=m)%*%matrix(rnorm(50*p),nc=p) 
  imiss = sample(mp,mp*missfrac, replace = FALSE)
  xna = x + matrix(rnorm(mp),nr = m,nc = p)
  xna[imiss]=NA
  I1 = row(xna)[!is.na(xna)]
  I2 = col(xna)[!is.na(xna)]
  Y = xna[!is.na(xna)]
  nsample =  (1-missfrac)*m*p
  obs = Y

# Gibbs sampler
k = 10 #rank of U,V
a = 1
b = 1/100
lambda = 1/(4*sigma^2)
Mstep = matrix(1 ,nr=m,nc=k)
Nstep = matrix(1,nr=p,nc=k)
gamma = rep(b/a,k)
X.Gibbs = matrix(data=0,nr=m,nc=p)
L2 = 2*lambda

Nmcmc = 200
burnin = 100
datalength = as.vector(1:nsample)
for (step in 1:(Nmcmc+burnin)){  
  # update M[i,j]
  for (i in 1:m) {
    seti = datalength[I1==i]
    seti = seti[!is.na(seti)]
    for (j in 1:k){
      Msteptrouj = Mstep[i,]
      Msteptrouj[j] = 0
      V = (1/gamma[j]) + L2*sum(Nstep[I2[seti],j]^2)
      D = sum(L2*(obs[seti] - Msteptrouj%*%t(Nstep[I2[seti],]) )*Nstep[I2[seti],j])
      Mstep[i,j] = rnorm(1,D/sqrt(V),1)/sqrt(V)
    }
  }
  
  # update N[i,j]
  for (i in 1:p) {
    seti = (1:nsample)[I2==i]
    for (j in 1:k){
      Nsteptrouj = Nstep[i,]
      Nsteptrouj[j] = 0
      V = (1/gamma[j]) + L2*sum(Mstep[I1[seti],j]^2)
      D = sum(L2*(obs[seti] - Nsteptrouj%*%t(Mstep[I1[seti],]))*Mstep[I1[seti],j])
      Nstep[i,j] = rnorm(1,D/sqrt(V),1)/sqrt(V)
    }
  }
  
  # update gamma
  for (j in 1:k) gamma[j] = 1/rgamma(1,a+(m+p)/2,b+(sum(Mstep[,j]^2)+sum(Nstep[,j]^2))/2)
  l1 <- list(X.Gibbs*(1-1/(step-burnin)), Mstep%*%t(Nstep)/(step-burnin) )
  res <- Reduce(`+`, lapply(l1, function(x)replace(x, is.na(x), 0)))
  X.Gibbs = res*NA^!Reduce(`+`, lapply(l1, function(x) !is.na(x)))
}
gibb[[ss]] = c(mean((X.Gibbs - x )^2),Matrix::rankMatrix(X.Gibbs,.1)[1],  mean((X.Gibbs - x )^2)/mean(x^2) ,mean((X.Gibbs[imiss] - x[imiss] )^2))

# MALA MC for MC  
a = 0
sig2 = 1
ystar = diag(nrow(xna))
Bm = matrix(data=0,nr=m,nc=p)
Iters = 100
burnin = 10
h =  1/p/m/400
lam = 1
M = X.Gibbs
eta = 2
for(s in 1:Iters){
  tam1 = glmnet::glmnet(M, ystar,
                        family = 'mgaussian',
                        thresh = 1e-3,
                        alpha = 0,lambda = lam^2,
                        standardize.response = F,
                        intercept = F)$beta
  tam = M - h*P_Omega(M-xna,imiss)/sig2/eta - 
    h*(m+p+2)*do.call(rbind, lapply(tam1,as.vector))+
    sqrt(2*h)*matrix(rnorm(mp),nr = m)
  
  logdet = determinant(lam^2*diag(m)+tam%*%t(tam))
  pro.tam = -sum(P_Omega(tam-xna,imiss)^2)/(eta*sig2) - 0.5*(p+m+2)*logdet$modulus*logdet$sign
  logdet = determinant(lam^2*diag(m)+M%*%t(M))
  pro.M = -sum(P_Omega(M-xna,imiss)^2)/(eta*sig2) - 0.5*(p+m+2)*logdet$modulus*logdet$sign
  
  tam2 = glmnet::glmnet(tam,ystar, family = 'mgaussian',
                        alpha = 0,lambda = lam^2,
                        thresh = 1e-3,
                        standardize.response = F,
                        intercept = F)$beta
  tran.m = -sum((M-tam-h*P_Omega(tam-xna,imiss)/sig2 - h*(p+m+2)*do.call(rbind, lapply(tam2,as.vector)) )^2)/(4*h)
  tran.tam = -sum((tam-M-h*P_Omega(M-xna,imiss)/sig2 - h*(p+m+2)*do.call(rbind, lapply(tam1,as.vector)) )^2)/(4*h)
  
  pro.trans = pro.tam+tran.m-pro.M-tran.tam
  if(log(runif(1)) < pro.trans){
    M = tam
    a = a+1
  } 
  if (s>burnin){
    Bm = Bm + M/(Iters-burnin)
  } 
}
mala[[ss]]= c(mean((Bm - x )^2),Matrix::rankMatrix(Bm,.1)[1],  mean((Bm - x )^2)/mean(x^2), mean((Bm[imiss] - x[imiss] )^2)  )

#Langevin MC for BRRR
h = 1/p/m/400
sig2 = 1
M = X.Gibbs
eta = 2
lam = 1
ystar = diag(nrow(xna))
X.lmc = matrix(data=0,nr=m,nc=p)
Iters = 100
burnin = 10
start_time <- Sys.time()
for(s in 1:Iters){
  tam = glmnet::glmnet(M, ystar,
                       family = 'mgaussian',
                       thresh = 1e-3,
                       alpha = 0,lambda = lam^2,
                       standardize.response = F,
                       intercept = F)$beta
  M = M - h*P_Omega(xna-M,imiss)/sig2/eta - 
    h*(m+p+2)*get.mat(tam) +
    sqrt(2*h)*matrix(rnorm(mp),nr = m)
  if (s>burnin) X.lmc = X.lmc + M/(Iters-burnin)
}
lmc[[ss]] = c(mean((X.lmc - x )^2),Matrix::rankMatrix(X.lmc ,.1)[1], mean((X.lmc - x )^2)/mean(x^2) ,mean((X.lmc[imiss] - x[imiss] )^2) )
print(ss)
}

save.image(file = '/data2/thetm/matrixCompletion/fastSampling/r5pp100miss08.rda')





colMeans(do.call(rbind,gibb) )
apply(do.call(rbind,gibb), 2, sd)
###
colMeans(do.call(rbind,lmc))
apply(do.call(rbind,lmc), 2, sd)
###
colMeans(do.call(rbind,mala))
apply(do.call(rbind,mala), 2, sd)



