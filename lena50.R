
get.mat = function(x){
  return(do.call(rbind, lapply(x,as.vector)) )
}
get.mat <- compiler::cmpfun(get.mat)
P_Omega = function(a,entri){
  a[entri] = 0
  return(a)
}
P_Omega <- compiler::cmpfun(P_Omega)

library(filling)
missfrac = 0.5
x = lena256
## transform 5% of entries into missing
gibb = lmc = mala = softimpt = list()
for(ss in 1:30){
  xna <- aux.rndmissing(lena256, x=missfrac)
  lamlam = seq(100,10)
  fill <- fill.SoftImpute(xna, lambdas= lamlam,
                          tol = 1e-4,maxiter = 1e3)
  rmse<-sapply(1:length(lamlam),function(s) mean((fill$X[,,s] - x)^2) ) 
  I1 = row(xna)[!is.na(xna)]
  I2 = col(xna)[!is.na(xna)]
  Y = xna[!is.na(xna)]
  m = nrow(xna)
  p = ncol(xna)
  nsample =  sum(!is.na(xna))
  obs = Y
  sigma = 1
  mp = m*p
  imiss = (1:mp)[is.na(xna)]
  
  # MALA MC for MC  
  a = 0
  sig2 = 1
  ystar = diag(nrow(xna))
  Bm = matrix(data=0,nr=m,nc=p)
  Iters = 200
  burnin = 10
  h =  1/p/m/50
  lam = 1
  M = fill$X[,,which.min(rmse)]
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
  #Langevin MC for BRRR
  h =  1/p/m/50
  sig2 = 1
  M = fill$X[,,which.min(rmse)]
  eta = 2
  lam = 1
  ystar = diag(nrow(xna))
  X.lmc = matrix(data=0,nr=m,nc=p)
  Iters = 200
  burnin = 10
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
    print(s)
    if (s>burnin) X.lmc = X.lmc + M/(Iters-burnin)
  }
 
  # Gibbs sampler
  I1 = row(xna)[!is.na(xna)]
  I2 = col(xna)[!is.na(xna)]
  Y = xna[!is.na(xna)]
  m = nrow(xna)
  p = ncol(xna)
  nsample =  sum(!is.na(xna))
  obs = Y
  sigma = 1
  
  k = 10 #rank of U,V
  a = 1
  b = 1/100
  lambda = 1/(4*sigma^2)
  Mstep = matrix(1 ,nr=m,nc=k)
  Nstep = matrix(1 ,nr=p,nc=k)
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
  lmc[[ss]] = c(mean((X.lmc - x )^2),Matrix::rankMatrix(X.lmc ,.1)[1], mean((X.lmc - x )^2)/mean(x^2) ,mean((X.lmc[imiss] - x[imiss] )^2) )
  mala[[ss]] = c(mean((Bm - x )^2),Matrix::rankMatrix(Bm,.1)[1],  mean((Bm - x )^2)/mean(x^2), mean((Bm[imiss] - x[imiss] )^2)  )
  softimpt[[s]] = c(mean((fill$X[,,which.min(rmse)] - x )^2),Matrix::rankMatrix(fill$X[,,which.min(rmse)],.1)[1],  mean((fill$X[,,which.min(rmse)] - x )^2)/mean(x^2), mean((fill$X[,,which.min(rmse)][imiss] - x[imiss] )^2)  )
}
save(gibb,lmc,mala,softimpt,file = '/data2/thetm/matrixCompletion/fastSampling/lena50out.rda')

colMeans(do.call(rbind,gibb) )
apply(do.call(rbind,gibb), 2, sd)
###
colMeans(do.call(rbind,lmc))
apply(do.call(rbind,lmc), 2, sd)
###
colMeans(do.call(rbind,mala))
apply(do.call(rbind,mala), 2, sd)


