# rm(list = ls())

setwd("C:/Users/klebe/Dropbox/PIBIC/Script")
source("gamlss_BP.r")
source("glmBP.r")
source("Estimation.r")
source("simu_BP_v11.r")

library(gamlss)
library(gamlss.dist)
library(gamlss.data)
library(splines)

set.seed(144)

n=50          # tamanho amostra
j=4           # indice para usar tamanhos maiores 
NR=10000      # Numero de replicas de MC
pw=3          # potencia variavel de teste
n1=n*j        # Quando n for 100, n1=n*j=50*2

beta=c(1.9,-1.5) #valores reais vetor beta
x0=rep(1,n)
x01=runif(n)
x02=runif(n)

x0=rep(x0,j)
x01=rep(x01,j)
x02=rep(x02,j)
X=cbind(x0,x01) #matriz de regressores de mu
eta1=X%*%beta   #preditor linear
mu = exp(eta1)
summary(mu)


# Precisão pequena, a média proxima de 5
# lambda=c(1.8,-0.4)
# z0=rep(1,n)
# z01=runif(n)
# z02=runif(n)
# 
# z0=rep(z0,j)
# z01=rep(z01,j)
# z02=rep(z02,j)
# Z=cbind(z0,z01)
# kappa=Z%*%lambda  #preditor linear
# phi = exp(kappa)
# summary(phi)

# Precisão grande, a média proxima de 60
lambda=c(3.5,1)
z0=rep(1,n)
z01=runif(n)
z02=runif(n)

z0=rep(z0,j)
z01=rep(z01,j)
z02=rep(z02,j)
Z=cbind(z0,z01)
kappa=Z%*%lambda  #preditor linear
phi = exp(kappa)
summary(phi)

Vmu=mu*(1+mu)
vary= Vmu/phi
summary(vary)

#Inicializando vetores para guardar estatisticas de teste
Etest1=rep(0,NR)
Etest2=rep(0,NR)
Etest3=rep(0,NR)
Etest4=rep(0,NR)

contC=0
contC1=0

# PRIMEIRO LAçO DE MONTE CARLO

for(i in 1:NR) {
  
  y=rBP(n1,mu,phi) 
  ystar=ifelse(y==1,0,log(y/(1+y))) 
  
  fit <- gamlss(y~(X-1),sigma.formula=~(Z-1),family = BP(mu.link = "log",sigma.link = "log"), 
                trace=FALSE) 
  
  
  if(fit$conv==FALSE){
    
    i=i-1
    #print("i - CONVERGENCIA")
    #print(i)
    #NREP = NREP+1
    contC=contC+1
    next
    
  }
  
  # ADIÇAO DAS NOVAS VARIAVEIS DE TESTE
  
  etatest=fit$mu.lp         #Preditor linear
  mutest=fitted.values(fit) #Valor ajustado
  x1=mutest^pw 
  x1=as.vector(x1)
  Xnew=cbind(X,x1)          #Nova matriz de regressores de mu
  
  
  fitnew <- gamlss(y~(Xnew-1),sigma.formula=~(Z-1),family = BP(mu.link = "log",sigma.link = "log"), 
                   trace=FALSE)
  
  if (fitnew$conv==FALSE){
    
    i=i-1
    contC1=contC1+1
    next
    
  }
  
  # FUNÇÃO ESCORE VEROSSIMILHANÇA RESTRITA E INVERSA DE FISHER RESTRITA
  
  mx1=2
  mx2=1
  mx3=3
  zerox=rep(0,mx2)
  betatil=c(fit$mu.coef, zerox) 
  etatil1=Xnew%*%betatil        
  mutil=exp(etatil1)
  
  ms1=2
  ms2=1
  lambdatil=c(fit$sigma.coef) 
  etatil2=Z%*%lambdatil
  phitil=exp(etatil2)
  
  Phitil=diag(as.vector(1+phitil))
  D1 <- diag(as.vector(exp(etatil1)))
  ystartil=ystar
  mustartil=digamma(mutil*(1+phitil))-digamma(mutil*(1+phitil)+phitil+2)
  escorebetatil=t(Xnew)%*%Phitil%*%D1%*%(ystartil-mustartil)
  
  # D2 <- diag(as.vector(exp(etatil2)))
  # ystar1 <- ifelse(y==1,0,log((y^mutil)/((1+y)^(1+mutil))))
  # gama <- digamma(mutil*(1+phitil)+phitil+2)-digamma(phitil+2)
  # mustar1 <- mutil*(mustartil-(gama/(mutil)))
  # escorelambdatil <- t(Z)%*%D2%*%(ystar1-mustar1)
  
  
  # INVERSA DE FISHER RESTRITA
  
  deltatil=1+as.vector(phitil)
  Deltatil=as.vector(deltatil)
  
  psi1til=trigamma(mutil*(1+phitil))
  psi2til=trigamma(mutil*(1+phitil)+phitil+2)
  wtil=psi1til-psi2til
  Wtil=as.vector(wtil)
  
  psi3til=trigamma(phitil+2)
  ctil=(mutil^2)*psi1til-((1+mutil)^2)*psi2til+psi3til
  Ctil=as.vector(ctil)
  
  dtil=psi2til+mutil*(psi2til-psi1til)
  Dtil=as.vector(dtil)
  
  e3til=(Deltatil^2)*Wtil*as.vector(exp(etatil1))^2
  e4til=Ctil*as.vector(exp(etatil2))^2
  e5til=(-Deltatil)*Dtil*as.vector(exp(etatil1))*as.vector(exp(etatil2))
  
  K.bbtil=t(Xnew)%*%diag(as.vector(e3til))%*%Xnew
  K.bltil=t(Xnew)%*%diag(as.vector(e5til))%*%Z
  K.lbtil=t(Z)%*%diag(as.vector(e5til))%*%Xnew
  K.lltil=t(Z)%*%diag(as.vector(e4til))%*%Z
  
  bloc.suptil=cbind(K.bbtil,K.bltil)
  bloc.inftil=cbind(K.lbtil,K.lltil)
  
  Fisher.bltil=rbind(bloc.suptil,bloc.inftil)
  
  Inver.Fisher.bltil <- try(solve(Fisher.bltil))
  
  
  
  # INVERSA DE FISHER IRRESTRITA
  
  muhat=fitnew$mu.fv
  phihat=fitnew$sigma.fv
  
  eta1hat=fitnew$mu.lp
  eta2hat=fitnew$sigma.lp
  
  betahat=fitnew$mu.coefficients
  lambdahat=fitnew$sigma.coefficients
  
  deltahat=1+phihat
  Delta=as.vector(deltahat)
  
  psi1=trigamma(muhat*(1+phihat))
  psi2=trigamma(muhat*(1+phihat)+phihat+2)
  w=psi1-psi2
  W=as.vector(w)
  
  psi3=trigamma(phihat+2)
  c=(muhat^2)*psi1-((1+muhat)^2)*psi2+psi3
  C=as.vector(c)
  
  d=psi2+muhat*(psi2-psi1)
  D=as.vector(d)
  
  e3= (Delta^2)*W*(as.vector(exp(eta1hat))^2)
  e4= C*(as.vector(exp(eta2hat))^2)
  e5= -Delta*D*as.vector(exp(eta1hat))*as.vector(exp(eta2hat))
  
  
  K.bb=t(Xnew)%*%diag(as.vector(e3))%*%Xnew
  K.bl=t(Xnew)%*%diag(as.vector(e5))%*%Z
  K.lb=t(Z)%*%diag(as.vector(e5))%*%Xnew
  K.ll=t(Z)%*%diag(as.vector(e4))%*%Z
  
  bloc.sup=cbind(K.bb,K.bl)
  bloc.inf=cbind(K.lb,K.ll)
  
  Fisher.bl=rbind(bloc.sup,bloc.inf)
  
  Inver.Fisher.bl <- try(solve(Fisher.bl))
  
  
  # ESTATISTICAS DOS TESTES
  
  # RAZAO DE VEROSSIMILHANÇA
  
  lhat<-logLik(fitnew)
  ltil<-logLik(fit)
  
  Etest1[i]=2*(lhat[1]-ltil[1])
  
  # ESCORE
  
  Etest2[i]=  t(escorebetatil[mx3])%*%(Inver.Fisher.bltil[mx3,mx3])%*%escorebetatil[mx3]
  
  
  # WALD
  
  Etest3[i]= t(betahat[mx3])%*%((Inver.Fisher.bl[mx3,mx3])^(-1))%*%betahat[mx3]
  
  
  # Gradiente
  
  Etest4[i]= t(escorebetatil[mx3])%*%betahat[mx3]
  
  qchisq90=qchisq(0.90,1)
  qchisq95=qchisq(0.95,1)
  qchisq99=qchisq(0.99,1)
  
  
  # FIM DO LAÇO DE MONTE CARLO
}

print(contC)
print(contC1)
print(max(phi)/min(phi))
print(max(vary)/min(vary))

cont90RV=ifelse(Etest1 > qchisq90,1,0)
result_tam90RV=((sum(cont90RV))/NR)*100

cont90E=ifelse(Etest2 > qchisq90,1,0)
result_tam90E=((sum(cont90E))/NR)*100

cont90W=ifelse(Etest3 > qchisq90,1,0)
result_tam90W=((sum(cont90W))/NR)*100

cont90G=ifelse(Etest4 > qchisq90,1,0)
result_tam90G=((sum(cont90G))/NR)*100


cont95RV=ifelse(Etest1 > qchisq95,1,0)
result_tam95RV=((sum(cont95RV))/NR)*100

cont95E=ifelse(Etest2 > qchisq95,1,0)
result_tam95E=((sum(cont95E))/NR)*100

cont95W=ifelse(Etest3 > qchisq95,1,0)
result_tam95W=((sum(cont95W))/NR)*100

cont95G=ifelse(Etest4 > qchisq95,1,0)
result_tam95G=((sum(cont95G))/NR)*100


cont99RV=ifelse(Etest1 > qchisq99,1,0)
result_tam99RV=((sum(cont99RV))/NR)*100

cont99E=ifelse(Etest2 > qchisq99,1,0)
result_tam99E=((sum(cont99E))/NR)*100

cont99W=ifelse(Etest3 > qchisq99,1,0)
result_tam99W=((sum(cont99W))/NR)*100

cont99G=ifelse(Etest4 > qchisq99,1,0)
result_tam99G=((sum(cont99G))/NR)*100


print("resultRV")
print(result_tam90RV)
print(result_tam95RV)
print(result_tam99RV)

print("resultE")
print(result_tam90E)
print(result_tam95E)
print(result_tam99E)

print("resultW")
print(result_tam90W)
print(result_tam95W)
print(result_tam99W)

print("resultG")
print(result_tam90G)
print(result_tam95G)
print(result_tam99G)

# QUANTIS ESTIMADOS

quantilRV=quantile(Etest1, probs=c(.90,.95,.99))
quantilE=quantile(Etest2, probs=c(.90,.95,.99))
quantilW=quantile(Etest3, probs=c(.90,.95,.99))
quantilG=quantile(Etest4, probs=c(.90,.95,.99))

print("quantilRV")
print(quantilRV)

print("quantilE")
print(quantilE)

print("quantilW")
print(quantilW)

print("quantilG")
print(quantilG)

d<-date()
print(d)

resultRV=cbind(result_tam90RV,result_tam95RV,result_tam99RV)
resultE=cbind(result_tam90E,result_tam95E,result_tam99E)
resultW=cbind(result_tam90W,result_tam95W,result_tam99W)
resultG=cbind(result_tam90G,result_tam95G,result_tam99G)

result100=rbind(resultRV,resultE,resultW,resultG) # resultRV
# write.table(result100, file = "Mu3EtatestTam200.txt", sep = " ")

qestimado100=rbind(quantilRV,quantilE,quantilW,quantilG) # quantilRV,quantilE,quantilW
# write.table(qestimado100, file = "qestimMu3Etatest200.txt", sep = " ")
