rm(list=ls(all=TRUE)) 
library(MASS)
library(mvtnorm)
library(Rcpp)
setwd('C:/Users/lcgd4/iCloudDrive/Documents/Mizzou/Dissertation 1 - plm/Programs/cpp')
#setwd('/Users/lichen/Dropbox/Project/cpp')
sourceCpp("cpart_2.cpp")

####SIMULATE THE DATA

Nsample=1000
n=100
h.can=seq(n^(-.6), n^(-.5), length=10)
ha=n^-.4
Func='sin'
Btsp=1
Iter=100
xi=rep(1,n)

#h=h[1]

ny=rpois(n,5)+1  # # of y obs
nx=rpois(n,5)+1  # # of x obs

#ny=rep(1,n)
#nx=rep(1,n)

ty=as.list(NA)
tx=as.list(NA)
corr.x=as.list(NA)
x=as.list(NA)
y=as.list(NA)
corr.e=as.list(NA)

K1=function(x,h){
  y=0.75*(1-(x/h)^2)/h*(abs(x)<h)
  #y=dnorm(x/h)/h
  return(y)
}

K2=function(t,s,h1,h2){
  #temp.t=1-(t/h1)^2
  #temp.s=1-(s/h2)^2
  #if(temp.t<0) temp.t=0
  #if(temp.s<0) temp.s=0
  #y=0.5625*temp.t*temp.s/h1/h2
  y=1/(2*pi)*exp(-((t/h1)^2+(s/h2)^2)/2)/h1/h2
  return(y)
}

U.ini=function(b,h){
  U.dot=0
  U=0
  for(i in 1:n){
    for(j in 1:length(ty[[i]])){
      for(k in 1:length(tx[[i]])){
        U.dot=U.dot-x[[i]][k]^2*K1(ty[[i]][j]-tx[[i]][k],h)
        U=U+x[[i]][k]*(y[[i]][j]-b*x[[i]][k])*K1(ty[[i]][j]-tx[[i]][k],h)
      }
    }
  }
  return(c(U,U.dot))
}

b.sample=rep(NA,Nsample)
a.sample=matrix(NA,nrow=Nsample,ncol=101)
se.sample=rep(NA,Nsample)
aSE.sample=matrix(NA,nrow=Nsample,ncol=3)
a1.sample=matrix(NA,nrow=Nsample,ncol=101)
aBtsp=matrix(0,ncol=3,nrow=Nsample)
sample=1
na.count=0

h=rep(0, Nsample)

start=proc.time()

repeat{
for(i in 1:n){
  ty[[i]]=as.matrix(runif(ny[i]))    # y obs time
  tx[[i]]=as.matrix(runif(nx[i]))   # x obs time
  t.temp=rbind(tx[[i]],ty[[i]])
  n.temp=nx[i]+ny[i]
  corr=exp(-abs(rep(1,n.temp)%*%t(t.temp)-t.temp%*%t(rep(1,n.temp))))
  corr.e=2^(-abs(rep(1,n.temp)%*%t(t.temp)-t.temp%*%t(rep(1,n.temp))))
  #corr.x[[i]]=exp(-abs(rep(1,nx[i])%*%t(tx[[i]])-tx[[i]]%*%t(rep(1,nx[i]))))
  x.temp=mvrnorm(1,rep(0,n.temp),corr)
  x[[i]]=as.matrix(x.temp[1:nx[i]])
  if(Func=='sin'){a=sin(t.temp*2*pi)}
  if(Func=='sqrt'){a=sqrt(t.temp)}
  if(Func=='lin'){a=.4*t.temp+.5}
  b=-2         # true value of beta
  y.temp=a+x.temp*b+as.matrix(mvrnorm(1,rep(0,n.temp),corr.e))
  y[[i]]=as.matrix(y.temp[-(1:nx[i])])
}

b.can=rep(0,length(h.can))
se.can=b.can

for(h.ind in 1:length(h.can)){
  
  b.ini=1      # initial value of beta estimate
  #count=0
  
  repeat{
    U.temp=U.ini(b.ini, h.can[h.ind])
    U=U.temp[1]
    U.dot=U.temp[2]
    b.temp=b.ini-U/U.dot
    #count=count+1
    if(abs(b.temp-b.ini)>0.0005) b.ini=b.temp
    else break
  }
  b.est=b.ini      # initial value of beta estimate
  #count=0
  repeat{
    U.temp=u_b(b.est,tx,ty,x,y,h.can[h.ind],ha,n, xi)
    #U.temp=U.b(b.est)
    U=U.temp[1]
    U.dot=U.temp[2]
    U.sq=U.temp[3]
    b.temp=b.est-U/U.dot
    #count=count+1
    if(abs(b.temp-b.est)>0.0005) b.est=b.temp
    else break
  }
  
  gr=sample(rep(c(1:2), length =n ))
  b1=b.ini
  repeat{
    U.temp=u_b(b1,tx[gr==1],ty[gr==1],x[gr==1],y[gr==1],h.can[h.ind],ha,sum(gr==1), xi)
    #U.temp=U.b(b.est)
    U=U.temp[1]
    U.dot=U.temp[2]
    U.sq=U.temp[3]
    b.temp=b1-U/U.dot
    #count=count+1
    if(abs(b.temp-b1)>0.0005) b1=b.temp
    else break
  }
  b2=b.ini
  repeat{
    U.temp=u_b(b2,tx[gr==2],ty[gr==2],x[gr==2],y[gr==2],h.can[h.ind],ha,sum(gr==2), xi)
    #U.temp=U.b(b.est)
    U=U.temp[1]
    U.dot=U.temp[2]
    U.sq=U.temp[3]
    b.temp=b2-U/U.dot
    #count=count+1
    if(abs(b.temp-b2)>0.0005) b2=b.temp
    else break
  }
  
    b.can[h.ind]=b.est
    #se.can[h.ind]=U.sq^.5/abs(U.dot)
    se.can[h.ind]=(b1-b2)^2/4
}

#h.reg=h.can
h.reg=h.can^2
C=lm(b.can~h.reg)$coef[2]
temp=C^2*h.can^4+se.can

h[sample]=h.can[which.min(temp)]
#h
#se.can
#C^2*h.can^4
#temp


b.est=b.ini      # initial value of beta estimate
#count=0
repeat{
  U.temp=u_b(b.est,tx,ty,x,y,h[sample],ha,n, xi)
  #U.temp=U.b(b.est)
  U=U.temp[1]
  U.dot=U.temp[2]
  U.sq=U.temp[3]
  b.temp=b.est-U/U.dot
  #count=count+1
  if(abs(b.temp-b.est)>0.0005) b.est=b.temp
  else break
}

if(is.na(b.est)==0){
  
  x.seq=seq(0,1,length.out = 101)
  a.seq=rep(0,101)
  temp.a=matrix(0,ncol=3,nrow=Iter)
  #a1.sample=aSE.sample
  a.theo=x.seq
  aSE=rep(0,3)
  ind=1
  
  if(0){
    for(i in 1:length(x.seq)){
      S_temp=S(x.seq[i],tx,ty,x,y,b.est,ha,n)
      a.seq[i]=(S_temp[3]*S_temp[4]-S_temp[2]*S_temp[5] )/(S_temp[1]*S_temp[3]-S_temp[2]*S_temp[2] )
      #a1=(S_temp[1]*S_temp[5]-S_temp[2]*S_temp[4] )/(S_temp[1]*S_temp[3]-S_temp[2]*S_temp[2] )
      #a1=(S_temp[4]-a.seq[i]*S_temp[1])/S_temp[2]
      #a1.sample[sample,i]=a1
      #if(0){
      if(i%%25==1 & i>1 & i<100){
        a1=2*pi*cos(2*pi*x.seq[i])
        V=aV(a.seq[i],a1,b.est,x.seq[i],tx,ty,x,y,ha,n)
        aSE[ind]=(S_temp[2]^2*V[2]-S_temp[3]^2*V[1])/(S_temp[2]^4-S_temp[1]^2*S_temp[3]^2)
        #a1.sample[sample,ind]=a1
        ind=ind+1
      }
      if(Func=='sin'){a.theo[i]=sin(x.seq[i]*2*pi)}
      if(Func=='sqrt'){a.theo[i]=sqrt(x.seq[i])}
      if(Func=='lin'){a.theo[i]=.4*x.seq[i]+.5}
      #}
    }
    #plot(x.seq,a.seq,type='l')
    #lines(x.seq,a.theo,lty=2)
    
  }
  
  
  
  b.sample[sample]=b.est
  a.sample[sample,]=a.seq
  se.sample[sample]=U.sq^.5/abs(U.dot)
  aSE.sample[sample,]=sqrt(aSE)
  
  if(Btsp!=1){
    for(iter in 1:Iter){
      Btsp.ind=sample(n,ceiling(n*Btsp),replace=T)
      tx.Btsp=tx[Btsp.ind]
      ty.Btsp=ty[Btsp.ind]
      x.Btsp=x[Btsp.ind]
      y.Btsp=y[Btsp.ind]
      for(i in 1:3){
        S_temp=S(x.seq[i*25+1],tx.Btsp,ty.Btsp,x.Btsp,y.Btsp,b.est,ha,ceiling(n*Btsp))
        temp.a[iter,i]=(S_temp[3]*S_temp[4]-S_temp[2]*S_temp[5] )/(S_temp[1]*S_temp[3]-S_temp[2]*S_temp[2] )
      }
    }
    aBtsp[sample,]=apply(temp.a,2,sd)
  }
  
  sample=sample+1
  if(sample %% 50 ==0) print(sample)
  if(sample>Nsample) break
}else{na.count=na.count+1}

}

proc.time()-start
