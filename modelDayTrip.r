
library(glmnet)
library(car)

dat = read.csv( "./allBoarderData.csv" )

# write.csv(file="allBoarderData.csv",dat5,row.names=F)

dat <- subset(  
  dat, (prov!="Canada")&(prov!="Nova Scotia")
)

dat <- transform( 
  dat, date=as.Date(date), 
  prov=as.factor(prov),
  month=as.factor(month),
  weekday=as.factor(weekday),
  covid=as.factor(covid)
)

#contrasts( dat$prov  ) <- contr.sum
contrasts( dat$month ) <- contr.sum
contrasts( dat$covid ) <- contr.treatment(n=3,base=3)

md = lm( 
  dayTrip~prov+oil+xchange+month+weekday+unemploy+covid,
  data=dat
)

md.a = lm( 
  log1p(dayTrip)~prov+oil+xchange+month+weekday+unemploy,
  data=dat,subset=(covid=='prepan')
)
md.b = lm( 
  log1p(dayTrip)~prov+oil+xchange+month+weekday+unemploy,
  data=dat,subset=(covid=='lockdown')
)
md.c = lm( 
  log1p(dayTrip)~prov+oil+xchange+month+weekday+unemploy,
  data=dat,subset=(covid=='postpan')
)

###############

md.log = lm( 
  log1p(dayTrip) ~ covid*prov + oil + covid*xchange +
    month + covid*weekday + covid*unemploy,
  data=dat
)

mm = model.matrix(md.log)

md.net = glmnet(
  x=mm, y=log1p(dat$dayTrip), family='gaussian'
)

# 10-fold cv for lmb
md.net.cv = cv.glmnet(
  x=mm, y=log1p(dat$dayTrip), family='gaussian'
)

cbind( coef(md.net.cv), c(0,coef(md.log)) )

pr  = predict(md.net.cv,newx=mm,s="lambda.min")
rr  = log1p(dat$dayTrip)- pr
lmb = md.net.cv$lambda.min;  # lmb minimizes MSE
plot(pr,rr)

yy.org  = log1p(dat$dayTrip)
md.base = glmnet( 
  x=mm, y=yy.org, family='gaussian', 
  lambda=lmb
)
cf.base = coef(md.base)
pr.org  = predict(md.base,newx=mm,s=lmb)
rr.org  = yy.org - pr.org

## Single Wild Bootstrap ##

nn = length(yy.org)
BB = 10000
pb = txtProgressBar(min=0,max=BB,style=3)
res= c() 
for(ii in 1:BB){
  eps = 2*rbinom(nn,1,0.5)-1
  yy.wb  = pr.org + eps*rr.org
  md.tmp = glmnet( 
    x=mm, y=yy.wb , family='gaussian', 
    lambda=lmb
  )
  res = cbind(res,coef(md.tmp))
  setTxtProgressBar(pb,ii)
}
close(pb)

## Double Wild Bootstrap ##

nn = length(yy.org)
BB1= 100
BB2= 100
pb = txtProgressBar(min=0,max=BB1,style=3)
res= list( up=c(), lw=c() ); 
for(ii in 1:BB1){
  eps = 2*rbinom(nn,1,0.5)-1
  yy.wb  = pr.org + eps*rr.org
  md.tmp = glmnet( 
    x=mm, y=yy.wb , family='gaussian', 
    lambda=lmb
  )
  pr.org2  = predict(md.tmp,newx=mm,s=lmb)
  rr.org2  = yy.org - pr.org2
  res2=c();
  for( jj in 1:BB2 ){
    eps2 = 2*rbinom(nn,1,0.5)-1
    yy.wb2  = pr.org2 + eps2*rr.org2
    md.tmp2 = glmnet( 
      x=mm, y=yy.wb2, family='gaussian', 
      lambda=lmb
    )
    res2 = cbind(res2,coef(md.tmp2))
  }
  res$up = cbind(res$up,apply(res2,1,quantile,prob=0.975))
  res$lw = cbind(res$lw,apply(res2,1,quantile,prob=0.025))
  setTxtProgressBar(pb,ii)
}
close(pb)
fin = cbind(  
  apply(res$up,1,median),
  apply(res$lw,1,median)
)



## Fast Double Wild Bootstrap ##

nn = length(yy.org)
BB= 10000
pb = txtProgressBar(min=0,max=BB1,style=3)
res.f = c() 
for(ii in 1:BB){
  eps = 2*rbinom(nn,1,0.5)-1
  yy.wb  = pr.org + eps*rr.org
  md.tmp = glmnet( 
    x=mm, y=yy.wb , family='gaussian', 
    lambda=lmb
  )
  pr.org2  = predict(md.tmp,newx=mm,s=lmb)
  rr.org2  = yy.org - pr.org2
  res2=c();
  for( jj in 1:BB2 ){
    eps2 = 2*rbinom(nn,1,0.5)-1
    yy.wb2  = pr.org2 + eps2*rr.org2
    md.tmp2 = glmnet( 
      x=mm, y=yy.wb2, family='gaussian', 
      lambda=lmb
    )
    res2 = cbind(res2,coef(md.tmp2))
  }
  res$up = cbind(res$up,apply(res2,1,quantile,prob=0.975))
  res$lw = cbind(res$lw,apply(res2,1,quantile,prob=0.025))
  setTxtProgressBar(pb,ii)
}
close(pb)
fin = cbind(  
  apply(res$up,1,median),
  apply(res$lw,1,median)
)

####
# New Code, Sept 2024
####

library(MASS)

md.log = lm( 
  log1p(dayTrip) ~ covid*prov + oil + covid*xchange +
    month + covid*weekday + covid*unemploy,
  data=dat
)

mm = model.matrix(md.log)
xx <- apply(mm,2,scale)[,-1];
yy <- scale(log1p(dat$dayTrip),scale=F)

# 10-fold cv for lmb

ELNET=1;

md.cv = cv.glmnet(
  x=xx, y=yy, family='gaussian',
  alpha=ELNET
)

lmb.rd = 20*md.cv$lambda.min;
md.rd  = glmnet( 
  x=xx, y=yy, family='gaussian',
  alpha=ELNET, lambda=lmb.rd
)
bet.org= drop( md.rd$beta )
bet.rd = bet.org;

if(ELNET==1){
  penMat = nn*lmb.rd*ginv(diag(abs(bet.org)))
} else {
  penMat = diag(lmb.rd,pp)
}
#bet.rd= solve(
#  crossprod(xx) + penMat,
#  crossprod(xx,yy)
#);

ssexp =
  ( t(bet.rd) %*%
    ( t(xx)%*%xx + penMat ) %*%
    (bet.rd)
  )[1]
yfit = xx%*%bet.rd
res  = drop(yy - yfit)
ssres= sum(res^2)

alf= 0.05
nn = nrow(xx);
pp = ncol(xx);

fthrsh <- qf( alf, pp, nn-pp,lower.tail=F )*(pp)*ssres/(nn-pp)

plot( yfit,res, xlab="fitted values", ylab="residuals" )
abline(h=0,col='red')

rrad <- function(n){
  2*rbinom(n,1,1/2)-1
}
const = exp(1)/sqrt(pi)+exp(1)/pi
cc = 8

BS = 20;

hatmat =
  crossprod( t(xx), 
    solve( crossprod(xx) + penMat, t(xx))
  )

out.db = c();
out.pv = c();

pb = txtProgressBar( min=0,max=BS,style=3 )
for( ii in 1:BS ){
  res.wb <- res*rrad(nn);
  yy.wb  <- scale(
    crossprod(t(xx),bet.rd) + res.wb,scale=F
  )
  md.wb  = glmnet( 
    x=xx, y=yy.wb, family='gaussian',
    alpha=ELNET, lambda=lmb.rd
  )
  bet.wb = drop( md.wb$beta )
  #bet.wb= solve(
  #  crossprod(xx) + penMat,
  #  crossprod(xx,yy.wb)
  #);

  res2.wb = res.wb*rrad(nn);
  yy2.wb  = scale(
    crossprod(t(xx),bet.wb) + res2.wb,scale=F
  );
  md2.wb  = glmnet( 
    x=xx, y=yy2.wb, family='gaussian',
    alpha=ELNET, lambda=lmb.rd
  )
  bet2.wb = drop( md2.wb$beta )
  #bet2.wb= solve(
  #  crossprod(xx) + penMat,
  #  crossprod(xx,yy2.wb)
  #);

  new.t  =
    crossprod(
      bet2.wb-bet.wb,
      crossprod( 
        crossprod(xx) + penMat,
        bet2.wb-bet.wb
      )
    )[1];
  out.db = c( out.db, new.t )
  out.pv = c( out.pv,
    exp( -new.t/( 4*cc*sqrt(
      sum(t(res^2*(hatmat^2))*res^2 )
    ) )  )
  )
  setTxtProgressBar(pb,ii)
}
close(pb)

bsthrsh = quantile( out.db, 1-alf )

ABS = 10;

out.pv <- out.pv[1:ABS];
mu  = mean(out.pv);
vr  = var(out.pv);
th1 = mu^2*( 1-mu )/vr - mu
th2 = (mu*( 1-mu )/vr - 1)*(1 - mu)
if( th1  < th2 ) # right skewed 
  astar = qbeta(alf/const,th1,th2,lower.tail=F);
if( th1 > th2 ) # left skewed
  astar = qbeta(alf/const,th1,th2,lower.tail=T);
#astar = qbeta( alf, th1, th2 )
cnthrsh = -log(astar)*4*cc*sqrt(
  sum(t(res^2*(hatmat^2))*res^2 )
  #sum(hatmat^2*res^2 )
);

c(ssexp,fthrsh,bsthrsh,cnthrsh)


#
#  Individual Contrasts
#

indx=c(24,44);

c.xch   <- rep(0,pp);
c.xch[11] <- 1; 
c.une   <- rep(0,pp);
c.xch[24] <- 1; 

c.vec = rep(0,pp);
c.vec[indx] <- 1;

t.thrsh <- qt( alf, nn-pp,lower.tail=F )*sqrt(
  t(c.vec)%*%solve(crossprod(xx))%*%c.vec*ssres/(nn-pp)
)

BS = 100

out.db2 = c();
out.pv2 = c();

pb = txtProgressBar( min=0,max=BS,style=3 )
for( ii in 1:BS ){
  res.wb <- res*rrad(nn);
  yy.wb  <- scale(
    crossprod(t(xx),bet.rd) + res.wb,scale=F
  )
  md.wb  = glmnet( 
    x=xx, y=yy.wb, family='gaussian',
    alpha=ELNET, lambda=lmb.rd
  )
  bet.wb = drop( md.wb$beta )

  res2.wb = res.wb*rrad(nn);
  yy2.wb  = scale(
    crossprod(t(xx),bet.wb) + res2.wb,scale=F
  );
  md2.wb  = glmnet( 
    x=xx, y=yy2.wb, family='gaussian',
    alpha=ELNET, lambda=lmb.rd
  )
  bet2.wb = drop( md2.wb$beta )

  new.t.xch  =
    crossprod(
      c.vec,
      bet2.wb-bet.wb
    )[1];
  vvec = crossprod(
    c.vec,
    solve(
      crossprod(xx) + penMat,
      t(xx)
    )
  )
  out.db2 = c( out.db2, new.t.xch )
  out.pv2 = c( out.pv2,
    exp( -new.t.xch/( 4*cc*sqrt(
      sum( (vvec*res)^2 ) 
    ) )  )
  )
  setTxtProgressBar(pb,ii)
}
close(pb)

bsthrsh2 = quantile( out.db2, 1-alf )

ABS = 10;

out.pv <- out.pv2[1:ABS];
mu  = mean(out.pv2);
vr  = var(out.pv2);
th1 = mu^2*( 1-mu )/vr - mu
th2 = (mu*( 1-mu )/vr - 1)*(1 - mu)
if( th1  < th2 ) # right skewed 
  astar = qbeta(alf/const,th1,th2,lower.tail=F);
if( th1 > th2 ) # left skewed
  astar = qbeta(alf/const,th1,th2,lower.tail=T);
#astar = qbeta( alf, th1, th2 )
cnthrsh2 = -log(astar)*4*cc*sqrt(
  sum( (vvec*res)^2 ) 
  #sum(hatmat^2*res^2 )
);

c(t.thrsh,bsthrsh2,cnthrsh2)




##### plots #####

dat.on <- subset(dat,dat$prov=='Ontario')
pts <- c(1,500,1000,1500,2000,2500)

plot(
  dat.on$dayTrip+1,log='y',
  col=dat.on$covid,las=1,pch=as.numeric(dat.on$covid),
  ylab='day-trips',xlab='',xaxt='n'
)
axis(
  side=1,at=pts,
  labels=dat.on$date[pts]
)
title("Day-Trippers from US to Ontario, Canada")

plot(
  dat.on$xchange,
  col=dat.on$covid,las=1,
  ylab='Value (USD)',xlab='',xaxt='n'
)
axis(
  side=1,at=pts,
  labels=dat.on$date[pts]
)
title("Value of Canadian Dollar")

plot(
  dat.on$unemploy,log='y',
  col=dat.on$covid,las=1,
  ylab='unemployment rate',xlab='',xaxt='n'
)
axis(
  side=1,at=pts,
  labels=dat.on$date[pts]
)
title("US Unemployment Rate")




