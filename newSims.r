
library(MASS)
library(glmnet)

##########
# Simulate
##########

#set.seed(413)

ELNET = 1;  # 1 = LASSO, 0 = Ridge
DO_IWLS=0;
#HETRO = T;
output <- c();
reps = 200

pb = txtProgressBar(style=3,min=0,max=reps)
for( rr in 1:reps){

#nn  = 100;
sig = 3;
pp  = 5;
alf = 0.05;
xx  = matrix( runif(nn*pp,0,5),nn,pp );
xx <- apply(xx,2,scale)
bet = runif( pp, 0, 3 );
if(HETRO){
  eps = rnorm( nn, 0, xx[,1]^2 );
} else {
  eps = rnorm( nn, 0, sig );
}
yy  = scale(xx%*%bet + eps,scale=F);

xx <- t(t(xx) - apply(xx,2,mean))

######
# fit Model
######

md.cv = cv.glmnet( xx, yy, nfolds=10, alpha=ELNET );
lmb.rd= md.cv$lambda.min;
#md.rd = glmnet( xx, yy, alpha=ELNET, lambda=lmb.rd );
  
bet.org = drop(md.cv$glmnet.fit$beta[,md.cv$index[1]]);
if(ELNET==1){
  penMat = lmb.rd*ginv(diag(abs(bet.org)))
} else {
  penMat = diag(lmb.rd,pp)
}
bet.rd= solve( 
  crossprod(xx) + penMat, 
  crossprod(xx,yy) 
);

#####
# Test fit: F
#####

ssexp = 
  ( t(bet.rd - bet) %*% 
    ( t(xx)%*%xx + penMat ) %*% 
    (bet.rd - bet)
  )[1]
yfit = xx%*%bet.rd
res  = (yy - yfit)[1:nn]
ssres= sum(res^2)

fstat  <- ( ssexp / (pp) ) / ( ssres/(nn-pp) )
#fpv <- pf( 1-alf, pp, nn-pp )
#fthrsh <- qf( alf, pp+1, nn-pp-1,lower.tail=F )*(pp+1)*ssres/(nn-pp-1)
fthrsh <- qf( alf, pp, nn-pp,lower.tail=F )*(pp)*ssres/(nn-pp)

#print(paste(fstat,fthrsh))

######
# Test fit: DB
######

rrad <- function(n){
  2*rbinom(n,1,1/2)-1
}
const = exp(1)/sqrt(pi)+exp(1)/pi 
cc = 8

BS = 200;

out.db = c();
out.pv = c();
out.iw = c();

hatmat = 
  crossprod( t(xx), 
    solve( crossprod(xx) + penMat, t(xx))
  )

for( ii in 1:BS ){
  res.wb <- res*rrad(nn);
  yy.wb  <- scale(
    crossprod(t(xx),bet.rd) + res.wb,scale=F
  )
  #md.wb  <- glmnet( xx, yy.wb, alpha=ELNET, lambda=lmb.rd );
  #bet.wb <- md.wb$beta;
  bet.wb= solve( 
    crossprod(xx) + penMat, 
    crossprod(xx,yy.wb) 
  );
  
  # yfit.wb = xx%*%bet.wb;
  # wrong # res.new = yy - yfit.wb;
  res2.wb = res.wb*rrad(nn);
  yy2.wb  = scale(
    crossprod(t(xx),bet.wb) + res2.wb,scale=F
  );
  #md2.wb  = glmnet( xx, yy2.wb, alpha=ELNET, lambda=lmb.rd );
  #bet2.wb = md2.wb$beta;
  bet2.wb= solve( 
    crossprod(xx) + penMat, 
    crossprod(xx,yy2.wb) 
  );

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

  if(DO_IWLS){
    iter = 1;
    iter.max = 100;
    crit = 1;
    bet.iwls = rep(0,pp);
    lmb.iwls = rep(1,pp);
    fit.iwls = rep(0,pp);
    WW   = diag(1,nn);
    while( (crit>0.0001) & (iter<iter.max) ){
      iter = iter+1;
      bet.iwls.old <- bet.iwls;
      zz <- crossprod( t(xx),bet.iwls ) + ( yy - fit.iwls );
      bet.iwls = drop(solve(
        crossprod(xx) + lmb.rd*lmb.iwls,
        crossprod(xx,zz)
      ));
      fit.iwls <- crossprod(t(xx),bet.iwls);
      lmb.iwls <- ginv( diag( abs( bet.iwls ) ) );
      crit = max( abs(bet.iwls-bet.iwls.old) );
    }
    penMat.iwls = lmb.rd*lmb.iwls;
    hat.iwls = 
      crossprod( t(xx), 
        solve( 
          crossprod(xx) + penMat.iwls, 
          t(xx)
        )
      )
    res.iwls = drop(yy - crossprod(t(xx),bet.iwls));
    res.iwls.wb <- res.iwls*rrad(nn);
    yy.iwls.wb  <- scale(
      crossprod(t(xx),bet.iwls) + res.iwls.wb,scale=F
    )
    bet.iwls.wb= solve( 
      crossprod(xx) + penMat.iwls, 
      crossprod(xx,yy.iwls.wb) 
    );
    
    res2.iwls.wb = res.iwls.wb*rrad(nn);
    yy2.iwls.wb  = scale(
      crossprod(t(xx),bet.iwls.wb) + res2.iwls.wb,scale=F
    );
    bet2.iwls.wb= solve( 
      crossprod(xx) + penMat.iwls, 
      crossprod(xx,yy2.iwls.wb) 
    );
 
    new.iwls.t  = 
      crossprod( 
        bet2.iwls.wb-bet.iwls.wb,  
        crossprod( 
          crossprod(xx) + penMat.iwls,
          bet2.iwls.wb-bet.iwls.wb 
        )
      )[1];
    out.iw = c( out.iw, 
      exp( -new.iwls.t/( 4*cc*sqrt( 
        sum(t(res.iwls^2*(hat.iwls^2))*res.iwls^2 )
      ) )  )
    )  
  }
}

bsthrsh = quantile( out.db, 1-alf )


######
# Test fit: ADB
######

BS = 10;

out.pv <- out.pv[1:BS];
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

if(DO_IWLS){
  out.iw <- out.iw[1:BS];
  mu  = mean(out.iw);
  vr  = var(out.iw);
  th1 = mu^2*( 1-mu )/vr - mu
  th2 = (mu*( 1-mu )/vr - 1)*(1 - mu)
  if( th1  < th2 ) # right skewed 
    astar = qbeta(alf/const,th1,th2,lower.tail=F);
  if( th1 > th2 ) # left skewed
    astar = qbeta(alf/const,th1,th2,lower.tail=T);
  cn.iwls = -log(astar)*4*cc*sqrt( 
    sum(t(res.iwls^2*(hat.iwls^2))*res.iwls^2 )
  );
} else {
  cn.iwls=0
}

output = rbind(
   output,
   c(ssexp,fthrsh,bsthrsh,cnthrsh,cn.iwls)
)

setTxtProgressBar(pb,rr)

}

close(pb)

getCoverage <- function( out, IWLS=T ){
  nc = ncol(out)
  #if(IWLS){
  #  cv = out[,1] < out[,2:5];
  #} else {
    cv = out[,1] < out[,2:4];
  #}
  return( 
    apply(cv,2,mean)
  )
}


