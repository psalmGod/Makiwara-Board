library(stats4)
library(matlib)
maki = read.table("makiwaraboard.txt",header = TRUE)
head(maki)
tail(maki)
maki$WoodType = factor(maki$WoodType,labels=c(2,1,3,4))
maki$BoardType = factor(maki$BoardType,labels=c(1,2))
X1 = maki$WoodType
X2 = maki$BoardType
Y = maki$Deflection
# Building model without interactions
loglik1 = function(beta0,beta2,beta3,beta4,alpha2,mu,sigma) {
  R = Y - I(X1==2)*beta2 - I(X1==3)*beta3 - I(X1==4)*beta4 - I(X2==2)*alpha2 - beta0
  R = suppressWarnings(dnorm(R, mu, sigma))
  -sum(log(R))
}
m = mle(loglik1, start=list(beta0 = 100, beta2= 6.5,beta3=-17,beta4=-21,alpha2=-37, sigma=55,mu=0))
summary(m)
# Building model with interactions
loglik2 = function(beta0,beta2,beta3,beta4,alpha2,gam22, gam32, gam42, mu,sigma) {
  R = Y - I(X1==2)*beta2 - I(X1==3)*beta3 - I(X1==4)*beta4 - I(X2==2)*alpha2 - beta0 - (I(X1==2)*I(X2==2))*gam22 -(I(X1==3)*I(X2==2))*gam32 - (I(X1==4)*I(X2==2))*gam42
  R = suppressWarnings(dnorm(R, mu, sigma))
  -sum(log(R))
}
m_inter = mle(loglik2, start=list(beta0 = 100, beta2= 6.5,beta3=-17,beta4=-21,alpha2=-37, gam22=0, gam32=0, gam42=0, sigma=55,mu=0))
summary(m_inter)
# 2a 
pred_42 = (coef(m)["beta0"] + coef(m)["beta4"] + coef(m)["alpha2"])
pred_42
# 2b
# Testing for significance of interaction effect
K = rbind(matrix(c(0,0,0,0,0,1,0,0),nrow=1,ncol=8),matrix(c(0,0,0,0,0,0,1,0),nrow=1,ncol=8),matrix(c(0,0,0,0,0,0,0,1),nrow=1,ncol=8))
K # K matrix
# Wald Statistic
W = ((((coef(m_inter)["sigma"]**2)*(336-8)) - ((coef(m)["sigma"]**2))*(336-5))/3)/(coef(m_inter)["sigma"]**2)
W
# 2c
diff = (coef(m_inter)["beta0"] + coef(m_inter)["beta4"] + coef(m_inter)["alpha2"] + coef(m_inter)["gam42"])-coef(m_inter)["beta0"]
X = cbind(1,X1,X2)
V = (t(X)%*%X)
V= inv(V)
NewX = matrix(c(1,2,1))
s_err = (coef(m_inter)["sigma"])*sqrt(2 + t(NewX)%*%V%*%(NewX))
t_value = diff/s_err
t_value
pt(t_value,df=328,lower.tail = FALSE)# Accept Null hypothesis as p_value > 0.05

