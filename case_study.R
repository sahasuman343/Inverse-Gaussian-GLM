#Load packages
if(!require(SuppDists)) install.packages("SuppDists")
#library(SuppDists)
source("main.R")
set.seed(123)

#simulate Data
n=20
x1=rnorm(n,29,8.39)#rbinom(n,1,0.65)
x2=rgamma(n,3,4.78)
u=1/(16.56+4.56*x1+7.89*x2)
y=rinvGauss(n,u,1)
#X=as.matrix(cbind(rep(1,n),x1,x2))
#print(u)
#print(y)

#estimate parameters
#beta=IWRLS(y,X,100)
#print("Estimated betas:")
#print(beta)

#True value vs estimated value

#print("True mean vs estimated mean of the distribution of the response variable")

#print("SSE")
#e=u-update_mu(beta)
#print(sum(e^2))

x=cbind(x1,x2)
X=(cbind(rep(1,n),x))
#IWRLS(y,X,100)

model=GI_fit(y,x)
model$beta
confint(model)
#summary(fit)
plot(res(fit=model))
sum(res(fit=model)^2)

dev(model)


