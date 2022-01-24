


# Trying out Maximum Likelihood

library(tidyverse)
library(faraway)
library(foreign)
library(lme4)
library(nlme)


library(rstanarm)
library(arm)
library(bayesplot)
library(janitor)
library(metR)
library(rethinking)
library(modelsummary)

# alcohol cross sectional
d <- read.spss("data/_lecture alcoholcrosssectional.sav", 
               to.data.frame = TRUE) %>% 
              mutate(alcohol=as.numeric(alcuse>0))



ggplot(d, aes(x=age, y=alcohol, color=factor(coa))) + 
  geom_jitter(height=.1) + 
  geom_smooth(method="glm", method.args=list(family=binomial)) + 
  scale_y_continuous(breaks=c(0,1))

m1 <- glm(alcohol~age, data=d, family="binomial") 
m2 <- glm(alcohol~age + coa, data=d, family="binomial")
m3 <- glm(alcohol~age + coa + coa*age, data=d, family="binomial")
modelsummary(list(m1,m2,m3), stars=TRUE)

# Create likelihood heat map

set.seed(111)
n <- 100
x <- rnorm(n, mean=3, sd=5)
mean(x)
sd(x)

# likelihood: given the parameters of a normal distribution, 
# sigma and mu, what is the probability of observing 
# a certain outcome, x
# p(x|mu,sigma) = dnorm(x, mean=mu, sd=sigma )
lik <- function(x, mu, sigma){
  return(sum(log(dnorm(x, mean=mu, sd=sigma))))  
}

df <- expand.grid(mu = seq(2,4, b=.1), 
                  sigma = seq(4,6, b=.1))
df$loglik <- rep(NA, nrow(df))
for(i in 1:nrow(df)){
  mu    <- df$mu[i]
  sigma <- df$sigma[i]
  df$loglik[i] <- lik(x,mu,sigma)
}

# subtract the log likelihood from its maximum value
qplot(-log(max(df$loglik) - df$loglik + 1))
df$n.loglik <- round(-log(max(df$loglik)-df$loglik+1), digits=2)

df %>%  
  ggplot(aes(x=mu, y=sigma, 
             fill=n.loglik)) +
  geom_tile()


df %>%  
  ggplot(aes(x=mu, y=sigma, 
             z=n.loglik)) +
  geom_contour() + 
  geom_text_contour() 

b <- c(seq(-310,-308,by=.2))
df %>%  
  ggplot(aes(x=mu, y=sigma, 
             z=loglik)) +
  geom_contour(breaks=b) + 
  geom_text_contour(breaks=b) 


# geom_vline(xintercept=coef(glm(y~x, 
  #                                family="binomial"))[1]) + 
  # geom_hline(yintercept=coef(glm(y~x, 
  #                                family="binomial"))[2]) 






y.link <- 2*x+.5
y.link.scaled <- scale(y.link)
qplot(y.link)
qplot(y.link.scaled)

p.y <- ilogit(y.link)
p.y.scaled <- ilogit(y.link.scaled)

qplot(p.y)
qplot(p.y.scaled)

y <- rbinom(n,1,p.y)
y <- rbinom(n,1,p.y.scaled)

qplot(x,y) +
  geom_smooth(method="glm", 
              method.args=list(family="binomial"))


lik <- function(y, b, X){
# p(x|mu,sigma) = dbinom(y, size=1,p=ilogit(X %*% b) )
  return(sum(log(dbinom(y, size=1,p=ilogit(X %*% b) ))))  
}

b0 <- coef(glm(y~x, family="binomial"))[1]
b1 <- coef(glm(y~x, family="binomial"))[2]
max.lik <- lik(y, c(b0,b1), cbind(1,x))

b0.seq <- seq(b0-0.5, b0+0.5, length.out=50)
b1.seq <- seq(b1-0.5, b1+0.5, length.out=50)


df <- expand.grid(b0=b0.seq, b1=b1.seq)

df$lik <- rep(NA, dim(df)[1])
for(i in seq(dim(df)[1])){
  df$lik[i] <- lik(y, as.numeric(df[i,1:2]), cbind(1,x))
}

max.lik <- lik(y, c(b0,b1), cbind(1,x))
breaks <- round(seq(max.lik-.5,max.lik+.5,b=.05),digits=2)

df %>%  
  ggplot(aes(x=b0, y=b1, z=lik)) +
  geom_contour(breaks=breaks) + 
  geom_text_contour(breaks=breaks) + 
  geom_vline(xintercept=coef(glm(y~x, 
                                 family="binomial"))[1]) + 
  geom_hline(yintercept=coef(glm(y~x, 
                                 family="binomial"))[2]) 



#  Day 2
#  Dupuytren case study

d <- read.spss("data/practical2_dupuytren.sav", 
               to.data.frame = TRUE) 

ml.mdl <- glm(cbind(Count, TOTAL-Count)~Age, d, family=binomial)
summary(ml.mdl)

# Likelihood Ration Test (LRT) is a Chi-square test, 
# Comparing our model (residuals) to the null model
# In this case the deviance of our model is 982.23, 
# The deviance of the null is 5823.45
# 5823.45 - 982.23 = 4841.22 on 1 df (2df-1df=1df)
# Significance means that there is a significant diff
# B/w our model and the null
# Howeveer, our model's Residual deviance of 982.23 on 65
# DoF is significantly different from the saturated model
# Which suggests a bad fit

# dispersion
chisq.pearson <- sum(residuals(ml.mdl, type="pearson")^2)
disp.pearson <- chisq.pearson/ml.mdl$df.residual
disp.deviance <- ml.mdl$deviance/ml.mdl$df.residual

summary(ml.mdl, dispersion=disp.pearson)
summary(ml.mdl, dispersion=disp.deviance)


library(rethinking) 
data(UCBadmit)
d <-UCBadmit
d$gid <-ifelse(d$applicant.gender=="male",1L,2L)
d$gid_1 <- d$gid-1
dat <-list(A=d$admit,N=d$applications,gid=d$gid,gid_1=d$gid_1)
m12.1a <-ulam(
  alist(
    A ~dbinom(N,pbar),
    logit(pbar) <-b0+b1*gid_1,
    b0 ~dnorm(0,0.5),
    b1 ~dnorm(0,0.5)
  ), data=dat,chains=4, iter=2000)
precis( m12.1a,depth=2)




m12.1b <-ulam(
  alist(
    A ~dbetabinom(N,pbar,theta),
    logit(pbar) <-a[gid],
    a[gid] ~dnorm(0,1.5),
    # transpars> theta<<-phi+2.0,
    theta ~dexp(1)
  ), data=dat,chains=4, iter=2000)

precis( m12.1b,depth=2)

post <-extract.samples(m12.1)
post$da <-post$a[,1]-post$a[,2]
precis( post,depth=2)

mdl.admit <- glm(cbind(admit,reject)~gid_1, d, family=binomial)
summary(mdl.admit)
pchisq(783.61, df=10, lower.tail = FALSE)

# Add the department to the model and we can see 
# no significant in the chi-sq likelihood ratio test 
# relative to the saturated model 
mdl.admit.dept <- glm(cbind(admit,reject)~gid+dept, d, family=binomial)
summary(mdl.admit.dept)
pchisq(10.204,df=5, lower.tail = FALSE)


# check the dispersion level - pearson and deviance
chisq.pearson <- 
  sum(residuals(mdl.admit, type="pearson")^2)/mdl.admit$df.residual
chisq.deviance <- 
  mdl.admit$deviance/mdl.admit$df.residual

summary(mdl.admit, dispersion = chisq.pearson)
summary(mdl.admit, dispersion = chisq.deviance)

qmdl.admit <- 
  glm(cbind(admit,reject)~gid, d, family=quasibinomial)
summary(qmdl.admit)



#
#
# 

library(rethinking) 
data(reedfrogs)
d <-reedfrogs
str(d)
d$tank <-1:nrow(d)
rf.mdl <- glm(cbind(surv, density-surv) ~ factor(tank),  
              data=d, family=binomial)
summary(rf.mdl)

# make the tank cluster variable
d$tank <-1:nrow(d)
dat <-list(
  S =d$surv,
  N =d$density,
  tank =d$tank)

# approximate posterior
m13.1 <-ulam(
  alist(
    S ~dbinom(N,p),
    logit(p) <-a[tank],
    a[tank] ~dnorm(0,1.5)
  ), data=dat,chains=4,log_lik=TRUE)
precis(m13.1, depth=2)



m13.2 <-ulam(
  alist(
    S ~dbinom(N,p),
    logit(p) <-a[tank],
    a[tank] ~dnorm(abar,sigma), 
    abar ~ dnorm(0,1.5), 
    sigma ~ dexp(1)
  ), data=dat,chains=4,log_lik=TRUE)
precis(m13.2, depth=2)
compare( m13.1,m13.2)
post <-extract.samples(m13.2)

#  Day 3
#  _lecture Left-Submandibular-Glands
# Estimating fixed effects AND TWO standard deviation: 
# On sigma associated with the random effect and the other
# associated with the residual.
glands <- read.spss("data/_lecture mean_left_glands.sav", 
                    to.data.frame = TRUE) 

glands$Subject <- factor(glands$Subject)
ggplot(glands, aes(x=Oncologist, y=meanvolume, color=Subject)) + 
  geom_point() + geom_line()
  

glands$fct.onc <- factor(glands$Oncologist) %>% 
  fct_relevel("5")
glands$fct.sub <- factor(glands$Subject)


lm1 <- lm(meanvolume ~ fct.sub + fct.onc, data=glands)
summary(lm1)

glands  %>%  group_by(Subject) %>% 
  summarize(mn.s=mean(meanvolume), sd.s=sd(meanvolume)) %>% 
  summarise(mn=mean(mn.s), sd=sd(mn.s))
sd(glands$meanvolume)

library(lmerTest)
mod.lmer <- lme4::lmer(meanvolume ~ 1+ fct.onc1 + (1|fct.sub), 
                           data=glands, REML=TRUE)
mod.lmer <- lmerTest::lmer(meanvolume ~ 1 + fct.onc + 
                             (1|fct.sub), 
                 data=glands, REML=TRUE)

summary(mod.lmer)

anova(mod.lmer)
confint(mod.lmer)

dat <-list(
  V =glands$meanvolume,
  S =glands$Subject,
  O =glands$Oncologist)

m13.0 <-ulam(
  alist(
    V ~ dnorm(mu,sigma),
    mu <-bs[S]+ bo[O],
    bo[O] ~ dnorm(0, 1000),
    bs[S] ~dnorm(abar,sigma), 
    abar ~ dnorm(10000,1000), 
    sigma ~ dnorm(0,1000)
  ), data=dat,chains=4,log_lik=TRUE)
precis(m13.0, depth=2)








mod.glands <- lmer(meanvolume ~ 1+ fct.onc + (1|fct.sub),  
                           data=glands, REML=T)

summary(mod.glands)
2882188/(2882188+1117886  )
var(glands$meanvolume)
4429701
1117886 + 2882188
var(c(10007.5, 1160.0, -899.2, 1808.3, 1215.0))
mean(c(9140.7151, 12187.95076))

mod.lme <- lme(Volume ~ 1+ fct.onc, 
                random= ~1|fct.sub, 
                data=glands, 
                method="REML")
summary(mod.lme)




data(pulp, package="faraway")
d <- pulp %>% mutate(operator = )

# op <- options(contrasts=c("contr.sum", "contr.poly"))
lm(bright~operator, d) %>% summary()
lmod <- aov(bright ~ operator, pulp)
summary(lmod)

coef(lmod)
mod.lmer <- lmerTest::lmer(bright ~ 1 + (1|operator), 
                           data=d, REML=TRUE)
summary(mod.lmer)
d$op <- as.numeric(factor(d$operator))
  
mdl.q <- quap(
  alist(
    bright ~ dnorm(mu, sigma),
    mu  <-  mu0 + mu.op[op],
    mu0 ~ dnorm(60, 10),
    mu.op[op]  ~ dnorm(0, 1), 
    sigma ~ dunif(0,50)
  ), data=d)
precis(mdl.q, depth=2)

DF <- as.data.frame(UCBAdmissions)
## Now 'DF' is a data frame with a grid of the factors and the counts
## in variable 'Freq'.
DF
## Nice for taking margins ...
xtabs( ~ Gender + Admit, DF)


df <- data.frame(x1=rnorm(100), x2=rnorm(100)) %>% mutate(y=.3*x1+.7*x2+rnorm(100))
mdl <- lm(y~., data=df) 
summary(mdl)
cor(predict(mdl), df$y)
