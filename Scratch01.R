


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

pois.mdl <- glm(Count ~ offset(log(TOTAL)) + Age, family=poisson, data=d)
summary(pois.mdl)
sum(residuals(pois.mdl, type="pearson")^2)/pois.mdl$df.residual
pois.mdl$deviance/pois.mdl$df.residual
mean(d$Count/d$TOTAL)
var(d$Count/d$TOTAL)


library(MASS) # for Insurance dataset
n <- 10000
x <- rpois(n,lambda=runif(n,3,5))
x <- rpois(n,lambda=4)

x <- rbinom(n, 10, p=.5)
x <- rbinom(n, 10, p=runif(n,.3,.7))

10*.5*.5
var(x)
mean(x)


# modelling the claim rate, with exposure as a weight
# use quasipoisson family to stop glm complaining about nonintegral response
glm(Claims/Holders ~ District + Group + Age,
    family=quasipoisson, data=Insurance, weights=Holders) %>% summary()
glm(Claims ~ District + Group + Age + offset(log(Holders)),
    family=poisson, data=Insurance) %>% summary()
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
library(haven)
write_sav(d, "data/practical2_UCBadmit.sav")
d <- read_sav("data/practical2_UCBadmit.sav")


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






mtcars %>%
  group_by(cyl, am) %>%
  summarise(mpg = mean(mpg), .groups = "drop") %>%
  spread(am, mpg)


d %>% mutate(rate=admit/applications) %>% 
  spread(applicant.gender, rate)


d %>% mutate(rate=admit/applications) %>% 
  pivot_wider(names_from=applicant.gender, values_from = rate)

xtabs(~dept+applicant.gender, data=d)


x <- rbinom(10000, 10, 0.40)
mean(x)
sd(x)
sqrt(10*.4*.6)

# Lab 03
# 
library(foreign)
library(tidyverse)
library(lme4)
library(sjstats)

d <- read.spss("data/_lecture alcoholpp.sav", 
               to.data.frame = TRUE) 

ggplot(d, 
       aes(y=alcuse, x=age_14, color=factor(coa))) + 
  geom_smooth(method="lm", se=FALSE, 
              aes(group=factor(id)),  
              linetype="longdash") +  
  geom_smooth(method="lm", se=FALSE, size=3) + 
  theme(legend.position = "none") + theme_minimal() 



ggplot(d %>% filter(coa==1), 
       aes(y=alcuse, x=age_14)) + 
  geom_smooth(method="lm", se=FALSE, 
              aes(group=factor(id)), size=.5) +  
  geom_smooth(method="lm", se=FALSE, color="black") + 
  theme(legend.position = "none")


# Model A
# 
# One intercept for everyone 
# One random variable intercept
# 
# i is the person (level 2), j is the time point (level 1)
# Y_ij = b_0i + e_ij, e_ij \sim N(0, sigma_r)
# b_0i = beta_00 + a_0i, a_0i \sim N(0, tau_0)
# estimate: beta_00, sigma_s and sigma_r
# ICC = tau_0^2/(tau_0^2 + sigma_r^2)

ma <- lmer(alcuse ~ 1 + ( 1 | id ), data=d, REML=TRUE)
summary(ma)
icc(ma)
performance::icc(ma)
# conclusion: 0.5731/(0.5731 + 0.5617 ) = 50.5%
# An estimated 51% of alcohol use is 
# Variation attributable to differences between subjects 



# Model B
# i is the person (level 2), j is the time point (level 1)
# Y_ij = b_0i + b_1i X_time + e_ij, where e_ij \sim N(0,sigma_r^2)
# b_0i = b_00 + a_0i, where a_0i \sim N(0,tau_0^2)
# b_1i = b_10 + a_1i, where a_1i \sim N(0,tau_1^2)
# estimate fixed effects: b_00 and b_10
# estimate random effects: tau_0, tau_1

mb <- lmer(alcuse ~ 1 + age_14 + ( 1 + age_14 | id ), 
           data=d, REML=TRUE)
summary(mb)
modelsummary(list(ma,mb), stars=TRUE)
performance::icc(mb)

# estimate fixed effects: b_00=0.65130 and b_10=0.27065    
# estimate random effects: 
# tau_0 = 0.6355, tau_1 = 0.1552, sigma_r^2 = 0.3373   


# Model C
# i is the person (level 2), j is the time point (level 1)
# Y_ij = b_0i + b_1i X_time  + b_2i X_coa + e_ij, where e_ij \sim N(0,sigma_r^2)
# b_0i = b_00 + a_0i, where a_0i \sim N(0, tau_0^2)
# b_1i = b_10 + a_1i, where a_1i \sim N(0, tau_1^2)
# b_2i = b_20 + a_2i, where a_2i \sim N(0, tau_2^2)
# estimate fixed effects: b_00, b_10 and b_20
# estimate random effects: tau_0, tau_1 and tau_2

mc <- lmer(alcuse ~ 1 + age_14 + coa + age_14*coa + ( 1 + age_14 | id ), 
           data=d, REML=TRUE)
summary(mc)


mc.1 <- lmer(alcuse ~ 1 + age_14 + coa + ( 1 + age_14 | id ), 
           data=d, REML=TRUE)
summary(mc.1)

2*logLik(mc)- 
2*logLik(mc.1)
pchisq(631.923 - 629.78, df=1, lower.tail = FALSE)
# Doesnt look significant difference with and 
# without interaction. In other words, there is no evidence
# that the coa has a different impact on different ages. 


md <- lm(alcuse ~ coa + age_14 + age_14*coa, d)
summary(md)
modelsummary(list(md, mc), stars=TRUE)

x <- rnorm(10)
y <- rnorm(10) + x
m1 <- lm(y~x) 
m2 <- glm(y~x) 

confint(m1)
confint(m2)




# Day 4
# 


d <- read.spss("data/_lecture alcoholpp.sav", 
               to.data.frame = TRUE) %>%
  mutate(wave=factor(age_14))
                                                                                
library(nlme)
mdl <- gls(alcuse ~ wave+ coa + wave*coa,
    data = d,
    correlation = corSymm(form= ~1|id),
    weights =varIdent(form= ~1| wave))
summary(mdl)




library(rethinking)
library(lme4)
data("reedfrogs")
d <- reedfrogs
my.var <- function(x) return(var(x)*(length(x)-1)/length(x))


nj <- 7
ni <- 100
j <- 1:nj
mu_j <- rnorm(nj, 100, 20)


mu_ij <- rep(mu_j,ni) 
ij <- factor(rep(j,ni))
y_ij <- rnorm(ni*nj, mu_ij, 10)




data.frame(y_ij, ij) %>% group_by(ij) %>% 
  summarise(xbar=mean(y_ij), sd=sd(y_ij)) 

data.frame(y_ij, ij) %>% group_by(ij) %>% 
  summarise(Y_.j=mean(y_ij), sd=sd(y_ij)) %>% rename(group=ij) %>% 
  kbl() %>%
  kable_paper("hover", full_width = F) 

data.frame(y_ij, ij) %>% group_by(ij) %>% 
  summarise(xbar=mean(y_ij), sd=sd(y_ij)) %>% pull(xbar)->xbar




# The variance between the groups
var(rep(xbar,ni))
sum((rep(xbar,ni)-mean(rep(xbar,ni)))^2)/(ni*nj-1)

var(xbar)
sum((xbar-mean(xbar))^2)/(nj-1)

ni*sum((xbar-mean(xbar))^2)/(ni*nj-1)

# The variance within each group
var(y_ij - rep(xbar,ni))
sum((y_ij - rep(xbar,ni))^2)/(ni*nj-1)



# The total variance
var(y_ij)
sum((y_ij - mean(y_ij))^2)/(ni*nj-1)




# sigma_s2
var(xbar)


mdl.re <- lmer(y_ij ~  1 + (1|ij))
summary(mdl.re)

y.hat <- predict(mdl.re)
var(y.hat)
table(y.hat)


a  <-3.5 # average morning wait time
b <-(-1) # average difference afternoon wait time
sigma_a <-1# std dev in intercepts
sigma_b <-0.5#std dev in slopes
rho <-(-0.7)#correlation between intercepts and slopes





d <- read.spss("data/_lecture Dupuytren.sav", 
               to.data.frame = TRUE) %>% 
  mutate(Study=factor(Study))



library(gee)
mdl <- gee(cbind(Count, TOTAL-Count)~Age, id=Study, 
           corstr = "exchangeable", 
           scale.fix=TRUE, 
           family=binomial, 
           data=d)
summary(mdl)
