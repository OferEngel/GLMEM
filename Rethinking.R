
library(tidyverse)
library(rstanarm)
library(arm)
library(bayesplot)
library(janitor)
library(metR)
library(rethinking)
library(modelsummary)



y <-c(-1,1)
set.seed(11)
m9.2 <-quap(
  alist(
    y ~ dnorm(mu,sigma),
    mu <- alpha, 
    alpha ~ dnorm(0,1000), 
    sigma ~ dexp(.01)), data = data.frame(y=y))
precis(m9.2)

m9.2 <-ulam(
  alist(
    y ~ dnorm(mu,sigma),
    mu <- alpha, 
    alpha ~ dnorm(1,10), 
    sigma ~ dexp(1)), data = list(y=y), chains=3, cores=3)

traceplot(m9.2)
precis(m9.2)
pr <- extract.prior(m9.2)
ggplot()  + geom_density(aes(x=pr$alpha))

po <- extract.samples(m9.2)
ggplot()  + geom_density(aes(x=po$alpha))





# Chapter 11
# 
# 



data("chimpanzees")
d <- chimpanzees
glm(pulled_left~prosoc_left+ condition + prosoc_left*condition, 
    data=d, family=binomial) %>% 
  summary()
d$treatment <-1+d$prosoc_left+2*d$condition
xtabs( ~treatment+prosoc_left+condition,d)
xtabs( ~pulled_left+prosoc_left+condition,d)



m11.1  <-quap(
  alist(
    pulled_left ~dbinom(1,p),
    logit(p) <- a ,
    a ~dnorm(0,.5)
  ) ,data=d)
set.seed(1999) 11.5
prior <-extract.prior(m11.1,n=1e4)
qplot(prior$a)
qplot(invlogit(prior$a))



m11.2  <-quap(
  alist(
    pulled_left ~dbinom(1,p),
    logit(p) <- a + b[treatment],
    a ~ dnorm(0,1.5),
    b[treatment] ~ dnorm(0,.5)
  ) ,data=d)
set.seed(1999) 
prior <-extract.prior(m11.2,n=1e4)
p <-sapply(1:4,function(k)inv_logit(prior$a+prior$b[,k]))
qplot(abs(p[,1]-p[,2]))
mean( abs(p[,1]-p[,2]))
qplot(prior$a)
qplot(invlogit(prior$a))

d$treatment

dat_list <-list(
  pulled_left =d$pulled_left,
  actor =d$actor,
  treatment =as.integer(d$treatment))


m11.4 <-ulam(
  alist(
    pulled_left ~dbinom(1,p),
    logit(p) <-a[actor]+b[treatment],
    a[actor] ~dnorm(0,1.5),
    b[treatment] ~dnorm(0,0.5)
  ) ,data=dat_list,chains=4,log_lik=TRUE, cores=4)
precis( m11.4,depth=2)

post <- extract.samples(m11.4)
p_left <-inv_logit(post$a)
plot( precis(as.data.frame(p_left)),xlim=c(0,1))


labs  <-c("R/N","L/N","R/P","L/P")
plot( precis(m11.4,depth=2,pars="b"),labels=labs)




data(UCBadmit)
d <-UCBadmit
d %>% mutate(pct=admit/applications) %>% 
  group_by(applicant.gender) %>% 
  summarize(admit=sum(admit), applications=sum(applications)) %>% 
  mutate(pct=admit/applications)


d %>% mutate(pct=admit/applications) %>% 
  group_by(dept, applicant.gender) %>% 
  summarize(admit=sum(admit), applications=sum(applications)) %>% 
  mutate(pct=admit/applications)

data.frame(g1=d$applicant.gender, 
           g2=as.numeric(d$applicant.gender))
d %>% mutate(gid=ifelse(applicant.gender=="male",1,2), 
             did=as.numeric(dept)) %>% 
  dplyr::select(admit, applications, gid, did) -> 
  d2


dat_list <- list(
  admit=d2$admit, 
  applications = d2$applications, 
  gid = d2$gid, 
  did = d2$did)

m11.7 <- ulam(
  alist(
    admit ~dbinom(applications,p), 
    logit(p) <-a[gid], 
    a[gid] ~ dnorm(0,1.5)
  ), data=dat_list, chains = 4, cores=4)


precis( m11.7,depth=2)

post <- extract.samples(m11.7)
diff.a <- post$a[,1] - post$a[,2]
diff.p <- invlogit(post$a[,1]) - invlogit(post$a[,2])
precis( list(diff.a=diff.a,diff.p=diff.p), hist=FALSE)




m11.8 <- ulam(
  alist(
    admit ~dbinom(applications,p), 
    logit(p) <-a[gid] + b[did], 
    a[gid] ~ dnorm(0,1.5),
    b[did] ~ dnorm(0,1.5)
  ), data=dat_list, chains = 4, cores=4)


precis( m11.8,depth=2)
post <- extract.samples(m11.8)
diff.p <- invlogit(post$a[,1]) - invlogit(post$a[,2])
diff.a <- (post$a[,1]) - (post$a[,2])
precis( list(diff.a=diff.a,diff.p=diff.p), hist=FALSE)
postcheck(m11.8)




# Chapter 12
# 

data(UCBadmit)
d <-UCBadmit
d %>% mutate(gid = ifelse(applicant.gender=="male",1,2)) %>% 
  dplyr::select(applications, gid, admit) -> d1

dat <- list(N=d1$applications, A=d1$admit, gid=d1$gid)

m12.1 <- ulam(
  alist(
    A ~ dbetabinom(N, pbar, theta), 
    logit(pbar) <- b[gid], 
    b[gid] ~ dnorm(0, 1.5), 
    transpars> theta<<-phi+2.0,
    phi ~ dexp(1)
  ), data=dat, chains=4, cores=4, iter=2000)

post <-extract.samples(m12.1)
post$da <-post$b[,1]-post$b[,2]
precis( post,depth=2, hist=FALSE)





library(rethinking) 
data(UCBadmit)
d <-UCBadmit
d$gid <-ifelse(d$applicant.gender=="male",1L,2L)
dat <-list(A=d$admit,N=d$applications,gid=d$gid)
m12.1 <-ulam(
  alist(
    A ~dbetabinom(N,pbar,theta),
    logit(pbar) <-a[gid],
    a[gid] ~dnorm(0,1.5),
    transpars> theta<<-phi+2.0,
    phi ~dexp(1)
  ), data=dat,chains=4)
precis(m12.1, depth=2)
post <-extract.samples(m12.1)
post$da <-post$a[,1]-post$a[,2]
precis( post,depth=2, hist=FALSE)
pairs(m12.1)


library(rethinking) 
data(UCBadmit)
d <-UCBadmit
d$gid <-ifelse(d$applicant.gender=="male",1L,2L)
dat <-list(A=d$admit,N=d$applications,gid=d$gid)
m12.1 <-ulam(
  alist(
    A ~dbetabinom(N,pbar,theta),
    logit(pbar) <-a[gid],
    a[gid] ~dnorm(0,1.5),
    transpars> theta<<-phi+2.0,
    phi ~dexp(1)
  ), data=dat,chains=4)

post <-extract.samples(m12.1)
post$da <-post$a[,1]-post$a[,2]
precis( post,depth=2, hist=FALSE)
precis(m12.1, depth=2)

data(orings, package="faraway")
o <- orings
o$N <- 6
glm.mdl  <- glm(cbind(damage,N-damage)~temp, data=o, family=binomial) 
glm.mdl %>% summary()

(pearson.chisq <- sum((residuals(glm.mdl))^2))

pchisq(38.898, df=22, lower.tail = FALSE)
pchisq(16.912, df=21, lower.tail = FALSE)  # Close to thingy

olist <- list(
  N=o$N, 
  T=o$temp, 
  D=o$damage)

omdl <- ulam(alist(
  D ~ dbinom(N, pbar), 
  logit(pbar) <- b0+b1*T, 
  b0 ~ dnorm(0,7), 
  b1 ~ dnorm(0, 1.5)
), data=olist, chains=4, cores=4)
precis(omdl)



library(rethinking)
data(Kline)
d <-Kline
d$P <-standardize(log(d$population))
d$contact_id <-ifelse(d$contact=="high",2L,1L)
dat2 <-list(
  T =d$total_tools,
  P =d$population,
  cid =d$contact_id)
m12.2 <-ulam(
  alist(
    T ~dgampois(lambda,phi),
    lambda <-exp(a[cid])*P^b[cid]/g,
    a[cid] ~dnorm(1,1),
    b[cid] ~dexp(1),
    g ~dexp(1),
    phi ~dexp(1)
  ), data=dat2,chains=4,log_lik=TRUE)

precis(m12.2, depth=2)



df <- data.frame(x1=rnorm(100), x2=rnorm(100)) 
df$y <- .3*df$x1+.6*df$x2+ rnorm(100)
mdl <- lm(y~x1+x2, data=df)
summary(mdl)
cor(predict(mdl), df$y)^2

data(mammalsleep, package="faraway")
mammalsleep$pdr <- with(mammalsleep, dream/sleep)
data_m <- list(
  P = mammalsleep$pdr, 
  LB = log(mammalsleep$body), 
  LL = log(mammalsleep$lifespan)
)


mdl <- ulam(
  
)
