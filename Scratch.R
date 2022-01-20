
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


set.seed(123)
n <- 100
x <- rnorm(n)
p.y <- ilogit(2*x+.5)
qplot(2*x+.5)
qplot(p.y)

y <- rbinom(n,1,p.y)
qplot(x,y) +
  geom_smooth(method="glm", 
              method.args=list(family="binomial"))


lik <- function(y, b, X){
# p(x|mu,sigma) = dbinom(y, size=1,p=ilogit(X %*% b) )
  return(sum(log(dbinom(y, size=1,p=ilogit(X %*% b) ))))  
}

b0 <- coef(glm(y~x, family="binomial"))[1]
b1 <- coef(glm(y~x, family="binomial"))[2]
max.lik <- lik(y, c(b0,b1), X)

b0.seq <- seq(b0-0.5, b0+0.5, length.out=50)
b1.seq <- seq(b1-0.5, b1+0.5, length.out=50)


df <- expand.grid(b0=b0.seq, b1=b1.seq)

df$lik <- rep(NA, dim(df)[1])
for(i in seq(dim(df)[1])){
  df$lik[i] <- lik(y, as.numeric(df[i,1:2]), X)
}

max.lik <- lik(y, c(b0,b1), X)
breaks <- round(seq(max.lik-.5,max.lik+.5,b=.05),digits=2)

df %>%  
  ggplot(aes(x=b0, y=b1, z=lik)) +
  geom_contour(breaks=breaks) + 
  geom_text_contour(breaks=breaks) + 
  geom_vline(xintercept=coef(glm(y~x, 
                                 family="binomial"))[1]) + 
  geom_hline(yintercept=coef(glm(y~x, 
                                 family="binomial"))[2]) 



#  Day 3
#  _lecture Left-Submandibular-Glands

glands <- read.spss("data/_lecture mean_left_glands.sav", 
                    to.data.frame = TRUE) 
glands <- read.spss("data/_lecture Left-Submandibular-Glands.sav", 
                    to.data.frame = TRUE) 

glands$meanvolume <- glands$Volume

mdl <- quap(
  alist(
    meanvolume ~ dnorm(mu.v, s.r), 
    mu.v ~  dnorm(10664, 5000),
    s.r ~  dexp(2)), data=glands)


precis(mdl, depth=2)


glands$fct.onc <- as.character(glands$Oncologist) 
glands$fct.sub <- as.character(glands$Subject)

glands$fct.onc <- factor(glands$Oncologist) %>% 
  fct_relevel("3")

glands$fct.sub <- factor(glands$Subject)%>% 
  fct_relevel("3")

lm1 <- lm(meanvolume ~ fct.sub + fct.onc, data=glands)
summary(lm1)

glands  %>%  group_by(Subject) %>% 
  summarize(mn=mean(meanvolume), sd=sd(meanvolume)) %>% 
  summarise(mn=mean(mn))

library(lmerTest)
mod.lmer <- lme4::lmer(meanvolume ~ 1+ fct.onc + (1|fct.sub), 
                           data=glands, REML=TRUE)
mod.lmer <- lmerTest::lmer(meanvolume ~ 1+ fct.onc + (1|fct.sub), 
                 data=glands, REML=TRUE)

summary(mod.lmer)
confint(mod.lmer)

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
d <- pulp
op <- options(contrasts=c("contr.sum", "contr.poly"))
lmod <- aov(bright ~ operator, pulp)
summary(lmod)
coef(lmod)
mod.lmer <- lmerTest::lmer(bright ~ 1 + (1|operator), 
                           data=d, REML=TRUE)
summary(mod.lmer)
