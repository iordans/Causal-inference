# Read the data
dvt <- read.csv("../pt_2010.csv")

# outcome variable:   tugdiff
# treatment variable: trt

#----------------------------------------------------------
#
# Step1: Model the treatment and the GPS
#
#----------------------------------------------------------

# Restrict to plausible average # of visits per week, e.g. up to 4
d4 <- dvt[dvt$trt<=4, ]

# (i) Select the variables correlated with the "Treatment"
var.trt <- c("adl_sev9", "AGE", ... ,"MUSCLE_CONDITIONS")

# (ii) Fit the treatment model
fla.temp <- paste(var.trt, sep="", collapse="+")
( fla.trt <- as.formula(paste("trt ~", fla.temp, collapse="+")) )

mtrt.all2 <- lm(fla.trt, data=d4)
(s.trt <- summary(mtrt.all2))

# (iii) GPS (propensity function) as normal density
R <- exp(-(s.trt$residuals/s.trt$sigma)^2/2)/sqrt(2*pi*s.trt$sigma^2)

#------------------------------------------------------------
#
# Step2: Check balance (with t-stats) for non-controlled vars
#
#------------------------------------------------------------
vars <- c("tugdiff","trt",var.trt)

# classify the variables as continuous, binary, ordinal
vars.c <- c("AGE", ...)
vars.o <- c("Pain_freq", ...)
vars.b <- vars[ -which(vars %in% c("tugdiff", "trt", vars.c, vars.o))]

tstat <- NULL
# (i) continuous
for (nm in vars.c) {
  if (nm=="adl_sev9")   {
    m  <- lm(d4[ ,nm] ~ trt, data=d4)   # ADL severity may be 0 but is not skewed
  }
  else  {
    m  <- lm(log(d4[ ,nm]) ~ trt, data=d4)
  }
  s <- coef(summary(m))
  tstat <- rbind(tstat, data.frame(var=nm,tstat=s[2, "t value"],dtype="c") )
}

# (ii) binomial
for (nm in vars.b) {
  m  <- glm(d4[ ,nm] ~ trt, family=binomial, data=d4)    # NOTE the glm()
  s <- coef(summary(m))
  tstat <- rbind(tstat, data.frame(var=nm,tstat=s[2, "z value"],dtype="b") )
  # Here t-stat is z-stat
}

# (iii) ordinal
for (nm in vars.o) {
  m  <- polr(as.factor(d4[ ,nm]) ~ trt, data=d4)  # NOTE the ordinal regression
  s <- coef(summary(m))
  tstat <- rbind(tstat, data.frame(var=nm,tstat=s[1, "t value"],dtype="o") )
}

#------------------------------------------------------------
#
# Step3: Adjust with the PS and try the balance again
#
#------------------------------------------------------------
tstat.a <- NULL
# (i) continuous  - shall we log() then? And the trt?
for (nm in vars.c) {
  if (nm=="adl_sev9")   {
    m  <- lm(d4[ ,nm] ~ trt + prop1, data=d4)
    # ADL sev could be 0 but is not skewed
  }
  else  {
    m  <- lm(log(d4[ ,nm]) ~ trt + prop1, data=d4)
  }
  s <- coef(summary(m))
  tstat.a <- rbind(tstat.a, data.frame(var=nm,tstat=s[2, "t value"],dtype="c") )
}

# (ii) binomial: warning for M1030_2
for (nm in vars.b) {
  m  <- glm(d4[ ,nm] ~ trt + prop1, family=binomial, data=d4)    # note the glm()
  s <- coef(summary(m))
  tstat.a <- rbind(tstat.a, data.frame(var=nm,tstat=s[2, "z value"],dtype="b") )
  # Here t-stat is z-stat
}

# (iii) ordinal: error for M1340 (surgical wound)
for (nm in vars.o) {
  m  <- polr(as.factor(d4[ ,nm]) ~ trt + prop1, data=d4)    # ordinal regression
  s <- coef(summary(m))
  tstat.a <- rbind(tstat.a, data.frame(var=nm,tstat=s[1, "t value"],dtype="o") )
  # Here t-stat is t-stat again
}


# "Plot the balance"
temp <- merge(tstat, tstat.a, by="var")
tstat.m <- data.frame(temp$var, dtype=temp$dtype.x, tstat.pre=temp$tstat.x, tstat.post=temp$tstat.y)

capt1 <- qqnorm(tstat.m$tstat.pre, main="Normal Q-Q plot \n unadjusted scores")
capt2 <- qqnorm(tstat.m$tstat.post, main="Normal Q-Q plot \n adjusted scores")
ymax  <- ceiling(max(tstat.m$tstat.pre,tstat.m$tstat.post))
ymin  <- floor(min(tstat.m$tstat.pre,tstat.m$tstat.post))
plot(capt1$x, capt1$y, xlim=c(-3,3),ylim=c(ymin,ymax)
        , xlab = "Theoretical Quantiles", ylab = "Sample Quantiles"
        , main="Unadjusted and adjusted t-values \n with 2 sd. band")
points(capt2$x, capt2$y, pch=16)

#add reference ranges
abline(h=2, lty=2); abline(h=-2, lty=2)

# add a legend
temp <- legend("topleft"
               , legend = c(" ", " ")
               , text.width = strwidth("1,000,000")
               , pch=c(1,16), col = c("black","black")
               , xjust = 1, yjust = 1)
text(temp$rect$left + temp$rect$w, temp$text$y,
     c("Unadjusted", "Adjusted"), pos=2)
     

#----------------------------------#
#
#  Step 4. DRF1: the GPS approach
#
#----------------------------------#
     
#--------------------------------#
#  The broken stick fit          #
#--------------------------------#

#  "Break" points suggested by DRF2 GAM model below
rhs <- function(x,c) ifelse(x>c, x-c, 0)
lhs <- function(x,c) ifelse(x<c, c-x, 0)

# Final broken (hockey) stick model
lm.p <- lm(tugdiff ~  rhs(trt, 1.25) + lhs(trt, 2.25) + R, data=d4)
s <- summary(lm.p)

#--------------------------------#
#  Bootstrap CI                  #
#--------------------------------#

# NOTE: The linear models and the function are all dependent on the dataset
dpr   <- d4
m.trt <- mtrt.all2  # model for 'trt'
m.drf <- lm.p       # DRF model

# Define the GPS function
GPS <- function(t, m=m.trt) {
  # Takes as paremeter an LM model
  m.s <- summary(m)
  exp(-(t-m$fitted.values)^2/(2*m.s$sigma^2))/sqrt(2*pi*m.s$sigma^2)
}

# 2-stage bootstrap: residuals sampling
drf.final <- function(tt, ds=dpr, m.t=m.trt, m.d=m.drf, replic=10) {
  # compute the reference curve - should be the mean of the bootstrap
  dpred <- data.frame(trt=rep(tt, nrow(ds)), R=GPS(tt, m.t))
  avg.c <- mean(predict(m.d, dpred))

  # prepare for the bootsrap
  ds.boot<- data.frame(ds, trt.i=ds$trt, tugdiff.i=ds$tugdiff, R.i=dpred$R)
  fb.trt <- formula(m.t)
  fb.trt <- update(fb.trt, trt.i ~ .)
  fb.out <- formula(m.d)
  fb.out <- update(fb.out, tugdiff.i ~ . - R + R.i)

  # Bootstrap for 95%CI
  numrows <- length(m.t$fitted.values)
  stat <- NULL
  for (i in 1:replic) {
     # 1st draw - for the treatment model
     ind <- sample(numrows, replace=TRUE)

     # sample residuals, add to predicted values
     ds.boot$trt.i <-  m.t$residuals[ind] + m.t$fitted.values

     # bootstrapped model of "trt" (m.i) => bootstrapped R estimation
     m.i <- lm(fb.trt, data=ds.boot, y=TRUE)
     summary(m.i)
     R.i <- GPS(ds.boot$trt, m.i)         # use the original 'trt' values!
     ds.boot$R.i <- R.i                   # pertubed 'R' (diff. model coeff.)

     # 2nd draw - for the outcome model
     ind2 <- sample(numrows, replace=TRUE)

     # sampled outcome
     ds.boot$tugdiff.i <-  m.d$residuals[ind2] + m.d$fitted.values

     # "Doubly bootstrapped" DRF model
     # The GPS ('R') and the outcome ('tugdiff') are peturbed (but not 'trt')
     m.d.i <- lm(fb.out, data= ds.boot)
     summary(m.d.i)

     # Use the last model to predict the outcome (on the original 'R' - denote R.i)
     dpred.i <- data.frame(trt=rep(tt, nrow(ds)), R.i=GPS(tt, m.i))
     stat <- c(stat, mean(predict(m.d.i, dpred.i)))
  }
  return(c(avg.c=avg.c, avg=mean(stat), std=sd(stat)))
}

# Plot
nmin <- min(c(y.drf["avg", ]-1.96*y.drf["std", ], y.drf["avg.c", ]))
nmax <- max(c(y.drf["avg", ]+1.96*y.drf["std", ], y.drf["avg.c", ]))
plot( x.trt, y.drf["avg.c", ], ylim=c(nmin, nmax)) #, type="n")
lines(x.trt, y.drf["avg", ])
lines(x.trt, y.drf["avg", ]-1.96*y.drf["std", ], lty=2, col="blue")
lines(x.trt, y.drf["avg", ]+1.96*y.drf["std", ], lty=2, col="blue")


# DRF2: using Simon Wood's library (for GAM)
library(mgcv)

t.gam1 <- gam( tugdiff ~ s(trt)+R, data=d4,  family=gaussian)
summary(t.gam1)
plot(t.gam1, residuals=F, xlab="Mean # visits per week", rug=F
, ylab="TUG change"
, main="Relative contribution of trt. to TUG change\n (Smooth transformation of treatment)")