library(Countr)
library(fitdistrplus)
library(KScorrect)
library(goftest)
library(ggplot2)
library(car)
library(poweRlaw)
library(reshape2)
library(flexsurv)
library(data.table)
library(scales)
library(msm)

load("/home/dmorina/Insync/2102177@uab.cat/OneDrive Biz/Docència/UAB/2025-2026/Primer Semestre/Taller de Modelització/Data/dstDEF.RData")
load("/home/dmorina/Insync/2102177@uab.cat/OneDrive Biz/Docència/UAB/2025-2026/Primer Semestre/Taller de Modelització/Data/diff.timesDEF.RData")
load("/home/dmorina/Insync/2102177@uab.cat/OneDrive Biz/Docència/UAB/2025-2026/Primer Semestre/Taller de Modelització/Data/Quebec.RData")
sunspots <- read.table("/home/dmorina/Insync/2102177@uab.cat/OneDrive Biz/Docència/UAB/2025-2026/Primer Semestre/Taller de Modelització/Data/SN_d_tot_V2.0.txt",
                       header=TRUE, sep=";")
sunspots$X <- NULL
sunspots$TOTAL[sunspots$TOTAL == -1] <- NA
sunspots <- sunspots[!is.na(sunspots$TOTAL), ]
sunspots$CDATE2 <- as.Date(paste0(sunspots$YEAR, "-", sunspots$MONTH, "-", sunspots$DAY))
dst$CDATE2 <- as.Date(dst$CDATE)
sunspots <- sunspots[sunspots$CDATE %in% dst$CDATE2, ]
dst2 <- merge(dst, sunspots[, c("CDATE2", "TOTAL")], by="CDATE2")
dst2 <- dst2[order(dst2$CDATE), ]
for (i in 1:40)
{
  assign(paste0("sunspots", i), vector())
  assign(paste0("storm", i), dst2[dst2[,paste0("STORM", i)] == 1, ])
  for (j in 2:eval(parse(text=paste0("(dim(storm",i,")[1])"))))
  {
    eval(parse(text=paste0("sunspots", i,"[", j,"]",
                           "<- storm", i, "$TOTAL[j]")))
  }
}

for (i in 1:40)
{
  eval(parse(text=paste0("length(sunspots",i,") <- length(diff.times$diff1)")))
}

rm(storm1, storm2, storm3, storm4, storm5, storm6, storm7, storm8, storm9, storm10, storm11, storm12,
   storm13, storm14, storm15, storm16, storm17, storm18, storm19, storm20, storm21, storm22, storm23,
   storm24, storm25, storm26, storm27, storm28, storm29, storm30, storm31, storm32, storm33, storm34,
   storm35, storm36, storm37, storm38, storm39, storm40, i, j)
sunspots<- data.frame(sunspots1,sunspots2,sunspots3,sunspots4,sunspots5,sunspots6,
                         sunspots7,sunspots8,sunspots9,sunspots10,sunspots11,sunspots12,
                         sunspots13,sunspots14,sunspots15,sunspots16,sunspots17,sunspots18,
                         sunspots19,sunspots20,sunspots21,sunspots22,sunspots23,sunspots24,
                         sunspots25,sunspots26,sunspots27,sunspots28,sunspots29,sunspots30,
                         sunspots31,sunspots32,sunspots33,sunspots34,sunspots35,sunspots36,
                         sunspots37,sunspots38,sunspots39,sunspots40)
diff.times <- cbind(diff.times, sunspots)

bi.test <- function(data)
{
  h <- vector()
  for (i in 3:(dim(data)[1]-2))
  {
    xpre  <- difftime(data$CDATE[i], data$CDATE[i-1], units="days")
    xpost <- difftime(data$CDATE[i+1], data$CDATE[i], units="days")
    if (xpre <= xpost)
    {
      x <- xpre
      y <- difftime(data$CDATE[i-1], data$CDATE[i-2], units="days")
    }else{
      x <- xpost
      y <- difftime(data$CDATE[i+2], data$CDATE[i+1], units="days")
    }        
    h[i] <- as.numeric(x)/as.numeric((x+y/2))
  }
  pval <- ad.test(h[!is.na(h)], "punif")$p.value
  return(list(h=h, pvalue=pval)) 
}

# Multiple plot function (used to generate figure 1)
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

### Test the Poisson assumption
bi.test(dst[dst$STORM25==1, ])

### Test if inter-occurrence times are Weibull or exponential for super-storms
LcKS(na.omit(diff.times$diff25), cdf="pweibull")$p.value
LcKS(na.omit(diff.times$diff25), cdf="pexp")$p.value

### Test if Dst tails follow a power law / exponential distribution
dst_tail <- dst[dst$DST<0, ]
m_pl <- conpl$new(-1*dst_tail$DST) ### power law
est <- estimate_xmin(m_pl)
m_pl$setXmin(est)
bootstrap_p(m_pl, threads = 4)$p
m_pl <- conexp$new(-1*dst_tail$DST) ### exponential
est <- estimate_xmin(m_pl)
m_pl$setXmin(est)
bootstrap_p(m_pl, threads = 4)$p

### Anual number of storms dispersion index (Dst < -250nT)
data <- dst[dst$STORM25==1, c(1, 2, 27)]
data$YEAR <- as.numeric(substr(data$CDATE,1,4))
for (i in 2:dim(data)[1])
{
  if (as.numeric(difftime(data$CDATE[i], data$CDATE[i-1])) < 2) data$STORM25[i] <- 0
}
eval(parse(text=paste0("df.YEAR <- aggregate(data$STORM25, by=list(data$YEAR), FUN=sum)")))
while(dim(df.YEAR)[1]!=(df.YEAR$Group.1[length(df.YEAR$Group.1)]-df.YEAR$Group.1[1]+1))
{
  for (i in 2:dim(df.YEAR)[1])
  {
    if ((df.YEAR$Group.1[i] != df.YEAR$Group.1[i-1]+1))
    {
      df.YEAR <- rbind(df.YEAR, c(df.YEAR$Group.1[i-1]+1, 0))
    }
    df.YEAR <- df.YEAR[order(df.YEAR$Group.1), ]
  }
}
while(df.YEAR$Group.1[length(df.YEAR$Group.1)] != 2017)
{
  df.YEAR <- rbind(df.YEAR, c(df.YEAR$Group.1[length(df.YEAR$Group.1)]+1, 0))
}
var(df.YEAR$x)/mean(df.YEAR$x)

### Figure 1
dst_march_1989 <- dst[substr(dst$CDATE, 1, 4)=="1989" & substr(dst$CDATE, 6, 7)=="03", ]
thresholds <- data.frame(x = seq(1:length(dst_march_1989$DST)), y = c(-50, -100, -250))
p1 <- ggplot(dst_march_1989, aes(x=seq(1:length(dst_march_1989$DST)), y=DST)) + geom_point() + geom_line() + xlab("Time (hours)") + ylab("Dst") + 
  geom_line(aes( x, y, linetype = factor(y)), thresholds) + scale_linetype_discrete(name="Magnitude", labels=c("Super storm", "Intense storm", "Moderate storm")) + 
  scale_x_discrete(limits=c(1:744), breaks=c(1,200,400,600))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position=c(0.8, 0.2)) + ggtitle("a. Dst index in March 1989")

quebec <- quebec[1:2000, ]
thresholds <- data.frame(x = seq(1:length(quebec$SYM.H)), y = c(50, -50, -100, -250))
p2 <- ggplot(quebec, aes(x=seq(1:length(quebec$SYM.H)), y=SYM.H)) + geom_point() + geom_line() + xlab("Time (minutes)") + ylab("SYM-H") + 
  geom_line(aes( x, y, linetype = factor(y)), thresholds) + scale_linetype_discrete(name="Magnitude", labels=c("Super storm", "Intense storm", "Moderate storm", "50 nT")) + 
  scale_x_discrete(limits=c(1:2000), breaks=c(1,500,1000,1500,2000))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position=c(0.2, 0.2)) + ggtitle("b. SYM-H around the March 1989 event")

multiplot(p1, p2, cols=2)

### Fit a Weibull distribution to inter-occurrence times
shape1 <- vector()
scale1 <- vector()
for (i in 1:26)
{
  eval(parse(text=paste0("shape1[", i,"]<-fitdist(as.numeric(na.omit(diff.times$diff", i+14,")), 'weibull')$estimate[1]")))
  eval(parse(text=paste0("scale1[", i,"]<-fitdist(as.numeric(na.omit(diff.times$diff", i+14,")), 'weibull')$estimate[2]")))
}
shape1[11]; scale1[11] # Weibull parameters corresponding to super-storms (Dst < -250nT)

### Model for shape and scale parameters (excluding inter-occurrence times below the 48h limit)
mdata <- melt(as.data.table(diff.times), measure = patterns("diff", "sunspots"))
for (i in 1:40)
{
  eval(parse(text=paste0("mdata$threshold[mdata$variable==", i, "] <- -10*i")))
}
mdata <- mdata[!is.na(mdata$value1), ]
### remove lower thresholds (independence is not clear) and fit the model
mdata <- mdata[mdata$threshold<=-150, ]

### Maximum (sunspots over 125)
mdata_max <- mdata[mdata$value2>125, ]
reg.exp     <- survreg(Surv(value1)~threshold + value2, data=mdata_max, dist="exp")
reg.weibull <- survreg(Surv(value1)~threshold + value2, data=mdata_max, dist="weibull")
reg.gamma   <- flexsurvreg(Surv(value1)~threshold + value2, data=mdata_max, dist="gamma")
AIC(reg.exp); AIC(reg.weibull); AIC(reg.gamma) ### AIC corresponding to the Weibull model is lower, so this is the preferred model

### Minimum (sunspots below 30)
mdata_min <- mdata[mdata$value2<30, ]
reg.exp     <- survreg(Surv(value1)~threshold + value2, data=mdata_min, dist="exp")
reg.weibull <- survreg(Surv(value1)~threshold + value2, data=mdata_min, dist="weibull")
reg.gamma   <- flexsurvreg(Surv(value1)~threshold + value2, data=mdata_min, dist="gamma")
AIC(reg.exp); AIC(reg.weibull); AIC(reg.gamma) ### AIC corresponding to the Weibull model is lower, so this is the preferred model

### Overall including sunspots
reg.exp     <- survreg(Surv(value1)~threshold + value2, data=mdata, dist="exp")
reg.weibull <- survreg(Surv(value1)~threshold + value2, data=mdata, dist="weibull")
reg.gamma   <- flexsurvreg(Surv(value1)~threshold + value2, data=mdata, dist="gamma")
AIC(reg.exp); AIC(reg.weibull); AIC(reg.gamma) ### AIC corresponding to the Weibull model is lower, so this is the preferred model

### Figure 3
par(mfrow=c(2, 3))
  plot(acf(na.omit(diff.times$diff5), plot=FALSE)[1:12],   ylim=c(-0.5, 0.5), main="a. Dst < -50 nT")
  plot(acf(na.omit(diff.times$diff10), plot=FALSE)[1:12],  ylim=c(-0.5, 0.5), main="b. Dst < -100 nT")
  plot(acf(na.omit(diff.times$diff15), plot=FALSE)[1:12],  ylim=c(-0.5, 0.5), main="c. Dst < -150 nT")
  plot(acf(na.omit(diff.times$diff20), plot=FALSE)[1:12],  ylim=c(-0.5, 0.5), main="d. Dst < -200 nT")
  plot(acf(na.omit(diff.times$diff25), plot=FALSE)[1:12],  ylim=c(-0.5, 0.5), main="e. Dst < -250 nT")
  plot(acf(na.omit(diff.times$diff30), plot=FALSE)[1:12],  ylim=c(-0.5, 0.5), main="f. Dst < -300 nT")
dev.off()

### Figure 4
par(mfrow=c(1, 2))
  plot(unique(mdata$threshold), log(shape1), xlab="Threshold", ylab="log(shape parameter)", main="a. Shape parameter")
  points(unique(mdata$threshold), rep(log(1/reg.weibull$scale), 26), type="l")
  plot(unique(mdata$threshold), log(scale1), xlab="Threshold", ylab="log(scale parameter)", main="b. Scale parameter")
  points(unique(mdata$threshold), coef(reg.weibull)[1]+coef(reg.weibull)[2]*unique(mdata$threshold)+coef(reg.weibull)[3]*unique(mdata$value2), type="l")
dev.off()

### Probability of a Carrington event (Dst < -850nT) in the next decade (bootstrap)
carringtonTime <- as.numeric(difftime(Sys.Date(), as.Date("1859-09-01")))
futurTime <- as.POSIXlt(Sys.time())
futurTime$year <- futurTime$year+10
futurTime <- as.numeric(difftime(as.Date(futurTime), as.Date("1859-09-01")))
pcarr <- vector()
for (i in 1:100000)
{
  boots.df_i  <- mdata[sample(nrow(mdata), dim(mdata)[1], replace=TRUE), ]
  reg.weibull <- survreg(Surv(value1)~threshold, data=boots.df_i, dist="weibull")
  shape       <- 1/reg.weibull$scale
  scale       <- exp(coef(reg.weibull)[1]+coef(reg.weibull)[2]*-850)
  pcarr[i]    <- ((1-exp(-(futurTime/scale)^(shape))) - 
                    (1-exp(-(carringtonTime/scale)^(shape))))/(exp(-(carringtonTime/scale)^(shape)))
}

median(pcarr); sd(pcarr); quantile(pcarr, 0.025); quantile(pcarr, 0.975)

### Probability of a Carrington event (Dst < -850nT) in the next decade (delta method)
reg.weibull <- survreg(Surv(value1)~threshold, data=mdata, dist="weibull")
carringtonTime <- as.numeric(difftime(Sys.Date(), as.Date("1859-09-01")))
futurTime <- as.POSIXlt(Sys.time())
futurTime$year <- futurTime$year+10
futurTime <- as.numeric(difftime(as.Date(futurTime), as.Date("1859-09-01")))
pr <- c(coef(reg.weibull), reg.weibull$scale)
names(pr) <- c("beta0", "beta1", "shape")
vcov.total <- as.data.frame(vcov(reg.weibull))
pcarr <- deltaMethod(pr, "((1-exp(-(futurTime/exp(beta0+beta1*-850))^(1/shape))) - 
                     (1-exp(-(carringtonTime/exp(beta0+beta1*-850))^(1/shape))))/(exp(-(carringtonTime/exp(beta0+beta1*-850))^(1/shape)))", vcov=as.matrix(vcov.total))$Estimate
se <- deltaMethod(pr, "((1-exp(-(futurTime/exp(beta0+beta1*-850))^(1/shape))) - 
                     (1-exp(-(carringtonTime/exp(beta0+beta1*-850))^(1/shape))))/(exp(-(carringtonTime/exp(beta0+beta1*-850))^(1/shape)))", vcov=as.matrix(vcov.total))$SE
pcarr; se ### probability and 95% confidence interval

### Probability of a Carrington event (Dst < -1760nT) in the next decade (bootstrap)
carringtonTime <- as.numeric(difftime(Sys.Date(), as.Date("1859-09-01")))
futurTime <- as.POSIXlt(Sys.time())
futurTime$year <- futurTime$year+10
futurTime <- as.numeric(difftime(as.Date(futurTime), as.Date("1859-09-01")))
pcarr <- vector()
for (i in 1:100000)
{
  boots.df_i  <- mdata[sample(nrow(mdata), dim(mdata)[1], replace=TRUE), ]
  reg.weibull <- survreg(Surv(value1)~threshold, data=boots.df_i, dist="weibull")
  shape       <- 1/reg.weibull$scale
  scale       <- exp(coef(reg.weibull)[1]+coef(reg.weibull)[2]*-1760)
  pcarr[i]    <- ((1-exp(-(futurTime/scale)^(shape))) - 
                    (1-exp(-(carringtonTime/scale)^(shape))))/(exp(-(carringtonTime/scale)^(shape)))
}

median(pcarr); sd(pcarr); quantile(pcarr, 0.025); quantile(pcarr, 0.975)

### Table 1
reg.weibull <- survreg(Surv(value)~threshold, data=mdata, dist="weibull")
# -100nT:
t1   <- 365
sh   <- 1/reg.weibull$scale
sc   <- exp(coef(reg.weibull)[1]+coef(reg.weibull)[2]*-100)
scal <- sc^(-sh) 
t1p  <- sum(dWeibullCount(0:100,sh,scal,method = c("series_mat"),time=t1)*seq(0,100,1))

# -200nT:
t1   <- 365
sh   <- 1/reg.weibull$scale
sc   <- exp(coef(reg.weibull)[1]+coef(reg.weibull)[2]*-200)
scal <- sc^(-sh) 
t2   <- sum(dWeibullCount(0:100,sh,scal,method = c("series_mat"),time=t1)*seq(0,100,1))

# -400nT:
t1   <- 365*10
sh   <- 1/reg.weibull$scale
sc   <- exp(coef(reg.weibull)[1]+coef(reg.weibull)[2]*-400)
scal <- sc^(-sh) 
t3   <- sum(dWeibullCount(0:100,sh,scal,method = c("series_mat"),time=t1)*seq(0,100,1))

# -800nT:
t1   <- 365*1000
sh   <- 1/reg.weibull$scale
sc   <- exp(coef(reg.weibull)[1]+coef(reg.weibull)[2]*-800)
scal <- sc^(-sh) 
t4   <- sum(dWeibullCount(0:100,sh,scal,method = c("series_mat"),time=t1)*seq(0,100,1))

# -1600nT:
t1   <- 365*1000000
sh   <- 1/reg.weibull$scale
sc   <- exp(coef(reg.weibull)[1]+coef(reg.weibull)[2]*-1600)
scal <- sc^(-sh) 
t5   <- sum(dWeibullCount(0:100,sh,scal,method = c("series_mat"),time=t1)*seq(0,100,1))

### Table S1
fitdist(as.numeric(na.omit(diff.times$diff15)), "weibull")$aic
fitdist(as.numeric(na.omit(diff.times$diff15)), "gamma")$aic
fitdist(as.numeric(na.omit(diff.times$diff15)), "exp")$aic

fitdist(as.numeric(na.omit(diff.times$diff20)), "weibull")$aic
fitdist(as.numeric(na.omit(diff.times$diff20)), "gamma", method="mme")$aic ### MLE does not converge
fitdist(as.numeric(na.omit(diff.times$diff20)), "exp")$aic

fitdist(as.numeric(na.omit(diff.times$diff25)), "weibull")$aic
fitdist(as.numeric(na.omit(diff.times$diff25)), "gamma")$aic ### MLE does not converge
fitdist(as.numeric(na.omit(diff.times$diff25)), "exp")$aic

fitdist(as.numeric(na.omit(diff.times$diff30)), "weibull")$aic
fitdist(as.numeric(na.omit(diff.times$diff30)), "gamma", method="mme")$aic  ### MLE does not converge
fitdist(as.numeric(na.omit(diff.times$diff30)), "exp", method="mme")$aic ### MLE does not converge

### Lilliefors corrected K-S test
LcKS(na.omit(diff.times$diff15), cdf="pweibull")$p.value
LcKS(na.omit(diff.times$diff15), cdf="pgamma")$p.value
LcKS(na.omit(diff.times$diff15), cdf="pexp")$p.value

LcKS(na.omit(diff.times$diff20), cdf="pweibull")$p.value
LcKS(na.omit(diff.times$diff20), cdf="pgamma")$p.value
LcKS(na.omit(diff.times$diff20), cdf="pexp")$p.value

LcKS(na.omit(diff.times$diff25), cdf="pweibull")$p.value
LcKS(na.omit(diff.times$diff25), cdf="pgamma")$p.value
LcKS(na.omit(diff.times$diff25), cdf="pexp")$p.value

LcKS(na.omit(diff.times$diff30), cdf="pweibull")$p.value
LcKS(na.omit(diff.times$diff30), cdf="pgamma")$p.value
LcKS(na.omit(diff.times$diff30), cdf="pexp")$p.value
