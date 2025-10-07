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

setwd("C:/Users/marre/Desktop/Modelling2025")
load("dstDEF.RData")
load("diff.timesDEF.RData")
load("Quebec.RData")
dst$CDATE2 <- as.Date(dst$CDATE)
dst2 <- dst[order(dst$CDATE), ]

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
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {  # Esta definint una funció anomenada multiplot
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

### Test if inter-occurrence times are Weibull or exponential for super-storms (si p < 0.05, podem concloure que no segueix aquella distribució )
                                                                                # (si p > 0.05, no podem concloure res)
LcKS(na.omit(diff.times$diff25), cdf="pweibull")$p.value
LcKS(na.omit(diff.times$diff25), cdf="pexp")$p.value

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
shape1[11]; scale1[11] # Weibull parameters corresponding to super-storms (Dst < -250nT) (i = 11 + 14 = 25)

### Model for shape and scale parameters (excluding inter-occurrence times below the 48h limit)
mdata <- melt(as.data.table(diff.times), measure = patterns("diff"))
for (i in 1:40)
{
  mdata$threshold[mdata$variable==paste0("diff", i)] <- -10*i
}
mdata <- mdata[!is.na(mdata$value), ]
### remove lower thresholds (independence is not clear) and fit the model
mdata <- mdata[mdata$threshold<=-150, ]

### Overall without sunspots
reg.exp     <- survreg(Surv(value)~threshold, data=mdata, dist="exp")
reg.weibull <- survreg(Surv(value)~threshold, data=mdata, dist="weibull")
reg.gamma   <- flexsurvreg(Surv(value)~threshold, data=mdata, dist="gamma")
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
  points(unique(mdata$threshold), coef(reg.weibull)[1]+coef(reg.weibull)[2]*unique(mdata$threshold), type="l")
dev.off()

### Probability of a Carrington event (Dst < -850nT) in the next decade (bootstrap)
carringtonTime <- as.numeric(difftime(Sys.Date(), as.Date("1859-09-01")))
futurTime <- as.POSIXlt(Sys.time())
futurTime$year <- futurTime$year+10
futurTime <- as.numeric(difftime(as.Date(futurTime), as.Date("1859-09-01")))
pcarr <- vector()
for (i in 1:1000)  # Fa subconjunts i recalcula amb els subconjunts (quan hi ha poques dades per calcular eficaçment els parameters, i fa la mitja despres)
{
  boots.df_i  <- mdata[sample(nrow(mdata), dim(mdata)[1], replace=TRUE), ]
  reg.weibull <- survreg(Surv(value)~threshold, data=boots.df_i, dist="weibull")
  shape       <- 1/reg.weibull$scale
  scale       <- exp(coef(reg.weibull)[1]+coef(reg.weibull)[2]*-850)
  pcarr[i]    <- ((1-exp(-(futurTime/scale)^(shape))) - 
                    (1-exp(-(carringtonTime/scale)^(shape))))/(exp(-(carringtonTime/scale)^(shape)))
}

median(pcarr); sd(pcarr); quantile(pcarr, 0.025); quantile(pcarr, 0.975)

### Probability of a Carrington event (Dst < -1760nT) in the next decade (bootstrap)
carringtonTime <- as.numeric(difftime(Sys.Date(), as.Date("1859-09-01")))
futurTime <- as.POSIXlt(Sys.time())
futurTime$year <- futurTime$year+10
futurTime <- as.numeric(difftime(as.Date(futurTime), as.Date("1859-09-01")))
pcarr <- vector()
for (i in 1:1000)
{
  boots.df_i  <- mdata[sample(nrow(mdata), dim(mdata)[1], replace=TRUE), ]
  reg.weibull <- survreg(Surv(value1)~threshold, data=boots.df_i, dist="weibull")
  shape       <- 1/reg.weibull$scale
  scale       <- exp(coef(reg.weibull)[1]+coef(reg.weibull)[2]*-1760)
  pcarr[i]    <- ((1-exp(-(futurTime/scale)^(shape))) - 
                    (1-exp(-(carringtonTime/scale)^(shape))))/(exp(-(carringtonTime/scale)^(shape)))
}

median(pcarr); sd(pcarr); quantile(pcarr, 0.025); quantile(pcarr, 0.975)

### Table 1 (del article, mirar l'article)
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

