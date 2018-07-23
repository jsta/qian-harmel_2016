base <- getwd()
  ### change this line to your working directory
dataDIR <- paste(base, "data", sep="/")

plotDIR <- paste(base, "manuscript", sep="/")
require(arm)
require(Matching)
require(lattice)
source("~/dropbox/rutil/util.r")

### MANAGER
load(paste(dataDIR, "manage.RData", sep="/"))

## processing data for the best phosphorus model:

LoadApr2007$AvgTPLoad[LoadApr2007$AvgTPLoad==0] <- 0.002 ## not using log(y+1)
LoadApr2007$AvgTNLoad[LoadApr2007$AvgTNLoad==0] <- 0.02 ## not using log(y+1)
LoadApr2007$yp <- log(LoadApr2007$AvgTPLoad)
LoadApr2007$xp <- log(LoadApr2007$Avg.P.Applied+1)  ## log(x+1) is fine
LoadApr2007$xr <- log(LoadApr2007$Avg.Runoff)
LoadApr2007$xs <- log(LoadApr2007$Avg.SoilL+1)
LoadApr2007$LU.App <- paste(LoadApr2007$ShortLU, LoadApr2007$F.AppMethod1)
LoadApr2007$F.AppMethod1 <- as.character(LoadApr2007$F.AppMethod1)
LoadApr2007$F.AppMethod1[LoadApr2007$F.AppMethod1=="Other"] <-"Surface Applied"
LoadApr2007$F.AppMethod1[LoadApr2007$F.AppMethod1=="Other2"] <-"Unknown"
LoadApr2007$F.AppMethod1 <- factor(LoadApr2007$F.AppMethod1)
LoadApr2007$CFarm <- 0
LoadApr2007$CFarm[LoadApr2007$CP1=="Contour Farming" | LoadApr2007$CP2=="Contour Farming" |
                  LoadApr2007$CP3=="Contour Farming" ]<-1
LoadApr2007$Buffer <- 0
LoadApr2007$Buffer[LoadApr2007$CP1=="Filter Strip" | LoadApr2007$CP2=="Filter Strip" |
                   LoadApr2007$CP3=="Filter Strip" |
                   LoadApr2007$CP1=="Riparian Buffer" | LoadApr2007$CP2=="Riparian Buffer"|
                   LoadApr2007$CP3=="Riparian Buffer"]<-1
LoadApr2007$Terrace <- 0
LoadApr2007$Terrace[LoadApr2007$CP1=="Terrace" | LoadApr2007$CP2=="Terrace" |
                    LoadApr2007$CP3=="Terrace" ]<-1
LoadApr2007$Waterway <- 0
LoadApr2007$Waterway[LoadApr2007$CP1=="Waterway" | LoadApr2007$CP2=="Waterway" |
                     LoadApr2007$CP3=="Waterway" ]<-1
LoadApr2007$AnyCP <- as.numeric(LoadApr2007$Waterway|LoadApr2007$Buffer|LoadApr2007$CFarm|
                                LoadApr2007$Terrace)
LoadApr2007$NumCP <- LoadApr2007$Waterway+LoadApr2007$Buffer+LoadApr2007$CFarm+
  LoadApr2007$Terrace
LoadApr2007$CP <- paste(LoadApr2007$Waterway, LoadApr2007$Buffer, LoadApr2007$CFarm,
                        LoadApr2007$Terrace)

## Figure 1 TP loading without matching

postscript(file=paste(plotDIR, "with_outCP.eps", sep="/"),
    height=5.5, width=5.5, horizontal=F)
par(mfcol=c(2,2), mar=c(2,3,1,1), mgp=c(1.75,0.125,0), las=1, tck=0.01)
plot(AvgTPLoad ~ factor(AnyCP), data=LoadApr2007, log="y", axes=F, xlab="",
     ylab="Average TP Load (kg/ha/yr)")
axis(1, at=c(1,2), labels=c("Without", "With"))
axis(2, at=c(0.01,0.05, 0.1,0.5,1,5,10), labels=c(0.01,0.05,0.1,0.5,1,5,10))
box()
plot(I(Avg.P.Applied+1) ~ factor(AnyCP), data=LoadApr2007, log="y", axes=F, xlab="",
     ylab="P Fertilizer Application (kg/ha/yr)")
axis(1, at=c(1,2), labels=c("Without", "With"))
axis(2)
box()
plot(AvgTNLoad ~ factor(AnyCP), data=LoadApr2007, log="y", axes=F, xlab="",
     ylab="Average TN Load (kg/ha/yr)")
axis(1, at=c(1,2), labels=c("Without", "With"))
axis(2, at=c(0.05, 0.1,0.5,1,5,10, 50), labels=c(0.05,0.1,0.5,1,5,10, 50))
box()
plot(I(Avg.N.Applied+1) ~ factor(AnyCP), data=LoadApr2007, log="y", axes=F, xlab="",
     ylab="N Fertilizer Application (kg/ha/yr)")
axis(1, at=c(1,2), labels=c("Without", "With"))
axis(2)
box()
dev.off()


## centering predictor
LoadApr2007$std.PApplied <- scale(log(LoadApr2007$Avg.P.Applied+1), scale=F)
LoadApr2007$std.NApplied <- scale(log(LoadApr2007$Avg.N.Applied+1), scale=F)
LoadApr2007$std.Runoff <- scale(log(LoadApr2007$Avg.Runoff), scale=F)
LoadApr2007$std.SoilL <- scale(log(LoadApr2007$Avg.SoilL+1), scale=F)
LoadApr2007$std.TPLoad <- scale(log(LoadApr2007$AvgTPLoad), scale=F)

## matching
clean.data <- LoadApr2007[,  c("AnyCP", "yp","std.PApplied","std.Runoff","std.SoilL",
                               "ShortLU","F.AppMethod1","Tillage","LU.App")]
clean.data <- na.omit(clean.data)
ps.m1 <- glm(AnyCP~std.PApplied*std.Runoff*std.SoilL*
             ShortLU*F.AppMethod1*Tillage, data=clean.data,
             family=binomial(link="logit"))
### a most complicated model, often overly complicated and not necessary

pscores1 <- predict(ps.m1, type="link")
matches1 <- matching(z=clean.data$AnyCP, score=pscores1)

matched1 <- clean.data[matches1$matched, ]
b.stats <- balance(clean.data, matched1, ps.m1)
#print(b.stats)
plot(b.stats)
by(matched1$yp, matched1$AnyCP, mean)
by(clean.data$yp, clean.data$AnyCP, mean)


## alternative methods -- may lead to different matched subset, hence, different effect
### 1. using package Matching: 
rr1 <- Match(Y=clean.data$yp, Tr=clean.data$AnyCP, X=ps.m1$fitted)
summary(rr1)
## a smaller effect

### 2. using additive model
ps.m2 <- glm(AnyCP~. , data=clean.data[,-2],
             family=binomial(link="logit"))
pscores2 <- predict(ps.m2, type="link")
matches2 <- matching(z=clean.data$AnyCP, score=pscores2)

matched2 <- clean.data[matches2$matched, ]
by(matched2$yp, matched2$AnyCP, mean)
### slightly larger effect
##### end of alternatives ###

## Figure 3 before and after matching
ctr <- attr(LoadApr2007$std.PApplied, "scaled:center")
postscript(file=paste(plotDIR, "psmatching.eps", sep="/"),
           width=5, height=2.5, horizontal=F)
par(mfrow=c(1,2),tck=0.01,las=1,mgp=c(1.5,0.125,0))
par(mar=c(2, 3, 2,0.25))
plot(std.PApplied~factor(AnyCP), data=clean.data, names=c("Without","With"),
     xlab=" ", ylab="P Applied (kg/ha/yr)", axes=F)
axis(2, at=log(c(0,5,10,50,100,250)+1)-ctr, labels=c(0,5,10,50,100,250))
axis(1, at=c(1,2), labels=c("Without","With"))
box()
title(main="before matching", cex=0.75)
par(mar=c(2,0.25,2,3))
plot(std.PApplied~factor(AnyCP), data=matched2, names=c("Without","With"),
     xlab=" ", ylab=" ", axes=F)
axis(1, at=c(1,2), labels=c("Without","With"))
axis(4, at=log(c(0,5,10,50,100,250)+1)-ctr, labels=c(0,5,10,50,100,250))
box()
title(main="after matching", cex=0.75)
dev.off()
##

## regression after matching
### t-test
Loadlm0.P <- lm(yp ~ AnyCP, data=matched2)
display(Loadlm0.P)

### with P.applied CP interaction
Loadlm1.P <- lm(yp~std.Runoff+ std.SoilL+std.PApplied*AnyCP, data=matched2)
display(Loadlm1.P)

plot(yp~std.PApplied, data=matched2)
plot(yp~std.PApplied, data=clean.data)

### removing interaction
Loadlm11.P <- lm(yp~std.Runoff+ std.SoilL+std.PApplied+AnyCP, data=matched2)
display(Loadlm11.P)
### removing P.applied
Loadlm12.P <- lm(yp~std.Runoff+ std.SoilL + AnyCP, data=matched2)
display(Loadlm12.P)


## without matching
Loadlm2.P <- lm(yp~std.Runoff+ std.SoilL+std.PApplied*AnyCP, data=clean.data)
display(Loadlm2.P)

### multilevel modeling
M10.1 <- lmer(yp ~ xp + xr + xs + (1+xp|ShortLU) +
            (xp-1|F.AppMethod1)+(xr-1|Tillage) + (xp-1|LU.App) + (1+xp|NumCP),
            data=LoadApr2007)


M10.2 <- lmer(yp ~ std.PApplied + std.Runoff+ std.SoilL+
              (1+std.PApplied|ShortLU) +
              (1+std.PApplied| F.AppMethod1)+
              (1+std.Runoff|Tillage) +(1+std.PApplied|LU.App),
              data=LoadApr2007)

par(mfrow=c(1,2), mgp =c(1.5,.5,0), mar=c(3, 5, 1, 0.25))
plot.lmer.ranef(M=M10.2, ranV=3, column=1, xlab="Intercept", ylab=" ", column.fixed=1, y.cex=0.7)
par(mar=c(3, 0.25, 1, 5))
plot.lmer.ranef(M=M10.2, ranV=3, column=2, xlab="Runoff", ylab=" ", yaxis=4, column.fixed=2, y.cex=0.7)


M10.3 <- lmer(yp ~ std.PApplied + std.Runoff+ std.SoilL+
              (1+std.PApplied|ShortLU) +
              (1+std.PApplied| F.AppMethod1)+
              (1+std.Runoff|Tillage) +
              (1+std.PApplied|LU.App)+(1+std.PApplied|NumCP),
            data=LoadApr2007)

M10.4 <- lmer(yp ~ std.PApplied + std.Runoff+ std.SoilL + AnyCP + std.PApplied:AnyCP +
              (1+std.PApplied+ AnyCP+std.PApplied:AnyCP|ShortLU) +
              ( AnyCP-1| F.AppMethod1)+(AnyCP-1+std.PApplied:AnyCP|Tillage) +
              (1|LU.App),
              data=LoadApr2007)

M10.5 <- lmer(yp ~ std.Runoff+ std.SoilL + std.PApplied * AnyCP+
              (AnyCP+std.PApplied : AnyCP-1| F.AppMethod1)+
              (1+AnyCP+std.PApplied : AnyCP|ShortLU),
              data=LoadApr2007)
## the best model

mm<-M10.5
display(mm)


fixed <- fixef(mm)
fixed.se <- se.fixef(mm)

ran <- ranef(mm)
ran.se <- se.ranef(mm)

random <- ran[[1]]
n.levels <- length(row.names(random))

plot.data1 <- data.frame(LU = rep(row.names(random), 4), lu.row =rep(1:n.levels, 4),
                        CP = rep(c(1,0), each=2*n.levels), x = rep(range(LoadApr2007$std.PApplied, na.rm=T), each=n.levels))
plot.data1$y <- (fixed[1] + random[plot.data1$lu.row, 1])+
  (fixed[4])* plot.data1$x +
  (fixed[5] + random[plot.data1$lu.row, 2]) * plot.data1$CP +
  (fixed[6] + random[plot.data1$lu.row, 3]) * plot.data1$x*plot.data1$CP
plot.data1$yp <- NA
plot.data1 <- plot.data1[,-2]

plot.data2 <- na.omit(data.frame(LU=LoadApr2007$ShortLU, CP=as.numeric(LoadApr2007$AnyCP),
                                 x=LoadApr2007$std.PApplied, yp=LoadApr2007$yp))
plot.data2$y<- NA

plot.data <- make.groups(plot.data1, plot.data2)


temp <- plot.data$CP
plot.data$CP1 <- "With"
plot.data$CP1[temp==0] <- "Without"

plot.data$CP1
temp <- plot.data$LU!="Other" & plot.data$LU!="Soybeans"
plot.data <- plot.data[temp,]
plot.data$LU <- ordered(as.character(plot.data$LU))

## figure 3
trellis.device(postscript, file=paste(plotDIR, "manageEff1.eps", sep="/"),
               width=4.5, height=5, horizontal=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.set(list(fontsize=list(text=8),
                 par.xlab.text=list(cex=1.25),
                     add.text=list(cex=1.25),
                     superpose.symbol=list(cex=1)))
key <- simpleKey(unique(plot.data$CP1), lines=T, points=F, space = "top", columns=2)
key$text$cex <- 1.25
xyplot(y ~ x|LU, data=plot.data, type="l",
    group=plot.data$CP1, key=key, xlab="centered log P Applied", ylab="Log TP Loading",
    panel=function(x,y,...){
        panel.xyplot(x,y,lwd=1.5,...)
        panel.grid()
    }, layout=c(3,3)
    )
dev.off()

## data points (not used in the paper)
trellis.device(postscript, file="manageEff12.eps", width=4.5, height=5, horizontal=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.set(list(fontsize=list(text=8),
                 par.xlab.text=list(cex=1.25),
                     add.text=list(cex=1.25),
                     superpose.symbol=list(cex=1)))
key <- simpleKey(unique(plot.data$CP1), lines=F, points=T, space = "top", columns=2)
key$text$cex <- 1.25
xyplot(yp ~ x|LU, data=plot.data,type="p",
    group=plot.data$CP1, key=key, xlab="log TP Applied", ylab="Log TP Loading",
#    panel=panel.superpose,
#    distribute.type=TRUE,
#    type=c("p","l")
    panel=function(x,y,...){
        panel.xyplot(x,y,lwd=1.5,...)
        panel.grid()
    }#,
    #scales=list(x=list(at=log(c(0, 10, 50, 100, 500, 1000)+1), labels=as.character(c(0, 10, 50, 100, 500, 1000))))
    )
dev.off()

## other models tried
M10.2 <- lmer(yp ~ std.PApplied + std.Runoff+ std.SoilL+std.PApplied*AnyCP+
            (1+std.PApplied| F.AppMethod1)+(1+std.Runoff|Tillage) +(1+std.PApplied|LU.App),
            data=LoadApr2007)
par(mfrow=c(1,2), mgp =c(1.5,.5,0), mar=c(3, 5, 1, 0.25))
plot.lmer.ranef(M=M10.2, ranV=1, column=1, xlab="Intercept", ylab=" ", column.fixed=1, y.cex=0.7)
par(mar=c(3, 0.25, 1, 5))
plot.lmer.ranef(M=M10.2, ranV=1, column=2, xlab="P Applied", ylab=" ", yaxis=4, column.fixed=2, y.cex=0.7)

M10.3 <- lmer2(yp ~ std.PApplied + std.Runoff+ std.SoilL+std.PApplied*AnyCP+
            (1+std.PApplied| F.AppMethod1)+(1+std.Runoff|Tillage) +(1|LU.App),
            data=LoadApr2007)
par(mfrow=c(1,2), mgp =c(1.5,.5,0), mar=c(3, 5, 1, 0.25))
plot.lmer.ranef(M=M10.3, ranV=3, column=1, xlab="Intercept", ylab=" ", column.fixed=1, y.cex=0.7)
par(mar=c(3, 0.25, 1, 5))
plot.lmer.ranef(M=M10.3, ranV=3, column=2, xlab="P Applied", ylab=" ", yaxis=4, column.fixed=2, y.cex=0.7)

M10.4 <- lmer(yp ~ std.PApplied + std.Runoff+ std.SoilL+std.PApplied*AnyCP+
            (std.PApplied-1| F.AppMethod1)+(1+std.Runoff|Tillage) +(1|LU.App),
            data=LoadApr2007, subset=xp>0)
par(mfrow=c(1,2), mgp =c(1.5,.5,0), mar=c(3, 5, 1, 0.25))
plot.lmer.ranef(M=M10.4, ranV=3, column=1, xlab="Intercept", ylab=" ", column.fixed=1, y.cex=0.7)
par(mar=c(3, 0.25, 1, 5))
plot.lmer.ranef(M=M10.4, ranV=3, column=2, xlab="P Applied", ylab=" ", yaxis=4, column.fixed=2, y.cex=0.7)

M10.5 <- lmer(yp ~ std.Runoff+ std.SoilL + std.PApplied * AnyCP+
            (std.PApplied+AnyCP| F.AppMethod1)+(1+std.Runoff+AnyCP|Tillage) +(1|LU.App),
            data=LoadApr2007, subset=xp>0)

M10.55 <- lmer(yp ~ std.PApplied + std.Runoff+ std.SoilL + std.PApplied * AnyCP+
            (std.PApplied| F.AppMethod1)+(1+std.Runoff|Tillage) +(1|LU.App),
            data=LoadApr2007)

par(mfrow=c(1,2), mgp =c(1.5,.5,0), mar=c(3, 5, 1, 0.25))
plot.lmer.ranef(M=M10.5, ranV=2, column=3, xlab="AnyCP", ylab=" ", column.fixed=5, y.cex=0.7)
par(mar=c(3, 0.25, 1, 5))
plot.lmer.ranef(M=M10.5, ranV=3, column=3, xlab="AnyCP", ylab=" ", yaxis=4, column.fixed=5, y.cex=0.7)

### using CP (and P applied > 0)
M10.6 <- lmer(yp ~ std.PApplied + std.Runoff+ std.SoilL+ (1+std.PApplied|CP)+
            (std.PApplied| F.AppMethod1)+(1+std.Runoff|Tillage) +(1|LU.App),
            data=LoadApr2007, subset=xp>0)

## Figure 2
postscript(file=paste(plotDIR, "mlindCPs.eps", sep="/"),
           height=4, width=4.5, horizontal=F)
par(mfrow=c(1,2), mgp =c(1.5,.5,0), mar=c(3, 5, 1, 0.25))
plot.lmer.ranef(M=M10.6, ranV=2, column=1, xlab="Intercept", ylab=" ", column.fixed=1, y.cex=0.7)
par(mar=c(3, 0.25, 1, 5))
plot.lmer.ranef(M=M10.6, ranV=2, column=2, xlab="log P Applied", ylab=" ", yaxis=4, column.fixed=2, y.cex=0.7)
dev.off()
##

## other exploratory plots
xyplot(yp~xp|ShortLU, data=LoadApr2007, subset=xp>0&CP=="0 0 0 0")
xyplot(yp~xp|ShortLU, data=LoadApr2007, subset=xp>0&CP!="0 0 0 0")
xyplot(yp~xp|ShortLU*CP, data=LoadApr2007, subset=xp)
xyplot(yp~xp|ShortLU*CP, data=LoadApr2007, subset=xp>0)
xyplot(yp~xp|ShortLU*CP, data=LoadApr2007)

dotplot(CP~yp|Tillage, data=LoadApr2007)
dotplot(CP~yp|Tillage, data=LoadApr2007, subset=xp>0)

