#6 March 2020
#TNS

#script to read in results of flow cytometry of isogenic mutants (for various purposes), fit titration curves and some interpretation from each experiment

setwd("~/bloom_j/computational_notebooks/tstarr/2019/NIH45-46_DMS")

library(ggplot2)
library(data.table)

#dir.create(file.path("results/isogenic_titrations"))

###################################
#4 March 2020 experiment: tested 10 single mutants to compare isogenic ddG (and E) measurements to those determined from library experiments. Also tested four double mutants for epistasis
dt_20200304 <- read.csv(file="data/isogenic_titrations/20200304.csv", stringsAsFactors=F)
dt_20200304$conc_nM <- dt_20200304$conc_M*1e9

#fit curves to each titration
#fit binding curves to estimate Kd for "mean bin" metric
fit_20200304_1 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                data=dt_20200304[dt_20200304$titration==1,],start=list(a=3,b=1,Kd=0.1))
fit_20200304_2 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                data=dt_20200304[dt_20200304$titration==2,],start=list(a=3,b=1,Kd=0.1))
fit_20200304_3 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                data=dt_20200304[dt_20200304$titration==3,],start=list(a=3,b=1,Kd=0.1))
fit_20200304_4 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                data=dt_20200304[dt_20200304$titration==4,],start=list(a=3,b=1,Kd=5))
fit_20200304_5 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                data=dt_20200304[dt_20200304$titration==5,],start=list(a=3,b=1,Kd=5))
fit_20200304_6 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                data=dt_20200304[dt_20200304$titration==6,],start=list(a=3,b=1,Kd=5))
fit_20200304_7 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                data=dt_20200304[dt_20200304$titration==7,],start=list(a=3,b=1,Kd=5))
fit_20200304_8 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                data=dt_20200304[dt_20200304$titration==8,],start=list(a=3,b=1,Kd=5))
fit_20200304_9 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                data=dt_20200304[dt_20200304$titration==9,],start=list(a=3,b=1,Kd=5))
fit_20200304_10 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                 data=dt_20200304[dt_20200304$titration==10,],start=list(a=3,b=1,Kd=5))
fit_20200304_11 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                 data=dt_20200304[dt_20200304$titration==11,],start=list(a=3,b=1,Kd=5))
fit_20200304_12 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                 data=dt_20200304[dt_20200304$titration==12,],start=list(a=3,b=1,Kd=5))
fit_20200304_13 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                 data=dt_20200304[dt_20200304$titration==13,],start=list(a=3,b=1,Kd=5))
fit_20200304_14 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                 data=dt_20200304[dt_20200304$titration==14,],start=list(a=3,b=1,Kd=5))
fit_20200304_15 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                 data=dt_20200304[dt_20200304$titration==15,],start=list(a=3,b=1,Kd=5))


######
pdf(file="results/isogenic_titrations/20200304.pdf",width=6,height=10,useDingbats=F)
#plot titration curves for mean bin
options(repr.plot.width=6, repr.plot.height=10, scipen=1)
par(mfrow=c(5,3),oma=c(4,4,0,0),mar=c(1,1,1,1))

ymax <- 4
ymin <- 1

#titration==1, fit_20200304_1
plot(dt_20200304[dt_20200304$titration==1,"conc_nM"],
     dt_20200304[dt_20200304$titration==1,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200304[dt_20200304$titration==1,"genotype"])),xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200304[dt_20200304$titration==1 & dt_20200304$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200304_1,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200304_1)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200304_1)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==2, fit_20200304_2
plot(dt_20200304[dt_20200304$titration==2,"conc_nM"],
     dt_20200304[dt_20200304$titration==2,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200304[dt_20200304$titration==2,"genotype"])),xaxt='n',yaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200304[dt_20200304$titration==2 & dt_20200304$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F); axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200304_2,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200304_2)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200304_2)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==3, fit_20200304_3
plot(dt_20200304[dt_20200304$titration==3,"conc_nM"],
     dt_20200304[dt_20200304$titration==3,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200304[dt_20200304$titration==3,"genotype"])),yaxt='n',xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200304[dt_20200304$titration==3 & dt_20200304$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F);axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200304_3,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200304_3)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200304_3)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==4, fit_20200304_4
plot(dt_20200304[dt_20200304$titration==4,"conc_nM"],
     dt_20200304[dt_20200304$titration==4,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200304[dt_20200304$titration==4,"genotype"])),xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200304[dt_20200304$titration==4 & dt_20200304$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200304_4,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200304_4)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200304_4)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==5, fit_20200304_5
plot(dt_20200304[dt_20200304$titration==5,"conc_nM"],
     dt_20200304[dt_20200304$titration==5,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200304[dt_20200304$titration==5,"genotype"])),xaxt='n',yaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200304[dt_20200304$titration==5 & dt_20200304$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F); axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200304_5,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200304_5)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200304_5)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==6, fit_20200304_6
plot(dt_20200304[dt_20200304$titration==6,"conc_nM"],
     dt_20200304[dt_20200304$titration==6,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200304[dt_20200304$titration==6,"genotype"])),yaxt='n',xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200304[dt_20200304$titration==6 & dt_20200304$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F);axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200304_6,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200304_6)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200304_6)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==7, fit_20200304_7
plot(dt_20200304[dt_20200304$titration==7,"conc_nM"],
     dt_20200304[dt_20200304$titration==7,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200304[dt_20200304$titration==7,"genotype"])),xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200304[dt_20200304$titration==7 & dt_20200304$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200304_7,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200304_7)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200304_7)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==8, fit_20200304_8
plot(dt_20200304[dt_20200304$titration==8,"conc_nM"],
     dt_20200304[dt_20200304$titration==8,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200304[dt_20200304$titration==8,"genotype"])),xaxt='n',yaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200304[dt_20200304$titration==8 & dt_20200304$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F); axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200304_8,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200304_8)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200304_8)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==9, fit_20200304_9
plot(dt_20200304[dt_20200304$titration==9,"conc_nM"],
     dt_20200304[dt_20200304$titration==9,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200304[dt_20200304$titration==9,"genotype"])),yaxt='n',xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200304[dt_20200304$titration==9 & dt_20200304$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F);axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200304_9,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200304_9)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200304_9)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==10, fit_20200304_10
plot(dt_20200304[dt_20200304$titration==10,"conc_nM"],
     dt_20200304[dt_20200304$titration==10,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200304[dt_20200304$titration==10,"genotype"])),yaxt='n',xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200304[dt_20200304$titration==10 & dt_20200304$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F);axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200304_10,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200304_10)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200304_10)$coefficients["Kd","Std. Error"],digits=1),"nM")))
#titration==11, fit_20200304_11
plot(dt_20200304[dt_20200304$titration==11,"conc_nM"],
     dt_20200304[dt_20200304$titration==11,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200304[dt_20200304$titration==11,"genotype"])),yaxt='n',xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200304[dt_20200304$titration==11 & dt_20200304$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F);axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200304_11,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200304_11)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200304_11)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==12, fit_20200304_12
plot(dt_20200304[dt_20200304$titration==12,"conc_nM"],
     dt_20200304[dt_20200304$titration==12,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200304[dt_20200304$titration==12,"genotype"])),yaxt='n',xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200304[dt_20200304$titration==12 & dt_20200304$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F);axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200304_12,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200304_12)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200304_12)$coefficients["Kd","Std. Error"],digits=1),"nM")))
#titration==13, fit_20200304_13
plot(dt_20200304[dt_20200304$titration==13,"conc_nM"],
     dt_20200304[dt_20200304$titration==13,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200304[dt_20200304$titration==13,"genotype"])),ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200304[dt_20200304$titration==13 & dt_20200304$conc_nM==0, "mean_bin"][1])
lines(seq(0.0001,3200,0.001),predict(fit_20200304_13,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200304_13)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200304_13)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==14, fit_20200304_14
plot(dt_20200304[dt_20200304$titration==14,"conc_nM"],
     dt_20200304[dt_20200304$titration==14,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200304[dt_20200304$titration==14,"genotype"])),yaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200304[dt_20200304$titration==14 & dt_20200304$conc_nM==0, "mean_bin"][1])
axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200304_14,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200304_14)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200304_14)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==15, fit_20200304_15
plot(dt_20200304[dt_20200304$titration==15,"conc_nM"],
     dt_20200304[dt_20200304$titration==15,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200304[dt_20200304$titration==15,"genotype"])),yaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200304[dt_20200304$titration==15 & dt_20200304$conc_nM==0, "mean_bin"][1])
axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200304_15,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200304_15)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200304_15)$coefficients["Kd","Std. Error"],digits=1),"nM")))
axis(side=2,labels=F)
mtext("ligand concentration (nM, log scale)", side=1, outer=TRUE, line=2)
mtext("mean bin of FITC+ cells", side=2, outer=TRUE, line=2)
dev.off()
###################################
#11 March 2020 experiment: tested same mutants, but note the change in order of variants. Included a higher concentration (3.2uM) which before was no good because of either lowered volume that I attempted or incomplete emptying of wells during washing
dt_20200311 <- read.csv(file="data/isogenic_titrations/20200311.csv", stringsAsFactors=F)
dt_20200311$conc_nM <- dt_20200311$conc_M*1e9

#fit curves to each titration
#fit binding curves to estimate Kd for "mean bin" metric
fit_20200311_1 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                      data=dt_20200311[dt_20200311$titration==1,],start=list(a=3,b=1,Kd=0.1))
fit_20200311_2 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                      data=dt_20200311[dt_20200311$titration==2,],start=list(a=3,b=1,Kd=0.1))
fit_20200311_3 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                      data=dt_20200311[dt_20200311$titration==3,],start=list(a=3,b=1,Kd=0.1))
fit_20200311_4 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                      data=dt_20200311[dt_20200311$titration==4,],start=list(a=3,b=1,Kd=5))
fit_20200311_5 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                      data=dt_20200311[dt_20200311$titration==5,],start=list(a=3,b=1,Kd=5))
fit_20200311_6 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                      data=dt_20200311[dt_20200311$titration==6,],start=list(a=3,b=1,Kd=5))
fit_20200311_7 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                      data=dt_20200311[dt_20200311$titration==7,],start=list(a=3,b=1,Kd=5))
fit_20200311_8 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                      data=dt_20200311[dt_20200311$titration==8,],start=list(a=3,b=1,Kd=5))
fit_20200311_9 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                      data=dt_20200311[dt_20200311$titration==9,],start=list(a=3,b=1,Kd=5))
fit_20200311_10 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                       data=dt_20200311[dt_20200311$titration==10,],start=list(a=3,b=1,Kd=5))
fit_20200311_11 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                       data=dt_20200311[dt_20200311$titration==11,],start=list(a=3,b=1,Kd=5))
fit_20200311_12 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                       data=dt_20200311[dt_20200311$titration==12,],start=list(a=3,b=1,Kd=5))
fit_20200311_13 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                       data=dt_20200311[dt_20200311$titration==13,],start=list(a=3,b=1,Kd=5))
fit_20200311_14 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                       data=dt_20200311[dt_20200311$titration==14,],start=list(a=3,b=1,Kd=5))
fit_20200311_15 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                       data=dt_20200311[dt_20200311$titration==15,],start=list(a=3,b=1,Kd=5))


######
pdf(file="results/isogenic_titrations/20200311.pdf",width=6,height=10,useDingbats=F)
#plot titration curves for mean bin
options(repr.plot.width=6, repr.plot.height=10, scipen=1)
par(mfrow=c(5,3),oma=c(4,4,0,0),mar=c(1,1,1,1))

ymax <- 4
ymin <- 1

#titration==1, fit_20200311_1
plot(dt_20200311[dt_20200311$titration==1,"conc_nM"],
     dt_20200311[dt_20200311$titration==1,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200311[dt_20200311$titration==1,"genotype"])),xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200311[dt_20200311$titration==1 & dt_20200311$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200311_1,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200311_1)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200311_1)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==2, fit_20200311_2
plot(dt_20200311[dt_20200311$titration==2,"conc_nM"],
     dt_20200311[dt_20200311$titration==2,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200311[dt_20200311$titration==2,"genotype"])),xaxt='n',yaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200311[dt_20200311$titration==2 & dt_20200311$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F); axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200311_2,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200311_2)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200311_2)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==3, fit_20200311_3
plot(dt_20200311[dt_20200311$titration==3,"conc_nM"],
     dt_20200311[dt_20200311$titration==3,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200311[dt_20200311$titration==3,"genotype"])),yaxt='n',xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200311[dt_20200311$titration==3 & dt_20200311$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F);axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200311_3,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200311_3)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200311_3)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==4, fit_20200311_4
plot(dt_20200311[dt_20200311$titration==4,"conc_nM"],
     dt_20200311[dt_20200311$titration==4,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200311[dt_20200311$titration==4,"genotype"])),xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200311[dt_20200311$titration==4 & dt_20200311$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200311_4,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200311_4)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200311_4)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==5, fit_20200311_5
plot(dt_20200311[dt_20200311$titration==5,"conc_nM"],
     dt_20200311[dt_20200311$titration==5,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200311[dt_20200311$titration==5,"genotype"])),xaxt='n',yaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200311[dt_20200311$titration==5 & dt_20200311$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F); axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200311_5,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200311_5)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200311_5)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==6, fit_20200311_6
plot(dt_20200311[dt_20200311$titration==6,"conc_nM"],
     dt_20200311[dt_20200311$titration==6,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200311[dt_20200311$titration==6,"genotype"])),yaxt='n',xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200311[dt_20200311$titration==6 & dt_20200311$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F);axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200311_6,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200311_6)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200311_6)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==7, fit_20200311_7
plot(dt_20200311[dt_20200311$titration==7,"conc_nM"],
     dt_20200311[dt_20200311$titration==7,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200311[dt_20200311$titration==7,"genotype"])),xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200311[dt_20200311$titration==7 & dt_20200311$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200311_7,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200311_7)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200311_7)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==8, fit_20200311_8
plot(dt_20200311[dt_20200311$titration==8,"conc_nM"],
     dt_20200311[dt_20200311$titration==8,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200311[dt_20200311$titration==8,"genotype"])),xaxt='n',yaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200311[dt_20200311$titration==8 & dt_20200311$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F); axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200311_8,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200311_8)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200311_8)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==9, fit_20200311_9
plot(dt_20200311[dt_20200311$titration==9,"conc_nM"],
     dt_20200311[dt_20200311$titration==9,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200311[dt_20200311$titration==9,"genotype"])),yaxt='n',xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200311[dt_20200311$titration==9 & dt_20200311$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F);axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200311_9,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200311_9)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200311_9)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==10, fit_20200311_10
plot(dt_20200311[dt_20200311$titration==10,"conc_nM"],
     dt_20200311[dt_20200311$titration==10,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200311[dt_20200311$titration==10,"genotype"])),yaxt='n',xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200311[dt_20200311$titration==10 & dt_20200311$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F);axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200311_10,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200311_10)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200311_10)$coefficients["Kd","Std. Error"],digits=1),"nM")))
#titration==11, fit_20200311_11
plot(dt_20200311[dt_20200311$titration==11,"conc_nM"],
     dt_20200311[dt_20200311$titration==11,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200311[dt_20200311$titration==11,"genotype"])),yaxt='n',xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200311[dt_20200311$titration==11 & dt_20200311$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F);axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200311_11,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200311_11)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200311_11)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==12, fit_20200311_12
plot(dt_20200311[dt_20200311$titration==12,"conc_nM"],
     dt_20200311[dt_20200311$titration==12,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200311[dt_20200311$titration==12,"genotype"])),yaxt='n',xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200311[dt_20200311$titration==12 & dt_20200311$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F);axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200311_12,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200311_12)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200311_12)$coefficients["Kd","Std. Error"],digits=1),"nM")))
#titration==13, fit_20200311_13
plot(dt_20200311[dt_20200311$titration==13,"conc_nM"],
     dt_20200311[dt_20200311$titration==13,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200311[dt_20200311$titration==13,"genotype"])),ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200311[dt_20200311$titration==13 & dt_20200311$conc_nM==0, "mean_bin"][1])
lines(seq(0.0001,3200,0.001),predict(fit_20200311_13,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200311_13)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200311_13)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==14, fit_20200311_14
plot(dt_20200311[dt_20200311$titration==14,"conc_nM"],
     dt_20200311[dt_20200311$titration==14,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200311[dt_20200311$titration==14,"genotype"])),yaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200311[dt_20200311$titration==14 & dt_20200311$conc_nM==0, "mean_bin"][1])
axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200311_14,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200311_14)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200311_14)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==15, fit_20200311_15
plot(dt_20200311[dt_20200311$titration==15,"conc_nM"],
     dt_20200311[dt_20200311$titration==15,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200311[dt_20200311$titration==15,"genotype"])),yaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200311[dt_20200311$titration==15 & dt_20200311$conc_nM==0, "mean_bin"][1])
axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200311_15,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200311_15)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200311_15)$coefficients["Kd","Std. Error"],digits=1),"nM")))
axis(side=2,labels=F)
mtext("ligand concentration (nM, log scale)", side=1, outer=TRUE, line=2)
mtext("mean bin of FITC+ cells", side=2, outer=TRUE, line=2)
dev.off()



###################################
#13 March 2020 experiment: tested same mutants with same conditions as 20200311 experiment. Note I mixed up some wells, but this should be appropriately switched in the indexing in the source data sheet
dt_20200313 <- read.csv(file="data/isogenic_titrations/20200313.csv", stringsAsFactors=F)
dt_20200313$conc_nM <- dt_20200313$conc_M*1e9

#fit curves to each titration
#fit binding curves to estimate Kd for "mean bin" metric
fit_20200313_1 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                      data=dt_20200313[dt_20200313$titration==1,],start=list(a=3,b=1,Kd=0.1))
fit_20200313_2 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                      data=dt_20200313[dt_20200313$titration==2,],start=list(a=3,b=1,Kd=0.1))
fit_20200313_3 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                      data=dt_20200313[dt_20200313$titration==3,],start=list(a=3,b=1,Kd=0.1))
fit_20200313_4 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                      data=dt_20200313[dt_20200313$titration==4,],start=list(a=3,b=1,Kd=5))
fit_20200313_5 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                      data=dt_20200313[dt_20200313$titration==5,],start=list(a=3,b=1,Kd=5))
fit_20200313_6 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                      data=dt_20200313[dt_20200313$titration==6,],start=list(a=3,b=1,Kd=5))
fit_20200313_7 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                      data=dt_20200313[dt_20200313$titration==7,],start=list(a=3,b=1,Kd=5))
fit_20200313_8 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                      data=dt_20200313[dt_20200313$titration==8,],start=list(a=3,b=1,Kd=5))
fit_20200313_9 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                      data=dt_20200313[dt_20200313$titration==9,],start=list(a=3,b=1,Kd=5))
fit_20200313_10 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                       data=dt_20200313[dt_20200313$titration==10,],start=list(a=3,b=1,Kd=5))
fit_20200313_11 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                       data=dt_20200313[dt_20200313$titration==11,],start=list(a=3,b=1,Kd=5))
fit_20200313_12 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                       data=dt_20200313[dt_20200313$titration==12,],start=list(a=3,b=1,Kd=5))
fit_20200313_13 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                       data=dt_20200313[dt_20200313$titration==13,],start=list(a=3,b=1,Kd=5))
fit_20200313_14 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                       data=dt_20200313[dt_20200313$titration==14,],start=list(a=3,b=1,Kd=5))
fit_20200313_15 <- nls(mean_bin ~ a*(conc_nM/(conc_nM+Kd))+b,
                       data=dt_20200313[dt_20200313$titration==15,],start=list(a=3,b=1,Kd=5))


######
pdf(file="results/isogenic_titrations/20200313.pdf",width=6,height=10,useDingbats=F)
#plot titration curves for mean bin
options(repr.plot.width=6, repr.plot.height=10, scipen=1)
par(mfrow=c(5,3),oma=c(4,4,0,0),mar=c(1,1,1,1))

ymax <- 4
ymin <- 1

#titration==1, fit_20200313_1
plot(dt_20200313[dt_20200313$titration==1,"conc_nM"],
     dt_20200313[dt_20200313$titration==1,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200313[dt_20200313$titration==1,"genotype"])),xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200313[dt_20200313$titration==1 & dt_20200313$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200313_1,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200313_1)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200313_1)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==2, fit_20200313_2
plot(dt_20200313[dt_20200313$titration==2,"conc_nM"],
     dt_20200313[dt_20200313$titration==2,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200313[dt_20200313$titration==2,"genotype"])),xaxt='n',yaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200313[dt_20200313$titration==2 & dt_20200313$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F); axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200313_2,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200313_2)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200313_2)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==3, fit_20200313_3
plot(dt_20200313[dt_20200313$titration==3,"conc_nM"],
     dt_20200313[dt_20200313$titration==3,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200313[dt_20200313$titration==3,"genotype"])),yaxt='n',xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200313[dt_20200313$titration==3 & dt_20200313$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F);axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200313_3,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200313_3)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200313_3)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==4, fit_20200313_4
plot(dt_20200313[dt_20200313$titration==4,"conc_nM"],
     dt_20200313[dt_20200313$titration==4,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200313[dt_20200313$titration==4,"genotype"])),xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200313[dt_20200313$titration==4 & dt_20200313$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200313_4,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200313_4)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200313_4)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==5, fit_20200313_5
plot(dt_20200313[dt_20200313$titration==5,"conc_nM"],
     dt_20200313[dt_20200313$titration==5,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200313[dt_20200313$titration==5,"genotype"])),xaxt='n',yaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200313[dt_20200313$titration==5 & dt_20200313$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F); axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200313_5,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200313_5)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200313_5)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==6, fit_20200313_6
plot(dt_20200313[dt_20200313$titration==6,"conc_nM"],
     dt_20200313[dt_20200313$titration==6,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200313[dt_20200313$titration==6,"genotype"])),yaxt='n',xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200313[dt_20200313$titration==6 & dt_20200313$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F);axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200313_6,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200313_6)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200313_6)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==7, fit_20200313_7
plot(dt_20200313[dt_20200313$titration==7,"conc_nM"],
     dt_20200313[dt_20200313$titration==7,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200313[dt_20200313$titration==7,"genotype"])),xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200313[dt_20200313$titration==7 & dt_20200313$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200313_7,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200313_7)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200313_7)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==8, fit_20200313_8
plot(dt_20200313[dt_20200313$titration==8,"conc_nM"],
     dt_20200313[dt_20200313$titration==8,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200313[dt_20200313$titration==8,"genotype"])),xaxt='n',yaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200313[dt_20200313$titration==8 & dt_20200313$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F); axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200313_8,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200313_8)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200313_8)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==9, fit_20200313_9
plot(dt_20200313[dt_20200313$titration==9,"conc_nM"],
     dt_20200313[dt_20200313$titration==9,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200313[dt_20200313$titration==9,"genotype"])),yaxt='n',xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200313[dt_20200313$titration==9 & dt_20200313$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F);axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200313_9,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200313_9)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200313_9)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==10, fit_20200313_10
plot(dt_20200313[dt_20200313$titration==10,"conc_nM"],
     dt_20200313[dt_20200313$titration==10,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200313[dt_20200313$titration==10,"genotype"])),yaxt='n',xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200313[dt_20200313$titration==10 & dt_20200313$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F);axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200313_10,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200313_10)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200313_10)$coefficients["Kd","Std. Error"],digits=1),"nM")))
#titration==11, fit_20200313_11
plot(dt_20200313[dt_20200313$titration==11,"conc_nM"],
     dt_20200313[dt_20200313$titration==11,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200313[dt_20200313$titration==11,"genotype"])),yaxt='n',xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200313[dt_20200313$titration==11 & dt_20200313$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F);axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200313_11,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200313_11)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200313_11)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==12, fit_20200313_12
plot(dt_20200313[dt_20200313$titration==12,"conc_nM"],
     dt_20200313[dt_20200313$titration==12,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200313[dt_20200313$titration==12,"genotype"])),yaxt='n',xaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200313[dt_20200313$titration==12 & dt_20200313$conc_nM==0, "mean_bin"][1])
axis(side=1,labels=F);axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200313_12,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200313_12)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200313_12)$coefficients["Kd","Std. Error"],digits=1),"nM")))
#titration==13, fit_20200313_13
plot(dt_20200313[dt_20200313$titration==13,"conc_nM"],
     dt_20200313[dt_20200313$titration==13,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200313[dt_20200313$titration==13,"genotype"])),ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200313[dt_20200313$titration==13 & dt_20200313$conc_nM==0, "mean_bin"][1])
lines(seq(0.0001,3200,0.001),predict(fit_20200313_13,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200313_13)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200313_13)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==14, fit_20200313_14
plot(dt_20200313[dt_20200313$titration==14,"conc_nM"],
     dt_20200313[dt_20200313$titration==14,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200313[dt_20200313$titration==14,"genotype"])),yaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200313[dt_20200313$titration==14 & dt_20200313$conc_nM==0, "mean_bin"][1])
axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200313_14,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200313_14)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200313_14)$coefficients["Kd","Std. Error"],digits=1),"nM")))

#titration==15, fit_20200313_15
plot(dt_20200313[dt_20200313$titration==15,"conc_nM"],
     dt_20200313[dt_20200313$titration==15,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200313[dt_20200313$titration==15,"genotype"])),yaxt='n',ylim=c(ymin,ymax),xlim=c(0.000005,4000))
points(0.000005, dt_20200313[dt_20200313$titration==15 & dt_20200313$conc_nM==0, "mean_bin"][1])
axis(side=2,labels=F)
lines(seq(0.0001,3200,0.001),predict(fit_20200313_15,list(conc_nM=seq(0.0001,3200,0.001))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("Kd",format(summary(fit_20200313_15)$coefficients["Kd","Estimate"],digits=2),
                      "+/-",format(summary(fit_20200313_15)$coefficients["Kd","Std. Error"],digits=1),"nM")))
axis(side=2,labels=F)
mtext("ligand concentration (nM, log scale)", side=1, outer=TRUE, line=2)
mtext("mean bin of FITC+ cells", side=2, outer=TRUE, line=2)
dev.off()
############
#summary: calculate mean and SEM ddG for each variant from the three experiments
triplicate <- data.frame(genotype=unique(dt_20200304$genotype), Kd_1=NA, Kd_2=NA, Kd_3=NA)
triplicate[triplicate$genotype=="NIH45-46",c("Kd_1","Kd_2","Kd_3")] <- c(summary(fit_20200304_1)$coefficients["Kd","Estimate"],summary(fit_20200311_1)$coefficients["Kd","Estimate"],summary(fit_20200313_1)$coefficients["Kd","Estimate"])
triplicate[triplicate$genotype=="G54hW",c("Kd_1","Kd_2","Kd_3")] <- c(summary(fit_20200304_2)$coefficients["Kd","Estimate"],summary(fit_20200311_2)$coefficients["Kd","Estimate"],summary(fit_20200313_2)$coefficients["Kd","Estimate"])
triplicate[triplicate$genotype=="P33hY",c("Kd_1","Kd_2","Kd_3")] <- c(summary(fit_20200304_3)$coefficients["Kd","Estimate"],summary(fit_20200311_3)$coefficients["Kd","Estimate"],summary(fit_20200313_3)$coefficients["Kd","Estimate"])
triplicate[triplicate$genotype=="N35hF",c("Kd_1","Kd_2","Kd_3")] <- c(summary(fit_20200304_4)$coefficients["Kd","Estimate"],summary(fit_20200311_8)$coefficients["Kd","Estimate"],summary(fit_20200313_8)$coefficients["Kd","Estimate"])        
triplicate[triplicate$genotype=="N35hH",c("Kd_1","Kd_2","Kd_3")] <- c(summary(fit_20200304_5)$coefficients["Kd","Estimate"],summary(fit_20200311_9)$coefficients["Kd","Estimate"],summary(fit_20200313_9)$coefficients["Kd","Estimate"])        
triplicate[triplicate$genotype=="I37hV",c("Kd_1","Kd_2","Kd_3")] <- c(summary(fit_20200304_6)$coefficients["Kd","Estimate"],summary(fit_20200311_4)$coefficients["Kd","Estimate"],summary(fit_20200313_4)$coefficients["Kd","Estimate"])        
triplicate[triplicate$genotype=="I37hW",c("Kd_1","Kd_2","Kd_3")] <- c(summary(fit_20200304_7)$coefficients["Kd","Estimate"],summary(fit_20200311_5)$coefficients["Kd","Estimate"],summary(fit_20200313_5)$coefficients["Kd","Estimate"])        
triplicate[triplicate$genotype=="N100EhG",c("Kd_1","Kd_2","Kd_3")] <- c(summary(fit_20200304_8)$coefficients["Kd","Estimate"],summary(fit_20200311_10)$coefficients["Kd","Estimate"],summary(fit_20200313_10)$coefficients["Kd","Estimate"])        
triplicate[triplicate$genotype=="N100EhS",c("Kd_1","Kd_2","Kd_3")] <- c(summary(fit_20200304_9)$coefficients["Kd","Estimate"],summary(fit_20200311_6)$coefficients["Kd","Estimate"],summary(fit_20200313_6)$coefficients["Kd","Estimate"])
triplicate[triplicate$genotype=="N100EhY",c("Kd_1","Kd_2","Kd_3")] <- c(summary(fit_20200304_10)$coefficients["Kd","Estimate"],summary(fit_20200311_11)$coefficients["Kd","Estimate"],summary(fit_20200313_11)$coefficients["Kd","Estimate"])        
triplicate[triplicate$genotype=="F98lV",c("Kd_1","Kd_2","Kd_3")] <- c(summary(fit_20200304_11)$coefficients["Kd","Estimate"],summary(fit_20200311_7)$coefficients["Kd","Estimate"],summary(fit_20200313_7)$coefficients["Kd","Estimate"])        
triplicate[triplicate$genotype=="N35hH/N100EhS",c("Kd_1","Kd_2","Kd_3")] <- c(summary(fit_20200304_12)$coefficients["Kd","Estimate"],summary(fit_20200311_12)$coefficients["Kd","Estimate"],summary(fit_20200313_12)$coefficients["Kd","Estimate"])        
triplicate[triplicate$genotype=="N35hH/N100EhG",c("Kd_1","Kd_2","Kd_3")] <- c(summary(fit_20200304_13)$coefficients["Kd","Estimate"],summary(fit_20200311_13)$coefficients["Kd","Estimate"],summary(fit_20200313_13)$coefficients["Kd","Estimate"])        
triplicate[triplicate$genotype=="N35hF/N100EhS",c("Kd_1","Kd_2","Kd_3")] <- c(summary(fit_20200304_14)$coefficients["Kd","Estimate"],summary(fit_20200311_14)$coefficients["Kd","Estimate"],summary(fit_20200313_14)$coefficients["Kd","Estimate"])
triplicate[triplicate$genotype=="I37hW/F98lV",c("Kd_1","Kd_2","Kd_3")] <- c(summary(fit_20200304_15)$coefficients["Kd","Estimate"],summary(fit_20200311_15)$coefficients["Kd","Estimate"],summary(fit_20200313_15)$coefficients["Kd","Estimate"])

triplicate$dG_1 <- (1.987204*10^-3)*(296.15)*log(triplicate$Kd_1*1e-9)
triplicate$dG_2 <- (1.987204*10^-3)*(296.15)*log(triplicate$Kd_2*1e-9)
triplicate$dG_3 <- (1.987204*10^-3)*(296.15)*log(triplicate$Kd_3*1e-9)

triplicate$ddG_1 <- triplicate$dG_1-triplicate[triplicate$genotype=="NIH45-46","dG_1"]
triplicate$ddG_2 <- triplicate$dG_2-triplicate[triplicate$genotype=="NIH45-46","dG_2"]
triplicate$ddG_3 <- triplicate$dG_3-triplicate[triplicate$genotype=="NIH45-46","dG_3"]

plot(triplicate$ddG_1,triplicate$ddG_2,pch=16)
plot(triplicate$ddG_1,triplicate$ddG_3,pch=16)
plot(triplicate$ddG_2,triplicate$ddG_3,pch=16)

for(i in 1:nrow(triplicate)){
        triplicate$mean_ddG[i] <- mean(c(triplicate$ddG_1[i],triplicate$ddG_2[i],triplicate$ddG_3[i]))
        triplicate$SEM_ddG[i] <- sd(c(triplicate$ddG_1[i],triplicate$ddG_2[i],triplicate$ddG_3[i]))/sqrt(3)
}

#for double mutants, compare actual to additive ddG, make cycle plots
triplicate$additive_ddG <- NA
triplicate[triplicate$genotype=="N35hH/N100EhS","additive_ddG"] <- triplicate[triplicate$genotype=="N35hH","mean_ddG"]+triplicate[triplicate$genotype=="N100EhS","mean_ddG"]
triplicate[triplicate$genotype=="N35hH/N100EhG","additive_ddG"] <- triplicate[triplicate$genotype=="N35hH","mean_ddG"]+triplicate[triplicate$genotype=="N100EhG","mean_ddG"]
triplicate[triplicate$genotype=="N35hF/N100EhS","additive_ddG"] <- triplicate[triplicate$genotype=="N35hF","mean_ddG"]+triplicate[triplicate$genotype=="N100EhS","mean_ddG"]
triplicate[triplicate$genotype=="I37hW/F98lV","additive_ddG"] <- triplicate[triplicate$genotype=="I37hW","mean_ddG"]+triplicate[triplicate$genotype=="F98lV","mean_ddG"]

triplicate$additive_ddG_SEM <- NA
triplicate[triplicate$genotype=="N35hH/N100EhS","additive_ddG_SEM"] <- sqrt(triplicate[triplicate$genotype=="N35hH","SEM_ddG"]^2+triplicate[triplicate$genotype=="N100EhS","SEM_ddG"]^2)
triplicate[triplicate$genotype=="N35hH/N100EhG","additive_ddG_SEM"] <- sqrt(triplicate[triplicate$genotype=="N35hH","SEM_ddG"]^2+triplicate[triplicate$genotype=="N100EhG","SEM_ddG"]^2)
triplicate[triplicate$genotype=="N35hF/N100EhS","additive_ddG_SEM"] <- sqrt(triplicate[triplicate$genotype=="N35hF","SEM_ddG"]^2+triplicate[triplicate$genotype=="N100EhS","SEM_ddG"]^2)
triplicate[triplicate$genotype=="I37hW/F98lV","additive_ddG_SEM"] <- sqrt(triplicate[triplicate$genotype=="I37hW","SEM_ddG"]^2+triplicate[triplicate$genotype=="F98lV","SEM_ddG"]^2)

pdf(file="results/isogenic_titrations/double-mutant-cycles_ddG.pdf",width=6,height=7,useDingbats=F)
#options(repr.plot.width=6, repr.plot.height=10, scipen=1)
par(mfrow=c(2,2))
#N35hH/N100EhS cycle
plot(NA,NA,xlim=c(-0.25,2.25),ylim=c(6,-0.5),ylab="ddG (kcal/mol)",xlab="number mutations",main="N35hH/N100EhS")
points(0,triplicate[triplicate$genotype=="NIH45-46","mean_ddG"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="N35hH","mean_ddG"]-1.96*triplicate[triplicate$genotype=="N35hH","SEM_ddG"],x1=1,y1=triplicate[triplicate$genotype=="N35hH","mean_ddG"]+1.96*triplicate[triplicate$genotype=="N35hH","SEM_ddG"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="N35hH","mean_ddG"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="N100EhS","mean_ddG"]-1.96*triplicate[triplicate$genotype=="N100EhS","SEM_ddG"],x1=1,y1=triplicate[triplicate$genotype=="N100EhS","mean_ddG"]+1.96*triplicate[triplicate$genotype=="N100EhS","SEM_ddG"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="N100EhS","mean_ddG"],pch=16)
arrows(x0=2,y0=triplicate[triplicate$genotype=="N35hH/N100EhS","mean_ddG"]-1.96*triplicate[triplicate$genotype=="N35hH/N100EhS","SEM_ddG"],x1=2,y1=triplicate[triplicate$genotype=="N35hH/N100EhS","mean_ddG"]+1.96*triplicate[triplicate$genotype=="N35hH/N100EhS","SEM_ddG"],angle=90,code=3,length=0.05)
points(2,triplicate[triplicate$genotype=="N35hH/N100EhS","mean_ddG"],pch=16)
arrows(x0=2,y0=triplicate[triplicate$genotype=="N35hH/N100EhS","additive_ddG"]-1.96*triplicate[triplicate$genotype=="N35hH/N100EhS","additive_ddG_SEM"],x1=2,y1=triplicate[triplicate$genotype=="N35hH/N100EhS","additive_ddG"]+1.96*triplicate[triplicate$genotype=="N35hH/N100EhS","additive_ddG_SEM"],angle=90,code=3,length=0.05,col="gray80")
points(2,triplicate[triplicate$genotype=="N35hH/N100EhS","additive_ddG"],pch=16,col="gray80")

#N35hH/N100EhG cycle
plot(NA,NA,xlim=c(-0.25,2.25),ylim=c(6,-0.5),ylab="ddG (kcal/mol)",xlab="number mutations",main="N35hH/N100EhG")
points(0,triplicate[triplicate$genotype=="NIH45-46","mean_ddG"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="N35hH","mean_ddG"]-1.96*triplicate[triplicate$genotype=="N35hH","SEM_ddG"],x1=1,y1=triplicate[triplicate$genotype=="N35hH","mean_ddG"]+1.96*triplicate[triplicate$genotype=="N35hH","SEM_ddG"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="N35hH","mean_ddG"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="N100EhG","mean_ddG"]-1.96*triplicate[triplicate$genotype=="N100EhG","SEM_ddG"],x1=1,y1=triplicate[triplicate$genotype=="N100EhG","mean_ddG"]+1.96*triplicate[triplicate$genotype=="N100EhG","SEM_ddG"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="N100EhG","mean_ddG"],pch=16)
arrows(x0=2,y0=triplicate[triplicate$genotype=="N35hH/N100EhG","mean_ddG"]-1.96*triplicate[triplicate$genotype=="N35hH/N100EhG","SEM_ddG"],x1=2,y1=triplicate[triplicate$genotype=="N35hH/N100EhG","mean_ddG"]+1.96*triplicate[triplicate$genotype=="N35hH/N100EhG","SEM_ddG"],angle=90,code=3,length=0.05)
points(2,triplicate[triplicate$genotype=="N35hH/N100EhG","mean_ddG"],pch=16)
arrows(x0=2,y0=triplicate[triplicate$genotype=="N35hH/N100EhG","additive_ddG"]-1.96*triplicate[triplicate$genotype=="N35hH/N100EhG","additive_ddG_SEM"],x1=2,y1=triplicate[triplicate$genotype=="N35hH/N100EhG","additive_ddG"]+1.96*triplicate[triplicate$genotype=="N35hH/N100EhG","additive_ddG_SEM"],angle=90,code=3,length=0.05,col="gray80")
points(2,triplicate[triplicate$genotype=="N35hH/N100EhG","additive_ddG"],pch=16,col="gray80")

#N35hF/N100EhS cycle
plot(NA,NA,xlim=c(-0.25,2.25),ylim=c(6,-0.5),ylab="ddG (kcal/mol)",xlab="number mutations",main="N35hF/N100EhS")
points(0,triplicate[triplicate$genotype=="NIH45-46","mean_ddG"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="N35hF","mean_ddG"]-1.96*triplicate[triplicate$genotype=="N35hF","SEM_ddG"],x1=1,y1=triplicate[triplicate$genotype=="N35hF","mean_ddG"]+1.96*triplicate[triplicate$genotype=="N35hF","SEM_ddG"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="N35hF","mean_ddG"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="N100EhS","mean_ddG"]-1.96*triplicate[triplicate$genotype=="N100EhS","SEM_ddG"],x1=1,y1=triplicate[triplicate$genotype=="N100EhS","mean_ddG"]+1.96*triplicate[triplicate$genotype=="N100EhS","SEM_ddG"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="N100EhS","mean_ddG"],pch=16)
arrows(x0=2,y0=triplicate[triplicate$genotype=="N35hF/N100EhS","mean_ddG"]-1.96*triplicate[triplicate$genotype=="N35hF/N100EhS","SEM_ddG"],x1=2,y1=triplicate[triplicate$genotype=="N35hF/N100EhS","mean_ddG"]+1.96*triplicate[triplicate$genotype=="N35hF/N100EhS","SEM_ddG"],angle=90,code=3,length=0.05)
points(2,triplicate[triplicate$genotype=="N35hF/N100EhS","mean_ddG"],pch=16)
arrows(x0=2,y0=triplicate[triplicate$genotype=="N35hF/N100EhS","additive_ddG"]-1.96*triplicate[triplicate$genotype=="N35hF/N100EhS","additive_ddG_SEM"],x1=2,y1=triplicate[triplicate$genotype=="N35hF/N100EhS","additive_ddG"]+1.96*triplicate[triplicate$genotype=="N35hF/N100EhS","additive_ddG_SEM"],angle=90,code=3,length=0.05,col="gray80")
points(2,triplicate[triplicate$genotype=="N35hF/N100EhS","additive_ddG"],pch=16,col="gray80")

#I37hW/F98lV cycle
plot(NA,NA,xlim=c(-0.25,2.25),ylim=c(3,-0.5),ylab="ddG (kcal/mol)",xlab="number mutations",main="I37hW/F98lV")
points(0,triplicate[triplicate$genotype=="NIH45-46","mean_ddG"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="I37hW","mean_ddG"]-1.96*triplicate[triplicate$genotype=="I37hW","SEM_ddG"],x1=1,y1=triplicate[triplicate$genotype=="I37hW","mean_ddG"]+1.96*triplicate[triplicate$genotype=="I37hW","SEM_ddG"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="I37hW","mean_ddG"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="F98lV","mean_ddG"]-1.96*triplicate[triplicate$genotype=="F98lV","SEM_ddG"],x1=1,y1=triplicate[triplicate$genotype=="F98lV","mean_ddG"]+1.96*triplicate[triplicate$genotype=="F98lV","SEM_ddG"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="F98lV","mean_ddG"],pch=16)
arrows(x0=2,y0=triplicate[triplicate$genotype=="I37hW/F98lV","mean_ddG"]-1.96*triplicate[triplicate$genotype=="I37hW/F98lV","SEM_ddG"],x1=2,y1=triplicate[triplicate$genotype=="I37hW/F98lV","mean_ddG"]+1.96*triplicate[triplicate$genotype=="I37hW/F98lV","SEM_ddG"],angle=90,code=3,length=0.05)
points(2,triplicate[triplicate$genotype=="I37hW/F98lV","mean_ddG"],pch=16)
arrows(x0=2,y0=triplicate[triplicate$genotype=="I37hW/F98lV","additive_ddG"]-1.96*triplicate[triplicate$genotype=="I37hW/F98lV","additive_ddG_SEM"],x1=2,y1=triplicate[triplicate$genotype=="I37hW/F98lV","additive_ddG"]+1.96*triplicate[triplicate$genotype=="I37hW/F98lV","additive_ddG_SEM"],angle=90,code=3,length=0.05,col="gray80")
points(2,triplicate[triplicate$genotype=="I37hW/F98lV","additive_ddG"],pch=16,col="gray80")
dev.off()

#look at expression phenotype (and, whether there is epistasis in the mutational effects on expression in the double mutants?)
triplicate[,c("FITCpos_1","FITCpos_2","FITCpos_3","meanFITC_1","meanFITC_2","meanFITC_3")] <- NA
triplicate[triplicate$genotype=="NIH45-46",c("FITCpos_1","FITCpos_2","FITCpos_3")] <- c(mean(dt_20200304[dt_20200304$titration==1,"FITCpos"]),mean(dt_20200311[dt_20200311$titration==1,"FITCpos"]),mean(dt_20200313[dt_20200313$titration==1,"FITCpos"]))
triplicate[triplicate$genotype=="NIH45-46",c("meanFITC_1","meanFITC_2","meanFITC_3")] <- c(mean(log(dt_20200304[dt_20200304$titration==1,"geomean_FITC"])),mean(log(dt_20200311[dt_20200311$titration==1,"geomean_FITC"])),mean(log(dt_20200313[dt_20200313$titration==1,"geomean_FITC"])))        

triplicate[triplicate$genotype=="G54hW",c("FITCpos_1","FITCpos_2","FITCpos_3")] <- c(mean(dt_20200304[dt_20200304$titration==2,"FITCpos"]),mean(dt_20200311[dt_20200311$titration==2,"FITCpos"]),mean(dt_20200313[dt_20200313$titration==2,"FITCpos"]))
triplicate[triplicate$genotype=="G54hW",c("meanFITC_1","meanFITC_2","meanFITC_3")] <- c(mean(log(dt_20200304[dt_20200304$titration==2,"geomean_FITC"])),mean(log(dt_20200311[dt_20200311$titration==2,"geomean_FITC"])),mean(log(dt_20200313[dt_20200313$titration==2,"geomean_FITC"])))        

triplicate[triplicate$genotype=="P33hY",c("FITCpos_1","FITCpos_2","FITCpos_3")] <- c(mean(dt_20200304[dt_20200304$titration==3,"FITCpos"]),mean(dt_20200311[dt_20200311$titration==3,"FITCpos"]),mean(dt_20200313[dt_20200313$titration==3,"FITCpos"]))
triplicate[triplicate$genotype=="P33hY",c("meanFITC_1","meanFITC_2","meanFITC_3")] <- c(mean(log(dt_20200304[dt_20200304$titration==3,"geomean_FITC"])),mean(log(dt_20200311[dt_20200311$titration==3,"geomean_FITC"])),mean(log(dt_20200313[dt_20200313$titration==3,"geomean_FITC"])))        

triplicate[triplicate$genotype=="N35hF",c("FITCpos_1","FITCpos_2","FITCpos_3")] <- c(mean(dt_20200304[dt_20200304$titration==4,"FITCpos"]),mean(dt_20200311[dt_20200311$titration==8,"FITCpos"]),mean(dt_20200313[dt_20200313$titration==8,"FITCpos"]))
triplicate[triplicate$genotype=="N35hF",c("meanFITC_1","meanFITC_2","meanFITC_3")] <- c(mean(log(dt_20200304[dt_20200304$titration==4,"geomean_FITC"])),mean(log(dt_20200311[dt_20200311$titration==8,"geomean_FITC"])),mean(log(dt_20200313[dt_20200313$titration==8,"geomean_FITC"])))        

triplicate[triplicate$genotype=="N35hH",c("FITCpos_1","FITCpos_2","FITCpos_3")] <- c(mean(dt_20200304[dt_20200304$titration==5,"FITCpos"]),mean(dt_20200311[dt_20200311$titration==9,"FITCpos"]),mean(dt_20200313[dt_20200313$titration==9,"FITCpos"]))
triplicate[triplicate$genotype=="N35hH",c("meanFITC_1","meanFITC_2","meanFITC_3")] <- c(mean(log(dt_20200304[dt_20200304$titration==5,"geomean_FITC"])),mean(log(dt_20200311[dt_20200311$titration==9,"geomean_FITC"])),mean(log(dt_20200313[dt_20200313$titration==9,"geomean_FITC"])))        

triplicate[triplicate$genotype=="I37hV",c("FITCpos_1","FITCpos_2","FITCpos_3")] <- c(mean(dt_20200304[dt_20200304$titration==6,"FITCpos"]),mean(dt_20200311[dt_20200311$titration==4,"FITCpos"]),mean(dt_20200313[dt_20200313$titration==4,"FITCpos"]))
triplicate[triplicate$genotype=="I37hV",c("meanFITC_1","meanFITC_2","meanFITC_3")] <- c(mean(log(dt_20200304[dt_20200304$titration==6,"geomean_FITC"])),mean(log(dt_20200311[dt_20200311$titration==4,"geomean_FITC"])),mean(log(dt_20200313[dt_20200313$titration==4,"geomean_FITC"])))

triplicate[triplicate$genotype=="I37hW",c("FITCpos_1","FITCpos_2","FITCpos_3")] <- c(mean(dt_20200304[dt_20200304$titration==7,"FITCpos"]),mean(dt_20200311[dt_20200311$titration==5,"FITCpos"]),mean(dt_20200313[dt_20200313$titration==5,"FITCpos"]))
triplicate[triplicate$genotype=="I37hW",c("meanFITC_1","meanFITC_2","meanFITC_3")] <- c(mean(log(dt_20200304[dt_20200304$titration==7,"geomean_FITC"])),mean(log(dt_20200311[dt_20200311$titration==5,"geomean_FITC"])),mean(log(dt_20200313[dt_20200313$titration==5,"geomean_FITC"])))

triplicate[triplicate$genotype=="N100EhG",c("FITCpos_1","FITCpos_2","FITCpos_3")] <- c(mean(dt_20200304[dt_20200304$titration==8,"FITCpos"]),mean(dt_20200311[dt_20200311$titration==10,"FITCpos"]),mean(dt_20200313[dt_20200313$titration==10,"FITCpos"]))
triplicate[triplicate$genotype=="N100EhG",c("meanFITC_1","meanFITC_2","meanFITC_3")] <- c(mean(log(dt_20200304[dt_20200304$titration==8,"geomean_FITC"])),mean(log(dt_20200311[dt_20200311$titration==10,"geomean_FITC"])),mean(log(dt_20200313[dt_20200313$titration==10,"geomean_FITC"])))

triplicate[triplicate$genotype=="N100EhS",c("FITCpos_1","FITCpos_2","FITCpos_3")] <- c(mean(dt_20200304[dt_20200304$titration==9,"FITCpos"]),mean(dt_20200311[dt_20200311$titration==6,"FITCpos"]),mean(dt_20200313[dt_20200313$titration==6,"FITCpos"]))
triplicate[triplicate$genotype=="N100EhS",c("meanFITC_1","meanFITC_2","meanFITC_3")] <- c(mean(log(dt_20200304[dt_20200304$titration==9,"geomean_FITC"])),mean(log(dt_20200311[dt_20200311$titration==6,"geomean_FITC"])),mean(log(dt_20200313[dt_20200313$titration==6,"geomean_FITC"])))

triplicate[triplicate$genotype=="N100EhY",c("FITCpos_1","FITCpos_2","FITCpos_3")] <- c(mean(dt_20200304[dt_20200304$titration==10,"FITCpos"]),mean(dt_20200311[dt_20200311$titration==11,"FITCpos"]),mean(dt_20200313[dt_20200313$titration==11,"FITCpos"]))
triplicate[triplicate$genotype=="N100EhY",c("meanFITC_1","meanFITC_2","meanFITC_3")] <- c(mean(log(dt_20200304[dt_20200304$titration==10,"geomean_FITC"])),mean(log(dt_20200311[dt_20200311$titration==11,"geomean_FITC"])),mean(log(dt_20200313[dt_20200313$titration==11,"geomean_FITC"])))

triplicate[triplicate$genotype=="F98lV",c("FITCpos_1","FITCpos_2","FITCpos_3")] <- c(mean(dt_20200304[dt_20200304$titration==11,"FITCpos"]),mean(dt_20200311[dt_20200311$titration==7,"FITCpos"]),mean(dt_20200313[dt_20200313$titration==7,"FITCpos"]))
triplicate[triplicate$genotype=="F98lV",c("meanFITC_1","meanFITC_2","meanFITC_3")] <- c(mean(log(dt_20200304[dt_20200304$titration==11,"geomean_FITC"])),mean(log(dt_20200311[dt_20200311$titration==7,"geomean_FITC"])),mean(log(dt_20200313[dt_20200313$titration==7,"geomean_FITC"])))

triplicate[triplicate$genotype=="N35hH/N100EhS",c("FITCpos_1","FITCpos_2","FITCpos_3")] <- c(mean(dt_20200304[dt_20200304$titration==12,"FITCpos"]),mean(dt_20200311[dt_20200311$titration==12,"FITCpos"]),mean(dt_20200313[dt_20200313$titration==12,"FITCpos"]))
triplicate[triplicate$genotype=="N35hH/N100EhS",c("meanFITC_1","meanFITC_2","meanFITC_3")] <- c(mean(log(dt_20200304[dt_20200304$titration==12,"geomean_FITC"])),mean(log(dt_20200311[dt_20200311$titration==12,"geomean_FITC"])),mean(log(dt_20200313[dt_20200313$titration==12,"geomean_FITC"])))

triplicate[triplicate$genotype=="N35hH/N100EhG",c("FITCpos_1","FITCpos_2","FITCpos_3")] <- c(mean(dt_20200304[dt_20200304$titration==13,"FITCpos"]),mean(dt_20200311[dt_20200311$titration==13,"FITCpos"]),mean(dt_20200313[dt_20200313$titration==13,"FITCpos"]))
triplicate[triplicate$genotype=="N35hH/N100EhG",c("meanFITC_1","meanFITC_2","meanFITC_3")] <- c(mean(log(dt_20200304[dt_20200304$titration==13,"geomean_FITC"])),mean(log(dt_20200311[dt_20200311$titration==13,"geomean_FITC"])),mean(log(dt_20200313[dt_20200313$titration==13,"geomean_FITC"])))

triplicate[triplicate$genotype=="N35hF/N100EhS",c("FITCpos_1","FITCpos_2","FITCpos_3")] <- c(mean(dt_20200304[dt_20200304$titration==14,"FITCpos"]),mean(dt_20200311[dt_20200311$titration==14,"FITCpos"]),mean(dt_20200313[dt_20200313$titration==14,"FITCpos"]))
triplicate[triplicate$genotype=="N35hF/N100EhS",c("meanFITC_1","meanFITC_2","meanFITC_3")] <- c(mean(log(dt_20200304[dt_20200304$titration==14,"geomean_FITC"])),mean(log(dt_20200311[dt_20200311$titration==14,"geomean_FITC"])),mean(log(dt_20200313[dt_20200313$titration==14,"geomean_FITC"])))

triplicate[triplicate$genotype=="I37hW/F98lV",c("FITCpos_1","FITCpos_2","FITCpos_3")] <- c(mean(dt_20200304[dt_20200304$titration==15,"FITCpos"]),mean(dt_20200311[dt_20200311$titration==15,"FITCpos"]),mean(dt_20200313[dt_20200313$titration==15,"FITCpos"]))
triplicate[triplicate$genotype=="I37hW/F98lV",c("meanFITC_1","meanFITC_2","meanFITC_3")] <- c(mean(log(dt_20200304[dt_20200304$titration==15,"geomean_FITC"])),mean(log(dt_20200311[dt_20200311$titration==15,"geomean_FITC"])),mean(log(dt_20200313[dt_20200313$titration==15,"geomean_FITC"])))

for(i in 1:nrow(triplicate)){
        triplicate$mean_FITCpos[i] <- mean(c(triplicate$FITCpos_1[i],triplicate$FITCpos_2[i],triplicate$FITCpos_3[i]))
        triplicate$SEM_FITCpos[i] <- sd(c(triplicate$FITCpos_1[i],triplicate$FITCpos_2[i],triplicate$FITCpos_3[i]))/sqrt(3)
        triplicate$mean_FITC[i] <- mean(c(triplicate$meanFITC_1[i],triplicate$meanFITC_2[i],triplicate$meanFITC_3[i]))
        triplicate$SEM_FITC[i] <- sd(c(triplicate$meanFITC_1[i],triplicate$meanFITC_2[i],triplicate$meanFITC_3[i]))/sqrt(3)
}

pdf(file="results/isogenic_titrations/double-mutant-cycles_FITCpos.pdf",width=6,height=7,useDingbats=F)
#options(repr.plot.width=6, repr.plot.height=10, scipen=1)
par(mfrow=c(2,2))
#N35hH/N100EhS cycle
plot(NA,NA,xlim=c(-0.25,2.25),ylim=c(0,100),ylab="% FITC-positive",xlab="number mutations",main="N35hH/N100EhS")
arrows(x0=0,y0=triplicate[triplicate$genotype=="NIH45-46","mean_FITCpos"]-1.96*triplicate[triplicate$genotype=="NIH45-46","SEM_FITCpos"],x1=0,y1=triplicate[triplicate$genotype=="NIH45-46","mean_FITCpos"]+1.96*triplicate[triplicate$genotype=="NIH45-46","SEM_FITCpos"],angle=90,code=3,length=0.05)
points(0,triplicate[triplicate$genotype=="NIH45-46","mean_FITCpos"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="N35hH","mean_FITCpos"]-1.96*triplicate[triplicate$genotype=="N35hH","SEM_FITCpos"],x1=1,y1=triplicate[triplicate$genotype=="N35hH","mean_FITCpos"]+1.96*triplicate[triplicate$genotype=="N35hH","SEM_FITCpos"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="N35hH","mean_FITCpos"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="N100EhS","mean_FITCpos"]-1.96*triplicate[triplicate$genotype=="N100EhS","SEM_FITCpos"],x1=1,y1=triplicate[triplicate$genotype=="N100EhS","mean_FITCpos"]+1.96*triplicate[triplicate$genotype=="N100EhS","SEM_FITCpos"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="N100EhS","mean_FITCpos"],pch=16)
arrows(x0=2,y0=triplicate[triplicate$genotype=="N35hH/N100EhS","mean_FITCpos"]-1.96*triplicate[triplicate$genotype=="N35hH/N100EhS","SEM_FITCpos"],x1=2,y1=triplicate[triplicate$genotype=="N35hH/N100EhS","mean_FITCpos"]+1.96*triplicate[triplicate$genotype=="N35hH/N100EhS","SEM_FITCpos"],angle=90,code=3,length=0.05)
points(2,triplicate[triplicate$genotype=="N35hH/N100EhS","mean_FITCpos"],pch=16)
# arrows(x0=2,y0=triplicate[triplicate$genotype=="N35hH/N100EhS","additive_FITCpos"]-1.96*triplicate[triplicate$genotype=="N35hH/N100EhS","additive_FITCpos_SEM"],x1=2,y1=triplicate[triplicate$genotype=="N35hH/N100EhS","additive_FITCpos"]+1.96*triplicate[triplicate$genotype=="N35hH/N100EhS","additive_FITCpos_SEM"],angle=90,code=3,length=0.05,col="gray80")
# points(2,triplicate[triplicate$genotype=="N35hH/N100EhS","additive_FITCpos"],pch=16,col="gray80")

#N35hH/N100EhG cycle
plot(NA,NA,xlim=c(-0.25,2.25),ylim=c(0,100),ylab="% FITC-positive",xlab="number mutations",main="N35hH/N100EhG")
arrows(x0=0,y0=triplicate[triplicate$genotype=="NIH45-46","mean_FITCpos"]-1.96*triplicate[triplicate$genotype=="NIH45-46","SEM_FITCpos"],x1=0,y1=triplicate[triplicate$genotype=="NIH45-46","mean_FITCpos"]+1.96*triplicate[triplicate$genotype=="NIH45-46","SEM_FITCpos"],angle=90,code=3,length=0.05)
points(0,triplicate[triplicate$genotype=="NIH45-46","mean_FITCpos"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="N35hH","mean_FITCpos"]-1.96*triplicate[triplicate$genotype=="N35hH","SEM_FITCpos"],x1=1,y1=triplicate[triplicate$genotype=="N35hH","mean_FITCpos"]+1.96*triplicate[triplicate$genotype=="N35hH","SEM_FITCpos"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="N35hH","mean_FITCpos"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="N100EhG","mean_FITCpos"]-1.96*triplicate[triplicate$genotype=="N100EhG","SEM_FITCpos"],x1=1,y1=triplicate[triplicate$genotype=="N100EhG","mean_FITCpos"]+1.96*triplicate[triplicate$genotype=="N100EhG","SEM_FITCpos"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="N100EhG","mean_FITCpos"],pch=16)
arrows(x0=2,y0=triplicate[triplicate$genotype=="N35hH/N100EhG","mean_FITCpos"]-1.96*triplicate[triplicate$genotype=="N35hH/N100EhG","SEM_FITCpos"],x1=2,y1=triplicate[triplicate$genotype=="N35hH/N100EhG","mean_FITCpos"]+1.96*triplicate[triplicate$genotype=="N35hH/N100EhG","SEM_FITCpos"],angle=90,code=3,length=0.05)
points(2,triplicate[triplicate$genotype=="N35hH/N100EhG","mean_FITCpos"],pch=16)
# arrows(x0=2,y0=triplicate[triplicate$genotype=="N35hH/N100EhG","additive_FITCpos"]-1.96*triplicate[triplicate$genotype=="N35hH/N100EhG","additive_FITCpos_SEM"],x1=2,y1=triplicate[triplicate$genotype=="N35hH/N100EhG","additive_FITCpos"]+1.96*triplicate[triplicate$genotype=="N35hH/N100EhG","additive_FITCpos_SEM"],angle=90,code=3,length=0.05,col="gray80")
# points(2,triplicate[triplicate$genotype=="N35hH/N100EhG","additive_FITCpos"],pch=16,col="gray80")

#N35hF/N100EhS cycle
plot(NA,NA,xlim=c(-0.25,2.25),ylim=c(0,100),ylab="% FITC-positive",xlab="number mutations",main="N35hF/N100EhS")
arrows(x0=0,y0=triplicate[triplicate$genotype=="NIH45-46","mean_FITCpos"]-1.96*triplicate[triplicate$genotype=="NIH45-46","SEM_FITCpos"],x1=0,y1=triplicate[triplicate$genotype=="NIH45-46","mean_FITCpos"]+1.96*triplicate[triplicate$genotype=="NIH45-46","SEM_FITCpos"],angle=90,code=3,length=0.05)
points(0,triplicate[triplicate$genotype=="NIH45-46","mean_FITCpos"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="N35hF","mean_FITCpos"]-1.96*triplicate[triplicate$genotype=="N35hF","SEM_FITCpos"],x1=1,y1=triplicate[triplicate$genotype=="N35hF","mean_FITCpos"]+1.96*triplicate[triplicate$genotype=="N35hF","SEM_FITCpos"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="N35hF","mean_FITCpos"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="N100EhS","mean_FITCpos"]-1.96*triplicate[triplicate$genotype=="N100EhS","SEM_FITCpos"],x1=1,y1=triplicate[triplicate$genotype=="N100EhS","mean_FITCpos"]+1.96*triplicate[triplicate$genotype=="N100EhS","SEM_FITCpos"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="N100EhS","mean_FITCpos"],pch=16)
arrows(x0=2,y0=triplicate[triplicate$genotype=="N35hF/N100EhS","mean_FITCpos"]-1.96*triplicate[triplicate$genotype=="N35hF/N100EhS","SEM_FITCpos"],x1=2,y1=triplicate[triplicate$genotype=="N35hF/N100EhS","mean_FITCpos"]+1.96*triplicate[triplicate$genotype=="N35hF/N100EhS","SEM_FITCpos"],angle=90,code=3,length=0.05)
points(2,triplicate[triplicate$genotype=="N35hF/N100EhS","mean_FITCpos"],pch=16)
# arrows(x0=2,y0=triplicate[triplicate$genotype=="N35hF/N100EhS","additive_FITCpos"]-1.96*triplicate[triplicate$genotype=="N35hF/N100EhS","additive_FITCpos_SEM"],x1=2,y1=triplicate[triplicate$genotype=="N35hF/N100EhS","additive_FITCpos"]+1.96*triplicate[triplicate$genotype=="N35hF/N100EhS","additive_FITCpos_SEM"],angle=90,code=3,length=0.05,col="gray80")
# points(2,triplicate[triplicate$genotype=="N35hF/N100EhS","additive_FITCpos"],pch=16,col="gray80")

#I37hW/F98lV cycle
plot(NA,NA,xlim=c(-0.25,2.25),ylim=c(0,100),ylab="% FITC-positive",xlab="number mutations",main="I37hW/F98lV")
arrows(x0=0,y0=triplicate[triplicate$genotype=="NIH45-46","mean_FITCpos"]-1.96*triplicate[triplicate$genotype=="NIH45-46","SEM_FITCpos"],x1=0,y1=triplicate[triplicate$genotype=="NIH45-46","mean_FITCpos"]+1.96*triplicate[triplicate$genotype=="NIH45-46","SEM_FITCpos"],angle=90,code=3,length=0.05)
points(0,triplicate[triplicate$genotype=="NIH45-46","mean_FITCpos"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="I37hW","mean_FITCpos"]-1.96*triplicate[triplicate$genotype=="I37hW","SEM_FITCpos"],x1=1,y1=triplicate[triplicate$genotype=="I37hW","mean_FITCpos"]+1.96*triplicate[triplicate$genotype=="I37hW","SEM_FITCpos"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="I37hW","mean_FITCpos"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="F98lV","mean_FITCpos"]-1.96*triplicate[triplicate$genotype=="F98lV","SEM_FITCpos"],x1=1,y1=triplicate[triplicate$genotype=="F98lV","mean_FITCpos"]+1.96*triplicate[triplicate$genotype=="F98lV","SEM_FITCpos"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="F98lV","mean_FITCpos"],pch=16)
arrows(x0=2,y0=triplicate[triplicate$genotype=="I37hW/F98lV","mean_FITCpos"]-1.96*triplicate[triplicate$genotype=="I37hW/F98lV","SEM_FITCpos"],x1=2,y1=triplicate[triplicate$genotype=="I37hW/F98lV","mean_FITCpos"]+1.96*triplicate[triplicate$genotype=="I37hW/F98lV","SEM_FITCpos"],angle=90,code=3,length=0.05)
points(2,triplicate[triplicate$genotype=="I37hW/F98lV","mean_FITCpos"],pch=16)
# arrows(x0=2,y0=triplicate[triplicate$genotype=="I37hW/F98lV","additive_FITCpos"]-1.96*triplicate[triplicate$genotype=="I37hW/F98lV","additive_FITCpos_SEM"],x1=2,y1=triplicate[triplicate$genotype=="I37hW/F98lV","additive_FITCpos"]+1.96*triplicate[triplicate$genotype=="I37hW/F98lV","additive_FITCpos_SEM"],angle=90,code=3,length=0.05,col="gray80")
# points(2,triplicate[triplicate$genotype=="I37hW/F98lV","additive_FITCpos"],pch=16,col="gray80")
dev.off()

pdf(file="results/isogenic_titrations/double-mutant-cycles_FITC-meanF.pdf",width=6,height=7,useDingbats=F)
#options(repr.plot.width=6, repr.plot.height=10, scipen=1)
par(mfrow=c(2,2))
#N35hH/N100EhS cycle
plot(NA,NA,xlim=c(-0.25,2.25),ylim=c(7,8.5),ylab="FITC mean fluorescence (a.u.)",xlab="number mutations",main="N35hH/N100EhS")
arrows(x0=0,y0=triplicate[triplicate$genotype=="NIH45-46","mean_FITC"]-1.96*triplicate[triplicate$genotype=="NIH45-46","SEM_FITC"],x1=0,y1=triplicate[triplicate$genotype=="NIH45-46","mean_FITC"]+1.96*triplicate[triplicate$genotype=="NIH45-46","SEM_FITC"],angle=90,code=3,length=0.05)
points(0,triplicate[triplicate$genotype=="NIH45-46","mean_FITC"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="N35hH","mean_FITC"]-1.96*triplicate[triplicate$genotype=="N35hH","SEM_FITC"],x1=1,y1=triplicate[triplicate$genotype=="N35hH","mean_FITC"]+1.96*triplicate[triplicate$genotype=="N35hH","SEM_FITC"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="N35hH","mean_FITC"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="N100EhS","mean_FITC"]-1.96*triplicate[triplicate$genotype=="N100EhS","SEM_FITC"],x1=1,y1=triplicate[triplicate$genotype=="N100EhS","mean_FITC"]+1.96*triplicate[triplicate$genotype=="N100EhS","SEM_FITC"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="N100EhS","mean_FITC"],pch=16)
arrows(x0=2,y0=triplicate[triplicate$genotype=="N35hH/N100EhS","mean_FITC"]-1.96*triplicate[triplicate$genotype=="N35hH/N100EhS","SEM_FITC"],x1=2,y1=triplicate[triplicate$genotype=="N35hH/N100EhS","mean_FITC"]+1.96*triplicate[triplicate$genotype=="N35hH/N100EhS","SEM_FITC"],angle=90,code=3,length=0.05)
points(2,triplicate[triplicate$genotype=="N35hH/N100EhS","mean_FITC"],pch=16)
# arrows(x0=2,y0=triplicate[triplicate$genotype=="N35hH/N100EhS","additive_FITC"]-1.96*triplicate[triplicate$genotype=="N35hH/N100EhS","additive_FITC_SEM"],x1=2,y1=triplicate[triplicate$genotype=="N35hH/N100EhS","additive_FITC"]+1.96*triplicate[triplicate$genotype=="N35hH/N100EhS","additive_FITC_SEM"],angle=90,code=3,length=0.05,col="gray80")
# points(2,triplicate[triplicate$genotype=="N35hH/N100EhS","additive_FITC"],pch=16,col="gray80")

#N35hH/N100EhG cycle
plot(NA,NA,xlim=c(-0.25,2.25),ylim=c(7,8.5),ylab="FITC mean fluorescence (a.u.)",xlab="number mutations",main="N35hH/N100EhG")
arrows(x0=0,y0=triplicate[triplicate$genotype=="NIH45-46","mean_FITC"]-1.96*triplicate[triplicate$genotype=="NIH45-46","SEM_FITC"],x1=0,y1=triplicate[triplicate$genotype=="NIH45-46","mean_FITC"]+1.96*triplicate[triplicate$genotype=="NIH45-46","SEM_FITC"],angle=90,code=3,length=0.05)
points(0,triplicate[triplicate$genotype=="NIH45-46","mean_FITC"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="N35hH","mean_FITC"]-1.96*triplicate[triplicate$genotype=="N35hH","SEM_FITC"],x1=1,y1=triplicate[triplicate$genotype=="N35hH","mean_FITC"]+1.96*triplicate[triplicate$genotype=="N35hH","SEM_FITC"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="N35hH","mean_FITC"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="N100EhG","mean_FITC"]-1.96*triplicate[triplicate$genotype=="N100EhG","SEM_FITC"],x1=1,y1=triplicate[triplicate$genotype=="N100EhG","mean_FITC"]+1.96*triplicate[triplicate$genotype=="N100EhG","SEM_FITC"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="N100EhG","mean_FITC"],pch=16)
arrows(x0=2,y0=triplicate[triplicate$genotype=="N35hH/N100EhG","mean_FITC"]-1.96*triplicate[triplicate$genotype=="N35hH/N100EhG","SEM_FITC"],x1=2,y1=triplicate[triplicate$genotype=="N35hH/N100EhG","mean_FITC"]+1.96*triplicate[triplicate$genotype=="N35hH/N100EhG","SEM_FITC"],angle=90,code=3,length=0.05)
points(2,triplicate[triplicate$genotype=="N35hH/N100EhG","mean_FITC"],pch=16)
# arrows(x0=2,y0=triplicate[triplicate$genotype=="N35hH/N100EhG","additive_FITC"]-1.96*triplicate[triplicate$genotype=="N35hH/N100EhG","additive_FITC_SEM"],x1=2,y1=triplicate[triplicate$genotype=="N35hH/N100EhG","additive_FITC"]+1.96*triplicate[triplicate$genotype=="N35hH/N100EhG","additive_FITC_SEM"],angle=90,code=3,length=0.05,col="gray80")
# points(2,triplicate[triplicate$genotype=="N35hH/N100EhG","additive_FITC"],pch=16,col="gray80")

#N35hF/N100EhS cycle
plot(NA,NA,xlim=c(-0.25,2.25),ylim=c(7,8.5),ylab="FITC mean fluorescence (a.u.)",xlab="number mutations",main="N35hF/N100EhS")
arrows(x0=0,y0=triplicate[triplicate$genotype=="NIH45-46","mean_FITC"]-1.96*triplicate[triplicate$genotype=="NIH45-46","SEM_FITC"],x1=0,y1=triplicate[triplicate$genotype=="NIH45-46","mean_FITC"]+1.96*triplicate[triplicate$genotype=="NIH45-46","SEM_FITC"],angle=90,code=3,length=0.05)
points(0,triplicate[triplicate$genotype=="NIH45-46","mean_FITC"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="N35hF","mean_FITC"]-1.96*triplicate[triplicate$genotype=="N35hF","SEM_FITC"],x1=1,y1=triplicate[triplicate$genotype=="N35hF","mean_FITC"]+1.96*triplicate[triplicate$genotype=="N35hF","SEM_FITC"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="N35hF","mean_FITC"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="N100EhS","mean_FITC"]-1.96*triplicate[triplicate$genotype=="N100EhS","SEM_FITC"],x1=1,y1=triplicate[triplicate$genotype=="N100EhS","mean_FITC"]+1.96*triplicate[triplicate$genotype=="N100EhS","SEM_FITC"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="N100EhS","mean_FITC"],pch=16)
arrows(x0=2,y0=triplicate[triplicate$genotype=="N35hF/N100EhS","mean_FITC"]-1.96*triplicate[triplicate$genotype=="N35hF/N100EhS","SEM_FITC"],x1=2,y1=triplicate[triplicate$genotype=="N35hF/N100EhS","mean_FITC"]+1.96*triplicate[triplicate$genotype=="N35hF/N100EhS","SEM_FITC"],angle=90,code=3,length=0.05)
points(2,triplicate[triplicate$genotype=="N35hF/N100EhS","mean_FITC"],pch=16)
# arrows(x0=2,y0=triplicate[triplicate$genotype=="N35hF/N100EhS","additive_FITC"]-1.96*triplicate[triplicate$genotype=="N35hF/N100EhS","additive_FITC_SEM"],x1=2,y1=triplicate[triplicate$genotype=="N35hF/N100EhS","additive_FITC"]+1.96*triplicate[triplicate$genotype=="N35hF/N100EhS","additive_FITC_SEM"],angle=90,code=3,length=0.05,col="gray80")
# points(2,triplicate[triplicate$genotype=="N35hF/N100EhS","additive_FITC"],pch=16,col="gray80")

#I37hW/F98lV cycle
plot(NA,NA,xlim=c(-0.25,2.25),ylim=c(7,8.5),ylab="FITC mean fluorescence (a.u.)",xlab="number mutations",main="I37hW/F98lV")
arrows(x0=0,y0=triplicate[triplicate$genotype=="NIH45-46","mean_FITC"]-1.96*triplicate[triplicate$genotype=="NIH45-46","SEM_FITC"],x1=0,y1=triplicate[triplicate$genotype=="NIH45-46","mean_FITC"]+1.96*triplicate[triplicate$genotype=="NIH45-46","SEM_FITC"],angle=90,code=3,length=0.05)
points(0,triplicate[triplicate$genotype=="NIH45-46","mean_FITC"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="I37hW","mean_FITC"]-1.96*triplicate[triplicate$genotype=="I37hW","SEM_FITC"],x1=1,y1=triplicate[triplicate$genotype=="I37hW","mean_FITC"]+1.96*triplicate[triplicate$genotype=="I37hW","SEM_FITC"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="I37hW","mean_FITC"],pch=16)
arrows(x0=1,y0=triplicate[triplicate$genotype=="F98lV","mean_FITC"]-1.96*triplicate[triplicate$genotype=="F98lV","SEM_FITC"],x1=1,y1=triplicate[triplicate$genotype=="F98lV","mean_FITC"]+1.96*triplicate[triplicate$genotype=="F98lV","SEM_FITC"],angle=90,code=3,length=0.05)
points(1,triplicate[triplicate$genotype=="F98lV","mean_FITC"],pch=16)
arrows(x0=2,y0=triplicate[triplicate$genotype=="I37hW/F98lV","mean_FITC"]-1.96*triplicate[triplicate$genotype=="I37hW/F98lV","SEM_FITC"],x1=2,y1=triplicate[triplicate$genotype=="I37hW/F98lV","mean_FITC"]+1.96*triplicate[triplicate$genotype=="I37hW/F98lV","SEM_FITC"],angle=90,code=3,length=0.05)
points(2,triplicate[triplicate$genotype=="I37hW/F98lV","mean_FITC"],pch=16)
# arrows(x0=2,y0=triplicate[triplicate$genotype=="I37hW/F98lV","additive_FITC"]-1.96*triplicate[triplicate$genotype=="I37hW/F98lV","additive_FITC_SEM"],x1=2,y1=triplicate[triplicate$genotype=="I37hW/F98lV","additive_FITC"]+1.96*triplicate[triplicate$genotype=="I37hW/F98lV","additive_FITC_SEM"],angle=90,code=3,length=0.05,col="gray80")
# points(2,triplicate[triplicate$genotype=="I37hW/F98lV","additive_FITC"],pch=16,col="gray80")
dev.off()


#compare to expression and binding measurements from the bulk Tite-seq experiment!
betas <- data.table(read.csv(file="results/analyze_mut_effects/betas.csv",stringsAsFactors = F))

triplicate[,c("bind_latent","bind_observed","bind_observed_A","bind_observed_B","expr_latent","expr_latent_A","expr_latent_B","expr_observed")] <- NA
triplicate[triplicate$genotype=="G54hW",c("bind_latent","bind_observed","bind_observed_A","bind_observed_B","expr_latent","expr_latent_A","expr_latent_B","expr_observed")] <- betas[mutation=="G55W",c("bind_latent","bind_observed","bind_observed_A","bind_observed_B","expr_latent","expr_latent_A","expr_latent_B","expr_observed")]
triplicate[triplicate$genotype=="P33hY",c("bind_latent","bind_observed","bind_observed_A","bind_observed_B","expr_latent","expr_latent_A","expr_latent_B","expr_observed")] <- betas[mutation=="P33Y",c("bind_latent","bind_observed","bind_observed_A","bind_observed_B","expr_latent","expr_latent_A","expr_latent_B","expr_observed")]
triplicate[triplicate$genotype=="N35hF",c("bind_latent","bind_observed","bind_observed_A","bind_observed_B","expr_latent","expr_latent_A","expr_latent_B","expr_observed")] <- betas[mutation=="N35F",c("bind_latent","bind_observed","bind_observed_A","bind_observed_B","expr_latent","expr_latent_A","expr_latent_B","expr_observed")]
triplicate[triplicate$genotype=="N35hH",c("bind_latent","bind_observed","bind_observed_A","bind_observed_B","expr_latent","expr_latent_A","expr_latent_B","expr_observed")] <- betas[mutation=="N35H",c("bind_latent","bind_observed","bind_observed_A","bind_observed_B","expr_latent","expr_latent_A","expr_latent_B","expr_observed")]
triplicate[triplicate$genotype=="I37hV",c("bind_latent","bind_observed","bind_observed_A","bind_observed_B","expr_latent","expr_latent_A","expr_latent_B","expr_observed")] <- betas[mutation=="I37V",c("bind_latent","bind_observed","bind_observed_A","bind_observed_B","expr_latent","expr_latent_A","expr_latent_B","expr_observed")]
triplicate[triplicate$genotype=="I37hW",c("bind_latent","bind_observed","bind_observed_A","bind_observed_B","expr_latent","expr_latent_A","expr_latent_B","expr_observed")] <- betas[mutation=="I37W",c("bind_latent","bind_observed","bind_observed_A","bind_observed_B","expr_latent","expr_latent_A","expr_latent_B","expr_observed")]
triplicate[triplicate$genotype=="N100EhG",c("bind_latent","bind_observed","bind_observed_A","bind_observed_B","expr_latent","expr_latent_A","expr_latent_B","expr_observed")] <- betas[mutation=="N109G",c("bind_latent","bind_observed","bind_observed_A","bind_observed_B","expr_latent","expr_latent_A","expr_latent_B","expr_observed")]
triplicate[triplicate$genotype=="N100EhS",c("bind_latent","bind_observed","bind_observed_A","bind_observed_B","expr_latent","expr_latent_A","expr_latent_B","expr_observed")] <- betas[mutation=="N109S",c("bind_latent","bind_observed","bind_observed_A","bind_observed_B","expr_latent","expr_latent_A","expr_latent_B","expr_observed")]
triplicate[triplicate$genotype=="N100EhY",c("bind_latent","bind_observed","bind_observed_A","bind_observed_B","expr_latent","expr_latent_A","expr_latent_B","expr_observed")] <- betas[mutation=="N109Y",c("bind_latent","bind_observed","bind_observed_A","bind_observed_B","expr_latent","expr_latent_A","expr_latent_B","expr_observed")]
triplicate[triplicate$genotype=="F98lV",c("bind_latent","bind_observed","bind_observed_A","bind_observed_B","expr_latent","expr_latent_A","expr_latent_B","expr_observed")] <- betas[mutation=="F232V",c("bind_latent","bind_observed","bind_observed_A","bind_observed_B","expr_latent","expr_latent_A","expr_latent_B","expr_observed")]

triplicate$delta_FITC <- triplicate$mean_FITC-triplicate[triplicate$genotype=="NIH45-46","mean_FITC"]

pdf(file="results/isogenic_titrations/Titeseq-isogenic_ddG_correlation.pdf",width=3.5,height=4,useDingbats=F,bg="transparent")
plot(triplicate$mean_ddG, triplicate$bind_observed,pch=16,xlab="isogenic titration, ddG (kcal/mol)",ylab="bulk Tite-seq, ddG (kcal/mol)",ylim=c(-0.5,6))
arrows(x0=triplicate$mean_ddG, y0=triplicate$bind_observed_A,x1=triplicate$mean_ddG, y1=triplicate$bind_observed_B,angle=90,code=3,length=0.04)
arrows(x0=triplicate$mean_ddG-triplicate$SEM_ddG, y0=triplicate$bind_observed,x1=triplicate$mean_ddG+triplicate$SEM_ddG, y1=triplicate$bind_observed,angle=90,code=3,length=0.04)
model <- lm(triplicate$bind_observed~triplicate$mean_ddG);abline(model,lty=2,col="red");summary(model)
legend("topleft",legend=paste("R-squared =",round(summary(model)$r.squared,digits=3),"\nslope =",round(summary(model)$coefficients[2,"Estimate"],digits=2)),bty="n")
#abline(a=0,b=1)
dev.off()

plot(triplicate$mean_ddG, triplicate$bind_latent,pch=16,xlab="isogenic titration, ddG (kcal/mol)",ylab="bulk Tite-seq, ddG (kcal/mol, latent scale)")
abline(a=0,b=1)

pdf(file="results/isogenic_titrations/Sortseq-isogenic_dE_correlation.pdf",width=3.5,height=4,useDingbats=F,bg="transparent")
plot(triplicate$delta_FITC, triplicate$expr_latent,pch=16,xlab="isogenic titration, dE (a.u.)",ylab="bulk Sort-seq, dE (a.u.)",ylim=c(-0.7,0.3),xlim=c(-0.7,0.3))
arrows(x0=triplicate$delta_FITC, y0=triplicate$expr_latent_A,x1=triplicate$delta_FITC, y1=triplicate$expr_latent_B,angle=90,code=3,length=0.04)
arrows(x0=triplicate$delta_FITC-triplicate$SEM_FITC, y0=triplicate$expr_latent,x1=triplicate$delta_FITC+triplicate$SEM_FITC, y1=triplicate$expr_latent,angle=90,code=3,length=0.04)
model <- lm(triplicate$expr_latent~triplicate$delta_FITC);abline(model,lty=2,col="red");summary(model)
legend("topleft",legend=paste("R-squared =",round(summary(model)$r.squared,digits=3),"\nslope =",round(summary(model)$coefficients[2,"Estimate"],digits=2)),bty="n")
#abline(a=0,b=1)
dev.off()

plot(triplicate$delta_FITC, triplicate$expr_observed,pch=16,xlab="isogenic titration, dE (a.u.)",ylab="bulk Sort-seq, dE (a.u.)",ylim=c(-1.5,0.3),xlim=c(-0.7,0.3))
model <- lm(triplicate$expr_observed~triplicate$delta_FITC);abline(model,lty=2,col="red");summary(model)
legend("topleft",legend=paste("R-squared =",round(summary(model)$r.squared,digits=3),"\nslope =",round(summary(model)$coefficients[2,"Estimate"],digits=2)),bty="n")
#abline(a=0,b=1)
