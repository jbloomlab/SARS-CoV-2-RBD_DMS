#6 May 2020
#TNS

#script to read in results of flow cytometry of isogenic mutants (for various purposes), fit titration curves and some interpretation from each experiment

#library(ggplot2)
#library(data.table)

###################################
#1 May 2020 experiment: triplicate measurements of the SARS-CoV-2 RBD binding to human ACE2
dt_20200501 <- read.csv(file="20200501.csv", stringsAsFactors=F)

#fit curves to each titration, and the three pooled measurements
fit_20200501_1 <- nls(mean_bin ~ a*(conc_M/(conc_M+Kd))+b,
                data=dt_20200501[dt_20200501$titration==1,],start=list(a=3,b=1,Kd=5*10^-11))
fit_20200501_2 <- nls(mean_bin ~ a*(conc_M/(conc_M+Kd))+b,
                data=dt_20200501[dt_20200501$titration==2,],start=list(a=3,b=1,Kd=5*10^-11))
fit_20200501_3 <- nls(mean_bin ~ a*(conc_M/(conc_M+Kd))+b,
                data=dt_20200501[dt_20200501$titration==3,],start=list(a=3,b=1,Kd=5*10^-11))
fit_20200501_pooled <- nls(mean_bin ~ a*(conc_M/(conc_M+Kd))+b,
                data=dt_20200501,start=list(a=3,b=1,Kd=5*10^-11))


#plot results
pdf(file="20200501_individual.pdf",width=3,height=8,useDingbats=F)
#plot three separate fits with their Kd
options(scipen=1)
par(mfrow=c(3,1),oma=c(4,4,0,0),mar=c(1,1,1,1))

ymax <- 4
ymin <- 1

#titration==1, fit_20200501_1
plot(dt_20200501[dt_20200501$titration==1,"conc_M"],
     dt_20200501[dt_20200501$titration==1,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200501[dt_20200501$titration==1,"genotype"])),xaxt='n',ylim=c(ymin,ymax),xlim=c(1e-14,1e-7))
points(1e-14, dt_20200501[dt_20200501$titration==1 & dt_20200501$conc_M==0, "mean_bin"][1])
axis(side=1,labels=F)
lines(seq(1e-13,1e-7,1e-14),predict(fit_20200501_1,list(conc_M=seq(1e-13,1e-7,1e-14))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("log10(Kd)",formatC(log10(summary(fit_20200501_1)$coefficients["Kd","Estimate"]),format='f',digits=2),
                      "+/-",formatC(0.434*summary(fit_20200501_1)$coefficients["Kd","Std. Error"]/summary(fit_20200501_1)$coefficients["Kd","Estimate"],format='f',digits=2))))

#titration==2, fit_20200501_2
plot(dt_20200501[dt_20200501$titration==2,"conc_M"],
     dt_20200501[dt_20200501$titration==2,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200501[dt_20200501$titration==2,"genotype"])),xaxt='n',ylim=c(ymin,ymax),xlim=c(1e-14,1e-7),ylab="mean bin")
points(1e-14, dt_20200501[dt_20200501$titration==2 & dt_20200501$conc_M==0, "mean_bin"][1])
axis(side=1,labels=F)
lines(seq(1e-13,1e-7,1e-14),predict(fit_20200501_2,list(conc_M=seq(1e-13,1e-7,1e-14))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("log10(Kd)",formatC(log10(summary(fit_20200501_2)$coefficients["Kd","Estimate"]),format='f',digits=2),
                      "+/-",formatC(0.434*summary(fit_20200501_2)$coefficients["Kd","Std. Error"]/summary(fit_20200501_2)$coefficients["Kd","Estimate"],format='f',digits=2))))

#titration==3, fit_20200501_3
plot(dt_20200501[dt_20200501$titration==3,"conc_M"],
     dt_20200501[dt_20200501$titration==3,"mean_bin"],
     pch=19,log="x",main=as.character(unique(dt_20200501[dt_20200501$titration==3,"genotype"])),ylim=c(ymin,ymax),xlim=c(1e-14,1e-7),ylab="mean bin")
points(1e-14, dt_20200501[dt_20200501$titration==3 & dt_20200501$conc_M==0, "mean_bin"][1])
lines(seq(1e-13,1e-7,1e-14),predict(fit_20200501_3,list(conc_M=seq(1e-13,1e-7,1e-14))),col="red")
legend("topleft",bty="n",cex=0.75,
       legend=c(paste("log10(Kd)",formatC(log10(summary(fit_20200501_3)$coefficients["Kd","Estimate"]),format='f',digits=2),
                      "+/-",formatC(0.434*summary(fit_20200501_3)$coefficients["Kd","Std. Error"]/summary(fit_20200501_3)$coefficients["Kd","Estimate"],format='f',digits=2))))

mtext("ligand concentration (M, log scale)", side=1, outer=TRUE, line=2)
mtext("mean bin of FITC+ cells", side=2, outer=TRUE, line=2)
dev.off()


