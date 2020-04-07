#4 November 2019
#TNS
#script to analyze sequencing counts from Tite-seq experiments to generate Kd estimates for each barcode

setwd("~/bloom_j/computational_notebooks/tstarr/2019/NIH45-46_DMS")

#import R libraries
library(yaml)
library(data.table)
library(tidyverse)
library(Hmisc)
#library(fitdistrplus)

#read the configuration file
config <- read_yaml("config.yaml")

#make output directory
if(!file.exists(config$Titeseq_Kds_dir)){
  dir.create(file.path(config$Titeseq_Kds_dir))
}

#read dataframe with list of barcode runs
barcode_runs <- read.csv(file=config$barcode_runs,stringsAsFactors=F); barcode_runs <- subset(barcode_runs, select=-c(R1))

#read file giving count of each barcode in each sample partition, sort by barcode
counts <- data.table(read.csv(file=config$variant_counts_file,stringsAsFactors=F)); counts <- counts[order(counts$library, counts$barcode),]

#eliminate rows from counts that are not part of a Tite-seq bin
counts <- subset(counts, paste(library,sample) %in% sapply(1:nrow(barcode_runs[!is.na(barcode_runs$Titeseq_bin),c("library","sample")]), function(x) paste(barcode_runs[!is.na(barcode_runs$Titeseq_bin),"library"][x],barcode_runs[!is.na(barcode_runs$Titeseq_bin),"sample"][x])) )

#if reads>cells for a sort bin, normalize counts by the read:cell ratio for each of the sort bins
for(i in as.numeric(row.names(barcode_runs[!is.na(barcode_runs$Titeseq_bin),]))){
  lib <- as.character(barcode_runs$library[i])
  bin <- as.character(barcode_runs$sample[i])
  if(is.na(barcode_runs$number_cells[i])){
    counts[library==lib & sample==bin,count.norm := as.numeric(count)] 
    print(paste("no cell count given for",lib,bin,", un-normalized"))
  }else{
    if(sum(counts[library==lib & sample==bin,"count"]) < barcode_runs$number_cells[i]){
      counts[library==lib & sample==bin,count.norm := as.numeric(count)]
      print(paste("reads < cells for",lib,bin,", un-normalized"))
    }else{
      #normalize by read:cell ratio
      ratio <- sum(counts[library==lib & sample==bin,"count"])/barcode_runs$number_cells[i]
      counts[library==lib & sample==bin,count.norm := as.numeric(count/ratio)]
    }
  }
}

#annotate each barcode as to whether it's wildtype, synonymous muts only, stop, nonsynonymous, >1 nonsynonymous mutations
counts[n_codon_substitutions==0, variant_class := "wildtype"]
counts[n_codon_substitutions > 0 & n_aa_substitutions==0, variant_class := "synonymous"]
counts[n_aa_substitutions>0 & grepl("*",aa_substitutions,fixed=T), variant_class := "stop"]
counts[n_aa_substitutions == 1 & !grepl("*",aa_substitutions,fixed=T), variant_class := "1 nonsynonymous"]
counts[n_aa_substitutions > 1 & !grepl("*",aa_substitutions,fixed=T), variant_class := ">1 nonsynonymous"]

#cast the counts data frame into wide format for A and B replicates
counts_A <- dcast(counts[library=="libA2",], barcode + variant_call_support + variant_class + aa_substitutions + n_aa_substitutions + n_codon_substitutions ~ library + sample, value.var="count.norm")
counts_B <- dcast(counts[library=="libB2",], barcode + variant_call_support + variant_class + aa_substitutions + n_aa_substitutions + n_codon_substitutions ~ library + sample, value.var="count.norm")

#create table giving names of samples and gp120 concentrations
bins_A <- data.frame(bin=unique(barcode_runs[barcode_runs$library=="libA2" & !is.na(barcode_runs$Titeseq_bin),"Titeseq_bin"]),conc=unique(barcode_runs[barcode_runs$library=="libA2" & !is.na(barcode_runs$Titeseq_bin),"Titeseq_sample_concentration"]))
bins_B <- data.frame(bin=unique(barcode_runs[barcode_runs$library=="libB2" & !is.na(barcode_runs$Titeseq_bin),"Titeseq_bin"]),conc=unique(barcode_runs[barcode_runs$library=="libB2" & !is.na(barcode_runs$Titeseq_bin),"Titeseq_sample_concentration"]))

#function that returns mean bin and sum of counts for four bins at a given concentration
calc.meanbin <- function(vec){return( list((vec[1]*1+vec[2]*2+vec[3]*3+vec[4]*4)/(vec[1]+vec[2]+vec[3]+vec[4]),(vec[1]+vec[2]+vec[3]+vec[4])) )}

#iterate through samples, compute mean.bin and total.count
for(i in 1:nrow(bins_A)){
  meanbin.out <- paste(bins_A[i,"bin"],"_meanbin",sep="")
  totalcount.out <- paste(bins_A[i,"bin"],"_totalcount",sep="")
  b1.in <- paste(bins_A[i,"bin"],"-b1",sep="")
  b2.in <- paste(bins_A[i,"bin"],"-b2",sep="")
  b3.in <- paste(bins_A[i,"bin"],"-b3",sep="")
  b4.in <- paste(bins_A[i,"bin"],"-b4",sep="")
  counts_A[,c(meanbin.out,totalcount.out) := calc.meanbin(c(get(b1.in), get(b2.in), get(b3.in), get(b4.in))),by=barcode]
}

for(i in 1:nrow(bins_B)){
  meanbin.out <- paste(bins_B[i,"bin"],"_meanbin",sep="")
  totalcount.out <- paste(bins_B[i,"bin"],"_totalcount",sep="")
  b1.in <- paste(bins_B[i,"bin"],"-b1",sep="")
  b2.in <- paste(bins_B[i,"bin"],"-b2",sep="")
  b3.in <- paste(bins_B[i,"bin"],"-b3",sep="")
  b4.in <- paste(bins_B[i,"bin"],"-b4",sep="")
  counts_B[,c(meanbin.out,totalcount.out) := calc.meanbin(c(get(b1.in), get(b2.in), get(b3.in), get(b4.in))),by=barcode]
}


# #correlation in total counts per barcode across the 16 different concentrations
# png(filename=paste(config$Titeseq_Kds_dir,"/pairwise-sample-count-correlations%d.png",sep=""),width=960,height=960)
# pairs(counts_A[1:10000,c("libA2_191010_s01_totalcount","libA2_191010_s02_totalcount","libA2_191010_s03_totalcount","libA2_191010_s04_totalcount","libA2_191010_s05_totalcount","libA2_191010_s06_totalcount","libA2_191010_s07_totalcount","libA2_191010_s08_totalcount","libA2_191010_s09_totalcount","libA2_191010_s10_totalcount","libA2_191010_s11_totalcount","libA2_191010_s12_totalcount","libA2_191010_s13_totalcount","libA2_191010_s14_totalcount","libA2_191010_s15_totalcount","libA2_191010_s16_totalcount")]+0.5,log="xy",pch=19,col="#00000020")
# pairs(counts_B[1:10000,c("libB2_200102_s01_totalcount","libB2_200102_s02_totalcount","libB2_200102_s03_totalcount","libB2_200102_s04_totalcount","libB2_200102_s05_totalcount","libB2_200102_s06_totalcount","libB2_200102_s07_totalcount","libB2_200102_s08_totalcount","libB2_200102_s09_totalcount","libB2_200102_s10_totalcount","libB2_200102_s11_totalcount","libB2_200102_s12_totalcount","libB2_200102_s13_totalcount","libB2_200102_s14_totalcount","libB2_200102_s15_totalcount","libB2_200102_s16_totalcount")]+0.5,log="xy",pch=19,col="#00000020")
# dev.off()
# 
# png(filename=paste(config$Titeseq_Kds_dir,"/pairwise-sample-bin-count-correlations%d.png",sep=""),width=3840,height=3840)
# pairs(counts_A[1:10000,c("libA2_191010_s01-b1","libA2_191010_s01-b2","libA2_191010_s01-b3","libA2_191010_s01-b4","libA2_191010_s02-b1","libA2_191010_s02-b2","libA2_191010_s02-b3","libA2_191010_s02-b4","libA2_191010_s03-b1","libA2_191010_s03-b2","libA2_191010_s03-b3","libA2_191010_s03-b4","libA2_191010_s04-b1","libA2_191010_s04-b2","libA2_191010_s04-b3","libA2_191010_s04-b4","libA2_191010_s05-b1","libA2_191010_s05-b2","libA2_191010_s05-b3","libA2_191010_s05-b4","libA2_191010_s06-b1","libA2_191010_s06-b2","libA2_191010_s06-b3","libA2_191010_s06-b4","libA2_191010_s07-b1","libA2_191010_s07-b2","libA2_191010_s07-b3","libA2_191010_s07-b4","libA2_191010_s08-b1","libA2_191010_s08-b2","libA2_191010_s08-b3","libA2_191010_s08-b4","libA2_191010_s09-b1","libA2_191010_s09-b2","libA2_191010_s09-b3","libA2_191010_s09-b4","libA2_191010_s10-b1","libA2_191010_s10-b2","libA2_191010_s10-b3","libA2_191010_s10-b4","libA2_191010_s11-b1","libA2_191010_s11-b2","libA2_191010_s11-b3","libA2_191010_s11-b4","libA2_191010_s12-b1","libA2_191010_s12-b2","libA2_191010_s12-b3","libA2_191010_s12-b4","libA2_191010_s13-b1","libA2_191010_s13-b2","libA2_191010_s13-b3","libA2_191010_s13-b4","libA2_191010_s14-b1","libA2_191010_s14-b2","libA2_191010_s14-b3","libA2_191010_s14-b4","libA2_191010_s15-b1","libA2_191010_s15-b2","libA2_191010_s15-b3","libA2_191010_s15-b4","libA2_191010_s16-b1","libA2_191010_s16-b2","libA2_191010_s16-b3","libA2_191010_s16-b4")]+0.5,log="xy",pch=19,col="#00000020")
# pairs(counts_B[1:10000,c("libB2_200102_s01-b1","libB2_200102_s01-b2","libB2_200102_s01-b3","libB2_200102_s01-b4","libB2_200102_s02-b1","libB2_200102_s02-b2","libB2_200102_s02-b3","libB2_200102_s02-b4","libB2_200102_s03-b1","libB2_200102_s03-b2","libB2_200102_s03-b3","libB2_200102_s03-b4","libB2_200102_s04-b1","libB2_200102_s04-b2","libB2_200102_s04-b3","libB2_200102_s04-b4","libB2_200102_s05-b1","libB2_200102_s05-b2","libB2_200102_s05-b3","libB2_200102_s05-b4","libB2_200102_s06-b1","libB2_200102_s06-b2","libB2_200102_s06-b3","libB2_200102_s06-b4","libB2_200102_s07-b1","libB2_200102_s07-b2","libB2_200102_s07-b3","libB2_200102_s07-b4","libB2_200102_s08-b1","libB2_200102_s08-b2","libB2_200102_s08-b3","libB2_200102_s08-b4","libB2_200102_s09-b1","libB2_200102_s09-b2","libB2_200102_s09-b3","libB2_200102_s09-b4","libB2_200102_s10-b1","libB2_200102_s10-b2","libB2_200102_s10-b3","libB2_200102_s10-b4","libB2_200102_s11-b1","libB2_200102_s11-b2","libB2_200102_s11-b3","libB2_200102_s11-b4","libB2_200102_s12-b1","libB2_200102_s12-b2","libB2_200102_s12-b3","libB2_200102_s12-b4","libB2_200102_s13-b1","libB2_200102_s13-b2","libB2_200102_s13-b3","libB2_200102_s13-b4","libB2_200102_s14-b1","libB2_200102_s14-b2","libB2_200102_s14-b3","libB2_200102_s14-b4","libB2_200102_s15-b1","libB2_200102_s15-b2","libB2_200102_s15-b3","libB2_200102_s15-b4","libB2_200102_s16-b1","libB2_200102_s16-b2","libB2_200102_s16-b3","libB2_200102_s16-b4")]+0.5,log="xy",pch=19,col="#00000020")
# dev.off()

#see libB2 has much stronger right tail of highly represented variants.... this is not good! causes the mean totalcount for the rest of the variants to go down by ~10fold relative to libA2 because so many of the reads went toward these guys! Why did this happen?!
hist(log10(counts_A$libA2_191010_s01_totalcount+0.5))
hist(log10(counts_B$libB2_200102_s01_totalcount+0.5),add=T,breaks=30,col="#7f7f7f50") #better plot for this below on average count across the 16 samples

#want to fit titration curve as weighted least squares -- need an estimate of variance. Will use repeated observations of wildtype across different count depths at the s09 concentration (~Kd) to generate these estimated variances?
wt_A <- counts_A[variant_class %in% c("synonymous","wildtype") ,]
hist(wt_A[libA2_191010_s09_totalcount<2,libA2_191010_s09_meanbin],xlim=c(1,4),breaks=20)
#make bins of wildtype observations
n.breaks_A <- 30
wt.bins_A <- data.frame(bin=1:n.breaks_A)
breaks_A <- cut2(wt_A$libA2_191010_s09_totalcount,m=250,g=n.breaks_A,onlycuts=T)#; breaks_A[length(breaks_A)] <- breaks_A[length(breaks_A)]+500000
#pool wt observations that fall within the bin breaks, compute sd
for(i in 1:nrow(wt.bins_A)){
  wt.bins_A$range.cells[i] <- list(c(breaks_A[i],breaks_A[i+1]))
  data <- wt_A[libA2_191010_s09_totalcount >= wt.bins_A$range.cells[i][[1]][[1]] & libA2_191010_s09_totalcount < wt.bins_A$range.cells[i][[1]][[2]],]
  wt.bins_A$mean.meanbin[i] <- mean(data$libA2_191010_s09_meanbin,na.rm=T)
  wt.bins_A$sd.meanbin[i] <- sd(data$libA2_191010_s09_meanbin,na.rm=T)
  wt.bins_A$median.cells[i] <- median(data$libA2_191010_s09_totalcount,na.rm=T)
}
#look at relationship between variance and cell counts; fit curve to estimate variance from cell count
y_A <- (wt.bins_A$sd.meanbin[-1])^2
x_A <- wt.bins_A$median.cells[-1]
plot(x_A,y_A,xlab="number cells",ylab="variance")
plot(log(x_A),log(y_A),xlab="log(number cells)",ylab="log(variance)")
wt.fit_A <- lm(log(y_A) ~ log(x_A));summary(wt.fit_A);abline(wt.fit_A)

#repeat replicate B
wt_B <- counts_B[variant_class %in% c("synonymous","wildtype") ,]
hist(wt_B[libB2_200102_s09_totalcount<2,libB2_200102_s09_meanbin],xlim=c(1,4))
#make bins of wildtype observations
n.breaks_B <- 30
wt.bins_B <- data.frame(bin=1:n.breaks_B)
breaks_B <- cut2(wt_B$libB2_200102_s09_totalcount,g=n.breaks_B,onlycuts=T)#; breaks_B[length(breaks_B)] <- breaks_B[length(breaks_B)]+500000
#pool wt observations that fall within the bin breaks, compute sd
for(i in 1:nrow(wt.bins_B)){
  wt.bins_B$range.cells[i] <- list(c(breaks_B[i],breaks_B[i+1]))
  data <- wt_B[libB2_200102_s09_totalcount >= wt.bins_B$range.cells[i][[1]][[1]] & libB2_200102_s09_totalcount < wt.bins_B$range.cells[i][[1]][[2]],]
  wt.bins_B$mean.meanbin[i] <- mean(data$libB2_200102_s09_meanbin,na.rm=T)
  wt.bins_B$sd.meanbin[i] <- sd(data$libB2_200102_s09_meanbin,na.rm=T)
  wt.bins_B$median.cells[i] <- median(data$libB2_200102_s09_totalcount,na.rm=T)
}
#look at relationship between variance and cell counts; fit curve to estimate variance from cell count
y_B <- (wt.bins_B$sd.meanbin[-1])^2
x_B <- wt.bins_B$median.cells[-1]
#compare to rep A plot
#pdf(file=paste(config$Titeseq_Kds_dir,"/wt-variance_v_number-cells_both-libs.pdf",sep=""),width=4,height=4.5,useDingbats=F,bg="transparent")
plot(x_A,y_A,xlab="number cells",ylab="variance",ylim=c(0,0.6),pch=19,col="blue")
points(x_B[!is.na(y_B) & y_B>0],y_B[!is.na(y_B) & y_B>0],xlab="number cells",ylab="variance",pch=19,col="green")
#dev.off()
plot(log(x_B),log(y_B),xlab="log(number cells)",ylab="log(variance)")
wt.fit_B <- lm(log(y_B[!is.na(y_B) & y_B>0]) ~ log(x_B[!is.na(y_B) & y_B>0]));summary(wt.fit_B);abline(wt.fit_B)

#function that approximates variance from cell count, and not allow var to be less than our observed var at the highest count bin
est.var_A <- function(x){
  var <- exp(as.numeric(wt.fit_A$coefficients[1]) + as.numeric(wt.fit_A$coefficients[2]) * log(x))
  if(var < min(wt.bins_A$sd.meanbin[-1]^2,na.rm=T)){return(min(wt.bins_A$sd.meanbin[-1]^2,na.rm=T))}else{return(var)}
}

est.var_B <- function(x){
  var <- exp(as.numeric(wt.fit_B$coefficients[1]) + as.numeric(wt.fit_B$coefficients[2]) * log(x))
  if(var < min(wt.bins_B$sd.meanbin[-c(1,2)]^2,na.rm=T)){return(min(wt.bins_B$sd.meanbin[-c(1,2)]^2,na.rm=T))}else{return(var)}
}


#also for doing QC, output columns giving the average number of cells, number of samples for which there were fewer than "cutoff" cells
counts_A[,libA2_191010_avgcount := mean(c(libA2_191010_s01_totalcount,libA2_191010_s02_totalcount,libA2_191010_s03_totalcount,libA2_191010_s04_totalcount,libA2_191010_s05_totalcount,libA2_191010_s06_totalcount,libA2_191010_s07_totalcount,libA2_191010_s08_totalcount,libA2_191010_s09_totalcount,libA2_191010_s10_totalcount,libA2_191010_s11_totalcount,libA2_191010_s12_totalcount,libA2_191010_s13_totalcount,libA2_191010_s14_totalcount,libA2_191010_s15_totalcount,libA2_191010_s16_totalcount)),by=barcode]
counts_A[,libA2_191010_min.cell.filtered := sum(c(libA2_191010_s01_totalcount,libA2_191010_s02_totalcount,libA2_191010_s03_totalcount,libA2_191010_s04_totalcount,libA2_191010_s05_totalcount,libA2_191010_s06_totalcount,libA2_191010_s07_totalcount,libA2_191010_s08_totalcount,libA2_191010_s09_totalcount,libA2_191010_s10_totalcount,libA2_191010_s11_totalcount,libA2_191010_s12_totalcount,libA2_191010_s13_totalcount,libA2_191010_s14_totalcount,libA2_191010_s15_totalcount,libA2_191010_s16_totalcount)<2),by=barcode]

counts_B[,libB2_200102_avgcount := mean(c(libB2_200102_s01_totalcount,libB2_200102_s02_totalcount,libB2_200102_s03_totalcount,libB2_200102_s04_totalcount,libB2_200102_s05_totalcount,libB2_200102_s06_totalcount,libB2_200102_s07_totalcount,libB2_200102_s08_totalcount,libB2_200102_s09_totalcount,libB2_200102_s10_totalcount,libB2_200102_s11_totalcount,libB2_200102_s12_totalcount,libB2_200102_s13_totalcount,libB2_200102_s14_totalcount,libB2_200102_s15_totalcount,libB2_200102_s16_totalcount)),by=barcode]
counts_B[,libB2_200102_min.cell.filtered := sum(c(libB2_200102_s01_totalcount,libB2_200102_s02_totalcount,libB2_200102_s03_totalcount,libB2_200102_s04_totalcount,libB2_200102_s05_totalcount,libB2_200102_s06_totalcount,libB2_200102_s07_totalcount,libB2_200102_s08_totalcount,libB2_200102_s09_totalcount,libB2_200102_s10_totalcount,libB2_200102_s11_totalcount,libB2_200102_s12_totalcount,libB2_200102_s13_totalcount,libB2_200102_s14_totalcount,libB2_200102_s15_totalcount,libB2_200102_s16_totalcount)<1),by=barcode]

#pdf(file=paste(config$Titeseq_Kds_dir,"/rank-ordered-cell-counts.pdf",sep=""),width=4,height=4.5,useDingbats=F,bg="transparent")
plot(counts_B$libB2_200102_avgcount[order(counts_B$libB2_200102_avgcount)],log="y",type="l",col="green",ylab="cell counts, Tite-seq sort (green = libB, blue=libA)",xlab="variant, rank-ordered",lwd=3)
points(counts_A$libA2_191010_avgcount[order(counts_A$libA2_191010_avgcount)],log="y",type="l",col="blue",lwd=3)
#dev.off()

#function that fits a nls regression to the titration series, including weights from the counts and an option to filter below certain thresholds
fit.titration_A <- function(y.vals,x.vals,count.vals,min.cfu=2,min.means=0.6,min.average=5,Kd.start=5e-11,a.start=3,a.lower=2,a.upper=3,b.start=1,b.lower=1,b.upper=1.5){
  indices <- count.vals>min.cfu
  y <- y.vals[indices]
  x <- x.vals[indices]
  w <- 1/sapply(count.vals[indices], est.var_A)
  if((length(y) < min.means*length(y.vals)) | (mean(count.vals[indices],na.rm=T) < min.average)){ #return NAs if >min.means % of concentrations have below min.cfu counts or if the average count across all used concentrations is below min.average
    return(list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA)))
  }else{
    fit <- nls(y ~ a*(x/(x+Kd))+b,
               start=list(a=a.start,b=b.start,Kd=Kd.start),
               lower=list(a=a.lower,b=b.lower,Kd=min(x.vals[x.vals>0])/100),
               upper=list(a=a.upper,b=b.upper,Kd=max(x.vals[x.vals>0])*100),
               weights=w,algorithm="port")  
    return(list(as.numeric(summary(fit)$coefficients["Kd","Estimate"]),as.numeric(summary(fit)$coefficients["Kd","Std. Error"]),as.numeric(summary(fit)$coefficients["a","Estimate"]),as.numeric(summary(fit)$coefficients["b","Estimate"]),as.numeric(summary(fit)$sigma),list(fit)))
  }
}

#fit titration to 191010 Titeseq data for each barcode
counts_A[,c("Kd","Kd_SE","response","baseline","RSE","fit") := tryCatch(fit.titration_A(y.vals=c(libA2_191010_s01_meanbin,
                                                                                                 libA2_191010_s02_meanbin,
                                                                                                 libA2_191010_s03_meanbin,
                                                                                                 libA2_191010_s04_meanbin,
                                                                                                 libA2_191010_s05_meanbin,
                                                                                                 libA2_191010_s06_meanbin,
                                                                                                 libA2_191010_s07_meanbin,
                                                                                                 libA2_191010_s08_meanbin,
                                                                                                 libA2_191010_s09_meanbin,
                                                                                                 libA2_191010_s10_meanbin,
                                                                                                 libA2_191010_s11_meanbin,
                                                                                                 libA2_191010_s12_meanbin,
                                                                                                 libA2_191010_s13_meanbin,
                                                                                                 libA2_191010_s14_meanbin,
                                                                                                 libA2_191010_s15_meanbin,
                                                                                                 libA2_191010_s16_meanbin),
                                                                                        x.vals=bins_A$conc,
                                                                                        count.vals=c(libA2_191010_s01_totalcount,
                                                                                                 libA2_191010_s02_totalcount,
                                                                                                 libA2_191010_s03_totalcount,
                                                                                                 libA2_191010_s04_totalcount,
                                                                                                 libA2_191010_s05_totalcount,
                                                                                                 libA2_191010_s06_totalcount,
                                                                                                 libA2_191010_s07_totalcount,
                                                                                                 libA2_191010_s08_totalcount,
                                                                                                 libA2_191010_s09_totalcount,
                                                                                                 libA2_191010_s10_totalcount,
                                                                                                 libA2_191010_s11_totalcount,
                                                                                                 libA2_191010_s12_totalcount,
                                                                                                 libA2_191010_s13_totalcount,
                                                                                                 libA2_191010_s14_totalcount,
                                                                                                 libA2_191010_s15_totalcount,
                                                                                                 libA2_191010_s16_totalcount)),
                                                                          error=function(e){return(list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA)))}),by=barcode]

#save R object with all curve fits for interactive work, since I'll be filtering the curve fits in the main data frame.
counts_A.all <- copy(counts_A)
save(counts_A,file=paste(config$Titeseq_Kds_dir,"/dt.temp_A.Rda",sep=""))
#load(file=paste(config$Titeseq_Kds_dir,"/dt.temp_A.Rda",sep=""))

#repeat replicate B
fit.titration_B <- function(y.vals,x.vals,count.vals,min.cfu=1,min.means=0.6,min.average=3,Kd.start=5e-11,a.start=3,a.lower=2,a.upper=3,b.start=1,b.lower=1,b.upper=1.5){
  indices <- count.vals>min.cfu
  y <- y.vals[indices]
  x <- x.vals[indices]
  w <- 1/sapply(count.vals[indices], est.var_B)
  if((length(y) < min.means*length(y.vals)) | (mean(count.vals[indices],na.rm=T) < min.average)){ #return NAs if >min.means % of concentrations have below min.cfu counts or if the average count across all used concentrations is below min.average
    return(list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA)))
  }else{
    fit <- nls(y ~ a*(x/(x+Kd))+b,
               start=list(a=a.start,b=b.start,Kd=Kd.start),
               lower=list(a=a.lower,b=b.lower,Kd=min(x.vals[x.vals>0])/100),
               upper=list(a=a.upper,b=b.upper,Kd=max(x.vals[x.vals>0])*100),
               weights=w,algorithm="port")  
    return(list(as.numeric(summary(fit)$coefficients["Kd","Estimate"]),as.numeric(summary(fit)$coefficients["Kd","Std. Error"]),as.numeric(summary(fit)$coefficients["a","Estimate"]),as.numeric(summary(fit)$coefficients["b","Estimate"]),as.numeric(summary(fit)$sigma),list(fit)))
  }
}

#fit titration to 200102 Titeseq data for each barcode
counts_B[,c("Kd","Kd_SE","response","baseline","RSE","fit") := tryCatch(fit.titration_B(y.vals=c(libB2_200102_s01_meanbin,
                                                                                                 libB2_200102_s02_meanbin,
                                                                                                 libB2_200102_s03_meanbin,
                                                                                                 libB2_200102_s04_meanbin,
                                                                                                 libB2_200102_s05_meanbin,
                                                                                                 libB2_200102_s06_meanbin,
                                                                                                 libB2_200102_s07_meanbin,
                                                                                                 libB2_200102_s08_meanbin,
                                                                                                 libB2_200102_s09_meanbin,
                                                                                                 libB2_200102_s10_meanbin,
                                                                                                 libB2_200102_s11_meanbin,
                                                                                                 libB2_200102_s12_meanbin,
                                                                                                 libB2_200102_s13_meanbin,
                                                                                                 libB2_200102_s14_meanbin,
                                                                                                 libB2_200102_s15_meanbin,
                                                                                                 libB2_200102_s16_meanbin),
                                                                                        x.vals=bins_B$conc,
                                                                                        count.vals=c(libB2_200102_s01_totalcount,
                                                                                                     libB2_200102_s02_totalcount,
                                                                                                     libB2_200102_s03_totalcount,
                                                                                                     libB2_200102_s04_totalcount,
                                                                                                     libB2_200102_s05_totalcount,
                                                                                                     libB2_200102_s06_totalcount,
                                                                                                     libB2_200102_s07_totalcount,
                                                                                                     libB2_200102_s08_totalcount,
                                                                                                     libB2_200102_s09_totalcount,
                                                                                                     libB2_200102_s10_totalcount,
                                                                                                     libB2_200102_s11_totalcount,
                                                                                                     libB2_200102_s12_totalcount,
                                                                                                     libB2_200102_s13_totalcount,
                                                                                                     libB2_200102_s14_totalcount,
                                                                                                     libB2_200102_s15_totalcount,
                                                                                                     libB2_200102_s16_totalcount)),
                                                                          error=function(e){return(list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA)))}),by=barcode]

#save R object with all curve fits for interactive work, since I'll be filtering the curve fits in the main data frame.
counts_B.all <- copy(counts_B)
save(counts_B,file=paste(config$Titeseq_Kds_dir,"/dt.temp_B.Rda",sep=""))
#load(file=paste(config$Titeseq_Kds_dir,"/dt.temp_B.Rda",sep=""))

#make function that allows me to plot a titration for any given row
plot.titration.A <- function(row){
  y.vals <- c();for(bin in bins_A$bin){y.vals <- c(y.vals,paste(bin,"_meanbin",sep=""))};y.vals <- unlist(counts_A[row,y.vals,with=F])
  x.vals <- bins_A$conc
  count.vals <- c();for(bin in bins_A$bin){count.vals <- c(count.vals,paste(bin,"_totalcount",sep=""))};count.vals <- unlist(counts_A[row,count.vals,with=F])
  fit <- counts_A[row,fit[[1]]]
  plot(x.vals[count.vals>2],y.vals[count.vals>2],xlab="[gp120]",ylab="mean bin",log="x",ylim=c(1,4),xlim=c(1e-13,1e-6),pch=19,main=counts_A[row,aa_substitutions])
  lines(x.vals,predict(fit,newdata=list(x=x.vals)))
  counts_A[row,.(barcode,variant_class,aa_substitutions,libA2_191010_avgcount,libA2_191010_min.cell.filtered,Kd,Kd_SE,baseline,response,RSE)]
}

plot.titration.B <- function(row){
  y.vals <- c();for(bin in bins_B$bin){y.vals <- c(y.vals,paste(bin,"_meanbin",sep=""))};y.vals <- unlist(counts_B[row,y.vals,with=F])
  x.vals <- bins_B$conc
  count.vals <- c();for(bin in bins_B$bin){count.vals <- c(count.vals,paste(bin,"_totalcount",sep=""))};count.vals <- unlist(counts_B[row,count.vals,with=F])
  fit <- counts_B[row,fit[[1]]]
  plot(x.vals[count.vals>1],y.vals[count.vals>1],xlab="[gp120]",ylab="mean bin",log="x",ylim=c(1,4),xlim=c(1e-13,1e-6),pch=19,main=counts_B[row,aa_substitutions])
  lines(x.vals,predict(fit,newdata=list(x=x.vals)))
  counts_B[row,.(barcode,variant_class,aa_substitutions,libB2_200102_avgcount,libB2_200102_min.cell.filtered,Kd,Kd_SE,baseline,response,RSE)]
}


############################################################################################
##QC/sanity checks for curve fits
sum(!is.na(counts_A$Kd))/nrow(counts_A) #have Kds for 74.9% of my barcodes in rep A (117662)
sum(!is.na(counts_B$Kd))/nrow(counts_B) #have Kds for 50.5% of my barcodes in rep B (90202)
#look at average counts, number missing for variants with successful versus unsuccessfully inferred Kd
hist(log10(counts_A[!is.na(Kd),libA2_191010_avgcount]+0.5),breaks=20,xlim=c(0,5))
hist(log10(counts_A[is.na(Kd),libA2_191010_avgcount]+0.5),add=T,col="red")

hist(counts_A[!is.na(Kd),libA2_191010_min.cell.filtered],breaks=17,xlim=c(0,16))
hist(counts_A[is.na(Kd),libA2_191010_min.cell.filtered],breaks=17,add=T,col="red")
#can see that the vast majority of the NAs come from my minimum read cutoffs, meaning there weren't too many curves that failed to be fit.
hist(log10(counts_B[!is.na(Kd),libB2_200102_avgcount]+0.5),breaks=20,xlim=c(0,5))
hist(log10(counts_B[is.na(Kd),libB2_200102_avgcount]+0.5),add=T,col="red")

hist(counts_B[!is.na(Kd),libB2_200102_min.cell.filtered],breaks=17,xlim=c(0,16))
hist(counts_B[is.na(Kd),libB2_200102_min.cell.filtered],breaks=17,add=T,col="red")

#looking through those with high average counts (e.g. >50) and yet NA Kds: seems like pretty much all curves with high Kds, that also either have low saturation point or just not saturated. They all have considerable noise at the higher concentrations so it's hard to tell whether it's plateaued at ~2 or not. These could be saved by allowing a to be smaller, but there's only arouhnd one hundred so I'm ok just losing them
which(counts_A$libA2_191010_avgcount > 50 & is.na(counts_A$Kd))
plot.titration.A(which(counts_A$libA2_191010_avgcount > 50 & is.na(counts_A$Kd))[5])

which(counts_B$libB2_200102_avgcount > 50 & is.na(counts_B$Kd))
plot.titration.B(which(counts_B$libB2_200102_avgcount > 50 & is.na(counts_B$Kd))[2])

# row <- 25050; test.fit <- fit.titration(y.vals=counts[row,c(libA2_191010_s01_meanbin,libA2_191010_s02_meanbin,libA2_191010_s03_meanbin,libA2_191010_s04_meanbin,libA2_191010_s05_meanbin,libA2_191010_s06_meanbin,libA2_191010_s07_meanbin,libA2_191010_s08_meanbin,libA2_191010_s09_meanbin,libA2_191010_s10_meanbin,libA2_191010_s11_meanbin,libA2_191010_s12_meanbin,libA2_191010_s13_meanbin,libA2_191010_s14_meanbin,libA2_191010_s15_meanbin,libA2_191010_s16_meanbin)],
#                          x.vals=bins$conc,
#                          count.vals=counts[row,c(libA2_191010_s01_totalcount,libA2_191010_s02_totalcount,libA2_191010_s03_totalcount,libA2_191010_s04_totalcount,libA2_191010_s05_totalcount,libA2_191010_s06_totalcount,libA2_191010_s07_totalcount,libA2_191010_s08_totalcount,libA2_191010_s09_totalcount,libA2_191010_s10_totalcount,libA2_191010_s11_totalcount,libA2_191010_s12_totalcount,libA2_191010_s13_totalcount,libA2_191010_s14_totalcount,libA2_191010_s15_totalcount,libA2_191010_s16_totalcount)],
#                          a.lower=1.5)
# I could fit more of these by reducing the minimum response amount, but I'm not sure I believe these plateaus
#this means to me that I don't have major classes of "failed" fit that I need to go back to and figure out. Can leave the NAs as NA and just spot check the successfully inferred curves

#next, look at Kds
hist(log10(counts_A$Kd),breaks=30)
hist(log10(counts_B$Kd),breaks=40,col="#7f7f7f50",add=T)


#censored at maximum Kd
plot.titration.A(which(counts_A$Kd==7.27e-05)[3]) #all flatlines
plot.titration.B(which(counts_B$Kd==7.14e-5)[13]) #all flatlines

#put indicator on variants whose value is this censored max, for easy filtering downstream if desired
counts_A$max_cens <- 0
counts_A[Kd==7.27e-05,max_cens := 1]
counts_B$max_cens <- 0
counts_B[Kd==7.14e-5,max_cens := 1]


#within the 10-6 mode -- is there differentiation within this mode?
plot.titration.A(which(counts_A$Kd > 1e-6 & counts_A$Kd < 1.5e-6)[2])#typically have ~2 upturned points at end
plot.titration.A(which(counts_A$Kd > 1e-5 & counts_A$Kd < 1.5e-5)[7])#typically have 1 point at the end that is literally just above the baseline. 

plot.titration.B(which(counts_B$Kd > 1e-6 & counts_B$Kd < 1.5e-6)[5])#typically have ~2 upturned points at end (some appear to be bad fits but have large RSEs, shoudl be filtered)
plot.titration.B(which(counts_B$Kd > 1e-5 & counts_B$Kd < 1.5e-5)[4])#typically have 1 point at the end that is literally just above the baseline. 

#wonder if it's worth distinguishing these "max Kd" ones from these ~1e-5 ones. They are more or less the same, might be worth collabsing everything above e.g. 5e-6
#in the middle continuous range (1e-7)
plot.titration.A(which(counts_A$Kd > 1e-7 & counts_A$Kd < 1.5e-7)[4])#have not plateaued, but sampled several points into the response so curves look ok. But definitely extrapolating
plot.titration.B(which(counts_B$Kd > 1e-7 & counts_B$Kd < 1.5e-7)[4])
#in the middle continuous range (1e-8)
plot.titration.A(which(counts_A$Kd > 1e-8 & counts_A$Kd < 1.5e-8)[1])#have just plateuaed
plot.titration.B(which(counts_B$Kd > 1e-8 & counts_B$Kd < 1.5e-8)[4])#have just plateuaed
#at the minimum of the range

#low Kds (was getting censored at minimum Kd with clear errors, had to increase number of required points to fit curve as these were coming from curves with large numbers of missing points)
which(counts_A$Kd<10^-11.5)
plot.titration.A(which(counts_A$Kd < 1e-11)[6]) #have more missing points than average, wonder if believable?
which(counts_B$Kd<10^-11.5)
plot.titration.B(which(counts_B$Kd<10^-11.5)[2])


#next, look at distributions of baseline, spot check the boundary cases and see if I need to relax these limits
hist(counts_A$baseline) #substantial number falling at the enforced minimum baseline of 1.0. What do these look like?
hist(counts_B$baseline)
plot.titration.A(which(counts_A$baseline==1.0)[4]) #some of these (but not all) also have the enforced minimum response of 2.0
plot.titration.A(which(counts_A$baseline==1.5)[1])#spot checks don't really raise concerns for me
plot.titration.B(which(counts_B$baseline==1.0)[7]) #some of these (but not all) also have the enforced minimum response of 2.0
plot.titration.B(which(counts_B$baseline==1.5)[5])


#next look at response
hist(counts_A$response)
hist(counts_B$response)
plot.titration.A(which(counts_A$response==2.0)[4]) #many have Kd above the highest concentration point
plot.titration.A(which(counts_A$response==3)[5]) #some of these are all curving up right at the end of the concentration so thye hit the max response because they're just guessing where the plateau is
plot.titration.B(which(counts_B$response==2.0)[6]) #many have Kd above the highest concentration point
plot.titration.B(which(counts_B$response==3)[2]) #some of these are all curving up right at the end of the concentration so thye hit the max response because they're just guessing where the plateau is


#next, RSE. 
hist(counts_A$RSE,main="",xlab="RSE")
hist(counts_B$RSE,main="",xlab="RSE",col="#7f7f7f50",add=T)
#see that the high-RSE are weighted toward lower average cell counts, which makes sense
plot(log10(counts_A$libA2_191010_avgcount+0.5),counts_A$RSE,pch=19,col="#00000005")
plot(log10(counts_B$libB2_200102_avgcount+0.5),counts_B$RSE,pch=19,col="#00000005")
#RSE versus Kd
plot(log10(counts_A$Kd),counts_A$RSE,pch=19,col="#00000005")
plot(log10(counts_B$Kd),counts_B$RSE,pch=19,col="#00000005")

hist(log10(counts_A[RSE > quantile(counts_A$RSE,probs=0.95,na.rm=T),Kd]),freq=F,col="red",breaks=30)
hist(log10(counts_A[RSE < quantile(counts_A$RSE,probs=0.95,na.rm=T),Kd]),freq=F,add=T,col="blue",breaks=30) #high RSE is weighted toward high/intermediate Kds

hist(log10(counts_B[RSE > quantile(counts_B$RSE,probs=0.95,na.rm=T),Kd]),freq=F,col="red",breaks=30)
hist(log10(counts_B[RSE < quantile(counts_B$RSE,probs=0.95,na.rm=T),Kd]),freq=F,add=T,col="blue",breaks=30) #high RSE is weighted toward high/intermediate Kds

#what are the curves with very low RSE?
plot.titration.A(which(counts_A$RSE<0.05)[8]) #all flatlines/all bin1 (no binding) or very weak binding -- get max Kd possible
plot.titration.B(which(counts_B$RSE<0.05)[6]) #all flatlines/all bin1 (no binding) or very weak binding -- get max Kd possible
plot.titration(which(counts$RSE<0.2 & counts$RSE>0.1)[2]) #look great
#intermediate RSE
plot.titration.A(which(counts_A$RSE<0.55 & counts_A$RSE > 0.45)[9]) #all reasonable fits
plot.titration.B(which(counts_B$RSE<0.55 & counts_B$RSE > 0.45)[8]) #all reasonable fits; notice they have lower avg cell count than equivalently well fits from rep A
#med-high RSE
plot.titration.A(which(counts_A$RSE<1.05 & counts_A$RSE > 0.95)[9]) #ok fits
plot.titration.B(which(counts_B$RSE<1.05 & counts_B$RSE > 0.95)[6]) #ok fits
#high RSE
plot.titration.A(which(counts_A$RSE > 1.5)[1]) #v bad fits/noisy points
plot.titration.B(which(counts_B$RSE > 1.5)[7]) #v bad fits/noisy points (honestly, maybe not that bad? or not as bad as A)
cutoff_A <- quantile(counts_A$RSE,probs=0.95,na.rm=T)
cutoff_B <- quantile(counts_B$RSE,probs=0.95,na.rm=T)
#around proposed cutoff
plot.titration.A(which(counts_A$RSE< cutoff_A + 0.02 & counts_A$RSE > cutoff_A - 0.02)[4]) #ok fits, we'll go with them
plot.titration.B(which(counts_B$RSE< cutoff_B + 0.02 & counts_B$RSE > cutoff_B - 0.02)[5]) #ok fits, we'll go with them

#replace as NA values where RSE > cutoff
counts_A[RSE > cutoff_A,c("Kd", "Kd_SE", "response", "baseline", "RSE", "fit") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))]
sum(!is.na(counts_A$Kd))/nrow(counts_A) #have Kds for 71.2% of my barcodes (111,778)
counts_B[RSE > cutoff_B,c("Kd", "Kd_SE", "response", "baseline", "RSE", "fit") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))]
sum(!is.na(counts_B$Kd))/nrow(counts_B) #have Kds for 48.0% of my barcodes (85,691)

#convert from Kd to ddG. 
#dG = -RTln(Kd) or (1.987204*10^-3)*(296.15)*ln(Kd)
counts_A[,dG := (1.987204*10^-3)*(296.15)*log(Kd),by=barcode]
counts_B[,dG := (1.987204*10^-3)*(296.15)*log(Kd),by=barcode]
#dG_SE = RT*Kd_SE/Kd
counts_A[,dG_SE := (1.987204*10^-3)*(296.15)*Kd_SE/Kd,by=barcode]
counts_B[,dG_SE := (1.987204*10^-3)*(296.15)*Kd_SE/Kd,by=barcode]
#Kd_SE estimates are very large for Kds that are at the boundary. Make these NA? Or some very large maximum? (e.g. 4.2, says this is an uncertainty of three orders of magnitude, very large!). Global epistasis can't have NA variances, so I'll just put a very high maximum
hist(counts_A$dG_SE[counts_A$dG_SE<5],main="",xlab="SE of dG")
hist(counts_B$dG_SE[counts_B$dG_SE<5],main="",xlab="SE of dG")
plot.titration.A(which(counts_A$dG_SE > 5)[20])
plot.titration.B(which(counts_B$dG_SE > 5)[8]) #all with Kds above highest concentration, as far as I can see
counts_A[dG_SE>2.8,dG_SE:=2.8]
counts_B[dG_SE>2.8,dG_SE:=2.8]
hist(counts_A$dG_SE)
hist(counts_B$dG_SE,col="#7f7f7f50",add=T)
plot(counts_A$dG,counts_A$dG_SE,pch=19,col="#00000010") #clear bias in the dG_SE measurements. Will try global epistasis with and without error estimates. Will think if there's anything else I can do...
plot(counts_B$dG,counts_B$dG_SE,pch=19,col="#00000010") #clear bias in the dG_SE measurements. Will try global epistasis with and without error estimates. Will think if there's anything else I can do...
plot(counts_A$RSE,counts_A$dG_SE,pch=19,col="#00000010") #could consider using something based on RSE as a correlate of variability?
plot(counts_B$RSE,counts_B$dG_SE,pch=19,col="#00000010") #could consider using something based on RSE as a correlate of variability?
#further reduce the max error?. I think having these large errors gives the global epistasis too much "flexibility" -- I'd rather have it extrapolate above the max Kd because of its global curve fit, not because it thinks it can just put points there because of their inherent experimental variance
counts_A[dG_SE>1.4,dG_SE:=1.4]
counts_B[dG_SE>1.4,dG_SE:=1.4]
hist(counts_A$dG_SE)
hist(counts_A[max_cens==0,dG_SE])
#add another indicator variable for whether the dG_SE was greater than 1.4
counts_A$max_SE <- 0
counts_A[dG_SE==1.4,max_SE := 1]
counts_B$max_SE <- 0
counts_B[dG_SE==1.4,max_SE := 1]
table(counts_A[!is.na(Kd),max_cens],counts_A[!is.na(Kd),max_SE])
table(counts_B[!is.na(Kd),max_cens],counts_B[!is.na(Kd),max_SE])

#ddG = dG_mut - dG_wt
dG_wt_A <- mean(counts_A[variant_class %in% c("synonymous","wildtype"),dG],na.rm=T)
dG_wt_B <- mean(counts_B[variant_class %in% c("synonymous","wildtype"),dG],na.rm=T)
counts_A[,ddG := dG - dG_wt_A,by=barcode]
counts_B[,ddG := dG - dG_wt_B,by=barcode]
#ddG_SE = sqrt(dG_SE^2 + dG_SE_wt^2)
dG_SE_wt_A <- sd(counts_A[variant_class %in% c("synonymous","wildtype"),dG],na.rm=T)/sqrt(sum(!is.na(counts_A[variant_class %in% c("synonymous","wildtype"),dG])))
dG_SE_wt_B <- sd(counts_B[variant_class %in% c("synonymous","wildtype"),dG],na.rm=T)/sqrt(sum(!is.na(counts_B[variant_class %in% c("synonymous","wildtype"),dG])))
counts_A[,ddG_SE := sqrt(dG_SE^2 + dG_SE_wt_A^2),by=barcode]
counts_B[,ddG_SE := sqrt(dG_SE^2 + dG_SE_wt_B^2),by=barcode]

#histogram of ddG
#pdf(file=paste(config$Titeseq_Kds_dir,"/hist_libA_ddG-DFE.pdf",sep=""),width=5,height=4.5,useDingbats=F,bg="transparent")
hist(counts_A[variant_class %in% (c("1 nonsynonymous",">1 nonsynonymous")),ddG],col="gray40",xlim=c(-1.5,9),breaks=50,main="",xlab="ddG (kcal/mol)")
hist(counts_A[variant_class %in% (c("synonymous","wildtype")),ddG],col="blue",add=T,breaks=50)
#dev.off()
#make, stacking >1 nonsynon, 1 nonsyn, and wildtype distributions on top of one another
#pdf(file=paste(config$Titeseq_Kds_dir,"/hist_libA_ddG-DFE_stacked-single-multiple-wildtype.pdf",sep=""),width=5,height=4.5,useDingbats=F,bg="transparent")
hist(counts_A[variant_class %in% (c("1 nonsynonymous",">1 nonsynonymous","synonymous","wildtype")),ddG],col="gray40",xlim=c(-1.5,9),breaks=50,main="",xlab="ddG (kcal/mol)")
hist(counts_A[variant_class %in% (c("1 nonsynonymous","synonymous","wildtype")),ddG],col="gray80",xlim=c(-1.5,9),breaks=50,main="",xlab="ddG (kcal/mol)",add=T)
hist(counts_A[variant_class %in% (c("synonymous","wildtype")),ddG],col="#0000ff",add=T,breaks=50)
#dev.off()
#hist(rep(counts_A[variant_class %in% (c("stop")),ddG],25),col="red",add=T,breaks=50) #repeat 15 times to be able to see these low counts
#rep B
#pdf(file=paste(config$Titeseq_Kds_dir,"/hist_libB_ddG-DFE.pdf",sep=""),width=5,height=4.5,useDingbats=F,bg="transparent")
hist(counts_B[variant_class %in% (c("1 nonsynonymous",">1 nonsynonymous")),ddG],col="gray40",xlim=c(-1.5,9),breaks=50,main="",xlab="ddG (kcal/mol)")
hist(counts_B[variant_class %in% (c("synonymous","wildtype")),ddG],col="blue",add=T,breaks=50)
#hist(rep(counts_B[variant_class %in% (c("stop")),ddG],50),col="red",add=T,breaks=50) #repeat 50 times to be able to see these low counts
#dev.off()

#same histogram but for Kd
#pdf(file=paste(config$Titeseq_Kds_dir,"/hist_libA_Kd-DFE.pdf",sep=""),width=5,height=4.5,useDingbats=F,bg="transparent")
hist(log10(counts_A[variant_class %in% (c("1 nonsynonymous",">1 nonsynonymous")),Kd]),col="gray40",breaks=50,main="",xlab="Kd (M)",xlim=c(-12,-4))
hist(log10(counts_A[variant_class %in% (c("synonymous","wildtype")),Kd]),col="blue",add=T,breaks=40)
#dev.off()
#repB
#pdf(file=paste(config$Titeseq_Kds_dir,"/hist_libB_Kd-DFE.pdf",sep=""),width=5,height=4.5,useDingbats=F,bg="transparent")
hist(log10(counts_B[variant_class %in% (c("1 nonsynonymous",">1 nonsynonymous")),Kd]),col="gray40",breaks=50,main="",xlab="Kd (M)",xlim=c(-12,-4))
hist(log10(counts_B[variant_class %in% (c("synonymous","wildtype")),Kd]),col="blue",add=T,breaks=30)
#dev.off()

#stop variants seem a bit left-shifted even relative to the other library genotypes.
print(counts_A[variant_class == "stop" & !is.na(ddG),.(aa_substitutions,variant_call_support,libA2_191010_avgcount,ddG)],topn=100)
#seems like there are several factors: stops with high variant call support and high sequencing depth but ~neutral binding seem to be C-terminal biased stops. High depth and variant call support but N-term often severely decreased affinity. Other neutral stops often have low variant call support, low sequencing depth, or both
print(counts_B[variant_class == "stop" & !is.na(ddG),.(aa_substitutions,variant_call_support,libB2_200102_avgcount,ddG)],topn=100)

#make NA estimates for stop variants
counts_A[variant_class == "stop",c("Kd", "Kd_SE", "response", "baseline", "RSE", "fit", "dG", "dG_SE", "ddG", "ddG_SE") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))]
counts_B[variant_class == "stop",c("Kd", "Kd_SE", "response", "baseline", "RSE", "fit", "dG", "dG_SE", "ddG", "ddG_SE") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))]

#save csv with relevant information for next steps, and save the entire .Rda data frame for counts in case I want to load and look at fits, etc.
counts_A[,library:="libA2"]
counts_A[,sample:="191010"]
counts_B[,library:="libB2"]
counts_B[,sample:="200102"]

#want to show whether the max_cens and max_SE filtered variants are taking away many single mutants?
table(counts_A[!is.na(ddG),variant_class])
table(counts_A[!is.na(ddG) & max_cens==0,variant_class])
table(counts_A[!is.na(ddG) & max_SE==0,variant_class])

table(counts_B[!is.na(ddG),variant_class])
table(counts_B[!is.na(ddG) & max_cens==0,variant_class])
table(counts_B[!is.na(ddG) & max_SE==0,variant_class])


write.csv(rbind(counts_A[,.(library, sample, barcode, variant_call_support, average_counts=libA2_191010_avgcount, dG, dG_SE, ddG, ddG_SE, Kd, response, baseline, RSE, max_cens, max_SE, variant_class,aa_substitutions, n_aa_substitutions)],counts_B[,.(library, sample, barcode, variant_call_support, average_counts=libB2_200102_avgcount, dG, dG_SE, ddG, ddG_SE, Kd, response, baseline, RSE, max_cens, max_SE, variant_class, aa_substitutions, n_aa_substitutions)]),file=config$Titeseq_Kds_file)
save(counts_A,file=paste(config$Titeseq_Kds_dir,"/dt.temp.A.Rda",sep=""))
save(counts_B,file=paste(config$Titeseq_Kds_dir,"/dt.temp.B.Rda",sep=""))

###load(file=paste(config$Titeseq_Kds_dir,"/dt.temp.Rda",sep=""))
