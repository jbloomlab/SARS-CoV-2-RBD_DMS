#21 November 2019
#TNS
#script to analyze sequencing counts from Sort-seq experiments for surface expression

setwd("~/bloom_j/computational_notebooks/tstarr/2019/NIH45-46_DMS")

#import R libraries
library(yaml)
library(data.table)
library(tidyverse)
library(Hmisc)
library(fitdistrplus)

#read the configuration file
config <- read_yaml("config.yaml")

#make output files/directories
#make output directory
if(!file.exists(config$expression_sortseq_dir)){
  dir.create(file.path(config$expression_sortseq_dir))
}

#read dataframe with list of barcode runs
barcode_runs <- read.csv(file=config$barcode_runs,stringsAsFactors=F); barcode_runs <- subset(barcode_runs, select=-c(R1))

#read file giving count of each barcode in each sample partition, sort by barcode
counts <- data.table(read.csv(file=config$variant_counts_file,stringsAsFactors=F)); counts <- counts[order(counts$barcode),]

#eliminate rows from counts that are not part of an expression sort-seq bin
counts <- subset(counts, paste(library,sample) %in% sapply(1:nrow(barcode_runs[!is.na(barcode_runs$Sortseq_bin),c("library","sample")]), function(x) paste(barcode_runs[!is.na(barcode_runs$Sortseq_bin),"library"][x],barcode_runs[!is.na(barcode_runs$Sortseq_bin),"sample"][x])) )

#for each FITC bin in an experiment, normalize the read counts to the observed ratio of cell recovery among bins
for(i in as.numeric(row.names(barcode_runs[!is.na(barcode_runs$Sortseq_bin),]))){
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

#try out two metrics for functional scores from the Sortseq counts -- simplest is just weighted meanF by multiplying counts and mean fluor of cells in a bin (note, for bin1, this includes many cells that did not out-grow! definitely problematic. Won't effect ML calculation)
#more complicated is ML fit taking into account known fluorescence boundaries of the bin/interval censoring of the per-cell observations.
#mean bin
# calc.simplemean <- function(b1,b2,b3,b4,m1=1,m2=2,m3=3,m4=4){ #b1-4 gives cell observation counts in bins 1-4
#   return(sum(b1*m1,b2*m2,b3*m3,b4*m4)/sum(b1,b2,b3,b4))
# }
# counts_A[,simple_meanF := calc.simplemean(libA2_191113_FITC_bin1,libA2_191113_FITC_bin2,libA2_191113_FITC_bin3,libA2_191113_FITC_bin4,log(196),log(1505),log(4594),log(13135)),by=barcode]
# counts_B[,simple_meanF := calc.simplemean(libB2_191215_FITCbin1,libB2_191215_FITCbin2,libB2_191215_FITCbin3,libB2_191215_FITCbin4,log(250),log(2551),log(6380),log(16180)),by=barcode]
# 
# #the mean log-fluor for bin1 is misleading, because the vast majority of bin1 cells don't grow out (lost plasmid). This is evident in the colony counts. Therefore this bin's meanF gets dragged down which has a disproportionate pull down of meanF for variants with counts in bin1 that actually did grow out. I want to try altering this to make it more ~log-linear with the other three
# c(log(1505),log(4594),log(13135))
# b1.logF.est <- log(1505) - mean(c(log(13135)-log(4594), log(4594)-log(1505)))
# 
# counts[,simple_meanF2 := calc.simplemean(libA2_191113_FITC_bin1,libA2_191113_FITC_bin2,libA2_191113_FITC_bin3,libA2_191113_FITC_bin4,b1.logF.est,log(1505),log(4594),log(13135)),by=barcode]
# 
# hist(counts[variant_class %in% (c("1 nonsynonymous",">1 nonsynonymous")),simple_meanF2],col="gray40",breaks=50,main="",xlab="simple meanF")
# hist(counts[variant_class %in% (c("synonymous","wildtype")),simple_meanF2],col="blue",add=T,breaks=50)
# hist(counts[variant_class %in% (c("stop")),simple_meanF2],col="red",add=T,breaks=50)

#ML meanF
calc.MLmean <- function(b1,b2,b3,b4,min.b1,min.b2,min.b3,min.b4,max.b4,min.count=1){ #b1-4 gives cell observation counts in bins 1-4; remaining arguments give fluorescence boundaries of the respective bins; min.count gives minimum number of total observations needed across bins in order to calculate meanF (default 1)
  data <- data.frame(left=c(rep(min.b1,round(b1)),rep(min.b2,round(b2)),rep(min.b3,round(b3)),rep(min.b4,round(b4))),right=c(rep(min.b2,round(b1)),rep(min.b3,round(b2)),rep(min.b4,round(b3)),rep(max.b4,round(b4))))
  if(nrow(unique(data))>1 & nrow(data)>min.count){
    fit <- fitdistcens(data,"norm")
    return(list(as.numeric(summary(fit)$estimate["mean"]),as.numeric(summary(fit)$estimate["sd"])))
  } else {
    return(list(as.numeric(NA),as.numeric(NA)))
  }
}

# row <- 4
# counts[row,]
# data <- data.frame(left=c(rep(log(100),round(counts[row,libA2_191113_FITC_bin1]*2)),rep(log(827.5),round(counts[row,libA2_191113_FITC_bin2]*2)),rep(log(2832.5),round(counts[row,libA2_191113_FITC_bin3]*2)),rep(log(7418.5),round(counts[row,libA2_191113_FITC_bin4]*2))),
#                    right=c(rep(log(827.5),round(counts[row,libA2_191113_FITC_bin1]*2)),rep(log(2832.5),round(counts[row,libA2_191113_FITC_bin2]*2)),rep(log(7418.5),round(counts[row,libA2_191113_FITC_bin3]*2)),rep(log(200000),round(counts[row,libA2_191113_FITC_bin4]*2))))
# fit <- fitdistcens(data,"norm")
# summary(fit)

#multiply each count value by 50 to minimize rounding errors due to count normalization.
counts_A[,c("ML_meanF","ML_sdF") := tryCatch(calc.MLmean(b1=libA2_191113_FITC_bin1*50,
                                                       b2=libA2_191113_FITC_bin2*50,
                                                       b3=libA2_191113_FITC_bin3*50,
                                                       b4=libA2_191113_FITC_bin4*50,
                                                       min.b1=log(20),
                                                       min.b2=log(827.5),
                                                       min.b3=log(2832.5),
                                                       min.b4=log(7418.5),
                                                       max.b4=log(200000)),
                                           error=function(e){return(list(as.numeric(NA),as.numeric(NA)))}),by=barcode]
counts_A[,libA2_191113_count := sum(libA2_191113_FITC_bin1,libA2_191113_FITC_bin2,libA2_191113_FITC_bin3,libA2_191113_FITC_bin4),by=barcode]

counts_B[,c("ML_meanF","ML_sdF") := tryCatch(calc.MLmean(b1=libB2_191215_FITCbin1*50,
                                                         b2=libB2_191215_FITCbin2*50,
                                                         b3=libB2_191215_FITCbin3*50,
                                                         b4=libB2_191215_FITCbin4*50,
                                                         min.b1=log(20),
                                                         min.b2=log(1418.5),
                                                         min.b3=log(4178.5),
                                                         min.b4=log(9766.5),
                                                         max.b4=log(200000)),
                                             error=function(e){return(list(as.numeric(NA),as.numeric(NA)))}),by=barcode]
counts_B[,libB2_191215_count := sum(libB2_191215_FITCbin1,libB2_191215_FITCbin2,libB2_191215_FITCbin3,libB2_191215_FITCbin4),by=barcode]


#save R object with the fits since the ML meanF took >1hr to calculate
save(counts_A,file=paste(config$expression_sortseq_dir,"/dt.temp.A.Rda",sep=""))
save(counts_B,file=paste(config$expression_sortseq_dir,"/dt.temp.B.Rda",sep=""))

#load(file=paste(config$expression_sortseq_dir,"/dt.temp.A.Rda",sep=""))
#load(file=paste(config$expression_sortseq_dir,"/dt.temp.B.Rda",sep=""))

hist(counts_A[variant_class %in% (c("1 nonsynonymous",">1 nonsynonymous")),ML_meanF],col="gray40",breaks=50,main="",xlab="ML meanF")
hist(counts_A[variant_class %in% (c("synonymous","wildtype")),ML_meanF],col="#0000ff40",add=T,breaks=50)
hist(counts_A[variant_class %in% (c("stop")),ML_meanF],col="#ff000040",add=T,breaks=50)

hist(counts_B[variant_class %in% (c("1 nonsynonymous",">1 nonsynonymous")),ML_meanF],col="gray40",breaks=50,main="",xlab="ML meanF")
hist(counts_B[variant_class %in% (c("synonymous","wildtype")),ML_meanF],col="#0000ff40",add=T,breaks=50)
hist(counts_B[variant_class %in% (c("stop")),ML_meanF],col="#ff000040",add=T,breaks=50)

hist(log10(counts_A$libA2_191113_count+0.5),xlim=c(-0.5,5))
hist(log10(counts_B$libB2_191215_count+0.5),col="#7f7f7f50",add=T)

plot(counts_B$libB2_191215_count[order(counts_B$libB2_191215_count)],log="y",type="l",col="green",ylab="cell counts, FITC sort (green = libB, blue=libA)",xlab="variant, rank-ordered")
points(counts_A$libA2_191113_count[order(counts_A$libA2_191113_count)],log="y",type="l",col="blue")

# #look at NA ML_meanF -- make a combined meanF that uses simple meanF(2) for NA ML_meanFs?
# counts[,combined_meanF := ML_meanF]
# counts[!is.na(simple_meanF) & is.na(ML_meanF),combined_meanF := simple_meanF,by=barcode]

#estimate variance as a function of cell counts using repeat WT and STOP distributions
#I notice that a small number WT variants, even with high counts, have null-like meanFs. I don't think these are the "experimental noise" I want to capture in this analysis, but rather are genuinely unexpressing variants -- probably with a mutation outside of the PacBio sequencing or some other phenomenon. I will exclude these for the purpose of estimating standard error of the meanF estimates
wt_A <- counts_A[variant_class %in% c("synonymous","wildtype") & ML_meanF>7.5,]
stop_A <- counts_A[variant_class %in% c("stop") ,]
wt_B <- counts_B[variant_class %in% c("synonymous","wildtype") & ML_meanF>8,]
stop_B <- counts_B[variant_class %in% c("stop") ,]

#make bins of observations as a function of cell count
n.breaks.wt_A <- 30
n.breaks.stop_A <- 30
wt.bins_A <- data.frame(bin=1:n.breaks.wt_A)
stop.bins_A <- data.frame(bin=1:n.breaks.stop_A)
breaks.wt_A <- cut2(wt_A$libA2_191113_count,m=250,g=n.breaks.wt_A,onlycuts=T)
breaks.stop_A <- cut2(stop_A$libA2_191113_count,m=250,g=n.breaks.stop_A,onlycuts=T)
#pool observations that fall within the bin breaks, compute sd on different meanF types
for(i in 1:nrow(wt.bins_A)){
  wt.bins_A$range.cells[i] <- list(c(breaks.wt_A[i],breaks.wt_A[i+1]))
  data <- wt_A[libA2_191113_count >= wt.bins_A$range.cells[i][[1]][[1]] & libA2_191113_count < wt.bins_A$range.cells[i][[1]][[2]],]
  wt.bins_A$median.cells[i] <- median(data$libA2_191113_count,na.rm=T)
  #wt.bins_A$mean.simple_meanF[i] <- mean(data$simple_meanF,na.rm=T)
  #wt.bins_A$sd.simple_meanF[i] <- sd(data$simple_meanF,na.rm=T)
  wt.bins_A$mean.ML_meanF[i] <- mean(data$ML_meanF,na.rm=T)
  wt.bins_A$sd.ML_meanF[i] <- sd(data$ML_meanF,na.rm=T)
  #wt.bins_A$mean.combined_meanF[i] <- mean(data$combined_meanF,na.rm=T)
  #wt.bins_A$sd.combined_meanF[i] <- sd(data$combined_meanF,na.rm=T)
}
for(i in 1:nrow(stop.bins_A)){
  stop.bins_A$range.cells[i] <- list(c(breaks.stop_A[i],breaks.stop_A[i+1]))
  data <- stop_A[libA2_191113_count >= stop.bins_A$range.cells[i][[1]][[1]] & libA2_191113_count < stop.bins_A$range.cells[i][[1]][[2]],]
  stop.bins_A$median.cells[i] <- median(data$libA2_191113_count,na.rm=T)
  #stop.bins_A$mean.simple_meanF[i] <- mean(data$simple_meanF,na.rm=T)
  #stop.bins_A$sd.simple_meanF[i] <- sd(data$simple_meanF,na.rm=T)
  stop.bins_A$mean.ML_meanF[i] <- mean(data$ML_meanF,na.rm=T)
  stop.bins_A$sd.ML_meanF[i] <- sd(data$ML_meanF,na.rm=T)
  #stop.bins_A$mean.combined_meanF[i] <- mean(data$combined_meanF,na.rm=T)
  #stop.bins_A$sd.combined_meanF[i] <- sd(data$combined_meanF,na.rm=T)
}
#look at relationship between variance and cell counts; fit curve to estimate variance from cell count
#ML meanF
y1_A <- (wt.bins_A$sd.ML_meanF)^2;x1_A <- wt.bins_A$median.cells
y2_A <- (stop.bins_A$sd.ML_meanF)^2;x2_A <- stop.bins_A$median.cells
plot(x1_A,y1_A,xlab="number cells",ylab="variance",main="ML meanF",pch=19,col="blue",ylim=c(0,1));points(x2_A,y2_A,pch=19,col="red")
plot(log(x1_A),log(y1_A),xlab="log(number cells)",ylab="log(variance)",main="ML meanF",pch=19,col="blue");points(log(x2_A),log(y2_A),pch=19,col="red")
fit_ML_meanF_A <- lm(log(c(y1_A,y2_A)) ~ log(c(x1_A,x2_A)));summary(fit_ML_meanF_A);abline(fit_ML_meanF_A)

#repB -- noticed for repB that the stops are much more bimodal, and thus have an artificially high variance -- therefore, I will use just the WT sampling distribution to estimate variance
#make bins of observations as a function of cell count
n.breaks.wt_B <- 30
n.breaks.stop_B <- 30
wt.bins_B <- data.frame(bin=1:n.breaks.wt_B)
stop.bins_B <- data.frame(bin=1:n.breaks.stop_B)
breaks.wt_B <- cut2(wt_B$libB2_191215_count,m=250,g=n.breaks.wt_B,onlycuts=T)
breaks.stop_B <- cut2(stop_B$libB2_191215_count,m=250,g=n.breaks.stop_B,onlycuts=T)
#pool observations that fall within the bin breaks, compute sd on different meanF types
for(i in 1:nrow(wt.bins_B)){
  wt.bins_B$range.cells[i] <- list(c(breaks.wt_B[i],breaks.wt_B[i+1]))
  data <- wt_B[libB2_191215_count >= wt.bins_B$range.cells[i][[1]][[1]] & libB2_191215_count < wt.bins_B$range.cells[i][[1]][[2]],]
  wt.bins_B$median.cells[i] <- median(data$libB2_191215_count,na.rm=T)
  #wt.bins_B$mean.simple_meanF[i] <- mean(data$simple_meanF,na.rm=T)
  #wt.bins_B$sd.simple_meanF[i] <- sd(data$simple_meanF,na.rm=T)
  wt.bins_B$mean.ML_meanF[i] <- mean(data$ML_meanF,na.rm=T)
  wt.bins_B$sd.ML_meanF[i] <- sd(data$ML_meanF,na.rm=T)
  #wt.bins_B$mean.combined_meanF[i] <- mean(data$combined_meanF,na.rm=T)
  #wt.bins_B$sd.combined_meanF[i] <- sd(data$combined_meanF,na.rm=T)
}
for(i in 1:nrow(stop.bins_B)){
  stop.bins_B$range.cells[i] <- list(c(breaks.stop_B[i],breaks.stop_B[i+1]))
  data <- stop_B[libB2_191215_count >= stop.bins_B$range.cells[i][[1]][[1]] & libB2_191215_count < stop.bins_B$range.cells[i][[1]][[2]],]
  stop.bins_B$median.cells[i] <- median(data$libB2_191215_count,na.rm=T)
  #stop.bins_B$mean.simple_meanF[i] <- mean(data$simple_meanF,na.rm=T)
  #stop.bins_B$sd.simple_meanF[i] <- sd(data$simple_meanF,na.rm=T)
  stop.bins_B$mean.ML_meanF[i] <- mean(data$ML_meanF,na.rm=T)
  stop.bins_B$sd.ML_meanF[i] <- sd(data$ML_meanF,na.rm=T)
  #stop.bins_B$mean.combined_meanF[i] <- mean(data$combined_meanF,na.rm=T)
  #stop.bins_B$sd.combined_meanF[i] <- sd(data$combined_meanF,na.rm=T)
}
#look at relationship between variance and cell counts; fit curve to estimate variance from cell count
#ML meanF
y1_B <- (wt.bins_B$sd.ML_meanF)^2;x1_B <- wt.bins_B$median.cells
y2_B <- (stop.bins_B$sd.ML_meanF[-1])^2;x2_B <- stop.bins_B$median.cells[-1]
plot(x1_B,y1_B,xlab="number cells",ylab="variance",main="ML meanF",pch=19,col="blue",ylim=c(0,1));points(x2_B,y2_B,pch=19,col="red")
plot(log(x1_B),log(y1_B),xlab="log(number cells)",ylab="log(variance)",main="ML meanF",pch=19,col="blue",ylim=c(-3.5,-0.5));points(log(x2_B),log(y2_B),pch=19,col="red")
fit_ML_meanF_B <- lm(log(c(y1_B)) ~ log(c(x1_B)));summary(fit_ML_meanF_B);abline(fit_ML_meanF_B)


#assign estimated variance values for the different meanF estimates
est.var <- function(count,fit){
  return(exp(as.numeric(fit$coefficients[1]) + as.numeric(fit$coefficients[2]) * log(count)))
}
counts_A[,var_ML_meanF := est.var(libA2_191113_count,fit_ML_meanF_A),by=barcode]
counts_B[,var_ML_meanF := est.var(libB2_191215_count,fit_ML_meanF_B),by=barcode]

#filter for minimum five reads
min <- 5
hist(counts_A[libA2_191113_count>min & variant_class %in% (c("1 nonsynonymous",">1 nonsynonymous")),ML_meanF],col="gray40",breaks=50,main="",xlab="ML mean fluorescence")
hist(counts_A[libA2_191113_count>min & variant_class %in% (c("synonymous","wildtype")),ML_meanF],col="#0000ff40",add=T,breaks=50)
hist(counts_A[libA2_191113_count>min & variant_class %in% (c("stop")),ML_meanF],col="#ff000040",add=T,breaks=50)

hist(counts_B[libB2_191215_count>min & variant_class %in% (c("1 nonsynonymous",">1 nonsynonymous")),ML_meanF],col="gray40",breaks=50,main="",xlab="ML mean fluorescence")
hist(counts_B[libB2_191215_count>min & variant_class %in% (c("synonymous","wildtype")),ML_meanF],col="#0000ff40",add=T,breaks=50)
hist(counts_B[libB2_191215_count>min & variant_class %in% (c("stop")),ML_meanF],col="#ff000040",add=T,breaks=50)

counts.filtered_A <- copy(counts_A)
counts.filtered_B <- copy(counts_B)

counts.filtered_A[libA2_191113_count < 5, c("ML_meanF") := NA,by=barcode]
counts.filtered_B[libB2_191215_count < 5, c("ML_meanF") := NA,by=barcode]

sum(!is.na(counts.filtered_A$ML_meanF))/nrow(counts.filtered_A) #ML meanF for 72.5% of my barcodes (113820)
sum(!is.na(counts.filtered_B$ML_meanF))/nrow(counts.filtered_B) #ML meanF for 87.3% of my barcodes (155793)

#save mean fluorescence table for reading into global epistasis
counts.filtered_A[,library:="libA2"]
counts.filtered_A[,sample:="191113"]
counts.filtered_B[,library:="libB2"]
counts.filtered_B[,sample:="191215"]

write.csv(rbind(counts.filtered_A[,.(library, sample, barcode, variant_call_support, average_counts=libA2_191113_count, ML_meanF, var_ML_meanF, variant_class, aa_substitutions, n_aa_substitutions)],counts.filtered_B[,.(library, sample, barcode, variant_call_support, average_counts=libB2_191215_count, ML_meanF, var_ML_meanF, variant_class, aa_substitutions, n_aa_substitutions)]),file=config$expression_sortseq_file)
