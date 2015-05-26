setwd("/Users/andrewtrlica/Desktop/IPNI meta/")

data <- read.csv("SAS data_5-5-2015-1pm.csv", header=TRUE) ### read in an clean up data
data$yld <- as.numeric(as.character(data$yld))
data$yld.sd <- as.numeric(as.character(data$yld.sd))
#range(data$yld[which(!is.na(data$yld))])
#unique(data$Source)

urea.dat <- data[which(data$Source=="urea"),] ## restrict to only studies looking at urea
urea.dat <- urea.dat[which(urea.dat$yld>0),] ## studies must have reported yield

urea.dat$Author <- as.character(urea.dat$Author) ## more data clean up
urea.dat$Inhibitor <- as.character(urea.dat$Inhibitor)
length(unique(urea.dat$Author))
#unique(urea.dat$Inhibitor)
urea.dat <- urea.dat[which(urea.dat$Inhibitor!="S.R"),] ## eliminated, too few datapoints
#dim(urea.dat[which(urea.dat$Inhibitor=="NBPT"),])[1]

dum <- unique(urea.dat$Author) ## identify studies that did not include inhibitor
names <- NULL
for(i in 1:length(dum)){
  #i=1
  set <- urea.dat[which(urea.dat$Author==dum[i]),]
  hup <- length(which(set$Inhibitor!="none"))
  if(hup==0){
    names <- c(names,set$Author)
  }
}
names <- unique(names) ## studies that did not include any inhibitor
rm(set)
rm(hup)
can <- NULL
for (i in 1:length(names)){ ## eliminate studies not including inhibitor
  yip <- which(urea.dat$Author==names[i])
  can <- c(can, yip)
}
urea.dat <- urea.dat[-can,] 
rm(can)
rm(names)
rm(dum)
rm(i)
rm(yip)

#tapply(urea.dat$yld, INDEX= urea.dat$Author, FUN=length)

write.csv(urea.dat, file="urea.dat.csv") ## have a look at culled study entries


names <- unique(urea.dat$Author) ### list of kept individual studies
params <- colnames(urea.dat)[c(5,7,8,9,10,11,12,13,14,18,19,20,22,26,31,32)] ### experimental parameters to match between datapoints
#length(params)
#colnames(urea.dat)
#params
urea.sum <- matrix(NA,0,25) ## set up summary table
colnames(urea.sum) <- c("study","source","inhibitor","mean.ctrl","sd.ctrl","n.ctrl",
                        "mean.trt","sd.trt","n.trt", params)

for(k in 1:length(names)){ ## divide everything by study first
  #k=1
  study <- subset(urea.dat, Author == names[k])
  track <- matrix(NA, dim(study)[1], dim(study)[1])
  test <- NULL
  trts <- unique(study$Inhibitor)
  inhibs <- trts[which(trts!="none")]
  
  for(g in 1:length(inhibs)){
    #g=1
    set <- study[which(study$Inhibitor=="none" | study$Inhibitor==inhibs[g]),]
    
  for(l in 1:dim(set)[1]){
    for(d in 1:dim(set)[1]){ ## identify how many exp. parameters match between each line
      #d=3
      #l=1
      test <- sum(set[l, params] == set[d, params])
      track[l,d] <- test
      }
    }
  
  truck <- upper.tri(track, diag=TRUE) ## use only upper half of matching score matrix (only unique comb. of entries with identical params)
  for(e in 1:dim(track)[1]){
    #e=4treats <- as.character(unique(urea.res$inhibitor))
    sub <- set[track[e,]==16 & truck[e,]==TRUE,]
    if(dim(sub)[1]>1){ ##for all sets of identical experiments with potential for none/inhibitor comparison
      mean.ctrl <- sum(sub[which(sub$Inhibitor=="none"),"yld"])/length(sub[which(sub$Inhibitor=="none"),"yld"])
      mean.trt <- sum(sub[which(sub$Inhibitor!="none"),"yld"])/length(sub[which(sub$Inhibitor!="none"),"yld"])
      n.ctrl <- sum(sub[which(sub$Inhibitor=="none"),"yld.n"])
      n.trt <- sum(sub[which(sub$Inhibitor!="none"),"yld.n"])
      var.ctrl <- sum((sub[which(sub$Inhibitor=="none"),"yld.n"]-1)*(sub[which(sub$Inhibitor=="none"),"yld.sd"]^2))/sum(sub[which(sub$Inhibitor=="none"),"yld.n"]-1) ### used pooled sd
      sd.ctrl <- sqrt(var.ctrl)
      var.trt <- sum((sub[which(sub$Inhibitor!="none"),"yld.n"]-1)*(sub[which(sub$Inhibitor!="none"),"yld.sd"]^2))/sum(sub[which(sub$Inhibitor!="none"),"yld.n"]-1) ### used pooled sd
      sd.trt <- sqrt(var.trt)

      ### package the results in single output line
      pkg <- cbind(sub$Author[1], as.character(sub$Source[1]), inhibs[g], mean.ctrl, sd.ctrl, n.ctrl, mean.trt, sd.trt, n.trt)
      report <- cbind(pkg,sub[1,params])
      
      urea.sum <- rbind(urea.sum,report)
      }
    }
  }
}

### clean up results
urea.res <- urea.sum
colnames(urea.res)[1:3] <- c("study", "source", "inhibitor")
#class(urea.res$study)
urea.res$study <- as.character(urea.res$study)
urea.res <- urea.res[!is.na(urea.res$study),]
urea.res$mean.ctrl <- as.numeric(as.character(urea.res$mean.ctrl))
urea.res$sd.ctrl <- as.numeric(as.character(urea.res$sd.ctrl))
urea.res$n.ctrl <- as.numeric(as.character(urea.res$n.ctrl))
urea.res$mean.trt <- as.numeric(as.character(urea.res$mean.trt))
urea.res$sd.trt <- as.numeric(as.character(urea.res$sd.trt))
urea.res$n.trt <- as.numeric(as.character(urea.res$n.trt))

urea.res <- urea.res[!is.nan(urea.res$mean.ctrl)&!is.nan(urea.res$mean.trt),]
urea.res <- urea.res[!is.na(urea.res$sd.ctrl)&!is.na(urea.res$sd.trt),]
#range(urea.res$n.ctrl)
#range(urea.res$n.trt)
#sum(duplicated(urea.res))
#urea.res <- urea.res[!duplicated(urea.res),]
#range(urea.res$mean.ctrl)[1]*16

### perform calcs for mean effects and significance, summarize
treats <- as.character(unique(urea.res$inhibitor))
results <- data.frame(NULL)
for (b in 1:length(treats)){
  #b=2
  set <- urea.res[which(urea.res$inhibitor==treats[b]),]
  set$diff <- set$mean.trt-set$mean.ctrl
  set$diff.sd <- ((set$sd.trt^2)/set$n.trt)+((set$sd.ctrl^2)/set$n.ctrl) 
  
  set$Wi <- 1/(set$diff.sd^2)
  Q <- sum(set$Wi*(set$diff^2))-(sum(set$Wi*set$diff)^2)/sum(set$Wi)
  df <- dim(set)[1]-1
  C <- sum(set$Wi)-((sum(set$Wi^2))/sum(set$Wi))
  T2 <- (Q-df)/C
  
  set$Wistar <- 1/(set$diff.sd^2+T2)
  I2 <- (Q-df)/Q
  Var.M <- 1/(sum(set$Wistar))
  print(M.star <- sum(set$Wistar*set$diff)/(sum(set$Wistar)))
  print(Zstar <- M.star/sqrt(Var.M))
  print(pz <- 2*(1-pnorm(Zstar,mean = 0,sd = 1)))
  print(LL.mean <- M.star-(1.96*sqrt(Var.M)))
  print(UL.mean <- M.star+(1.96*sqrt(Var.M)))
  
  report <- data.frame(treats[b], dim(set)[1], cbind(M.star, Var.M, T2, I2, Q, C, Zstar, pz, LL.mean, UL.mean))
  results <- rbind(results,report)
}

colnames(results) <- c("inhibitor", "N", "Mean diff.", "Var Difference", "T2", "I2", "Q", "C", "Z*", "p(Z*)", "Lower 95%", "Upper 95%")

write.csv(results, file="urea.meta.results.csv")




