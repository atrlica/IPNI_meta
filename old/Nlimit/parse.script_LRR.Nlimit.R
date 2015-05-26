setwd("/Users/andrewtrlica/Desktop/IPNI meta/Nlimit")

data <- read.csv("../SAS data_5-5-2015-1pm.csv", header=TRUE) ### read in and clean up database
data$yld <- as.numeric(as.character(data$yld))
data$yld.sd <- as.numeric(as.character(data$yld.sd))

data <- data[which(data$ratetot<150),] ### select studies with lower N application rates
#range(data$yld[which(!is.na(data$yld))])
ferts <- c("urea", "UAN", "AA") ### parse database by source
for (f in 1:length(ferts)){
#f=2
fert.dat <- data[which(data$Source==ferts[f]),] ## restrict to only studies looking at single fert source
fert.dat <- fert.dat[which(fert.dat$yld>0),] ## studies must have reported yield

fert.dat$Author <- as.character(fert.dat$Author) ## more data clean up
fert.dat$Inhibitor <- as.character(fert.dat$Inhibitor)
#length(unique(fert.dat$Author))
#unique(fert.dat$Inhibitor)
#fert.dat <- fert.dat[which(fert.dat$Inhibitor!="S.R"),] ## eliminated, too few datapoints
#dim(fert.dat[which(fert.dat$Inhibitor=="NBPT"),])[1]

nope <- unique(fert.dat$Author) ## identify studies that did not include inhibitor
names <- NULL
for(i in 1:length(nope)){
  #i=1
  set <- fert.dat[which(fert.dat$Author==nope[i]),]
  hup <- length(which(set$Inhibitor!="none"))
  if(hup==0){
    names <- c(names,set$Author)
  }
}
names <- unique(names) ## studies that did not include any inhibitor
rm(set)
rm(hup)


can <- NULL ## eliminate studies not including inhibitor
for (i in 1:length(names)){ 
  yip <- which(fert.dat$Author==names[i])
  can <- c(can, yip)
}
fert.dat <- fert.dat[-can,] 
rm(can)
rm(names)
rm(nope)
rm(i)
rm(yip)

#tapply(urea.dat$yld, INDEX= urea.dat$Author, FUN=length)

#write.csv(fert.dat, file=paste(ferts[f],"culled-data.csv")) ## have a look at culled study entries
#ferts[2]

names <- unique(fert.dat$Author) ### list of kept individual studies
params <- colnames(fert.dat)[c(5,7,8,9,10,11,12,13,14,18,19,20,22,26,31,32)] ### experimental parameters to match between datapoints
#length(params)
#colnames(fert.dat)
#params
fert.sum <- matrix(NA,0,25) ## set up summary table
colnames(fert.sum) <- c("study","source","inhibitor","mean.ctrl","sd.ctrl","n.ctrl",
                        "mean.trt","sd.trt","n.trt", params)

for(k in 1:length(names)){ ## divide everything by study first
  #k=1
  study <- subset(fert.dat, Author == names[k])
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
  rm(test)

  
  truck <- upper.tri(track, diag=TRUE) ## use only upper half of matching score matrix (only unique comb. of entries with identical params)
  for(e in 1:dim(track)[1]){
    #e=4
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
      
      fert.sum <- rbind(fert.sum,report)
      }
    }
  }
}
rm(pkg)
rm(track)
rm(truck)
rm(inhibs)
rm(sub)
rm(study)
rm(set)
rm(mean.ctrl)
rm(mean.trt)
rm(n.ctrl)
rm(n.trt)
rm(sd.ctrl)
rm(sd.trt)
rm(var.ctrl)
rm(var.trt)
rm(params)
rm(names)


### clean up results
fert.res <- fert.sum
colnames(fert.res)[1:3] <- c("study", "source", "inhibitor")
#class(fert.res$study)
fert.res$study <- as.character(fert.res$study)
fert.res <- fert.res[!is.na(fert.res$study),]
fert.res$mean.ctrl <- as.numeric(as.character(fert.res$mean.ctrl))
fert.res$sd.ctrl <- as.numeric(as.character(fert.res$sd.ctrl))
fert.res$n.ctrl <- as.numeric(as.character(fert.res$n.ctrl))
fert.res$mean.trt <- as.numeric(as.character(fert.res$mean.trt))
fert.res$sd.trt <- as.numeric(as.character(fert.res$sd.trt))
fert.res$n.trt <- as.numeric(as.character(fert.res$n.trt))

fert.res <- fert.res[!is.nan(fert.res$mean.ctrl)&!is.nan(fert.res$mean.trt),]
fert.res <- fert.res[!is.na(fert.res$sd.ctrl)&!is.na(fert.res$sd.trt),]
#range(fert.res$n.ctrl)
#range(fert.res$n.trt)
#sum(duplicated(fert.res))
#fert.res <- fert.res[!duplicated(fert.res),]
#range(urea.res$mean.ctrl)[1]*16



### perform calcs for mean effects and significance, summarize
treats <- as.character(unique(fert.res$inhibitor))
results <- data.frame(NULL)
pairs.final <- data.frame(NULL)
for (b in 1:length(treats)){
  #b=1
  set <- fert.res[which(fert.res$inhibitor==treats[b]),] ### computing summary effect size for all of one inhibitor vs. one fert
  set$LRR <- log(set$mean.trt/set$mean.ctrl)  ### mean difference (D) as effect size
  set$LRR.var <- (set$sd.trt^2/(set$n.trt*(set$mean.trt^2)))+(set$sd.ctrl^2/(set$n.ctrl*(set$mean.ctrl^2)))  ### generalized var of LRR not assuming same underlying S2 -- see Hedges, L.V., Gurevich, J. and P.S. Curtis. 1999. Ecology 80(4):1150-1156.
  
  set$Wi <- 1/(set$LRR.var)
  Q <- sum(set$Wi*(set$LRR^2))-(sum(set$Wi*set$LRR)^2)/sum(set$Wi)
  df <- dim(set)[1]-1
  C <- sum(set$Wi)-((sum(set$Wi^2))/sum(set$Wi))
  T2 <- (Q-df)/C ## DerSimmonian, R. & N. Liard. 1986. Controlled Clinical Trials 7(3): 177-188.
  if(Q<df){
    T2 <- 0
  }
 
  set$Wistar <- 1/(set$LRR.var+T2)
  I2 <- (Q-df)/Q
  Var.M <- 1/(sum(set$Wistar))
  print(c(paste(ferts[f], treats[b], "Summary effect", M.star <- sum(set$Wistar*set$LRR)/(sum(set$Wistar)))))
  print(c(paste(ferts[f], treats[b], "Z", Zstar <- M.star/sqrt(Var.M))))
  print(c(paste(ferts[f], treats[b], "p(Z)", pz <- 2*(1-pnorm(abs(Zstar),mean = 0,sd = 1)))))
  print(c(paste(ferts[f], treats[b], "Lower 95%", LL.mean <- M.star-(1.96*sqrt(Var.M)))))
  print(c(paste(ferts[f], treats[b], "Upper 95%", UL.mean <- M.star+(1.96*sqrt(Var.M)))))
  
  report <- data.frame(treats[b], length(unique(set$study)), dim(set)[1], cbind(M.star, Var.M, T2, I2, Q, C, df, Zstar, pz, LL.mean, UL.mean))
  results <- rbind(results, report)
  chunk <- cbind(fert.res[which(fert.res$inhibitor==treats[b]),], set$LRR, set$LRR.var, set$Wistar)
  pairs.final <- rbind(pairs.final, chunk)
}

colnames(pairs.final)[26:28] <- c("effect size", "var. effect size", "R.E. weight")
colnames(results) <- c("inhibitor","studies N", "Effects N", "Summary effect", "Var Effect", "T2", "I2", "Q", "C", "df","Z*", "p(Z*)", "Lower 95%", "Upper 95%")
write.csv(results, file=paste(ferts[f],"meta.allsplit.LRR-results.Nlimit.csv"))
write.csv(pairs.final, file=paste(ferts[f],"pairs-LRR.Nlimit.csv"))

rm(T2)
rm(C)
rm(Q)
rm(I2)
rm(df)
rm(LL.mean)
rm(UL.mean)
rm(Zstar)
rm(M.star)
rm(Var.M)
rm(pz)
rm(treats)

}

rm(f)
rm(ferts)
rm(b)
rm(d)
rm(e)
rm(g)
rm(k)
rm(l)
rm(fert.dat)
rm(fert.res)
rm(fert.sum)
rm(report)
rm(results)
rm(trts)
rm(set)
rm(chunk)
rm(pairs.final)