setwd("/Users/andrewtrlica/Desktop/IPNI_meta/parsed")

data <- read.csv("../SAS data_5-5-2015-1pm.csv", header=TRUE) ### read in and clean up database
ded <- which(data$yld.sd==".")
data <- data[-ded,] ### remove studies without report of sd of mean yield
data$yld.sd <- as.numeric(as.character(data$yld.sd))
data$yld <- as.numeric(as.character(data$yld))
rm(ded)

Nlim <- c("off") ### set N limiter on/off (< 150 lbs N/ac)
if (Nlim=="on"){
  data <- data[which(data$ratetot<150),]
  setwd("./Nlimit")
}

ferts <- c("urea", "UAN", "AA") ### parse database by source
for (f in 1:length(ferts)){
#f=1
fert.dat <- data[which(data$Source==ferts[f]),] ## restrict to only studies looking at single fert source

fert.dat$Author <- as.character(fert.dat$Author) ## more data clean up
fert.dat$Inhibitor <- as.character(fert.dat$Inhibitor)

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
rm(i)

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

#tapply(fert.dat$yld, INDEX= fert.dat$Author, FUN=length)

if(Nlim=="on"){
  write.csv(fert.dat, file=paste(ferts[f],"split.culled-data.Nlim.csv")) ## have a look at culled study entries
}
else{
  write.csv(fert.dat, file=paste(ferts[f],"split.culled-data.csv")) ## have a look at culled study entries
}

####### Sort and pair Fert vs. Fert+inhibitor datapoints along matching experimental conditions inside each study #####
names <- unique(fert.dat$Author) ### list of kept individual studies
params <- colnames(fert.dat)[c(5,7,8,9,10,11,12,13,14,18,19,20,22,26,31,32)] ### experimental parameters to match between datapoints
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
    for(d in 1:dim(set)[1]){ ## identify how many exp. parameters match between each database entry
      #d=3
      #l=1
      test <- sum(set[l, params] == set[d, params])
      track[l,d] <- test
      }
    }
  rm(test)

  truck <- upper.tri(track, diag=TRUE) ## use only upper half of matching score matrix (only unique comb. of entries with identical params)
  
  for(e in 1:dim(track)[1]){  ### loop through each unique combination of observations with matching experimental conditions
    #e=4
    sub <- set[which(track[e,]==16 & truck[e,]==TRUE),] ### if sub has only 1 member this fucks up
    ctrl <- which(sub$Inhibitor=="none")
    trt <- which(sub$Inhibitor==paste(inhibs[g]))

    if(length(sub$yld[ctrl])>=1 & length(sub$yld[trt])>=1){ ##for all sets of identical experiments with possible none/inhibitor comparison
      mean.ctrl <- sum(sub$yld[ctrl]*sub$yld.n[ctrl])/sum(sub$yld.n[ctrl])
      mean.trt <- sum(sub$yld[trt]*sub$yld.n[ctrl])/sum(sub$yld.n[trt])
      n.ctrl <- sum(sub$yld.n[ctrl])
      n.trt <- sum(sub$yld.n[trt])
      var.ctrl <- sum((sub$yld.n[ctrl]-1)*(sub$yld.sd[ctrl]^2))/sum(sub$yld.n[ctrl]-1) ### used pooled sd
      sd.ctrl <- sqrt(var.ctrl)
      var.trt <- sum((sub$yld.n[trt]-1)*(sub$yld.sd[trt]^2))/sum(sub$yld.n[trt]-1) ### used pooled sd
      sd.trt <- sqrt(var.trt)
    
        ### package the results in single output line
        pkg <- cbind(sub$Author[1], as.character(sub$Source[1]), inhibs[g], mean.ctrl, sd.ctrl, n.ctrl, mean.trt, sd.trt, n.trt)
        report <- cbind(pkg,sub[1,params])
        
        fert.sum <- rbind(fert.sum,report)
        }
      }## close loop on pairing/pooling identical ctrl/trt
    }##close loop on inhibitors within author+fert combos
  }## close loop on authors within Fert categories

### clean-up garbage
rm(pkg)
rm(track)
rm(truck)
rm(sub)
rm(study)
rm(set)
rm(report)
rm(d)
rm(e)
rm(l)
rm(g)
rm(k)
rm(ctrl)
rm(trt)
rm(trts)
rm(inhibs)
rm(names)
rm(n.ctrl)
rm(n.trt)
rm(sd.trt)
rm(sd.ctrl)
rm(mean.trt)
rm(mean.ctrl)
rm(var.ctrl)
rm(var.trt)


### clean up results
fert.res <- fert.sum
colnames(fert.res)[1:3] <- c("study", "source", "inhibitor")
fert.res$study <- as.character(fert.res$study)
fert.res$mean.ctrl <- as.numeric(as.character(fert.res$mean.ctrl))
fert.res$sd.ctrl <- as.numeric(as.character(fert.res$sd.ctrl))
fert.res$n.ctrl <- as.numeric(as.character(fert.res$n.ctrl))
fert.res$mean.trt <- as.numeric(as.character(fert.res$mean.trt))
fert.res$sd.trt <- as.numeric(as.character(fert.res$sd.trt))
fert.res$n.trt <- as.numeric(as.character(fert.res$n.trt))



### perform calcs for mean effects and significance, summarize
treats <- as.character(unique(fert.res$inhibitor))
results <- data.frame(NULL)
pairs.final <- data.frame(NULL)
for (b in 1:length(treats)){
  #b=1
  set <- fert.res[which(fert.res$inhibitor==treats[b]),] ### computing summary effect size for all of one inhibitor vs. one fert
  set$diff <- set$mean.trt-set$mean.ctrl  ### mean difference (D) as effect size
  set$diff.var <- ((set$sd.trt^2)/set$n.trt)+((set$sd.ctrl^2)/set$n.ctrl) ### var D not assuming same underlying S2
  
  set$Wi <- 1/(set$diff.var)
  Q <- sum(set$Wi*(set$diff^2))-(sum(set$Wi*set$diff)^2)/sum(set$Wi)
  df <- dim(set)[1]-1
  C <- sum(set$Wi)-((sum(set$Wi^2))/sum(set$Wi))
  T2 <- (Q-df)/C ## DerSimmonian, R. & N. Liard. 1986. Controlled Clinical Trials 7(3): 177-188.
  if(Q<df){
    T2 <- 0
  }
  
  set$Wistar <- 1/(set$diff.var+T2)
  I2 <- (Q-df)/Q
  Var.M <- 1/(sum(set$Wistar))
  print(c(paste(ferts[f], treats[b], "Summary effect", M.star <- sum(set$Wistar*set$diff)/(sum(set$Wistar)))))
  print(c(paste(ferts[f], treats[b], "Z", Zstar <- M.star/sqrt(Var.M))))
  print(c(paste(ferts[f], treats[b], "p(Z)", pz <- 2*(1-pnorm(abs(Zstar),mean = 0,sd = 1)))))
  print(c(paste(ferts[f], treats[b], "Lower 95%", LL.mean <- M.star-(1.96*sqrt(Var.M)))))
  print(c(paste(ferts[f], treats[b], "Upper 95%", UL.mean <- M.star+(1.96*sqrt(Var.M)))))
  
  report <- data.frame(treats[b], length(unique(set$study)), dim(set)[1], cbind(M.star, Var.M, T2, I2, Q, C, df, Zstar, pz, LL.mean, UL.mean))
  results <- rbind(results,report)
  chunk <- cbind(fert.res[which(fert.res$inhibitor==treats[b]),], set$diff, set$diff.var, set$Wistar)
  pairs.final <- rbind(pairs.final, chunk)
}

colnames(pairs.final)[26:28] <- c("effect size", "var. effect size", "R.E. weight")
colnames(results) <- c("inhibitor","studies N", "Effects N", "Summary effect", "Var Effect", "T2", "I2", "Q", "C", "df","Z*", "p(Z*)", "Lower 95%", "Upper 95%")

if(Nlim=="on"){
  write.csv(results, file=paste(ferts[f],"meta.split.diff-results.Nlim.csv"))
  write.csv(pairs.final, file=paste(ferts[f],"pairs-diff.Nlim.csv"))
}
else{
  write.csv(results, file=paste(ferts[f],"meta.split.diff-results.csv"))
  write.csv(pairs.final, file=paste(ferts[f],"pairs-diff.csv")) 
}



rm(list="chunk","T2","C","Q","I2","df","LL.mean","UL.mean","Zstar","M.star","Var.M","pz","treats","b")


}## close loop on fert categories
rm(list="set","f","ferts","fert.dat","fert.res","fert.sum","report","results","pairs.final","params")

