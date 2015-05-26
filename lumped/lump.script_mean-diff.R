setwd("/Users/andrewtrlica/Desktop/IPNI meta/lumped")

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
#f=2
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
rm(list="set", "hup", "i")

can <- NULL ## eliminate studies not including inhibitor
  for (i in 1:length(names)){ 
    yip <- which(fert.dat$Author==names[i])
    can <- c(can, yip)
  }

fert.dat <- fert.dat[-can,] 
rm(list="can", "names", "nope", "i", "yip")


#tapply(fert.dat$yld, INDEX= fert.dat$Author, FUN=length)
auths <- unique(fert.dat$Author) ### compute composite mean/sd/n by inhibitor class within each study 
fert.res <- data.frame(NULL)
  for (a in 1:length(auths)){
    set <- fert.dat[which(fert.dat$Author==auths[a]),]
    pkg <- NULL
    auth.res <- data.frame(NULL)
    inhibs <- unique(set$Inhibitor)
    for (m in 1:length(inhibs)){
      #m=3
      #a=1
      sut <- set[which(set$Inhibitor==inhibs[m]),]
      mean.comp <- sum(sut$yld*sut$yld.n)/sum(sut$yld.n)
      sd.comp <- sqrt((sum((sut$yld.n-1)*(sut$yld.sd^2)))/(sum(sut$yld.n-1)))  
      n.comp <- sum(sut$yld.n)
      pkg <- cbind(auths[a], paste(ferts[f]), paste(inhibs[m]), mean.comp, sd.comp, n.comp)
      auth.res <- rbind(auth.res, pkg)
      }
    colnames(auth.res)[1:3] <- c("Study", "Fert", "Inhibitor")
    fert.res <- rbind(fert.res, auth.res) ### package up all composite means for each inhibitor in each study
    }
fert.res$mean.comp <- as.numeric(as.character(fert.res$mean.comp))
fert.res$sd.comp <- as.numeric(as.character(fert.res$sd.comp))
fert.res$n.comp <- as.numeric(as.character(fert.res$n.comp))

rm(list="auth.res","pkg","set","sut","mean.comp","n.comp","sd.comp") 

if(Nlim=="on"){
  write.csv(fert.res, file=paste(ferts[f],"lump.culled-data.Nlim.csv")) ## have a look at culled study entries
}
else{
  write.csv(fert.res, file=paste(ferts[f],"lump.culled-data.csv")) ## have a look at culled study entries
}

### make sure every study has at least one composite mean for a control and an inhibitor
keep <- NULL

for (q in 1:length(auths)){
  #q=25
  check <- fert.res[which(fert.res$Study==auths[q]),]
  has.inhib <- length(which(check$Inhibitor!="none"))
  has.ctrl <- length(which(check$Inhibitor=="none"))
  
  if(has.inhib >= 1 & has.ctrl >= 1){
    keep <- c(keep, auths[q])
  }
}
fert.res <- fert.res[fert.res$Study %in% keep,]
rm(list="keep","has.inhib","has.ctrl","q")

#### compute effect size (mean difference) and variance for each composited ctrl/trt combo within each study
auths <- unique(fert.res$Study)
effect <- data.frame(NULL)
  for (k in 1:length(auths)){
    #k=11
    sub <- fert.res[which(fert.res$Study==auths[k]),]
    inhibs <- unique(sub[which(sub$Inhibitor!="none"),"Inhibitor"])
    for(j in 1:length(inhibs)){
      #j=1
      ctrl <- which(sub$Inhibitor=="none")
      trt <- which(sub$Inhibitor==paste(inhibs[j]))
      mean.diff <- sub$mean.comp[trt]-sub$mean.comp[ctrl]
      diff.var <- ((sub$sd.comp[trt]^2)/sub$n.comp[trt])+((sub$sd.comp[ctrl]^2)/sub$n.comp[ctrl]) ### var D not assuming same underlying S2
      
      ### package the effect sizes/vars in single output line
      pkg <- cbind(as.character(sub$Study[1]), as.character(sub$Fert[1]), as.character(inhibs[j]), mean.diff, diff.var)
      effect <- rbind(effect, pkg)
    }
  }
colnames(effect)[1:3] <- c("Study", "Fert", "Inhibitor")
rm(list="mean.diff","diff.var","pkg","sub")

### break up effect sizes by inhibitor and across all qualifying studies compute R.E. models weights, summary effects/stats
inhibs <- unique(effect$Inhibitor)
results <- NULL
archive <- data.frame(NULL)
  for(b in 1:length(inhibs)){
    #b=4
    sab <- effect[which(effect$Inhibitor==inhibs[b]),]
    sab$mean.diff <- as.numeric(as.character(sab$mean.diff))
    sab$diff.var <- as.numeric(as.character(sab$diff.var))
    sab$Wi <- 1/(sab$diff.var)
    Q <- sum(sab$Wi*(sab$mean.diff^2))-(sum(sab$Wi*sab$mean.diff)^2)/sum(sab$Wi)
    df <- dim(sab)[1]-1
    C <- sum(sab$Wi)-((sum(sab$Wi^2))/sum(sab$Wi))
    T2 <- (Q-df)/C ## DerSimmonian, R. & N. Liard. 1986. Controlled Clinical Trials 7(3): 177-188.
      if(Q<df){
        T2 <- 0
      }  
    sab$Wistar <- 1/(sab$diff.var+T2)
    Var.M <- 1/(sum(sab$Wistar))
    if(dim(sab)[1]<=1){
      sab$Wistar <- rep(1, dim(sab)[1])
      Var.M <- 1/sab$Wi
    }    
    I2 <- (Q-df)/Q
    print(c(paste(ferts[f], inhibs[b], "Summary effect", M.star <- sum(sab$Wistar*sab$mean.diff)/(sum(sab$Wistar)))))
    print(c(paste(ferts[f], inhibs[b], "Z", Zstar <- M.star/sqrt(Var.M))))
    print(c(paste(ferts[f], inhibs[b], "p(Z)", pz <- 2*(1-pnorm(abs(Zstar),mean = 0,sd = 1)))))
    print(c(paste(ferts[f], inhibs[b], "Lower 95%", LL.mean <- M.star-(1.96*sqrt(Var.M)))))
    print(c(paste(ferts[f], inhibs[b], "Upper 95%", UL.mean <- M.star+(1.96*sqrt(Var.M)))))
    
    report <- data.frame(ferts[f], inhibs[b], length(unique(sab$Study)), cbind(M.star, Var.M, T2, I2, Q, C, df, Zstar, pz, LL.mean, UL.mean))
    results <- rbind(results,report)
    archive <- rbind(archive,sab)
  }## summary effect loop closure

colnames(results)[1:5] <- c("Fert", "Inhibitor", "N.Studies", "Summary effect", "Var. Sum. Effect")
colnames(archive)[6:7] <- c("Weight(Wi)", "R.E. weight")

if(Nlim=="on"){
  write.csv(results, file=paste(ferts[f],"meta.lump.diff-results.Nlim.csv", sep="."))
  write.csv(archive, file=paste(ferts[f],"effect.sizes.mean-diff.Nlim.csv", sep=".")) ## export tables of effect sizes
}
else{
  write.csv(results, file=paste(ferts[f],"meta.lump.diff-results.csv", sep="."))
  write.csv(archive, file=paste(ferts[f],"effect.sizes.mean-diff.csv", sep=".")) ## export tables of effect sizes
  }
}### fert loop closure


### clean-up garbage
kill <- ls()
kill <- kill[-which(kill=="data")]
rm(list=kill)
rm(kill)