setwd("/Users/andrewtrlica/Desktop/IPNI meta/")

data <- read.csv("SAS data_5-5-2015-1pm.csv", header=TRUE)
data$yld <- as.numeric(as.character(data$yld))
data$yld.sd <- as.numeric(as.character(data$yld.sd))
#range(data$yld[which(!is.na(data$yld))])
#unique(data$Source)

AA.dat <- data[which(data$Source=="AA"),]
AA.dat <- AA.dat[which(AA.dat$yld.sd!="."),]
AA.dat$Author <- as.character(AA.dat$Author)
AA.dat$Inhibitor <- as.character(AA.dat$Inhibitor)
#unique(AA.dat$Author)
#unique(AA.dat$Inhibitor)

dum <- unique(AA.dat$Author)
names <- NULL
for(i in 1:length(dum)){
  set <- AA.dat[which(AA.dat$Author==dum[i]),]
  hup <- length(which(set$Inhibitor=="nitrapyrin"))
  if(hup==0){
    names <- c(names,set$Author)
  }
}
names <- unique(names)

can <- NULL
for (i in 1:length(names)){
  yip <- which(AA.dat$Author==names[i])
  can <- c(can, yip)
}

AA.dat <- AA.dat[-can,]

tapply(AA.dat$yld,INDEX = AA.dat$Author, FUN=length)

write.csv(AA.dat, file="AA.dat.csv")


names <- unique(AA.dat$Author)
studs <- NULL
yld.diff <- NULL
yld.sd <- NULL

#i=1
for (i in 1:length(names)){
  set <- AA.dat[which(AA.dat$Author==names[i]),]
  if (dim(set)[1]==2){ ### simple case, only 2 lines for the study, only one comparison to make
    pos.ctrl <- which(set$Inhibitor=="none")
    pos.trt <- which(set$Inhibitor=="nitrapyrin")
    diff <- set$yld[pos.trt]-set$yld[pos.ctrl]
    sd.diff <- ((set$yld.sd[pos.trt]^2)/set$yld.n[pos.trt])+((set$yld.sd[pos.ctrl]^2)/set$yld.n[pos.ctrl])
    studs <- c(studs,names[i])
    yld.diff <- c(yld.diff,diff)
    yld.sd <- c(yld.sd,sd.diff)
  }
  
  else if(length(unique(set$ratetot)) > 1){ ## case in which different ratetots applied
    rate <- unique(set$ratetot)
    num.rates <- length(rate)
    
    for (j in 1:num.rates){
      set2 <- set[which(set$ratetot==rate[j]),]
      pos.ctrl <- which(set2$Inhibitor=="none")
      pos.trt <- which(set2$Inhibitor=="nitrapyrin")
      diff <- set2$yld[pos.trt]-set2$yld[pos.ctrl]
      sd.diff <- ((set2$yld.sd[pos.trt]^2)/set2$yld.n[pos.trt])+((set2$yld.sd[pos.ctrl]^2)/set2$yld.n[pos.ctrl])
      studs <- c(studs,paste(names[i], rate[j]))
      yld.diff <- c(yld.diff, diff)
      yld.sd <- c(yld.sd, sd.diff)
    }  
  }
  
  else{
    fall <- set[which(set$ratefall>0),]
    pos.ctrl <- which(fall$Inhibitor=="none")
    pos.trt <- which(fall$Inhibitor=="nitrapyrin")
    diff <- fall$yld[pos.trt]-fall$yld[pos.ctrl]
    sd.diff <- ((fall$yld.sd[pos.trt]^2)/fall$yld.n[pos.trt])+((fall$yld.sd[pos.ctrl]^2)/fall$yld.n[pos.ctrl])
    studs <- c(studs,paste(names[i], "fall"))
    yld.diff <- c(yld.diff, diff)
    yld.sd <- c(yld.sd, sd.diff)
    
    spring <- set[which(set$ratefall==0),]
    pos.ctrl <- which(spring$Inhibitor=="none")
    pos.trt <- which(spring$Inhibitor=="nitrapyrin")
    diff <- spring$yld[pos.trt]-spring$yld[pos.ctrl]
    sd.diff <- ((spring$yld.sd[pos.trt]^2)/spring$yld.n[pos.trt])+((spring$yld.sd[pos.ctrl]^2)/spring$yld.n[pos.ctrl])
    studs <- c(studs,paste(names[i], "spring"))
    yld.diff <- c(yld.diff, diff)
    yld.sd <- c(yld.sd, sd.diff)
  }    
}


AA.sum <- as.data.frame(cbind(studs,yld.diff,yld.sd))
AA.sum <- AA.sum[-c(3,5),]

AA.sum$yld.diff <- as.numeric(as.character(AA.sum$yld.diff))
AA.sum$yld.sd <- as.numeric(as.character(AA.sum$yld.sd))

AA.sum$Wi <- 1/(AA.sum$yld.sd^2)
Q <- sum(AA.sum$Wi*(AA.sum$yld.diff^2))-(sum(AA.sum$Wi*AA.sum$yld.diff)^2)/sum(AA.sum$Wi)
df <- dim(AA.sum)[1]-1
C <- sum(AA.sum$Wi)-((sum(AA.sum$Wi^2))/sum(AA.sum$Wi))
T2 <- (Q-df)/C
AA.sum$Wistar <- 1/(AA.sum$yld.sd^2+T2)
I2 <- (Q-df)/Q
Var.M <- 1/(sum(AA.sum$Wistar))
print(M.star <- sum(AA.sum$Wistar*AA.sum$yld.diff)/(sum(AA.sum$Wistar)))
Zstar <- M.star/sqrt(Var.M)
print(2*(1-pnorm(Zstar,mean = 0,sd = 1)))
print(LL.mean <- M.star-(1.96*sqrt(Var.M)))
print(UL.mean <- M.star+(1.96*sqrt(Var.M)))
