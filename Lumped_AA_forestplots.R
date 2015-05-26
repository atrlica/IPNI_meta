setwd("/Users/andrewtrlica/Desktop/IPNI meta/")


library(metafor)
AA.dat <- read.csv(file="lumped/AA.effect.sizes.mean-diff.csv")
AA.dat$Inhibitor <- as.character(AA.dat$Inhibitor)
AA.dat$Name <- paste(AA.dat$Study, AA.dat$Inhibitor) ## mark each study effect size with its Inhibitor

### this loop comes up with a scaling factor for the squares to apply within each Inhibitor group
inhib <- unique(AA.dat$Inhibitor)
shit <- data.frame(NULL)
for(j in 1:length(inhib)){
  gunk <- AA.dat[which(AA.dat$Inhibitor==inhib[j]),]
  total <- sum(gunk$Inhibitor==inhib[j])
  skim <- gunk[which(gunk$Inhibitor==inhib[j]),8]/sum(gunk[which(gunk$Inhibitor==inhib[j]),8])
  #place <- rank(skim)
  #piece <- 1.5/(length(place)-1)
  #rwt <- 2-((place-1)*piece)
  gunk$rewt <- 1+(2*skim)
  shit <- rbind(shit, gunk)
}

### do the rma meta analysis to get the summary effect polygon and basic frame
AA.shit <- rma(yi = shit$mean.diff, vi = shit$diff.var, 
                measure = "GEN", method="DL",
                slab = shit$Name)
## expand margins 
par(mar=c(4,4,1,2), font=1)

##### the forest plot needs to be specified which rows for placing the effect sizes; 
##### you create gaps in the rows to provide room for titles and summary polygons;
##### the rows count up from the bottom of the figure (presume the all-study summary is on row 1)
spc <- 3.5
br1 <- 3 #NBPT needs 2
br2 <- spc+3+br1 # NBPT+DCD, needs 4
br3 <- spc+br2 # thiosulfate, needs 1
br4 <- spc+br3 # S.R., needs 1
br5 <- spc+2+br4 # Nitrap, needs 3

forest(AA.shit, xlim=c(-8, 8), cex=.75, ylim=c(-1, 20),
       xlab="Mean difference", mlab="RE Model for All Studies", 
       #=c((br1-1):br1,(br2-3):br2,br3,br4,(br5-2):br5),
       psize=shit$rewt)
### set font expansion factor (as in forest() above) and use bold italic
### font and save original settings in object 'op'
#mtext("AA")
op <- par(cex=.75, font=4)

### add text for the subgroups
text(-8, cex = 1, c(br1+1,br2+1,br3+1,br4+1,br5+1), pos=4, c("NBPT",
                                                             "NBPT+DCD",
                                                             "Thiosulfate",
                                                             "S.R.",
                                                             "Nitrapyrin"))



#### add subgroup polygons for summary effects by doing subgroup meta analysis and adding diamonds in the row spaces
#res.sr <- rma(yi = mean.diff, vi = diff.var, 
#             measure = "GEN", method="DL",
#            data=shit, 
#           subset=(Inhibitor=="S.R"))
res.nbpt <- rma(yi = mean.diff, vi = diff.var, 
                measure = "GEN", method="DL",
                data=shit, 
                subset=(Inhibitor=="NBPT"))
res.nbdc <- rma(yi = mean.diff, vi = diff.var, 
                measure = "GEN", method="DL",
                data=shit, 
                subset=(Inhibitor=="NBPT+DCD"))
res.nitr <- rma(yi = mean.diff, vi = diff.var, 
                measure = "GEN", method="DL",
                data=shit, 
                subset=(Inhibitor=="nitrapyrin"))
#res.pcf <- rma(yi = mean.diff, vi = diff.var, 
#              measure = "GEN", method="DL",
#             data=shit, 
#            subset=(Inhibitor=="PCF"))

#res.nd <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR",
#            subset=(alloc=="random"), method="REML")
#res.n <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR",
#            subset=(alloc=="alternate"), method="REML")


#addpoly(res.sr, row=br5-5-1, cex=.75, mlab="RE Model for Subgroup")
addpoly(res.nbdc, row=br2-4, cex=.75, mlab="RE Model for Subgroup")
addpoly(res.nitr, row=br5-3, cex=.75, mlab="RE Model for Subgroup")
addpoly(res.nbpt, row=br1-2, cex=.75, mlab="RE Model for Subgroup")











