setwd("/Users/andrewtrlica/Desktop/IPNI meta/")


library(metafor)
urea.dat <- read.csv(file="lumped/urea.effect.sizes.mean-diff.csv")
urea.dat$Inhibitor <- as.character(urea.dat$Inhibitor)
urea.dat$Name <- paste(urea.dat$Study, urea.dat$Inhibitor) ## mark each study effect size with its Inhibitor

### this loop comes up with a scaling factor for the squares to apply within each Inhibitor group
inhib <- unique(urea.dat$Inhibitor)
shit <- data.frame(NULL)
for(j in 1:length(inhib)){
  gunk <- urea.dat[which(urea.dat$Inhibitor==inhib[j]),]
  total <- sum(gunk$Inhibitor==inhib[j])
  skim <- gunk[which(gunk$Inhibitor==inhib[j]),8]/sum(gunk[which(gunk$Inhibitor==inhib[j]),8])
  #place <- rank(skim)
  #piece <- 1.5/(length(place)-1)
  #rwt <- 2-((place-1)*piece)
  gunk$rewt <- 1+(1*skim)
  shit <- rbind(shit, gunk)
}

### do the rma meta analysis to get the summary effect polygon and basic frame
urea.shit <- rma(yi = shit$mean.diff, vi = shit$diff.var, 
               measure = "GEN", method="DL",
               slab = shit$Name)

##### the forest plot needs to be specified which rows for placing the effect sizes; 
##### you create gaps in the rows to provide room for titles and summary polygons;
##### the rows count up from the bottom of the figure (presume the all-study summary is on row 1)
spc <- 3
br1 <- 21 # PCF, needs 18
br2 <- spc+6+br1 # Nitrap, needs 5
br3 <- spc+2+br2 # NBPT, needs 1
br4 <- spc+9+br3 # NBPT+DCD, needs 8
br5 <- spc+6+br4 # S.R., needs 5


par(mar=c(4,4,1,2),  font=1)
forest(urea.shit, xlim=c(-8, 8), cex=.65, ylim=c(-1, 63), 
       rows=c((br1-17):br1,(br2-4):br2,br3,(br4-7):br4,(br5-4):br5),
       xlab="Mean difference", mlab="RE Model for All Studies", psize=shit$rewt)
text(0, 64, "Urea", font=2)



### set font expansion factor (as in forest() above) and use bold italic
### font and save original settings in object 'op'
#mtext("Urea")
par(cex=.75, font=4)

### add text for the subgroups
text(-8, cex = 1, c(br5+1.5,br4+1.5,br3+1.5,br2+1.5,br1+1.5), pos=4, 
     c("S.R.",
       "NBPT+DCD",
       "NBPT",
       "Nitrapyrin",
       "PCF"))



#### add subgroup polygons for summary effects by doing subgroup meta analysis and adding diamonds in the row spaces
res.sr <- rma(yi = mean.diff, vi = diff.var, 
             measure = "GEN", method="DL",
             data=shit, 
             subset=(Inhibitor=="S.R"))
#res.nbpt <- rma(yi = mean.diff, vi = diff.var, 
 #            measure = "GEN", method="DL",
  #           data=shit, 
   #          subset=(Inhibitor=="NBPT"))
res.nbdc <- rma(yi = mean.diff, vi = diff.var, 
             measure = "GEN", method="DL",
             data=shit, 
             subset=(Inhibitor=="NBPT+DCD"))
res.nitr <- rma(yi = mean.diff, vi = diff.var, 
             measure = "GEN", method="DL",
             data=shit, 
             subset=(Inhibitor=="nitrapyrin"))
res.pcf <- rma(yi = mean.diff, vi = diff.var, 
             measure = "GEN", method="DL",
             data=shit, 
             subset=(Inhibitor=="PCF"))

#res.nd <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR",
 #            subset=(alloc=="random"), method="REML")
#res.n <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR",
 #            subset=(alloc=="alternate"), method="REML")

par(cex=1.5, font = 2)
addpoly(res.sr, row=br5-5.5, mlab="RE Model for Subgroup")
addpoly(res.nbdc, row=br4-8.5, mlab="RE Model for Subgroup")
addpoly(res.nitr, row=br2-5.5,  mlab="RE Model for Subgroup")
addpoly(res.pcf, row=br1-18.5,  mlab="RE Model for Subgroup")
