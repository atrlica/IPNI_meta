setwd("/Users/andrewtrlica/Desktop/IPNI meta/")


library(metafor)
AA.dat <- read.csv(file = "lumped/AA.effect.sizes.mean-diff.csv")
urea.dat <- read.csv(file="lumped/urea.effect.sizes.mean-diff.csv")
UAN.dat <- read.csv(file="lumped/UAN.effect.sizes.mean-diff.csv")


#### the following bmeta analysis results and corresponding forest plots for lumped values -- use as approach to making other forest plots

### AA lumped effects
AA.dat$Study <- as.character(AA.dat$Study)
AA.meta <- rma(yi = AA.dat$mean.diff, vi = AA.dat$diff.var, 
              measure = "GEN", method="DL",
              slab = AA.dat$Study)
AA.meta
forest(AA.meta,refline=0, xlab="Mean difference (95%CI)", mlab="Summary Effect")
mtext("AA")


### urea lumped effect
urea.dat$Study <- as.character(urea.dat$Study)
urea.meta <- rma(yi = urea.dat$mean.diff, vi = urea.dat$diff.var, 
                 measure = "GEN", method="DL",
                 slab = urea.dat$Study)

urea.meta
forest(urea.meta, refline=0, xlab="Mean difference (95%CI)", mlab="Summary Effect")
mtext("Urea")

urea1 <- urea.dat[which(urea.dat$Inhibitor=="PCF"),]
urea1.m <- rma(yi = urea1$mean.diff, vi = urea1$diff.var, 
                            measure = "GEN", method="DL",
                            slab = urea1$Study)
forest(urea1.m, refline=0, xlab="Mean difference (95%CI)", mlab="Summary Effect", )
mtext("Urea + PCF")

urea.dat$Inhibitor <- as.character(urea.dat$Inhibitor)
urea.dat$Name <- paste(urea.dat$Study, urea.dat$Inhibitor) ## mark each study effect size with its Inhibitor
urea2.m <- rma(yi = urea.dat$mean.diff, vi = urea.dat$diff.var, 
               measure = "GEN", method="DL",
               slab = urea.dat$Name)
urea2.m
forest(urea2.m, refline=0, xlab="Mean difference (95% CI)", mlab="Summary Effect")


forest(urea2.m, xlim=c(-6, 6), cex=.75, ylim=c(-1, 63), rows=c(3:20,25:29,34,39:46,51:55),
       xlab="Mean difference", mlab="RE Model for All Studies", psize=.5)






urea.dat$Inhibitor <- as.character(urea.dat$Inhibitor)
urea.dat$Name <- paste(urea.dat$Study, urea.dat$Inhibitor) ## mark each study effect size with its Inhibitor



inhib <- unique(urea.dat$Inhibitor)
shit <- data.frame(NULL)
for(j in 1:length(inhib)){
  gunk <- urea.dat[which(urea.dat$Inhibitor==inhib[j]),]
  total <- sum(gunk$Inhibitor==inhib[j])
  skim <- gunk[which(gunk$Inhibitor==inhib[j]),8]/sum(gunk[which(gunk$Inhibitor==inhib[j]),8])
  #place <- rank(skim)
  #piece <- 1.5/(length(place)-1)
  #rwt <- 2-((place-1)*piece)
  gunk$rewt <- 1+(2*skim)
  shit <- rbind(shit, gunk)
}

urea.shit <- rma(yi = shit$mean.diff, vi = shit$diff.var, 
               measure = "GEN", method="DL",
               slab = shit$Name)


par(mar=c(4,4,1,2))

br1 <- 20 # PCF, needs 18
br2 <- 29 # Nitrap, needs 5
br3 <- 34 # NBPT, needs 1
br4 <- 46 # NBPT+DCD, needs 8
br5 <- 55 # S.R., needs 5



forest(urea.shit, xlim=c(-8, 8), cex=.7, ylim=c(-1, 63), 
       rows=c((br1-17):br1,(br2-4):br2,br3,(br4-7):br4,(br5-4):br5),
       xlab="Mean difference", mlab="RE Model for All Studies", psize=shit$rewt)
### set font expansion factor (as in forest() above) and use bold italic
### font and save original settings in object 'op'
mtext("Urea")
op <- par(cex=.75, font=4)

### add text for the subgroups
text(-8, cex = 1, c(56.5,47.5,35.5,30.5,21.5), pos=4, c("S.R.",
                               "NBPT+DCD",
                               "NBPT",
                               "Nitrapyrin",
                               "PCF"))
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


addpoly(res.sr, row=49.5, cex=.75, mlab="RE Model for Subgroup")
addpoly(res.nbdc, row=37.5, cex=.75, mlab="RE Model for Subgroup")
addpoly(res.nitr, row=23.5, cex=.75, mlab="RE Model for Subgroup")
addpoly(res.pcf, row=1.5, cex=.75, mlab="RE Model for Subgroup")


############
############
####### USE TO DO MULTIPLE SUBGROUPS
### to save as png file
png(filename="forest_plot_with_subgroups.png",
    res=95, width=680, height=680, type="cairo")

### decrease margins so the full space is used
par(mar=c(4,4,1,2))

### load BCG vaccine data
data(dat.bcg)

### fit random-effects model (use slab argument to define study labels)
res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR",
           slab=paste(author, year, sep=", "), method="REML")

### set up forest plot (with 2x2 table counts added; rows argument is used
### to specify exactly in which rows the outcomes will be plotted)
forest(res, xlim=c(-16, 6), at=log(c(.05, .25, 1, 4)), atransf=exp,
       ilab=cbind(dat.bcg$tpos, dat.bcg$tneg, dat.bcg$cpos, dat.bcg$cneg),
       ilab.xpos=c(-9.5,-8,-6,-4.5), cex=.75, ylim=c(-1, 27),
       order=order(dat.bcg$alloc), rows=c(3:4,9:15,20:23),
       xlab="Relative Risk", mlab="RE Model for All Studies", psize=1)

?forest
### set font expansion factor (as in forest() above) and use bold italic
### font and save original settings in object 'op'
op <- par(cex=.75, font=4)

### add text for the subgroups
text(-16, c(24,16,5), pos=4, c("Systematic Allocation",
                               "Random Allocation",
                               "Alternate Allocation"))
### switch to bold font
par(font=2)

### add column headings to the plot
text(c(-9.5,-8,-6,-4.5), 26, c("TB+", "TB-", "TB+", "TB-"))
text(c(-8.75,-5.25),     27, c("Vaccinated", "Control"))
text(-16,                26, "Author(s) and Year",     pos=4)
text(6,                  26, "Relative Risk [95% CI]", pos=2)

### set par back to the original settings
par(op)

### fit random-effects model in the three subgroups
res.s <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR",
             subset=(alloc=="systematic"), method="REML")
res.r <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR",
             subset=(alloc=="random"), method="REML")
res.a <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR",
             subset=(alloc=="alternate"), method="REML")

### add summary polygons for the three subgroups
addpoly(res.s, row=18.5, cex=.75, atransf=exp, mlab="RE Model for Subgroup")
addpoly(res.r, row= 7.5, cex=.75, atransf=exp, mlab="RE Model for Subgroup")
addpoly(res.a, row= 1.5, cex=.75, atransf=exp, mlab="RE Model for Subgroup")

dev.off()








### UAN lumped effects
UAN.dat$Study <- as.character(UAN.dat$Study)
UAN.meta <- rma(yi = UAN.dat$mean.diff, vi = UAN.dat$diff.var, 
                measure = "GEN", method="REML",
                slab = UAN.dat$Study)

UAN.meta
forest(UAN.meta, refline=0, xlab="Mean difference (95%CI)", mlab="Summary Effect")
mtext("UAN")


View(UAN.dat)




