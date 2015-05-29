setwd("/Users/andrewtrlica/Desktop/IPNI_meta/parsed/")
library(metafor)

###### AA funnel plots and trim and fill results
AA.eff <- read.csv("AA pairs-diff.csv")
View(AA.eff)
AA.eff$Name <- paste(AA.eff$study, AA.eff$X)

AA.res <- rma(yi=effect.size, vi=var..effect.size, data=AA.eff, measure = "GEN", method="DL",
              slab = Name )
AA.taf <- trimfill(AA.res)
AA.taf

funnel(AA.res, main="AA+nitrapyrin")
funnel(AA.taf, main="AA+nitrapyrin, trim and fill")

#addpoly(AA.res, row=30, cex=.75, mlab="RE Model for Subgroup")

###### UAN funnel plots and trim and fill results
UAN.eff <- read.csv("UAN pairs-diff.csv")
UAN.eff$Name <- paste(UAN.eff$study, UAN.eff$source, UAN.eff$inhibitor, UAN.eff$X)
UAN.res.NBPT <- rma(yi=effect.size, vi=var..effect.size, data=UAN.eff[which(UAN.eff$inhibitor=="NBPT"),], measure = "GEN", method="DL",
              slab = Name )

UAN.taf.NBPT <- trimfill(UAN.res.NBPT)
UAN.taf.NBPT
funnel(UAN.res.NBPT, main="UAN+NBPT")
funnel(UAN.taf.NBPT, main="UAN+NBPT, trim and fill")



###### Urea funnel plots and trim and fill results
urea.eff <- read.csv("urea pairs-diff.csv")
urea.eff$Name <- paste(urea.eff$study, urea.eff$source, urea.eff$inhibitor, urea.eff$X)
urea.res.NBPT <- rma(yi=effect.size, vi=var..effect.size, data=urea.eff[which(urea.eff$inhibitor=="NBPT"),], measure = "GEN", method="DL",
                    slab = Name )

urea.taf.NBPT <- trimfill(urea.res.NBPT)
urea.taf.NBPT
funnel(urea.res.NBPT, main="urea+NBPT")
funnel(urea.taf.NBPT, main="urea+NBPT, trim and fill")

urea.res.PCF <- rma(yi=effect.size, vi=var..effect.size, data=urea.eff[which(urea.eff$inhibitor=="PCF"),], measure = "GEN", method="DL",
                     slab = Name )

urea.taf.PCF <- trimfill(urea.res.PCF)
urea.taf.PCF
funnel(urea.res.NBPT, main="urea+PCF")
funnel(urea.taf.NBPT, main="urea+PCF, trim and fill")


