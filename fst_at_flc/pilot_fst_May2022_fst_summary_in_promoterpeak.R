# knitr::opts_knit$set(root.dir = "~/safedata/ath_evo/grenepilot")
setwd("~/safedata/ath_evo/grenepilot")
# setwd("~/grenepilot")


library(moiR)
library(ggplot2)
# library(ggpmisc)
library(cowplot)
theme_set(theme_cowplot())
library(RColorBrewer)
library(tidyverse)
library(data.table)
# devtools::install("~/safedata/genemaps/")
# library(genemaps)
# devtools::load_all(".")

# startflc=3167318
# endflc= 3185514

beforeflc=3167318
afterflc= 3185514
startflc= 3173382
endflc = 3179448

#######################################################################################################################
## Load Fst
# Read the table of all pairwise F_ST
# tab=read.csv(file = "data-raw/fst_pairs.csv", sep = "\t")
# tab=read.csv(file = "fst-width-1-FLC.csv", sep = "\t")
# tab=read.csv("analyses/fst-spence-hudson.csv", sep = "\t")
tab=read.csv("analyses/fst-spence-nei.csv", sep = "\t")
# head(tab)

# Make long format.
tab_long <-  tab %>%
  dplyr::select(-SNPS, -END) %>%
  rename(chr=CHROM, pos=START) %>%
  pivot_longer(cols=starts_with("S"), values_to='S', names_to='time')
# tab_long

tab_long$S[tab_long$S<0] <-0
#######################################################################################################################
## load sample info
# samples=read.csv("Table_IDseq.tsv", sep = "\t")
samples=read.csv("data-raw/Table_IDseq.tsv", sep = "\t")
samplerep<-c(0,0,0,1,1,1,2,2,2,3,3,3)

# head(samples)

samples$ID_Sequence<-paste0("S",samples$ID_Sequence)

timetranslator<-samplerep
names(timetranslator)<-paste0("S",1:12)

codesofinterest<-expand.grid(samples$ID_Sequence,samples$ID_Sequence) %>%
  dplyr::filter(Var1!=Var2) %>%
  dplyr::filter(Var1 %in% c("S1","S2","S3") | Var2 %in% c("S1","S2","S3")) %>%
  dplyr::mutate(codes=paste0(Var1,".",Var2)) %>%
  dplyr::mutate(comparison=ifelse(Var1 %in% c("S1","S2","S3"), Var2,Var1 )) %>%
  dplyr::mutate(realtime=  timetranslator[Var2]) %>%
  dplyr::filter(comparison>3)

#####
# more targeted analysis, select time 0 to all and create the mean Fst in promoter
#
# selcol = c("S4.S10", "S5.S11", "S6.S12")
selcol = c(
  "S1.S4","S1.S7","S1.S10",
  "S2.S5","S2.S8","S2.S11",
  "S3.S6","S3.S9","S3.S12"
)
# selcol =  "S2.S5"
# select timepoint 0 to all others
tmpfst<-subset(tab_long, time %in% selcol)
tmpfst$SNP<-paste0(tmpfst$chr, "_",tmpfst$pos)
# remove NAs
mfst<-na.omit(tmpfst)
# annotate promoter region broadly
mfst<- mfst %>%
  dplyr::mutate(promotor= pos>3179448-500  & pos < 3179448+2500 ) %>%
  dplyr::mutate(rep=substr(time,2,2))
# show there is a wilcox text difference with all of them
message("Computer Wilcox test and t test of Fst values within (TRUE) and outside (FALSE) promoter region")
wilcox.test(mfst$S ~ mfst$promotor)
t.test(mfst$S ~ mfst$promotor)
message("Computer the 95% percentile within and outside")
message("Within promoter: ", quantile(mfst$S[mfst$promotor==TRUE], probs = .95 ))
message("Outside promoter: ", quantile(mfst$S[mfst$promotor==FALSE], probs = .95 ))
