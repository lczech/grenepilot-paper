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
tab=read.csv(file = "fst-width-1-FLC.csv", sep = "\t")
tab=read.csv("analyses/fst-spence-hudson.csv", sep = "\t")
tab=read.csv("analyses/fst-spence-nei.csv", sep = "\t")
head(tab)

# Make long format.
tab_long <-  tab %>%
  dplyr::select(-SNPS, -END) %>%
  rename(chr=CHROM, pos=START) %>%
  pivot_longer(cols=starts_with("S"), values_to='S', names_to='time')
tab_long

tab_long$S[tab_long$S<0] <-0
#######################################################################################################################
## load sample info
samples=read.csv("Table_IDseq.tsv", sep = "\t")
samples=read.csv("data-raw/Table_IDseq.tsv", sep = "\t")
samplerep<-c(0,0,0,1,1,1,2,2,2,3,3,3)

head(samples)

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


#######################################################################################################################
## make plot of tiempoint 0 vs everything else
# samples=read.csv("data-raw/Table_IDseq.tsv", sep = "\t")
# head(samples)

selcol<-codesofinterest$codes

tab_long2<-merge(x=tab_long,y=codesofinterest,by.x='time', by.y='codes')

# Plot. Missing: vertical lines to mark the beginning and end of the FLC locus.
# Those need to be added at 3173382 and 3179448, but somehow that command got lost in my R history.
toplot<-subset(tab_long2, (time %in% selcol) & (realtime >0) )

p0<-
  ggplot()+
  geom_vline(xintercept=startflc, color='grey', lty='dotted')+
  geom_vline(xintercept=endflc, color='grey', lty='dotted')+
  geom_point(data=toplot,
             aes(x = pos, y=S,
                 color=factor(realtime))
  ) +
  # stat_summary(data=toplot, fun=max,geom='line',
  #           aes(x = round(pos/500)*500, y=S,
  #               # shape=factor(samplerep[comparison] ),
  #               color=factor(realtime)) ) +
  scale_color_manual(
    labels=c("T1", "T2", "T3"),
    values=wesanderson::wes_palettes$GrandBudapest2) +
  labs(y="Fst (seed-flo)",
       x='',
       # color='E&R plot\nreplicates',
       color='time sample',
       subtitle =  "Flowering Locus C (FLC)"#,
       # shape='time'
  ) +
  theme(plot.subtitle = element_text(hjust = 0.5))

fstseedsflowers<-p0
fstseedsflowers
save(file = "tmpobjects/fstseedsflowers.rda",fstseedsflowers)
# save(file = "figs/fst-hudson",fstseedsflowers)

################################################################################
## make plot of tiempoint 1 vs timepoint 3

# Somehow, R reshape tries to be smart and removes the initial "S" from the column names.
# Well, okay, not learning how to trick it into being less smart today, so we just roll with it,
# and select the columns we want without the initial "S"...
selcol = c("S4.S10", "S5.S11", "S6.S12")

# Plot. Missing: vertical lines to mark the beginning and end of the FLC locus.
# Those need to be added at 3173382 and 3179448, but somehow that command got lost in my R history.
p1<-ggplot(subset(tab_long, time %in% selcol))+
  aes(x = pos, y=S, color=time) +
  geom_vline(xintercept=startflc, color='grey', lty='dotted')+
  geom_vline(xintercept=endflc, color='grey', lty='dotted')+
  geom_point() +
  scale_color_manual(labels=c("R1", "R2", "R3"),
                     values=wesanderson::wes_palettes$GrandBudapest1) +
  labs(y="Fst (T1-T3)",
       x='',
       color='E&R plot\nreplicates',
       subtitle =  "") +
  theme(plot.subtitle = element_text(hjust = 0.5))
p1
fstt1t3<-p1

save(file = "tmpobjects/fstt1t3.rda",fstt1t3)

# p1+
#   geom_vline(xintercept=3179448-500, color='grey', lty='dotted')+
#   geom_vline(xintercept=3179448+2500, color='grey', lty='dotted')


#######################################################################################################################
## load gwa on flowering time or fitness
library(data.table)
# flo<-fread("~/safedata/natvar/phenotypes/1001_Consortium_Cell_2016_PID_27293186/gemma/FT10/output/1001_Consortium_Cell_2016_PID_27293186.assoc.txt")
# flo<-fread("~/safedata/natvar/phenotypes/1001_Consortium_Cell_2016_PID_27293186/1001/FT10/output/1001_Consortium_Cell_2016_PID_27293186.lmm.assoc.txt")
# flo<-fread("~/safedata/natvar/phenotypes/Exposito-Alonso_Nature_2019_PID_31462776/1001/rSeeds_thp/output/Exposito-Alonso_Nature_2019_PID_31462776.lmm.assoc.txt")
# flo<-fread("~/safedata/natvar/phenotypes/Exposito-Alonso_Nature_2019_PID_31462776/1001/rSeeds_thp/output/Exposito-Alonso_Nature_2019_PID_31462776norm.lm.assoc.txt")
# flo<-fread("~/safedata/natvar/phenotypes/Exposito-Alonso_Nature_2019_PID_31462776/gemma/rFitness_thi/output/Exposito-Alonso_Nature_2019_PID_31462776.lmm.assoc.txt")
# flo<-fread("~/safedata/natvar/phenotypes/Exposito-Alonso_Nature_2019_PID_31462776/gemma/rFitness_thp/output/Exposito-Alonso_Nature_2019_PID_31462776.lmm.assoc.txt")

flo<-fread("~/safedata/natvar/phenotypes/Exposito-Alonso_Nature_2019_PID_31462776/gemma/rSeeds_thp/output/Exposito-Alonso_Nature_2019_PID_31462776.lmm.assoc.txt")
# flo<-fread("~/grenepilot/gwas/rSeeds_thi/Exposito-Alonso_Nature_2019_PID_31462776.lmm.assoc.txt")
# flo<-fread("~/grenepilot/gwas/rSeeds_thp/Exposito-Alonso_Nature_2019_PID_31462776.lmm.assoc.txt")
# flo<-fread("~/grenepilot/gwas/rSeeds_thp/Exposito-Alonso_Nature_2019_PID_31462776.assoc.txt")
tmp<-dplyr::filter(flo, chr == fn(tab_long$chr)[1], ps>min(tab_long$pos), ps<max(tab_long$pos) )
head(tmp)

p2<-ggplot(tmp)+
  geom_vline(xintercept=startflc, color='grey', lty='dotted')+
  geom_vline(xintercept=endflc, color='grey', lty='dotted')+
  geom_point(aes(x = ps, y=-log10(p_score)), size=2, color='black')+
  geom_point(aes(x = ps, y=-log10(p_score)), size=1.5, color='white')+
  geom_point(aes(x = ps, y=-log10(p_score), color=abs(beta)),size=1)+
  scale_color_gradientn(colours = brewer.pal(n=9,name="Greys")) +
  labs(y="-log10 P", x='genomic position (bp)',color='effect in relative\nseed production') +
  theme(plot.subtitle = element_text(hjust = 0.5))
p2
# startflc=3167318
# endflc= 3185514
flcgwa<-p2
save(file = "tmpobjects/flcgwa.rda",flcgwa)

#######################################################################################################################
load("tmpobjects/flcgwa.rda")
load("tmpobjects/fstt1t3.rda")
load("tmpobjects/fstseedsflowers.rda")
load("tmpobjects//averagefreqchange.rda")
# save
# # combined<-plot_grid(p1,p2,ncol=1,align = 'hv')
# combined<-plot_grid(p0,p1,p2,ncol=1,align = 'hv')
# combined
# save_plot(
#   # "~/safedata/ath_evo/grenepilot//figs/preliminary-pilot_fst_floweringeffect.pdf",
#   "~/grenepilot/preliminary-pilot_fst_floweringeffect.pdf",
#   plot = combined,
#   base_width = 12,base_height = 6)

# combined<-plot_grid(p1,p2,ncol=1,align = 'hv')
combined<-plot_grid(
                    fstseedsflowers + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+xlim(beforeflc,afterflc)   ,
                    fstt1t3 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm") )+xlim(beforeflc,afterflc) ,
                    averagefreqchange + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +xlim(beforeflc,afterflc) ,
                    flcgwa + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +xlim(beforeflc,afterflc) ,
                    ncol=1,align = 'hv')
combined
save_plot(
  # "~/safedata/ath_evo/grenepilot//figs/preliminary-pilot_fst_floweringeffect.pdf",
  # "~/grenepilot/preliminary-pilot_fst_floweringeffect.pdf",
  "preliminary-pilot_fst_floweringeffect_withfounders.svg",
  plot = combined,
  base_width = 12,base_height = 6)

