#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}
fst_file=args[1]

#setwd("~/safedata/ath_evo/grenepilot")
# setwd("~/grenepilot")

fn <-
function (data.frame) 
{
    unlist(as.numeric(as.matrix(data.frame)))
}
# library(moiR)
library(ggplot2)
# library(ggpmisc)
library(cowplot)
# theme_set(theme_cowplot())
library(RColorBrewer)
library(tidyverse)
library(data.table)
# devtools::install("~/safedata/genemaps/")
# library(genemaps)
# devtools::load_all(".") 


## Load Fst


# Read the table of all pairwise F_ST
# tab=read.csv(file = "data-raw/fst_pairs.csv", sep = "\t")
#tab=read.csv(file = "fst-width-1-FLC.csv", sep = "\t")
tab=read.csv(file = fst_file, sep = "\t")
head(tab)

# Make long format.
tab_long <-  tab %>% 
  dplyr::select(-SNPS, -END) %>%
  rename(chr=CHROM, pos=START) %>%
  pivot_longer(cols=starts_with("S"), values_to='S', names_to='time')
tab_long

## load sample info
samples=read.csv("Table_IDseq.tsv", sep = "\t")
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
  geom_vline(xintercept=3173382, color='grey', lty='dotted')+
  geom_vline(xintercept=3179448, color='grey', lty='dotted')+
  geom_point(data=toplot,
             aes(x = pos, y=S, 
                 color=factor(realtime)) ) + 
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
  theme(plot.subtitle = element_text(hjust = 0.5))+
  ylim(0,1)

p0



## make plot of tiempoint 1 vs timepoint 3

# Somehow, R reshape tries to be smart and removes the initial "S" from the column names.
# Well, okay, not learning how to trick it into being less smart today, so we just roll with it,
# and select the columns we want without the initial "S"...
selcol = c("S4.S10", "S5.S11", "S6.S12")

# Plot. Missing: vertical lines to mark the beginning and end of the FLC locus.
# Those need to be added at 3173382 and 3179448, but somehow that command got lost in my R history.
p1<-ggplot(subset(tab_long, time %in% selcol))+
  aes(x = pos, y=S, color=time) + 
  geom_vline(xintercept=3173382, color='grey', lty='dotted')+
  geom_vline(xintercept=3179448, color='grey', lty='dotted')+
  geom_point() + 
  scale_color_manual(labels=c("R1", "R2", "R3"),
                     values=wesanderson::wes_palettes$GrandBudapest1) +
  labs(y="Fst (T1-T3)", 
       x='',
       color='E&R plot\nreplicates', 
       subtitle =  "") +
  theme(plot.subtitle = element_text(hjust = 0.5))+
  ylim(0,1)

p1



  
# combined<-plot_grid(p1,p2,ncol=1,align = 'hv')
combined<-plot_grid(p0,p1,ncol=1,align = 'hv')
combined
save_plot(
  # "~/safedata/ath_evo/grenepilot//figs/preliminary-pilot_fst_floweringeffect.pdf",
  paste0(fst_file,"-pilot_fst.png"),
  plot = combined,
  base_width = 12,base_height = 6)



