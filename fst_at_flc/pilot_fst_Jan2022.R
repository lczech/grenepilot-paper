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

startflc=3167318
endflc= 3185514
#######################################################################################################################
## Load Fst
# Read the table of all pairwise F_ST
# tab=read.csv(file = "data-raw/fst_pairs.csv", sep = "\t")
tab=read.csv(file = "fst-width-1-FLC.csv", sep = "\t")
head(tab)

# Make long format.
tab_long <-  tab %>%
  dplyr::select(-SNPS, -END) %>%
  rename(chr=CHROM, pos=START) %>%
  pivot_longer(cols=starts_with("S"), values_to='S', names_to='time')
tab_long

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
  theme(plot.subtitle = element_text(hjust = 0.5))

fstseedsflowers<-p0
save(file = "tmpobjects/fstseedsflowers.rda",fstseedsflowers)


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
  theme(plot.subtitle = element_text(hjust = 0.5))
p1
fstt1t3<-p1
save(file = "tmpobjects/fstt1t3.rda",fstt1t3)

p1+
  geom_vline(xintercept=3179448-500, color='grey', lty='dotted')+
  geom_vline(xintercept=3179448+2500, color='grey', lty='dotted')


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
  geom_vline(xintercept=3173382, color='grey', lty='dotted')+
  geom_vline(xintercept=3179448, color='grey', lty='dotted')+
  geom_point(aes(x = ps, y=-log10(p_score)), size=2, color='black')+
  geom_point(aes(x = ps, y=-log10(p_score)), size=1.5, color='white')+
  geom_point(aes(x = ps, y=-log10(p_score), color=abs(beta)),size=1)+
  scale_color_gradientn(colours = brewer.pal(n=9,name="Greys")) +
  labs(y="-log10 P", x='genomic position (bp)',color='effect in relative\nseed production') +
  theme(plot.subtitle = element_text(hjust = 0.5))
p2
startflc=3167318
endflc= 3185514
flcgwa<-p2
save(file = "tmpobjects/flcgwa.rda",flcgwa)

#######################################################################################################################
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
combined<-plot_grid(fstseedsflowers,fstt1t3,flcgwa,ncol=1,align = 'hv')
combined
save_plot(
  # "~/safedata/ath_evo/grenepilot//figs/preliminary-pilot_fst_floweringeffect.pdf",
  # "~/grenepilot/preliminary-pilot_fst_floweringeffect.pdf",
  "~/grenepilot/preliminary-pilot_fst_floweringeffect_withfounders.pdf",
  plot = combined,
  base_width = 12,base_height = 6)


#######################################################################################################################
## load gwa on flowering time or fitness
# Load all fitness
# codesnatpaper<-c("thp","thi","tlp","tli")
codesnatpaper<-c("thp","thi","mlp", "mli")
phenotypes<-c("rSeeds","rSurvival_fruit")
allcombs<-expand.grid(codesnatpaper, phenotypes)

### do the loading

tmpall<- lapply(1:nrow(allcombs), function(i){
  pheno<-allcombs$Var2[i]
  env<-allcombs$Var1[i]
  gwafile<-paste0("~/natvar/phenotypes/Exposito-Alonso_Nature_2019_PID_31462776/1001/",pheno,"_",env,
                  "/output/Exposito-Alonso_Nature_2019_PID_31462776.lmm.assoc.txt")
  print(gwafile)
  flo<-fread(gwafile)
  tmp<-dplyr::filter(flo, chr == fn(tab_long$chr)[1], ps>min(tab_long$pos), ps<max(tab_long$pos) )
  tmp$condition<- paste0(pheno,"_",env)
  return(tmp)
})
tmpall<- tmpall %>% do.call(rbind, .)
# save(file = "tmpobjects/flc_gwas_1001g_exposito-alonso.rda",tmpall)
# load(file = "tmpobjects/flc_gwas_1001g_exposito-alonso.rda")

# tmpall<- lapply(1:nrow(allcombs), function(i){
#   pheno<-allcombs$Var2[i]
#   env<-allcombs$Var1[i]
#   gwafile<-paste0("~/natvar/phenotypes/Exposito-Alonso_Nature_2019_PID_31462776/1001/",pheno,"_",env,
#                   "/output/Exposito-Alonso_Nature_2019_PID_31462776norm.lmm.assoc.txt")
#   print(gwafile)
#   flo<-fread(gwafile)
#   tmp<-dplyr::filter(flo, chr == fn(tab_long$chr)[1], ps>min(tab_long$pos), ps<max(tab_long$pos) )
#   tmp$condition<- paste0(pheno,"_",env)
#   return(tmp)
# })
# tmpall<- tmpall %>% do.call(rbind, .)
# save(file = "tmpobjects/flc_gwas_1001g_exposito-alonso_norm.rda",tmpall)
# load(file = "tmpobjects/flc_gwas_1001g_exposito-alonso_norm.rda")

tmpall<- lapply(1:nrow(allcombs), function(i){
  pheno<-allcombs$Var2[i]
  env<-allcombs$Var1[i]
  gwafile<-paste0("~/natvar/phenotypes/Exposito-Alonso_Nature_2019_PID_31462776/gemma/",pheno,"_",env,"/output/Exposito-Alonso_Nature_2019_PID_31462776.lmm.assoc.txt")
  print(gwafile)
  flo<-fread(gwafile)
  tmp<-dplyr::filter(flo, chr == fn(tab_long$chr)[1], ps>min(tab_long$pos), ps<max(tab_long$pos) )
  tmp$condition<- paste0(pheno,"_",env)
  return(tmp)
})
tmpall<- tmpall %>% do.call(rbind, .)
# save(file = "tmpobjects/flc_gwas_1001g_exposito-alonso_imputed.rda",tmpall)
# load(file = "tmpobjects/flc_gwas_1001g_exposito-alonso_imputed.rda")


# codesnatpaper<-c("thp","thi")
# phenotypes<-c("rFitness","rSeeds","rSurvival_fruit")
# phenotypes<-c("rSeeds","rSurvival_fruit")
# allcombs<-expand.grid(codesnatpaper, phenotypes)
# tmpall<- lapply(1:nrow(allcombs), function(i){
#   pheno<-allcombs$Var2[i]
#   env<-allcombs$Var1[i]
#   gwafile<-paste0("~/natvar/phenotypes/Exposito-Alonso_Nature_2019_PID_31462776/1001/",pheno,"_",env,
#                   "/output/Exposito-Alonso_Nature_2019_PID_31462776.lm.assoc.txt")
#   print(gwafile)
#   flo<-fread(gwafile)
#   tmp<-dplyr::filter(flo, chr == fn(tab_long$chr)[1], ps>min(tab_long$pos), ps<max(tab_long$pos) )
#   tmp$condition<- paste0(pheno,"_",env)
#   return(tmp)
# })
# tmpall<- tmpall %>% do.call(rbind, .)
# save(file = "tmpobjects/flc_gwas_1001g_exposito-alonso_nokinship.rda",tmpall)



phenotypes<-c("FT")
allcombs<-expand.grid(codesnatpaper, phenotypes)
tmpall<- lapply(1:nrow(allcombs), function(i){
  pheno<-allcombs$Var2[i]
  env<-allcombs$Var1[i]
  gwafile<-paste0("~/natvar/phenotypes/Exposito-Alonso_Nature_2019_PID_31462776/gemma/",pheno,"_",env,
                  "/output/Exposito-Alonso_Nature_2019_PID_31462776.lmm.assoc.txt")
  print(gwafile)
  flo<-fread(gwafile)
  tmp<-dplyr::filter(flo, chr == fn(tab_long$chr)[1], ps>min(tab_long$pos), ps<max(tab_long$pos) )
  tmp$condition<- paste0(pheno,"_",env)
  return(tmp)
})
tmpall<- tmpall %>% do.call(rbind, .)
save(file = "tmpobjects/flc_gwas_1001g_exposito-alonso_flowering.rda",tmpall)

### plot!
# px<-
  # ggplot(tmpall)+
    ggplot(ta)+
    geom_vline(xintercept=3173382, color='grey', lty='dotted')+
  geom_vline(xintercept=3179448, color='grey', lty='dotted')+
  geom_point(aes(x = ps, y=-log10(p_score)), size=2, color='black')+
  geom_point(aes(x = ps, y=-log10(p_score)), size=1.5, color='white')+
  geom_point(aes(x = ps, y=-log10(p_score), color=abs(beta)),size=1)+
  scale_color_gradientn(colours = brewer.pal(n=9,name="Greys")) +
  labs(y="-log10 P", x='genomic position (bp)',color='effect in relative\nseed production') +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  facet_wrap(~condition, ncol=1)
### plot!


# plot of seed gwa imputed with Madrid and Tuebinen
load(file = "tmpobjects/flc_gwas_1001g_exposito-alonso_imputed.rda")
pdf("figs/flc_gwas_1001g_exposito-alonso_imputed.pdf",height = 8,width = 12)
tmpall %>% dplyr::filter(grepl("Seed", condition)) %>%
ggplot(.)+
  geom_vline(xintercept=3173382, color='grey', lty='dotted')+
  geom_vline(xintercept=3179448, color='grey', lty='dotted')+
  geom_point(aes(x = ps, y=-log10(p_score)), size=2, color='black')+
  geom_point(aes(x = ps, y=-log10(p_score)), size=1.5, color='white')+
  geom_point(aes(x = ps, y=-log10(p_score), color=abs(beta)),size=1)+
  scale_color_gradientn(colours = brewer.pal(n=9,name="Greys")) +
  labs(y="-log10 P", x='genomic position (bp)',color='effect in relative\nseed production') +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  facet_wrap(~condition, ncol=1)
dev.off()

##############################################
# Odds ratio
tmpthp<-dplyr::filter(tmpall, grepl("thp", condition), grepl("Seed", condition))
tmpthp$SNP<-paste0(tmpthp$chr, "_",tmpthp$ps)
dim(tmpthp)
ggplot(tmpthp)+
  geom_vline(xintercept=3173382, color='grey', lty='dotted')+
  geom_vline(xintercept=3179448, color='grey', lty='dotted')+
  geom_point(aes(x = ps, y=-log10(p_score)), size=2, color='black')+
  geom_point(aes(x = ps, y=-log10(p_score)), size=1.5, color='white')+
  geom_point(aes(x = ps, y=-log10(p_score), color=abs(beta)),size=1)+
  scale_color_gradientn(colours = brewer.pal(n=9,name="Greys")) +
  labs(y="-log10 P", x='genomic position (bp)',color='effect in relative\nseed production') +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  geom_vline(xintercept=3179448-500, color='grey', lty='dotted')+
  geom_vline(xintercept=3179448+2500, color='grey', lty='dotted')


# more targeted analysis
# merge with FST
# selcol = c("S4.S10", "S5.S11", "S6.S12")
selcol = c(
           "S1.S4","S1.S7","S1.S10",
           "S2.S5","S2.S8","S2.S11",
           "S3.S6","S3.S9","S3.S12"
           )
# selcol =  "S2.S5"
tmpfst<-subset(tab_long, time %in% selcol)
tmpfst$SNP<-paste0(tmpfst$chr, "_",tmpfst$pos)
head(tmpfst)
mfst<-inner_join(x = tmpfst, y= tmpthp, by=c("SNP"))
dim(mfst)
mfst<-na.omit(mfst)
mfst<- mfst %>%
  dplyr::mutate(promotor= pos>3179448-500  & pos < 3179448+2500 ) %>%
  dplyr::mutate(rep=substr(time,2,2))
# show there is a wilcox text difference with all of them
wilcox.test(mfst$S ~ mfst$promotor)
t.test(mfst$S ~ mfst$promotor)
mfst %>% group_by(promotor) %>% dplyr::summarise(q75=quantile(S, probs = 0.95), mean= median(S))
boxplot(mfst$beta ~ mfst$promotor)
# show there is a replicated fashion
mfst %>%
    dplyr::group_by(rep) %>%
    dplyr::summarize(wilpval = wilcox.test(S ~ promotor)$p.value) ->
  wilcoxregiontest
write.csv(file = "tables/wilcoxregionpromoter_time0-123together.csv", wilcoxregiontest)

# what about the GWA?
tmpthp<-tmpthp %>%
  dplyr::mutate(promotor= ps>3179448-500  & ps < 3179448+2500 )
wilcox.test(tmpthp$beta ~ tmpthp$promotor)
t.test(tmpthp$beta ~ tmpthp$promotor)
boxplot(tmpthp$beta ~ tmpthp$promotor)
tmpthp %>% group_by(promotor) %>% dplyr::summarise(q75=quantile(beta, probs = 0.95), mean= median(beta))

# does the region -500  + 2500 have a peak in both?
# yes

