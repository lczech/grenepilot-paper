library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyr)
library(RColorBrewer)

# Read the table.
francois <- read.csv( "frequency.csv", sep="\t" )
#francois <- read.csv( "test.csv", sep="\t" )

# See if all is good
#head(francois)

# Key for samples:
# S1      SampleId1-0_S1
# S2      SampleId2-0_S2
# S3      SampleId3-0_S3
# S4      SampleId1-1_S4
# S5      SampleId2-1_S5
# S6      SampleId3-1_S6
# S7      SampleId1-2_S7
# S8      SampleId2-2_S8
# S9      SampleId3-2_S9
# S10     SampleId1-3_S10
# S11     SampleId2-3_S11
# S12     SampleId3-3_S12

# Add columns that are the differences between time 0 and time 1-3, for each replicate
francois$R1T1 <- francois$S4.FREQ  - francois$S1.FREQ
francois$R1T2 <- francois$S7.FREQ  - francois$S1.FREQ
francois$R1T3 <- francois$S10.FREQ - francois$S1.FREQ
francois$R2T1 <- francois$S5.FREQ  - francois$S2.FREQ
francois$R2T2 <- francois$S8.FREQ  - francois$S2.FREQ
francois$R2T3 <- francois$S11.FREQ - francois$S2.FREQ
francois$R3T1 <- francois$S6.FREQ  - francois$S3.FREQ
francois$R3T2 <- francois$S9.FREQ  - francois$S3.FREQ
francois$R3T3 <- francois$S12.FREQ - francois$S3.FREQ

p_R1T1 <- ggplot(francois, aes(x=S1.FREQ, y=R1T1)) +
    geom_hex() +
    # ggtitle("R1T1") +
    xlab("Seed AF") +
    ylab("Delta AF Flower - Seed") +
    xlim(0, 1) +
    ylim(-1, 1) +
    scale_fill_gradientn("# SNPS",colours = RColorBrewer::brewer.pal(9, "Greys")[-1], trans='log10')
    
p_R1T2 <- ggplot(francois, aes(x=S1.FREQ, y=R1T2)) +
    geom_hex() +
    # ggtitle("R1T2") +
    xlab("Seed AF") +
    ylab("Delta AF Flower - Seed") +
    xlim(0, 1) +
    ylim(-1, 1) +
    scale_fill_gradientn("# SNPS",colours = RColorBrewer::brewer.pal(9, "Greys")[-1], trans='log10')
    
p_R1T3 <- ggplot(francois, aes(x=S1.FREQ, y=R1T3)) +
    geom_hex() +
    # ggtitle("R1T3") +
    xlab("Seed AF") +
    ylab("Delta AF Flower - Seed") +
    xlim(0, 1) +
    ylim(-1, 1) +
    scale_fill_gradientn("# SNPS",colours = RColorBrewer::brewer.pal(9, "Greys")[-1], trans='log10')
    
p_R2T1 <- ggplot(francois, aes(x=S2.FREQ, y=R2T1)) +
    geom_hex() +
    # ggtitle("R2T1") +
    xlab("Seed AF") +
    ylab("Delta AF Flower - Seed") +
    xlim(0, 1) +
    ylim(-1, 1) +
    scale_fill_gradientn("# SNPS",colours = RColorBrewer::brewer.pal(9, "Greys")[-1], trans='log10')
    
p_R2T2 <- ggplot(francois, aes(x=S2.FREQ, y=R2T2)) +
    geom_hex() +
    # ggtitle("R2T2") +
    xlab("Seed AF") +
    ylab("Delta AF Flower - Seed") +
    xlim(0, 1) +
    ylim(-1, 1) +
    scale_fill_gradientn("# SNPS",colours = RColorBrewer::brewer.pal(9, "Greys")[-1], trans='log10')
    
p_R2T3 <- ggplot(francois, aes(x=S2.FREQ, y=R2T3)) +
    geom_hex() +
    # ggtitle("R2T3") +
    xlab("Seed AF") +
    ylab("Delta AF Flower - Seed") +
    xlim(0, 1) +
    ylim(-1, 1) +
    scale_fill_gradientn("# SNPS",colours = RColorBrewer::brewer.pal(9, "Greys")[-1], trans='log10')
    
p_R3T1 <- ggplot(francois, aes(x=S3.FREQ, y=R3T1)) +
    geom_hex() +
    # ggtitle("R3T1") +
    xlab("Seed AF") +
    ylab("Delta AF Flower - Seed") +
    xlim(0, 1) +
    ylim(-1, 1) +
    scale_fill_gradientn("# SNPS",colours = RColorBrewer::brewer.pal(9, "Greys")[-1], trans='log10')
    
p_R3T2 <- ggplot(francois, aes(x=S3.FREQ, y=R3T2)) +
    geom_hex() +
    # ggtitle("R3T2") +
    xlab("Seed AF") +
    ylab("Delta AF Flower - Seed") +
    xlim(0, 1) +
    ylim(-1, 1) +
    scale_fill_gradientn("# SNPS",colours = RColorBrewer::brewer.pal(9, "Greys")[-1], trans='log10')
    
p_R3T3 <- ggplot(francois, aes(x=S3.FREQ, y=R3T3)) +
    geom_hex() +
    # ggtitle("R3T3") +
    xlab("Seed AF") +
    ylab("Delta AF Flower - Seed") +
    xlim(0, 1) +
    ylim(-1, 1) +
    scale_fill_gradientn("# SNPS",colours = RColorBrewer::brewer.pal(9, "Greys")[-1], trans='log10')

        
ggsave("p_R1T1.png", bg="white", p_R1T1)
ggsave("p_R1T2.png", bg="white", p_R1T2)
ggsave("p_R1T3.png", bg="white", p_R1T3)
ggsave("p_R2T1.png", bg="white", p_R2T1)
ggsave("p_R2T2.png", bg="white", p_R2T2)
ggsave("p_R2T3.png", bg="white", p_R2T3)
ggsave("p_R3T1.png", bg="white", p_R3T1)
ggsave("p_R3T2.png", bg="white", p_R3T2)
ggsave("p_R3T3.png", bg="white", p_R3T3)

# Haha, the below code is... interesting! There must be better ways to do this in R... :-D

t1plot<-ggplot() +annotate(geom = "text", label='T=1', x=1,y=1, angle = 90) +theme(element_blank(),element_blank(),element_blank())
t2plot<-ggplot() +annotate(geom = "text", label='T=2', x=1,y=1, angle = 90) +theme(element_blank(),element_blank(),element_blank())
t3plot<-ggplot() +annotate(geom = "text", label='T=3', x=1,y=1, angle = 90) +theme(element_blank(),element_blank(),element_blank())
r0plot<-ggplot() +theme(element_blank(),element_blank(),element_blank())
r1plot<-ggplot() +annotate(geom = "text", label='Rep=1', x=1,y=1, angle = 0) +theme(element_blank(),element_blank(),element_blank())
r2plot<-ggplot() +annotate(geom = "text", label='Rep=2', x=1,y=1, angle = 0) +theme(element_blank(),element_blank(),element_blank())
r3plot<-ggplot() +annotate(geom = "text", label='Rep=3', x=1,y=1, angle = 0) +theme(element_blank(),element_blank(),element_blank())
columnt<-plot_grid(ncol=1, t1plot,t2plot,t3plot)
rowrep<-plot_grid(nrow=1, r0plot, plot_grid(r1plot,r2plot,r3plot,nrow=1), rel_widths = c(1,10))
coreplot<-plot_grid(p_R1T1, p_R1T2, p_R1T3, p_R2T1, p_R2T2, p_R2T3, p_R3T1, p_R3T2, p_R3T3, ncol = 3, nrow = 3)
tmp1<-plot_grid(ncol=2,nrow=1,rel_widths = c(1,10), plotlist=list(columnt,coreplot))
tmp2<-plot_grid(ncol=1,nrow=2,rel_heights = c(10,1),plotlist = list(tmp1, rowrep))

ggsave("delta.png", width=16, height=10, dpi=600, bg="white", tmp2)

warnings()
