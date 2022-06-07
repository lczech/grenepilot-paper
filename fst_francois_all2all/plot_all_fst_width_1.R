#Basic plottingi
library(tidyr)
library(data.table)
library(ggplot2)
library(cowplot)
#library(tidyverse)
theme_set(theme_cowplot())
suppressMessages(library(viridis))
library(plyr)
library(forcats)
library(qqman)
#library(grid)
#library(gridGraphics)
library(dplyr)

# Input data
print("reading")

bim <- fread( "/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/1001g/515g.bim" , sep=" ", header=FALSE )
colnames(bim)[which(names(bim) == "V2")] <- "CHROM_POS"

#infile="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas/fst_francois_all2all/fst-width-1.csv"
infile="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas_quick_for_paper_submission/grenepilot-paper/fst_francois_all2all/fst-width-1.csv"
#infile="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas/fst_francois_all2all/head.csv"
infile="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas_quick_for_paper_submission/grenepilot-paper/fst_francois_all2all/fst-width-1-S2-only.csv"
df <- fread( infile, sep="\t", header=TRUE )

# Need some conversions and shit for this to work...
df <- as.data.frame(df)
df <- df[df$CHROM != "chloroplast" & df$CHROM != "mitochondria", ]
df$CHROM <- as.numeric(as.character(df$CHROM))
df$CHROM_POS <- paste0( as.character(df$CHROM), "_", as.character(df$START) )

# add a snp column that has Xn at FLC
df$SNP <- ""
#df$SNP[df$CHROM == 5 & df$START >= 3173382 & df$START <= 3179448] <- "X"
df$SNP[df$CHROM == 5 & df$START >= 3167316 & df$START <= 3185514] <- "X"

# Check result
#print(head(df))
#print(tail(df))

# translation for samples, if needed
#S1  1_0 56  
#S4  1_1 80  
#S7  1_2 65  
#S10 1_3 19  
#S2  2_0 69  
#S5  2_1 160 
#S8  2_2 205 
#S11 2_3 50  
#S3  3_0 101 
#S6  3_1 200 
#S9  3_2 296 
#S12 3_3 97 

outdir="fst_width_1/"
outdir="fst_width_1_test/"

print("manhattan")
for(i in 5:(ncol(df)-2)) {

    #sample <- "S1.S2"
    sample <- colnames(df)[i]
    print(paste0("  At ", i, ": ", sample))

    # need a copy of the data, so that we can remove NA values without affecting all other columns...
    cp <- df[, "CHROM", drop=FALSE]
    cp$POS <- df$START
    cp$SNP <- df$SNP
    cp$FRQ <- df[, sample]
    cp$CHROM_POS <- df$CHROM_POS
    cp <- na.omit(cp)

    # just bona fide
    bf <- inner_join( cp, bim, by="CHROM_POS")

    if(TRUE){

    # If we want to save a png, we need to overwrite the output device manually first...
    #par(mar=c(1,1,1,1))
    png(paste0(outdir,sample,"-all.png"), width=2000, height=800)
    manhattan(cp, chr="CHROM", bp="POS", snp="SNP", p="FRQ", logp=FALSE, ylab="Fst", ylim = c(0, 1.05), highlight="X" )
    #manhattan(df, chr="CHROM", bp="START", snp=NULL, p=sample, logp=FALSE )
    dev.off()

    png(paste0(outdir,sample,"-bona-fide.png"), width=2000, height=800)
    manhattan(bf, chr="CHROM", bp="POS", snp="SNP", p="FRQ", logp=FALSE, ylab="Fst", ylim = c(0, 1.05), highlight="X" )
    #manhattan(df, chr="CHROM", bp="START", snp=NULL, p=sample, logp=FALSE )
    dev.off()

    }

    if(FALSE){

    ggplot(bf, aes(x=FRQ)) + geom_histogram(bins=100) + coord_flip() + xlab("") + ylab("") + xlim(-0.1, 1)
    ggsave(paste0(outdir,sample,"-hist.png"), bg="white")
    ggsave(paste0(outdir,sample,"-hist.svg"), bg="white")

    }
}

# Of course, we need special treatment to be able to save it as we want...
#p <- recordPlot()
#g <- grid.grabExpr(grid.echo(p))
#ggsave(paste0(sample,".png"), g, width=10, height=6, bg="white")

