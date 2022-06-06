# Basic plotting
library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
suppressMessages(library(viridis))

print("reading tables")

# Read the tables. Ours has proper tabs as delimiters,
# the plink ones just uses some arbitrary spaces between columns,
# so we need to set the sep accordingly.
# Naming: seeds (obviously), 1001g (One Thousand One G) = otog,
# and 231g (Two Three One G) = ttog
#seeds <- read.csv( "frequency.csv", sep="\t" )
#otog <- read.csv( "plink.frq", sep="" )
#ttog <- read.csv( "231g.frq", sep="" )

seeds <- fread( "frequency.csv", sep="\t" )
otog  <- fread( "plink.frq", sep="auto" )
ttog  <- fread( "231g.frq", sep="auto" )

# Test cases
#seeds <- read.csv( "head_seeds.csv", sep="\t" )
#otog <- read.csv( "head_plink.frq", sep="" )
#ttog <- read.csv( "head_231g.frq", sep="" )

# See if all is good
#head(seeds)
#head(otog)
#head(ttog)

print("computing frequencies and joined tables")

# The seeds freqs are unfolded, but the 1001g and 231g are not,
# so we havel to fold the seeds as well.
index <- seeds$TOTAL.FREQ > 0.5
seeds$TOTAL.FREQ[index] <- 1.0 - seeds$TOTAL.FREQ[index]
#seeds[seeds$TOTAL.FREQ > 0.5]$TOTAL.FREQ <- 1.0 - seeds[seeds$TOTAL.FREQ > 0.5]$TOTAL.FREQ

# Join seeds with the two others, using the SNP column to only select
# those rows that appear in both tables, omitting all other rows,
# and filter out all NAN rows.

# make merged tables of the data
joined_otog <- merge( x=seeds, y=otog, by.x="CHROM_POS", by.y="SNP", all=FALSE)
joined_otog <- na.omit(joined_otog)
joined_ttog <- merge( seeds, ttog, by.x="CHROM_POS", by.y="SNP", all=FALSE)
joined_ttog <- na.omit(joined_ttog)
joined_both <- merge( otog, ttog, by.x="SNP", by.y="SNP", all=FALSE )
joined_both <- na.omit(joined_both)

# consistent renaming of the important columns
colnames(joined_otog)[colnames(joined_otog) == 'TOTAL.FREQ'] <- 'MAF1'
colnames(joined_otog)[colnames(joined_otog) == 'MAF'] <- 'MAF2'
colnames(joined_ttog)[colnames(joined_ttog) == 'TOTAL.FREQ'] <- 'MAF1'
colnames(joined_ttog)[colnames(joined_ttog) == 'MAF'] <- 'MAF2'
colnames(joined_both)[colnames(joined_both) == 'MAF.x'] <- 'MAF1'
colnames(joined_both)[colnames(joined_both) == 'MAF.y'] <- 'MAF2'

# See if things are still good
#print(head(joined_otog))
#print(head(joined_ttog))
################################################################################
#### FILTER FOR QUALITY! Jan 4 2022 Moi
library(dplyr)

joined_otog$SNP <- paste0(joined_otog$CHROM, "_", joined_otog$POS)
joined_ttog$SNP <- paste0(joined_ttog$CHROM, "_", joined_ttog$POS)

# bona fide 515 genomes SNPs
s<-fread("515g.bim")
joined_otog_qual<-inner_join(joined_otog,s,by=c("SNP"="V2"))
joined_ttog_qual<-inner_join(joined_ttog,s,by=c("SNP"="V2"))

# minimum alternative count 3 or more
joined_otog_qual_mac <- dplyr::filter(joined_otog_qual, TOTAL.ALT_CNT >2)
joined_ttog_qual_mac <- dplyr::filter(joined_ttog_qual, TOTAL.ALT_CNT >2)

################################################################################
# -----------------------------------------------------------------
#     plot function
# -----------------------------------------------------------------

make_plots <- function( data, myxlab, myylab, target ) {

    print(target)

    mysd <- sd( data$MAF1 - data$MAF2)
    mypc <- cor( x=data$MAF1, y=data$MAF2 )

    # Plot dot
    ggplot(data, aes(x=MAF1, y=MAF2)) +
        geom_hline( yintercept=0, color="gray") +
        geom_vline( xintercept=0, color="gray") +
        geom_abline(slope = 1, colour="blue") +
        geom_point(alpha = 1/100, size=0.1) +
        xlab(myxlab) +
        ylab(myylab) +
        xlim(0, 0.5) +
        ylim(0, 0.5) +
        annotate("text", x = 0.02, y = Inf, hjust = 0, vjust=1, label = paste0(
            "SD: ", format(mysd, digits=3), ", r: ", format(mypc, digits=3))
        )

    ggsave(paste0(target, "_freqs.png"), bg="white")
    #ggsave("seeds_vs_1001g_freqs_large.png", width=16, height=15)

    # # Plot hex
    # ggplot(data, aes(x=MAF1, y=MAF2)) +
    #     geom_hex(bins=50) +
    #     #scale_fill_viridis(direction=-1) +
    #     scale_fill_viridis(direction=-1) +
    #     scale_color_viridis(direction=-1) +
    #     coord_fixed() +
    #     xlab(myxlab) +
    #     ylab(myylab) +
    #     xlim(0, 0.5) +
    #     ylim(0, 0.5) +
    #     annotate("text", x = 0.02, y = Inf, hjust = 0, vjust=1, label = paste0(
    #         "SD: ", format(mysd, digits=3), ", r: ", format(mypc, digits=3))
    #     )
    # 
    # ggsave(paste0(target, "_hex.png"), bg="white")
    # #ggsave("seeds_vs_1001g_freqs_large.png", width=16, height=15)

# 
#     # Plot hex
#     ggplot(data, aes(x=MAF1, y=MAF2)) +
#         geom_hex(bins=50) +
#         #scale_fill_viridis(direction=-1) +
#         scale_fill_viridis(direction=-1, trans="log10") +
#         scale_color_viridis(direction=-1, trans="log10") +
#         coord_fixed() +
#         xlab(myxlab) +
#         ylab(myylab) +
#         xlim(0, 0.5) +
#         ylim(0, 0.5) +
#         annotate("text", x = 0.02, y = Inf, hjust = 0, vjust=1, label = paste0(
#             "SD: ", format(mysd, digits=3), ", r: ", format(mypc, digits=3))
#         )
# 
#     ggsave(paste0(target, "_hex_log.png"), bg="white")
#     #ggsave("seeds_vs_1001g_freqs_large.png", width=16, height=15)

}

# -----------------------------------------------------------------
#     do plots
# -----------------------------------------------------------------

make_plots( joined_otog_qual_mac, "Seed MAF", "1001g MAF", "seeds_vs_1001g_quality515snps_2MAC" )
make_plots( joined_ttog_qual_mac, "Seed MAF", "231g MAF", "seeds_vs_231g_quality515snps_2MAC" )

