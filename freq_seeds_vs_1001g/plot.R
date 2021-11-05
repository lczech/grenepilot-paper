# Basic plotting
library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

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

print("computing frequencies and joined tabled")

# The seeds freqs are unfolded, but the 1001g and 231g are not,
# so we havel to fold the seeds as well.
index <- seeds$TOTAL.FREQ > 0.5
seeds$TOTAL.FREQ[index] <- 1.0 - seeds$TOTAL.FREQ[index]
#seeds[seeds$TOTAL.FREQ > 0.5]$TOTAL.FREQ <- 1.0 - seeds[seeds$TOTAL.FREQ > 0.5]$TOTAL.FREQ

# Join seeds with the two others, using the SNP column to only select
# those rows that appear in both tables, omitting all other rows,
# and filter out all NAN rows.

joined_otog <- merge( x=seeds, y=otog, by.x="CHROM_POS", by.y="SNP", all=FALSE)
joined_otog <- na.omit(joined_otog)

joined_ttog <- merge( seeds, ttog, by.x="CHROM_POS", by.y="SNP", all=FALSE)
joined_ttog <- na.omit(joined_ttog)

joined_both <- merge( otog, ttog, by.x="SNP", by.y="SNP", all=FALSE )
joined_both <- na.omit(joined_both)

sd_otog <- sd(joined_otog$TOTAL.FREQ - joined_otog$MAF)
pc_otog <- cor( x=joined_otog$TOTAL.FREQ, y=joined_otog$MAF )

sd_ttog <- sd(joined_ttog$TOTAL.FREQ - joined_ttog$MAF)
pc_ttog <- cor( x=joined_ttog$TOTAL.FREQ, y=joined_ttog$MAF )

sd_both <- sd(joined_both$MAF.x - joined_both$MAF.y)
pc_both <- cor( x=joined_both$MAF.x, y=joined_both$MAF.y )

# See if things are still good
#print(head(joined_otog))
#print(head(joined_ttog))

# -----------------------------------------------------------------
#     dot plots
# -----------------------------------------------------------------

print("plotting dots")

# Tests...
#plot(joined$"MAF.x", joined$"MAF.y", main="Seed vs 1001g AFs",  xlab="Seed MAF", ylab="1001g MAF", pch=19)
#ggplot(joined_otog, aes(x=MAF.x, y=MAF.y)) + geom_point(alpha = 1/100) #+ geom_smooth(method=lm)

# Plot otog
ggplot(joined_otog, aes(x=TOTAL.FREQ, y=MAF)) + 
    geom_hline( yintercept=0, color="gray") +
    geom_vline( xintercept=0, color="gray") +
    geom_abline(slope = 1, colour="blue") +
    geom_point(alpha = 1/100, size=0.1) + 
    xlab("Seed MAF") + 
    ylab("1001g MAF") +
    xlim(0, 0.5) +
    ylim(0, 0.5) +
    annotate("text", x = 0.02, y = Inf, hjust = 0, vjust=1, label = paste0("SD: ", format(sd_otog, digits=3), ", r: ", format(pc_otog, digits=3)))

#show()
ggsave("seeds_vs_1001g_freqs.png", bg="white")
#ggsave("seeds_vs_1001g_freqs_large.png", width=16, height=15)

# Plot ttog
ggplot(joined_ttog, aes(x=TOTAL.FREQ, y=MAF)) + 
    geom_hline( yintercept=0, color="gray") +
    geom_vline( xintercept=0, color="gray") +
    geom_point(alpha = 1/100, size=0.1) + 
    geom_abline(slope = 1, colour="blue") +
    xlab("Seed MAF") + 
    ylab("231g MAF")  +
    xlim(0, 0.5) +
    ylim(0, 0.5) +
    annotate("text", x = 0.02, y = Inf, hjust = 0, vjust=1, label = paste0("SD: ", format(sd_ttog, digits=3), ", r: ", format(pc_ttog, digits=3)))

ggsave("seeds_vs_231g_freqs.png", bg="white")
#ggsave("seeds_vs_231g_freqs_large.png", width=16, height=15)


# Plot both
ggplot(joined_both, aes(x=MAF.x, y=MAF.y)) +
    geom_hline( yintercept=0, color="gray") +
    geom_vline( xintercept=0, color="gray") +
    geom_point(alpha = 1/100, size=0.1) +
    geom_abline(slope = 1, colour="blue") +
    xlab("1001g MAF") +
    ylab("231g MAF")  +
    xlim(0, 0.5) +
    ylim(0, 0.5) +
    annotate("text", x = 0.02, y = Inf, hjust = 0, vjust=1, label = paste0("SD: ", format(sd_both, digits=3), ", r: ", format(pc_both, digits=3)))

ggsave("1001g_vs_231g_freqs.png", bg="white")
#ggsave("seeds_vs_231g_freqs_large.png", width=16, height=15)

# -----------------------------------------------------------------
#     scatter plots
# -----------------------------------------------------------------

print("plotting scatter")

# Plot otog
ggplot(joined_otog, aes(x=TOTAL.FREQ, y=MAF)) + 
    geom_hex(bins=50) +
    #scale_fill_viridis(direction=-1) +
    scale_fill_viridis(direction=-1) +
    scale_color_viridis(direction=-1) +
    coord_fixed() +
    xlab("Seed MAF") + 
    ylab("1001g MAF") +
    xlim(0, 0.5) +
    ylim(0, 0.5) +
    annotate("text", x = 0.02, y = Inf, hjust = 0, vjust=1, label = paste0("SD: ", format(sd_otog, digits=3), ", r: ", format(pc_otog, digits=3)))

#show()
ggsave("seeds_vs_1001g_hex.png", bg="white")
#ggsave("seeds_vs_1001g_freqs_large.png", width=16, height=15)

# Plot ttog
ggplot(joined_ttog, aes(x=TOTAL.FREQ, y=MAF)) + 
    geom_hex(bins=50) +
    #scale_fill_viridis(direction=-1) +
    scale_fill_viridis(direction=-1) +
    scale_color_viridis(direction=-1) +
    coord_fixed() +
    xlab("Seed MAF") + 
    ylab("231g MAF")  +
    xlim(0, 0.5) +
    ylim(0, 0.5) +
    annotate("text", x = 0.02, y = Inf, hjust = 0, vjust=1, label = paste0("SD: ", format(sd_ttog, digits=3), ", r: ", format(pc_ttog, digits=3)))

ggsave("seeds_vs_231g_hex.png", bg="white")
#ggsave("seeds_vs_231g_freqs_large.png", width=16, height=15)


# Plot both
ggplot(joined_both, aes(x=MAF.x, y=MAF.y)) +
    geom_hex(bins=50) +
    #scale_fill_viridis(direction=-1) +
    scale_fill_viridis(direction=-1) +
    scale_color_viridis(direction=-1) +
    coord_fixed() +
    xlab("1001g MAF") +
    ylab("231g MAF")  +
    xlim(0, 0.5) +
    ylim(0, 0.5) +
    annotate("text", x = 0.02, y = Inf, hjust = 0, vjust=1, label = paste0("SD: ", format(sd_both, digits=3), ", r: ", format(pc_both, digits=3)))

ggsave("1001g_vs_231g_hex.png", bg="white")
#ggsave("seeds_vs_231g_freqs_large.png", width=16, height=15)


