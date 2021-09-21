# Basic plotting
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# Read the tables. Ours has proper tabs as delimiters,
# the plink ones just uses some arbitrary spaces between columns,
# so we need to set the sep accordingly.
# Naming: seeds (obviously), 1001g (One Thousand One G) = otog,
# and 231g (Two Three One G) = ttog
seeds <- read.csv( "seed_freqs.csv", sep="\t" )
otog <- read.csv( "plink.frq", sep="" )
ttog <- read.csv( "231g.frq", sep="" )

# See if all is good
#head(seeds)
#head(otog)
#head(ttog)

# Join seeds with the two others, using the SNP column to only select
# those rows that appear in both tables, omitting all other rows,
# and filter out all NAN rows.
joined_otog <- merge( seeds, otog, by.x="SNP", by.y="SNP", all=FALSE)
joined_otog <- na.omit(joined_otog)
joined_ttog <- merge( seeds, ttog, by.x="SNP", by.y="SNP", all=FALSE)
joined_ttog <- na.omit(joined_ttog)

# See if things are still good
#head(joined_otog)
#head(joined_ttog)

# Tests...
#plot(joined$"MAF.x", joined$"MAF.y", main="Seed vs 1001g AFs",  xlab="Seed MAF", ylab="1001g MAF", pch=19)
#ggplot(joined_otog, aes(x=MAF.x, y=MAF.y)) + geom_point(alpha = 1/100) #+ geom_smooth(method=lm)

# Plot otog
ggplot(joined_otog, aes(x=MAF.x, y=MAF.y)) + geom_point(alpha = 1/100) + geom_abline(slope = 1, colour="blue") + xlab("Seed MAF") + ylab("1001g MAF")
#show()
# ggsave cannot do pixels. stupid. let's try to find a size that works.
ggsave("seeds_vs_1001g_freqs.png")
ggsave("seeds_vs_1001g_freqs_large.png", width=16, height=15)

# Plot ttog
ggplot(joined_ttog, aes(x=MAF.x, y=MAF.y)) + geom_point(alpha = 1/100) + geom_abline(slope = 1, colour="blue") + xlab("Seed MAF") + ylab("231g MAF")
ggsave("seeds_vs_231g_freqs.png")
ggsave("seeds_vs_231g_freqs_large.png", width=16, height=15)
