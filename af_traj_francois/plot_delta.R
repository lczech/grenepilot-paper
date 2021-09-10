# memex does not provide packages yet, so we have to install everything first...
# Also, of course, we need to do that locally, and provide CRAN, because of reasons.
# Update: switching to calc now, where this does not seem necessary.
#install.packages("ggplot2", lib="/Carnegie/DPB/Homes/Users/lczech/Rpackages", repos = "http://cran.us.r-project.org")
#install.packages("cowplot", lib="/Carnegie/DPB/Homes/Users/lczech/Rpackages", repos = "http://cran.us.r-project.org")
#install.packages("labeling", lib="/Carnegie/DPB/Homes/Users/lczech/Rpackages", repos = "http://cran.us.r-project.org")

# Basic plotting. Load libraries with path from above. Not needed any more on calc.
#library(labeling, lib.loc="/Carnegie/DPB/Homes/Users/lczech/Rpackages")
#library(ggplot2, lib.loc="/Carnegie/DPB/Homes/Users/lczech/Rpackages")
#library(cowplot, lib.loc="/Carnegie/DPB/Homes/Users/lczech/Rpackages")

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# Read the table.
francois <- read.csv( "all.csv", sep="," )

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
francois$R1T1 <- francois$S4  - francois$S1
francois$R1T2 <- francois$S7  - francois$S1
francois$R1T3 <- francois$S10 - francois$S1
francois$R2T1 <- francois$S5  - francois$S2
francois$R2T2 <- francois$S8  - francois$S2
francois$R2T3 <- francois$S11 - francois$S2
francois$R3T1 <- francois$S6  - francois$S3
francois$R3T2 <- francois$S9  - francois$S3
francois$R3T3 <- francois$S12 - francois$S3

# For testing
#francois <- head(francois, 10000)

# For testing: full plots with all titles and axes labelled
#p_R1T1 <- ggplot(francois, aes(x=S1, y=R1T1)) + geom_point(alpha = 1/100) + ggtitle("R1T1") + xlab("Seed AF") + ylab("Delta AF Seed - Flower")
#p_R1T2 <- ggplot(francois, aes(x=S1, y=R1T2)) + geom_point(alpha = 1/100) + ggtitle("R1T2") + xlab("Seed AF") + ylab("Delta AF Seed - Flower")
#p_R1T3 <- ggplot(francois, aes(x=S1, y=R1T3)) + geom_point(alpha = 1/100) + ggtitle("R1T3") + xlab("Seed AF") + ylab("Delta AF Seed - Flower")
#p_R2T1 <- ggplot(francois, aes(x=S2, y=R2T1)) + geom_point(alpha = 1/100) + ggtitle("R2T1") + xlab("Seed AF") + ylab("Delta AF Seed - Flower")
#p_R2T2 <- ggplot(francois, aes(x=S2, y=R2T2)) + geom_point(alpha = 1/100) + ggtitle("R2T2") + xlab("Seed AF") + ylab("Delta AF Seed - Flower")
#p_R2T3 <- ggplot(francois, aes(x=S2, y=R2T3)) + geom_point(alpha = 1/100) + ggtitle("R2T3") + xlab("Seed AF") + ylab("Delta AF Seed - Flower")
#p_R3T1 <- ggplot(francois, aes(x=S3, y=R3T1)) + geom_point(alpha = 1/100) + ggtitle("R3T1") + xlab("Seed AF") + ylab("Delta AF Seed - Flower")
#p_R3T2 <- ggplot(francois, aes(x=S3, y=R3T2)) + geom_point(alpha = 1/100) + ggtitle("R3T2") + xlab("Seed AF") + ylab("Delta AF Seed - Flower")
#p_R3T3 <- ggplot(francois, aes(x=S3, y=R3T3)) + geom_point(alpha = 1/100) + ggtitle("R3T3") + xlab("Seed AF") + ylab("Delta AF Seed - Flower")

# Could not figure out how to get proper row and column labels, so we misuse individual subplot labels for this purpose...
# No idea how cowplot wants us to do this properly...
# Also, doing loops over column names seems to be super tricky in R, as it confuses column names for variables,
# and so we do everyting manually here.... :-(
p_R1T1 <- ggplot(francois, aes(x=S1, y=R1T1)) + geom_point(alpha = 1/100, size=0.05) + xlab("") + ylab("Replicate 1") + geom_hline(yintercept = 0, color="gray") #+ ggtitle("R1T1")
p_R1T2 <- ggplot(francois, aes(x=S1, y=R1T2)) + geom_point(alpha = 1/100, size=0.05) + xlab("") + ylab("") + geom_hline(yintercept = 0, color="gray") #+ ggtitle("R1T2")
p_R1T3 <- ggplot(francois, aes(x=S1, y=R1T3)) + geom_point(alpha = 1/100, size=0.05) + xlab("") + ylab("") + geom_hline(yintercept = 0, color="gray") #+ ggtitle("R1T3")
p_R2T1 <- ggplot(francois, aes(x=S2, y=R2T1)) + geom_point(alpha = 1/100, size=0.05) + xlab("") + ylab("Replicate 2") + geom_hline(yintercept = 0, color="gray") #+ ggtitle("R2T1")
p_R2T2 <- ggplot(francois, aes(x=S2, y=R2T2)) + geom_point(alpha = 1/100, size=0.05) + xlab("") + ylab("") + geom_hline(yintercept = 0, color="gray") #+ ggtitle("R2T2")
p_R2T3 <- ggplot(francois, aes(x=S2, y=R2T3)) + geom_point(alpha = 1/100, size=0.05) + xlab("") + ylab("") + geom_hline(yintercept = 0, color="gray") #+ ggtitle("R2T3")
p_R3T1 <- ggplot(francois, aes(x=S3, y=R3T1)) + geom_point(alpha = 1/100, size=0.05) + xlab("Time 1") + ylab("Replicate 3")+ geom_hline(yintercept = 0, color="gray") #+ ggtitle("R3T1")
p_R3T2 <- ggplot(francois, aes(x=S3, y=R3T2)) + geom_point(alpha = 1/100, size=0.05) + xlab("Time 2") + ylab("")+ geom_hline(yintercept = 0, color="gray") #+ ggtitle("R3T2")
p_R3T3 <- ggplot(francois, aes(x=S3, y=R3T3)) + geom_point(alpha = 1/100, size=0.05) + xlab("Time 3") + ylab("")+ geom_hline(yintercept = 0, color="gray") #+ ggtitle("R3T3")

plot_grid(p_R1T1, p_R1T2, p_R1T3, p_R2T1, p_R2T2, p_R2T3, p_R3T1, p_R3T2, p_R3T3, ncol = 3, nrow = 3)

# Plot with different sizes and resolutions, as those plots look very different,
# and we kind of have to explore what the "correct" visualization is...

ggsave("delta_12x8_600_0.05.png", width=12, height=8, dpi=600)

#ggsave("delta_8x8_300.png", width=8, height=8, dpi=300)
#ggsave("delta_12x12_300.png", width=12, height=12, dpi=300)
#ggsave("delta_16x16_300.png", width=16, height=16, dpi=300)

#ggsave("delta_8x8_600.png", width=8, height=8, dpi=600)
#ggsave("delta_12x12_600.png", width=12, height=12, dpi=600)
#ggsave("delta_16x16_600.png", width=16, height=16, dpi=600)

#ggsave("delta_8x8_1200.png", width=8, height=8, dpi=1200)
#ggsave("delta_12x12_1200.png", width=12, height=12, dpi=1200)
#ggsave("delta_16x16_1200.png", width=16, height=16, dpi=1200)
