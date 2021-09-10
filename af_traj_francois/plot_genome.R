# Basic plotting
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

library(tidyverse)

# Read the table.
francois <- read.csv( "all.csv", sep="," )

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

# FLC: 	3173382 - 3179448
# Range around it: 3167316 - 3185514

# No idea how to read that syntax, but it seem to work
flc <- francois %>% filter(Chrom == 5, Pos >= 3167316, Pos <= 3185514)
head(flc)


1T1 <- ggplot(flc, aes(x=Pos, y=R1T1)) + geom_hline(yintercept = 0, color="gray") + geom_vline(xintercept = 3173382, color="gray") + geom_vline(xintercept = 3179448, color="gray")  + geom_point() + ylim(-1, 1)
p_R1T2 <- ggplot(flc, aes(x=Pos, y=R1T2)) + geom_hline(yintercept = 0, color="gray") + geom_vline(xintercept = 3173382, color="gray") + geom_vline(xintercept = 3179448, color="gray")  + geom_point() + ylim(-1, 1)
p_R1T3 <- ggplot(flc, aes(x=Pos, y=R1T3)) + geom_hline(yintercept = 0, color="gray") + geom_vline(xintercept = 3173382, color="gray") + geom_vline(xintercept = 3179448, color="gray")  + geom_point() + ylim(-1, 1)
p_R2T1 <- ggplot(flc, aes(x=Pos, y=R2T1)) + geom_hline(yintercept = 0, color="gray") + geom_vline(xintercept = 3173382, color="gray") + geom_vline(xintercept = 3179448, color="gray")  + geom_point() + ylim(-1, 1)
p_R2T2 <- ggplot(flc, aes(x=Pos, y=R2T2)) + geom_hline(yintercept = 0, color="gray") + geom_vline(xintercept = 3173382, color="gray") + geom_vline(xintercept = 3179448, color="gray")  + geom_point() + ylim(-1, 1)
p_R2T3 <- ggplot(flc, aes(x=Pos, y=R2T3)) + geom_hline(yintercept = 0, color="gray") + geom_vline(xintercept = 3173382, color="gray") + geom_vline(xintercept = 3179448, color="gray")  + geom_point() + ylim(-1, 1)
p_R3T1 <- ggplot(flc, aes(x=Pos, y=R3T1)) + geom_hline(yintercept = 0, color="gray") + geom_vline(xintercept = 3173382, color="gray") + geom_vline(xintercept = 3179448, color="gray")  + geom_point() + ylim(-1, 1)
p_R3T2 <- ggplot(flc, aes(x=Pos, y=R3T2)) + geom_hline(yintercept = 0, color="gray") + geom_vline(xintercept = 3173382, color="gray") + geom_vline(xintercept = 3179448, color="gray")  + geom_point() + ylim(-1, 1)
p_R3T3 <- ggplot(flc, aes(x=Pos, y=R3T3)) + geom_hline(yintercept = 0, color="gray") + geom_vline(xintercept = 3173382, color="gray") + geom_vline(xintercept = 3179448, color="gray")  + geom_point() + ylim(-1, 1)

plot_grid(p_R1T1, p_R1T2, p_R1T3, p_R2T1, p_R2T2, p_R2T3, p_R3T1, p_R3T2, p_R3T3, ncol = 3, nrow = 3)

#gsave("pos_grid.png", width = 20, height =12)
ggsave("grid_20x12_600.png", width = 20, height =12, dpi=600)
ggsave("grid_15x9_600.png", width = 15, height =9, dpi=600)
