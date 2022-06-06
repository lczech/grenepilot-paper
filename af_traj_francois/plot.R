library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyr)
library(RColorBrewer)

print("Reading data")

# Read the table. Kick out unneeded data.
francois <- read.csv( "frequency.csv", sep="\t" )
#francois <- read.csv( "test.csv", sep="\t" )
francois = francois[ ,!grepl("ALT_CNT",names(francois))]
francois = francois[ ,!grepl("REF_CNT",names(francois))]
#francois = francois[ ,!grepl("COV",names(francois))]

print("Transmogrifying data")

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

# We need to get tricky here. We only want certain combinations of columns
# (same as the above difference between samples from T0 to T1-3),
# which means, we have to manually specify which ones. Also, because I don't 
# know how to do this smartly in R, we do it manually...
# Also, so much copying going on here - so wasteful! Want to go back to C++!
data_freqs <- data.frame(
    Pair = c( 
        "R1T0.R1T1", "R1T0.R1T2", "R1T0.R1T3", 
        "R2T0.R2T1", "R2T0.R2T2", "R2T0.R2T3", 
        "R3T0.R3T1", "R3T0.R3T2", "R3T0.R3T3" 
    ),
    T0 = c(
        1 - francois[,"S1.FREQ"], 1 - francois[,"S1.FREQ"], 1 - francois[,"S1.FREQ"],
        1 - francois[,"S2.FREQ"], 1 - francois[,"S2.FREQ"], 1 - francois[,"S2.FREQ"],
        1 - francois[,"S3.FREQ"], 1 - francois[,"S3.FREQ"], 1 - francois[,"S3.FREQ"]
    ),
    TX = c(
        1 - francois[,"S4.FREQ"], 1 - francois[,"S7.FREQ"], 1 - francois[,"S10.FREQ"],
        1 - francois[,"S5.FREQ"], 1 - francois[,"S8.FREQ"], 1 - francois[,"S11.FREQ"],
        1 - francois[,"S6.FREQ"], 1 - francois[,"S9.FREQ"], 1 - francois[,"S12.FREQ"]
    ),
    Cov0 = c(
        francois[,"S1.COV"], francois[,"S1.COV"], francois[,"S1.COV"],
        francois[,"S2.COV"], francois[,"S2.COV"], francois[,"S2.COV"],
        francois[,"S3.COV"], francois[,"S3.COV"], francois[,"S3.COV"]
    ),
    CovX = c(
        francois[,"S4.COV"], francois[,"S7.COV"], francois[,"S10.COV"],
        francois[,"S5.COV"], francois[,"S8.COV"], francois[,"S11.COV"],
        francois[,"S6.COV"], francois[,"S9.COV"], francois[,"S12.COV"]
    )
)

data_freqs <- data_freqs[(data_freqs$Cov0 > 5 & data_freqs$CovX > 5),]

data_freqs$Diff <- data_freqs$TX - data_freqs$T0

data_freqs$T0[data_freqs$T0 == 0.0] <- NA
data_freqs$TX[data_freqs$TX == 0.0] <- NA
data_freqs$Diff[data_freqs$Diff == 0.0] <- NA
data_freqs$T0[data_freqs$T0 == 1.0] <- NA
data_freqs$TX[data_freqs$TX == 1.0] <- NA
data_freqs$Diff[data_freqs$Diff == 1.0] <- NA


data_freqs <- na.omit(data_freqs)

# Drop remaining unneeded data
francois = francois[ ,!grepl("FREQ",names(francois))]

print("francois")
print(head(francois))
print("data_freqs")
print(head(data_freqs))

francois <- na.omit(francois)

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

# Turn the table into R long format, using Sample as the column with the
# previous column names, and Freq for the frequency difference values.
#data <- gather(
#    francois, Sample, Freq, 
#    c(
#        `S1.FREQ`, `S2.FREQ`, `S3.FREQ`, `S4.FREQ`, `S5.FREQ`, `S6.FREQ`, 
#        `S7.FREQ`, `S8.FREQ`, `S9.FREQ`, `S10.FREQ`, `S11.FREQ`, `S12.FREQ` 
#    ), factor_key=TRUE
#)

data_diffs <- gather(
    francois, Sample, Freq,
    c(
        `R1T1`, `R1T2`, `R1T3`, `R2T1`, `R2T2`, `R2T3`, `R3T1`, `R3T2`, `R3T3`
    ), factor_key=TRUE
)

print("data_diffs")
print(head(data_diffs))

print("Plotting data")
#print(head(data))

# -------------------------------------------------------------------------
#     Histogram freq diff seeds and caller
# -------------------------------------------------------------------------

# Get standard deviation of the diff
stddev <- sd(data_diffs$Freq)

# Histogram of diff between the two
ggplot(data_diffs, aes(x=Freq)) +
    geom_histogram(bins=100) +
    xlab("Frequency Difference") +
    xlim(-1,1) +
    #labs(title=mytitle) +
    annotate("text", x = -1, y = Inf, hjust = 0, vjust=1, label = paste0("SD: ", format(stddev, digits=3)))

suppressMessages(ggsave("freq-diff-hist.png", bg="white"))

# -------------------------------------------------------------------------
#     Histograms of data freqs
# -------------------------------------------------------------------------

ggplot(data_freqs, aes(x=T0)) +
    geom_histogram(bins=100)
suppressMessages(ggsave("freq-hist-T0.png", bg="white"))

ggplot(data_freqs, aes(x=TX)) +
    geom_histogram(bins=100) 
suppressMessages(ggsave("freq-hist-TX.png", bg="white"))

ggplot(data_freqs, aes(x=Cov0)) +
    geom_histogram(bins=100)
suppressMessages(ggsave("freq-hist-Cov0.png", bg="white"))

ggplot(data_freqs, aes(x=CovX)) +
    geom_histogram(bins=100)
suppressMessages(ggsave("freq-hist-CovX.png", bg="white"))

# Alternative hist plot, this time with coverage filter
ggplot(data_freqs, aes(x=Diff)) +
    geom_histogram(bins=100) +
    xlim(-1,1) +
    xlab("Frequency Difference")

suppressMessages(ggsave("freq-diff-hist-cov5.png", bg="white"))

# -------------------------------------------------------------------------
#     Scatter plot seeds vs caller
# -------------------------------------------------------------------------

# Plot the scatter plot of the frequencies.
ggplot(data_freqs, aes(x=T0, y=TX)) +
    geom_vline( xintercept=0, color="gray") +
    geom_hline( yintercept=0, color="gray") +
    geom_point(alpha = 1/100, size=0.2) +
    # geom_point(alpha = 1/100, size=ds) +
    #geom_point(alpha = 1/100, size=ds, aes(color=COV)) +
    #scale_color_viridis(trans = "log10", direction=-1) +
    geom_abline(slope = 1, colour="blue") +
    coord_fixed() +
    xlab("Frequency at T0") +
    ylab("Frequency at Tx") +
    facet_wrap( ~Pair )
    #labs(title="Freq") +
    #annotate("text", x = 0.02, y = Inf, hjust = 0, vjust=1, label = paste0("SD: ", format(stddev, digits=3)))

# We save different sizes, as the pixel density changes.
suppressMessages(ggsave("freq-scatter.png", bg="white", width=10, height=10))
#ggsave(paste(f,"-large.png", sep=""), width=16, height=15)

# Same, but using grid
ggplot(data_freqs, aes(x=T0, y=TX)) +
    geom_bin2d() +
    geom_abline(slope = 1, colour="blue") +
    coord_fixed() +
    xlab("Frequency at T0") +
    ylab("Frequency at Tx") +
    scale_fill_gradientn("# SNPS",colours = RColorBrewer::brewer.pal(9, "Greys")[-1]) +
    facet_wrap( ~Pair )
    #labs(title="Freq") +
    #annotate("text", x = 0.02, y = Inf, hjust = 0, vjust=1, label = paste0("SD: ", format(stddev, digits=3)))

# We save different sizes, as the pixel density changes.
suppressMessages(ggsave("freq-scatter-bin-lin.png", bg="white", width=10, height=10))

# Same again, but log scaled color
ggplot(data_freqs, aes(x=T0, y=TX)) +
    geom_bin2d() +
    geom_abline(slope = 1, colour="blue") +
    coord_fixed() +
    xlab("Frequency at T0") +
    ylab("Frequency at Tx") +
    scale_fill_gradientn("# SNPS",colours = RColorBrewer::brewer.pal(9, "Greys")[-1], trans='log10') +
    facet_wrap( ~Pair )
    #labs(title="Freq") +
    #annotate("text", x = 0.02, y = Inf, hjust = 0, vjust=1, label = paste0("SD: ", format(stddev, digits=3)))

# We save different sizes, as the pixel density changes.
suppressMessages(ggsave("freq-scatter-bin-log.png", bg="white", width=10, height=10))

# -------------------------------------------------------------------------
#     Diff plots
# -------------------------------------------------------------------------

ggplot(data_freqs, aes(x=T0, y=Diff)) +
    geom_hex() +
    geom_hline( yintercept=0, color="blue") +
    xlab("Seed AF") +
    ylab("Delta AF Flower - Seed") +
    xlim(0, 1) +
    ylim(-1, 1) +
    #scale_fill_gradientn("# SNPS",colours = RColorBrewer::brewer.pal(9, "Greys")[-1]) +
    scale_fill_gradientn("# SNPS",colours = RColorBrewer::brewer.pal(9, "Greys")[-1], trans='log10') +
    facet_wrap( ~Pair )

suppressMessages(ggsave("freq-hex-log.png", bg="white", width=10, height=10))

ggplot(data_freqs, aes(x=T0, y=Diff)) +
    geom_hex() +
    geom_hline( yintercept=0, color="blue") +
    xlab("Seed AF") +
    ylab("Delta AF Flower - Seed") +
    xlim(0, 1) +
    ylim(-1, 1) +
    scale_fill_gradientn("# SNPS",colours = RColorBrewer::brewer.pal(9, "Greys")[-1]) +
    #scale_fill_gradientn("# SNPS",colours = RColorBrewer::brewer.pal(9, "Greys")[-1], trans='log10') +
    facet_wrap( ~Pair )

suppressMessages(ggsave("freq-hex-lin.png", bg="white", width=10, height=10))


warnings()
