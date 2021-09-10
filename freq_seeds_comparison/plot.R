# Basic plotting
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
suppressMessages(library(viridis))

# R is really weird concerning arrays/lists and index-based access...
# See https://stackoverflow.com/a/43534546/4184258

# Input files
vcf_freq_files <- c(
    "bcftools-c-filtered",
    "bcftools-m-filtered",
    "freebayes-kv-filtered",
    "freebayes-nkv-filtered",
    "haplotypecaller-kv-filtered",
    "haplotypecaller-nkv-filtered"
)
pileup_freq_file = "frequency-ath-evo-seeds-pileups-all-merged-samples.csv"

# Read all vcf freq tables.
# We do it all at once, in preparation for later, when we probably want to make
# one big facet plot with all of them in one figure.
vcf_freq_data <- list()
i <- 1
for (f in vcf_freq_files) {
    fn <- paste("frequency-ath-evo-seeds-", f, ".csv", sep="")
    vcf_freq_data[[i]] <- read.table(fn, sep="\t", header=TRUE)
    i <- i+1
}
head(vcf_freq_data)
head(vcf_freq_data[[1]])

# Read pileup freq table.
pileup_freq_data <- read.table(pileup_freq_file, sep="\t", header=TRUE)
head(pileup_freq_data)

# Go through all freq tables and plot their frequency vs the pileup.
i <- 1
for (f in vcf_freq_files) {

    # Make a joined table that only contains the rows of SNPs
    # that occur in both the pileup and the vcf.
    joined_data <- merge(
        pileup_freq_data, vcf_freq_data[[i]],
        by.x = "CHROM_POS", by.y = "CHROM_POS", all=FALSE
    )
    joined_data <- na.omit(joined_data)
    joined_data$COV <- joined_data$TOTAL.COV.x + joined_data$TOTAL.COV.y

    # Change ref freq to minor freq
    joined_data$TOTAL.FREQ.x <- 1.0 - joined_data$TOTAL.FREQ.x
    joined_data$TOTAL.FREQ.y <- 1.0 - joined_data$TOTAL.FREQ.y

    # allow to skip this part
    if(FALSE) {

    # Plot the scatter plot of the frequencies.
    # We want to use different dot sizes to see how this affects the plot.
    #for( ds in c(0.01, 0.03, 0.1, 0.3, 1.0)) { 
    for( ds in c()) {
        ggplot(joined_data, aes(x=TOTAL.FREQ.x, y=TOTAL.FREQ.y)) + 
            geom_point(alpha = 1/100, size=ds) +
            #geom_point(alpha = 1/100, size=ds, aes(color=COV)) + 
            #scale_color_viridis(trans = "log10", direction=-1) + 
            geom_abline(slope = 1, colour="blue") + 
            coord_fixed() +
            xlab("Pileup MAF") + 
            ylab("VCF MAF") + 
            labs(title=f)

        # We save different sizes, as the pixel density changes.
        ggsave(paste0(f,"-", ds, ".png"))
        #ggsave(paste(f,"-large.png", sep=""), width=16, height=15)
    }

    # geom hex produces ugly white dots in the output, 
    # so we have to force it to draw borders of the same color
    # around each hexagon... see https://github.com/tidyverse/ggplot2/issues/2157
    # and https://stackoverflow.com/questions/52006888/ggplot2-geom-hex-white-border
    ggplot(joined_data, aes(x=TOTAL.FREQ.x, y=TOTAL.FREQ.y)) +
        geom_hex(aes(colour = ..count..), bins=50) +
        #scale_fill_viridis(direction=-1) +
        scale_fill_viridis() +
        scale_color_viridis() +
        coord_fixed() +
        xlab("Pileup MAF") +
        ylab("VCF MAF") +
        labs(title=f)
    ggsave(paste0(f,"-hex.png"), width=9, height=7)

    ggplot(joined_data, aes(x=TOTAL.FREQ.x, y=TOTAL.FREQ.y)) +
        geom_hex(aes(colour = ..count..), bins=50) +
        #scale_fill_viridis(trans = "log10", direction=-1) +
        scale_fill_viridis(trans = "log10") +
        scale_color_viridis(trans="log10") +
        coord_fixed() +
        xlab("Pileup MAF") +
        ylab("VCF MAF") +
        labs(title=f)
    ggsave(paste0(f,"-hex-log.png"), width=9, height=7)

    }

    # Histogram of diff between the two
    joined_data$DIFF <- joined_data$TOTAL.FREQ.y - joined_data$TOTAL.FREQ.x
    ggplot(joined_data, aes(x=DIFF)) + 
        geom_histogram(bins=50) +
        xlab("Frequency Diff") +
        labs(title=f)

    ggsave(paste0(f,"-hist.png"), bg="white")

    i <- i+1
}

