# Basic plotting
library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
suppressMessages(library(viridis))

# R is really weird concerning arrays/lists and index-based access...
# See https://stackoverflow.com/a/43534546/4184258

print("Reading tables")

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

plot_titles <- c(
    "BCFtools (--consensus-caller)",
    "BCFtools (--multiallelic-caller)",
    "freebayes (with known variants)",
    "freebayes (without known variants)",
    "GATK HaplotypeCaller (with known variants)",
    "GATK HaplotypeCaller (without known variants)"
)

# Read all vcf freq tables.
# We do it all at once, in preparation for later, when we probably want to make
# one big facet plot with all of them in one figure.
vcf_freq_data <- list()
i <- 1
for (f in vcf_freq_files) {
    fn <- paste("frequency-ath-evo-seeds-", f, ".csv", sep="")
    #vcf_freq_data[[i]] <- read.table(fn, sep="\t", header=TRUE)
    vcf_freq_data[[i]] <- fread(fn, sep="\t", header=TRUE)

    i <- i+1
}
# head(vcf_freq_data)
# head(vcf_freq_data[[1]])

# Read pileup freq table.
#pileup_freq_data <- read.table(pileup_freq_file, sep="\t", header=TRUE)
pileup_freq_data <- fread(pileup_freq_file, sep="\t", header=TRUE)

################################################################################
#### FILTER FOR QUALITY! Jan 4 2022 Moi
# pileup_freq_data_backup<-pileup_freq_data

pileup_freq_data_mac <- dplyr::filter(pileup_freq_data, TOTAL.ALT_CNT >2)

s<-fread("515g.bim")
pileup_freq_data_mac_qual<-inner_join(pileup_freq_data_mac,s,by=c("CHROM_POS"="V2"))

pileup_freq_data<-pileup_freq_data_mac_qual

################################################################################
# head(pileup_freq_data)

# =================================================================================================
#     Plot
# =================================================================================================

make_plots <- function(target_dir, basename, joined_data, mytitle){
    dir.create(file.path(target_dir), showWarnings = FALSE)

    # Get standard deviation of the diff
    stddev <- sd(joined_data$DIFF)
    pcc <- cor(x=joined_data$TOTAL.FREQ.x, y=joined_data$TOTAL.FREQ.y)

    # allow to skip this part
    if(TRUE) {

        # -------------------------------------------------------------------------
        #     Scatter plot seeds vs caller
        # -------------------------------------------------------------------------

        # Plot the scatter plot of the frequencies.
        # We want to use different dot sizes to see how this affects the plot.
        #for( ds in c(0.01, 0.03, 0.1, 0.3, 1.0)) {
        # for( ds in c(0.1)) {
        # }
        ggplot(joined_data, aes(x=TOTAL.FREQ.x, y=TOTAL.FREQ.y)) +
            geom_vline( xintercept=0, color="gray") +
            geom_hline( yintercept=0, color="gray") +
            geom_point(alpha = 1/25, size=0.2) +
            # geom_point(alpha = 1/100, size=ds) +
            #geom_point(alpha = 1/100, size=ds, aes(color=COV)) +
            #scale_color_viridis(trans = "log10", direction=-1) +
            geom_abline(slope = 1, colour="blue") +
            coord_fixed() +
            xlab("Frequency (bam)") +
            ylab("Frequency (vcf)") +
            labs(title=mytitle) +
            annotate("text", x = 0.02, y = Inf, hjust = 0, vjust=1, label = paste0("SD: ", format(stddev, digits=3), ", r: ", format(pcc, digits=3)))

        # We save different sizes, as the pixel density changes.
        suppressMessages(ggsave(paste0(target_dir, basename, "-scatter.png"), bg="white"))
        #ggsave(paste(f,"-large.png", sep=""), width=16, height=15)

        # # -------------------------------------------------------------------------
        # #     Hex scatter plot seeds vs caller
        # # -------------------------------------------------------------------------
        # 
        # # geom hex produces ugly white dots in the output,
        # # so we have to force it to draw borders of the same color
        # # around each hexagon... see https://github.com/tidyverse/ggplot2/issues/2157
        # # and https://stackoverflow.com/questions/52006888/ggplot2-geom-hex-white-border
        # ggplot(joined_data, aes(x=TOTAL.FREQ.x, y=TOTAL.FREQ.y)) +
        #     #geom_hex(aes(colour = ..count..), bins=50) +
        #     geom_hex(bins=50) +
        #     #scale_fill_viridis(direction=-1) +
        #     scale_fill_viridis(direction=-1) +
        #     scale_color_viridis(direction=-1) +
        #     coord_fixed() +
        #     xlab("Frequency (bam)") +
        #     ylab("Frequency (vcf)") +
        #     labs(title=mytitle)  +
        #     annotate("text", x = 0.02, y = Inf, hjust = 0, vjust=1, label = paste0("SD: ", format(stddev, digits=3), ", r: ", format(pcc, digits=3)))
        # 
        # suppressMessages(ggsave(paste0(target_dir, basename,"-hex.png"), width=9, height=7, bg="white"))
        # 
        # ggplot(joined_data, aes(x=TOTAL.FREQ.x, y=TOTAL.FREQ.y)) +
        #     #geom_hex(aes(colour = ..count..), bins=50) +
        #     geom_hex(bins=50) +
        #     #scale_fill_viridis(trans = "log10", direction=-1) +
        #     scale_fill_viridis(trans = "log10", direction=-1) +
        #     scale_color_viridis(trans="log10", direction=-1) +
        #     coord_fixed() +
        #     xlab("Frequency (bam)") +
        #     ylab("Frequency (vcf)") +
        #     labs(title=mytitle)  +
        #     annotate("text", x = 0.02, y = Inf, hjust = 0, vjust=1, label = paste0("SD: ", format(stddev, digits=3), ", r: ", format(pcc, digits=3)))
        # 
        # suppressMessages(ggsave(paste0(target_dir, basename,"-hex-log.png"), width=9, height=7, bg="white"))

    }

    if(TRUE) {

        # -------------------------------------------------------------------------
        #     Histogram freq diff seeds and caller
        # -------------------------------------------------------------------------

        # Histogram of diff between the two
        ggplot(joined_data, aes(x=DIFF)) +
            geom_histogram(bins=50) +
            xlab("Frequency Difference") +
            labs(title=mytitle) +
            annotate("text", x = -1, y = Inf, hjust = 0, vjust=1, label = paste0("SD: ", format(stddev, digits=3), ", r: ", format(pcc, digits=3)))

        suppressMessages(ggsave(paste0(target_dir, basename,"-hist.png"), bg="white"))

        # Histogram of diff between the two, log scaled
        ggplot(joined_data, aes(x=DIFF)) +
            geom_histogram(bins=50) +
            scale_y_continuous(trans = "log10") +
            xlab("Frequency Difference") +
            labs(title=mytitle) +
            annotate("text", x = -1, y = Inf, hjust = 0, vjust=1, label = paste0("SD: ", format(stddev, digits=3), ", r: ", format(pcc, digits=3)))

        suppressMessages(ggsave(paste0(target_dir, basename,"-hist-log.png"), bg="white"))

    }

    if(TRUE) {

        # -------------------------------------------------------------------------
        #     Histogram frequency of caller
        # -------------------------------------------------------------------------

        # Histogram of frequencies
        ggplot(joined_data, aes(x=TOTAL.FREQ.y)) +
            geom_histogram(bins=50) +
            xlab("Frequency") +
            labs(title=mytitle)

        suppressMessages(ggsave(paste0(target_dir, basename,"-freq-hist.png"), bg="white"))
    }

    if(TRUE) {

        # -------------------------------------------------------------------------
        #     Histogram of seed read coverage
        # -------------------------------------------------------------------------

        # get the data that is not larger than 1000 cov, as we do not want these.
        cov_1000 <- joined_data[joined_data$SEED.COV <= 1000, ]

        # Histogram of frequencies
        ggplot(cov_1000, aes(x=SEED.COV)) +
            geom_histogram(bins=50) +
            xlab("Coverage") +
            labs(title=mytitle)

        suppressMessages(ggsave(paste0(target_dir, basename,"-cov-hist.png"), bg="white"))

    }
}

# =================================================================================================
#     Main loop over tables
# =================================================================================================

# Go through all freq tables and plot their frequency vs the pileup.
i <- 1
for (f in vcf_freq_files) {

    print(paste0("At: ", f))

    # Make a joined table that only contains the rows of SNPs
    # that occur in both the pileup and the vcf.
    joined_data <- merge(
        pileup_freq_data, vcf_freq_data[[i]],
        by.x = "CHROM_POS", by.y = "CHROM_POS", all=FALSE
    )
    joined_data <- na.omit(joined_data)
    joined_data$COV <- joined_data$TOTAL.COV.x + joined_data$TOTAL.COV.y
    joined_data$SEED.COV <- joined_data$TOTAL.COV.x

    # Change ref freq to minor freq
    joined_data$TOTAL.FREQ.x <- 1.0 - joined_data$TOTAL.FREQ.x
    joined_data$TOTAL.FREQ.y <- 1.0 - joined_data$TOTAL.FREQ.y

    joined_data$DIFF <- joined_data$TOTAL.FREQ.y - joined_data$TOTAL.FREQ.x

    # -------------------------------------------------------------------------
    #     All data
    # -------------------------------------------------------------------------

    print(paste0("    all: ", nrow(joined_data)))
    make_plots("all/", f, joined_data, plot_titles[i])

    # -------------------------------------------------------------------------
    #     COV >= 50 & COV <= 100
    # -------------------------------------------------------------------------

    joined_data_50_100 <- joined_data[joined_data$SEED.COV >= 50 & joined_data$SEED.COV <= 100, ]
    print(paste0("    coverage-50-100: ", nrow(joined_data_50_100)))
    make_plots("coverage-50-100/", f, joined_data_50_100, plot_titles[i])

    # -------------------------------------------------------------------------
    #     COV >= 100 & COV <= 250
    # -------------------------------------------------------------------------

    joined_data_100_250 <- joined_data[joined_data$SEED.COV >= 100 & joined_data$SEED.COV <= 250, ]
    print(paste0("    coverage-100-250: ", nrow(joined_data_100_250)))
    make_plots("coverage-100-250/", f, joined_data_100_250, plot_titles[i])

    # -------------------------------------------------------------------------
    #     COV >= 250 & COV <= 500
    # -------------------------------------------------------------------------

    joined_data_250_500 <- joined_data[joined_data$SEED.COV >= 250 & joined_data$SEED.COV <= 500, ]
    print(paste0("    coverage-250-500: ", nrow(joined_data_250_500)))
    make_plots("coverage-250-500/", f, joined_data_250_500, plot_titles[i])

    i <- i+1
}

warnings()
