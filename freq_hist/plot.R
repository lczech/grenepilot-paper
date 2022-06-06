# Basic plotting
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
suppressMessages(library(viridis))

make_plots <- function( data, file, mytitle ){

    # Plot a histogram.
    # We always use range 0-1 for comparability.
    # Naturally, this cannot be set directly in R, as of course setting the xlim min to 0
    # removes the leftmost bar from the plot. Totally makes sense. So, trick.
    ggplot(data, aes( x=FREQ )) +
        geom_histogram( bins=100 ) +
        xlab("Frequency") +
        ylab("Count") +
        #xlim( 0, 1 ) +
        coord_cartesian(xlim=c(0, 1)) +
        labs(title=mytitle) #+

    suppressMessages(ggsave(paste0(file,"_hist_1.0.png"), bg="white", width=12, height=4))

    # Also plot folded
    data <- na.omit(data)
    index <- data$FREQ > 0.5
    data$FREQ[index] <- 1.0 - data$FREQ[index]

    ggplot(data, aes( x=FREQ )) +
        geom_histogram( bins=100 ) +
        xlab("Frequency") +
        ylab("Count") +
        #xlim( 0, 1 ) +
        coord_cartesian(xlim=c(0, 0.5)) +
        labs(title=mytitle)

    suppressMessages(ggsave(paste0(file,"_hist_0.5_100.png"), bg="white", width=12, height=4))


    # Two times, because the bins is weird

    ggplot(data, aes( x=FREQ )) +
        geom_histogram( bins=50 ) +
        xlab("Frequency") +
        ylab("Count") +
        #xlim( 0, 1 ) +
        coord_cartesian(xlim=c(0, 0.5)) +
        labs(title=mytitle) 

    suppressMessages(ggsave(paste0(file,"_hist_0.5_50.png"), bg="white", width=12, height=4))

}

# ------------------------------------------
#     1001g
# ------------------------------------------

if(TRUE) {

print("1001g")
data <- read.csv( "../freq_seeds_vs_1001g/plink.frq", sep="" )
#data <- read.csv( "../freq_seeds_vs_1001g/head_plink.frq", sep="" )
colnames(data)[colnames(data) == 'MAF'] <- 'FREQ'
make_plots( data, "1001g", "1001g (VCF)" )

#stop()

# ------------------------------------------
#     231g
# ------------------------------------------

print("231g")
data <- read.csv( "../freq_seeds_vs_1001g/231g.frq", sep="" )
colnames(data)[colnames(data) == 'MAF'] <- 'FREQ'
make_plots( data, "231g", "231g (VCF)" )

}

# ------------------------------------------
#     Seeds HaplotypeCaller
# ------------------------------------------

if(TRUE) {

print("Seeds HaplotypeCaller")
data <- read.csv( "../freq_seeds_comparison/frequency-ath-evo-seeds-haplotypecaller-kv-genotyped.csv", sep="" )
colnames(data)[colnames(data) == 'TOTAL.FREQ'] <- 'FREQ'
data$FREQ =  1.0 - data$FREQ
make_plots( data, "seed_vcf", "Seeds (GATK HaplotypeCaller)" )

# ------------------------------------------
#     Seeds Pileup
# ------------------------------------------

print("Seeds Pileup")
data <- read.csv( "../freq_seeds_comparison/frequency-ath-evo-seeds-pileups-all-merged-units.csv", sep="" )
colnames(data)[colnames(data) == 'TOTAL.FREQ'] <- 'FREQ'
data$FREQ =  1.0 - data$FREQ
make_plots( data, "seed_pileup", "Seeds (Pileup)" )

}

# ------------------------------------------
#     Flowers HaplotypeCaller
# ------------------------------------------

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

# Adding up all non-seed sample ref and alt counts to compute combined freq

print("Flowers HaplotypeCaller")
data <- read.csv( "frequency-francois-vcf.csv", sep="" )
#data <- read.csv( "head-frequency-francois-vcf.csv", sep="" )
#head(data)

data$REF_SUM <- data$S4.REF_CNT + data$S5.REF_CNT + data$S6.REF_CNT + data$S7.REF_CNT + data$S8.REF_CNT + data$S9.REF_CNT + data$S10.REF_CNT + data$S11.REF_CNT + data$S12.REF_CNT
data$ALT_SUM <- data$S4.ALT_CNT + data$S5.ALT_CNT + data$S6.ALT_CNT + data$S7.ALT_CNT + data$S8.ALT_CNT + data$S9.ALT_CNT + data$S10.ALT_CNT + data$S11.ALT_CNT + data$S12.ALT_CNT
data$FREQ <- data$ALT_SUM / (data$REF_SUM + data$ALT_SUM)

make_plots( data, "flower_vcf", "Flowers (GATK HaplotypeCaller)" )

# ------------------------------------------
#     Flowers Pileup
# ------------------------------------------

print("Flowers Pileup")
data <- read.csv( "../af_traj_francois/frequency.csv", sep="" )

data$REF_SUM <- data$S4.REF_CNT + data$S5.REF_CNT + data$S6.REF_CNT + data$S7.REF_CNT + data$S8.REF_CNT + data$S9.REF_CNT + data$S10.REF_CNT + data$S11.REF_CNT + data$S12.REF_CNT
data$ALT_SUM <- data$S4.ALT_CNT + data$S5.ALT_CNT + data$S6.ALT_CNT + data$S7.ALT_CNT + data$S8.ALT_CNT + data$S9.ALT_CNT + data$S10.ALT_CNT + data$S11.ALT_CNT + data$S12.ALT_CNT
data$FREQ <- data$ALT_SUM / (data$REF_SUM + data$ALT_SUM)

make_plots( data, "flower_pileup", "Flowers (Pileup)" )



