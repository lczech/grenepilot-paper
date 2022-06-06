library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyr)
library(RColorBrewer)
library(data.table)
library(dplyr)

# -------------------------------------------------------------------------
#    Load data
# -------------------------------------------------------------------------

print("Load data")

# List of bona fide snps
bim <- fread( "/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/1001g/515g.bim" , sep=" ", header=FALSE )
colnames(bim)[which(names(bim) == "V2")] <- "CHROM_POS"

seeds <- fread("../freq_seeds_vs_1001g/frequency.csv", sep="\t")
#seeds <- read.csv( "../freq_seeds_vs_1001g/frequency.csv", sep="\t" )
#seeds <- read.csv( "../freq_seeds_vs_1001g/head_seeds.csv", sep="\t" )

# Replace the `TOTAL` in the column names in that table by `Seeds` to make it intuitive.
colnames(seeds) <- gsub("TOTAL", "Seeds", colnames(seeds))

flpools <- fread( "frequency-flowerpools-mpileup.csv", sep="\t" )
#flpools <- read.csv( "frequency-flowerpools-mpileup.csv", sep="\t" )
#flpools <- read.csv( "head-frequency-flowerpools-mpileup.csv", sep="\t" )

#        { "S1",  "Flowerpool1001und2" },
#        { "S2",  "Flowerpool50A" },
#        { "S3",  "Flowerpool50B" },
#        { "S4",  "Flowerpool25A" },
#        { "S5",  "Flowerpool25B" },
#        { "S6",  "Flowerpool10B2" },
#        { "S7",  "Flowerpool5B" },
#        { "S8",  "Plantpool100x" },
#        { "S9",  "Plantpool50a" },
#        { "S10", "Plantpool50b" },
#        { "S11", "Plantpool25a" },
#        { "S12", "Plantpool25b" },
#        { "S13", "Plantpool10b1" },
#        { "S14", "Plantpool10b2" },
#        { "S15", "Plantpool5x" },
#        { "S16", "PlantpoolB12345" }

smpsnums = c(
    "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", 
    "S9", "S10", "S11", "S12", "S13", "S14", "S15", "S16"
)

smpnames = c(
    "Flowerpool1001und2", "Flowerpool50A", "Flowerpool50B", "Flowerpool25A", 
    "Flowerpool25B", "Flowerpool10B2", "Flowerpool5B", "Plantpool100x", 
    "Plantpool50a", "Plantpool50b", "Plantpool25a", "Plantpool25b", 
    "Plantpool10b1", "Plantpool10b2", "Plantpool5x", "PlantpoolB12345"
)

# Replace column names in the flpools to have nice readable names.
# We need to do this backwards, to not replace S16 by S1 etc
for( i in length(smpsnums):1 ) {
    colnames(flpools) <- gsub(smpsnums[i], smpnames[i], colnames(flpools))
}

# Merge using chromosome and position
joined <- merge( seeds, flpools, by.x="CHROM_POS", by.y="CHROM_POS", all=FALSE)
#print(head(joined))

# Remove low cov and invariant sites using seeds.
# We later also remove low cov and invarants from the flower and leaf pools,
# but need to do this on a per sample basis, to just remove their particular low cov and invar rows.
joined <- joined[(joined$Seeds.COV >= 50 & joined$Seeds.REF_CNT > 0 & joined$Seeds.ALT_CNT > 0),]
#print(head(joined))
joined <- as.data.frame(joined)

# We use the initial letters of our two column names to figure out
# what they are: seed (S), flower (F), leaf (P, because that is how they are named...)
get_coltype <- function(str){
    init <- substring(str, 1,1)
    chrs <- c("S", "F", "P")
    strs <- c("Seeds", "Flowers", "Leaves")
    r <- strs[which(chrs == init)]
    return(r)
}

# We also need to translate from sample names to number of individuals.
get_indiv <- function(str) {
    # Simple list, same order as sample names above
    indivcnts = c(
        100, 50, 50, 25, 25, 10, 5, 100, 50, 50, 25, 25, 10, 10, 5, 12345
    )
    r <- indivcnts[which(smpnames == str)]
    return(r)
}

# -------------------------------------------------------------------------
#    Histograms of coverages
# -------------------------------------------------------------------------

print("Make coverage histograms")

make_cov_hist <- function( data, out ) {

    #print(paste0("  hist ", out))
    #print(head(data))

    # get the data that is not larger than 300 cov, as there is not much in there above that.
    cov_1000 <- data[data$COV <= 100, ]

    # Histogram of coverages
    ggplot(cov_1000, aes(x=COV)) +
        geom_histogram(bins=50) +
        xlab("Coverage") +
        xlim(0,100)
        #labs(title=mytitle)

    suppressMessages(ggsave(paste0(out, "-cov-hist.png"), bg="white"))
    #suppressMessages(ggsave(paste0(out, "-cov-hist.pdf"), bg="white"))
    #suppressMessages(ggsave(paste0(out, "-cov-hist.svg"), bg="white"))

}

#if(FALSE){
for( i in 1:16 ) {
    print(paste0( "  At ", smpnames[i] ))
    
    # Make a copy of the data, containing only the columns we are interested in,
    # so that we can remove the rows that we are not interested in.
    base <- joined[, "CHROM_POS", drop=FALSE]
    base$COV <- joined[, paste0(smpnames[i], ".COV")]
    #base$SEED_FREQ <- joined[, "TOTAL.FREQ.X"]
    #base$SMP_FREQ  <- joined[, paste0("S", i, ".FREQ")]

    # Now drop rows with low cov or invariants in the sample.
    #base <- base[( base$COV >= 25 & base$SMP_FREQ > 0 & base$SMP_FREQ < 1 ),]
    #print(head(base))

    # Full data
    make_cov_hist( base, paste0(smpnames[i], "-all") )

    # Subset of only 1001g bona fide snps
    bf <- inner_join( base, bim, by="CHROM_POS")
    make_cov_hist( bf, paste0(smpnames[i], "-bona-fide") )
}

# -------------------------------------------------------------------------
#    Compute correlations
# -------------------------------------------------------------------------

make_corr_plot <- function( df, suffix ) {

    #print(head(df))
    print("Compute correlations")

    # Ridiculously complicated way to create an empty dataframe in R with just
    # column names. Seems to be how this is supposed to be done though. Why, R?!
    rescolumns = c("Pair","Individuals","Coverage","Correlation") 
    result = data.frame(matrix(nrow = 0, ncol = length(rescolumns))) 
    colnames(result) = rescolumns

    compute_correlations <- function ( data, col1, col2 ) {
        print(paste0("  At ", col1, " ", col2))

        # Fucking R somehow loses track of what data it is,
        # so we need to manually convert back to a proper data frame again.
        # Both before and after conversion, it is a dataframe,
        # but it does not work without the conversion... who knows...
        #print(is.data.frame(data))
        data <- as.data.frame(data)
        #print(is.data.frame(data))

        #print(head(data))
        #print(head(data[, c(paste0(col1, ".COV"))]))

        # Make a copy of the data, containing only the columns we are interested in,
        # so that we can remove the rows that we are not interested in.
        base <- data[, "CHROM_POS", drop=FALSE]
        base$COV1  <- data[, c(paste0(col1, ".COV"))]
        base$COV2  <- data[, c(paste0(col2, ".COV"))]
        base$FREQ1 <- data[, c(paste0(col1, ".FREQ"))]
        base$FREQ2 <- data[, c(paste0(col2, ".FREQ"))]
        #print(head(base))

        #print("baaaaaase")
        #print(head(base))

        # Drop rows with low cov or invariant freqs, based on both cols.
        # If one of them is seeds, those have already been filtered, but that's okay.
        base <- base[( base$COV1 >= 25 & base$FREQ1 > 0 & base$FREQ1 < 1 ),]
        base <- base[( base$COV2 >= 25 & base$FREQ2 > 0 & base$FREQ2 < 1 ),]
        #print(head(base))

        #print("baaaaaase")
        #print(head(base))

        # Get lower range subsets of coverage and compute correlations.
        for( minc in seq(25,150,5) ) {
            copy <- data.frame(base)
            copy <- copy[( copy$COV1 >= minc & copy$COV2 >= minc),]
            copy <- na.omit(copy)
       
            #print(head(copy))
 
            #pc = cor( copy$FREQ1, copy$FREQ2, use="complete.obs" )
            pc = cor( copy$FREQ1, copy$FREQ2, use="na.or.complete" )
            #print(paste0( "    minc ", minc, " pc ", pc ))
        
            # Add the result, using a "simple" way to append to a df in R...
            # We here assume that col2 is always the flowers or leaves one, never the seeds.
            # We do not run seeds against seeds correlation, so when calling this, 
            # we just need to have seeds first, see below.
            # Also, R does not pass by reference, and so we need a trick (<<-) to be able
            # to change the result from within a function. There are other tricks, but this one works.
            result[nrow(result) + 1,] <<- list(paste0(get_coltype(col1), " vs ", get_coltype(col2)), get_indiv(col2), minc, pc)
        }
    }

    # Call all the combinations that we want.
    compute_correlations( df, "Seeds", "Plantpool10b2")
    compute_correlations( df, "Seeds", "Plantpool25b")
    compute_correlations( df, "Seeds", "Plantpool50a")
    compute_correlations( df, "Seeds", "Plantpool5x")
    compute_correlations( df, "Seeds", "Flowerpool10B2")
    compute_correlations( df, "Seeds", "Flowerpool25A")
    compute_correlations( df, "Seeds", "Flowerpool50B")
    compute_correlations( df, "Seeds", "Flowerpool5B")
    compute_correlations( df, "Plantpool10b1",    "Flowerpool10B2")
    compute_correlations( df, "Plantpool25a",     "Flowerpool25A")
    compute_correlations( df, "Plantpool50b",     "Flowerpool50A")
    compute_correlations( df, "Plantpool5x",      "Flowerpool5B")

    print(head(result))

    write.csv(result, paste0("correlations-", suffix, ".csv"))

    # -------------------------------------------------------------------------
    #    Plot correlations
    # -------------------------------------------------------------------------

    #ggplot(result, aes(x=Coverage, y=Correlation, group=as.factor(Individuals), color=as.factor(Pair))) +
    ggplot(result, aes(x=Coverage, y=Correlation, color=as.factor(Pair))) +
        geom_line(aes(linetype=as.factor(Individuals))) +
        geom_hline( yintercept=1.0, color="gray") +
        scale_linetype_manual(values=c("dotted", "dotdash", "dashed", "solid")) +
        ylab("Frequency Correlation") +
        labs(color="Comparison", linetype="Individuals")

    suppressMessages(ggsave(paste0("correlations-", suffix, ".png"), bg="white", width=9, height=6))
    suppressMessages(ggsave(paste0("correlations-", suffix, ".pdf"), bg="white", width=9, height=6))
    suppressMessages(ggsave(paste0("correlations-", suffix, ".svg"), bg="white", width=9, height=6))
}

#print(head(joined))
print("Plot all")

# Plot all
make_corr_plot( joined, "all" )

print("Plot bona fide")

# Subset for bona fide snps, and only plot those
#bf <- inner_join( joined, bim, by="CHROM_POS")
bf <- joined[joined$CHROM_POS %in% as.list(bim$CHROM_POS) ,]
make_corr_plot( bf, "bona-fide" )

