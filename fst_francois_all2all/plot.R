#Basic plottingi
library(tidyr)
library(ggplot2)
library(cowplot)
#library(tidyverse)
theme_set(theme_cowplot())
suppressMessages(library(viridis))
library(plyr)
library(forcats)
library("wesanderson")
#wes_palette("GrandBudapest1")

# =============================================================================
    print("Read data")
# =============================================================================

# Input data
infile="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas/fst_francois_all2all/fst-width-10k.csv"
#infile="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas/fst_francois_all2all/head.csv"
data = read.table( infile, sep="\t", header=TRUE )
#head(data)

#S1  1_0 56  
#S4  1_1 80  
#S7  1_2 65  
#S10 1_3 19  
#S2  2_0 69  
#S5  2_1 160 
#S8  2_2 205 
#S11 2_3 50  
#S3  3_0 101 
#S6  3_1 200 
#S9  3_2 296 
#S12 3_3 97 

# Turn it into the R long format, putting the previous column names into `Pair`,
# and the values into `FST`, taking all columns with S..S as needed.
data_long <- gather(data, Pair, FST, `S1.S2`:`S11.S12`, factor_key=TRUE)
#data_long <- gather(data, Pair, FST, `S1.S2`:`S11.S12`)


# Add a col with nice group names for those pairs that have nice groups.
# That is, pairs that belong to the same timepoint get a group Tx,
# and pairs that belong to the same replicate get a group Rx.

# All column names / factors
group_from <- c( 
    "S1.S2", "S1.S3", "S1.S4", "S1.S5", "S1.S6", "S1.S7", "S1.S8", "S1.S9", 
    "S1.S10", "S1.S11", "S1.S12", "S2.S3", "S2.S4", "S2.S5", "S2.S6", "S2.S7", 
    "S2.S8", "S2.S9", "S2.S10", "S2.S11", "S2.S12", "S3.S4", "S3.S5", "S3.S6", 
    "S3.S7", "S3.S8", "S3.S9", "S3.S10", "S3.S11", "S3.S12", "S4.S5", "S4.S6", 
    "S4.S7", "S4.S8", "S4.S9", "S4.S10", "S4.S11", "S4.S12", "S5.S6", "S5.S7", 
    "S5.S8", "S5.S9", "S5.S10", "S5.S11", "S5.S12", "S6.S7", "S6.S8", "S6.S9", 
    "S6.S10", "S6.S11", "S6.S12", "S7.S8", "S7.S9", "S7.S10", "S7.S11", "S7.S12", 
    "S8.S9", "S8.S10", "S8.S11", "S8.S12", "S9.S10", "S9.S11", "S9.S12", 
    "S10.S11", "S10.S12", "S11.S12"
)

# Some of them have a "nice" name. For example, samples S1, S2, and S3
# are all of time-point 0, so all pairs of these three get "T0" here;
# samples S1, S4, S7, S10 are all of the same replicate plot, so all
# pairs of these get "R1" here.
group_to <- c( 
    "T0", "T0", "R1", "-", "-", "R1", "-", "-", "R1", "-", "-", "T0", "-", "R2", 
    "-", "-", "R2", "-", "-", "R2", "-", "-", "-", "R3", "-", "-", "R3", "-", "-", 
    "R3", "T1", "T1", "R1", "-", "-", "R1", "-", "-", "T1", "-", "R2", "-", "-", 
    "R2", "-", "-", "-", "R3", "-", "-", "R3", "T2", "T2", "R1", "-", "-", "T2", 
    "-", "R2", "-", "-", "-", "R3", "T3", "T3", "T3"
)

data_long$Group <- mapvalues(data_long$Pair, 
          from=group_from,
          to=group_to
    )
data_long$Group <- factor(data_long$Group, levels = c("R1", "R2", "R3", "T0", "T1", "T2", "T3", "-"))

# Now also let's get some nice sample names.
# It seems there are no easy solutions to do all at once, so we do line by line...
data_long$Name <- data_long$Pair

# Need to do the two digit names first, to not accidentally replace the shorter ones...
data_long$Name <- sub( "S10", "T3R1", data_long$Name )
data_long$Name <- sub( "S11", "T3R2", data_long$Name )
data_long$Name <- sub( "S12", "T3R3", data_long$Name )
data_long$Name <- sub( "S1",  "T0R1", data_long$Name )
data_long$Name <- sub( "S2",  "T0R2", data_long$Name )
data_long$Name <- sub( "S3",  "T0R3", data_long$Name )
data_long$Name <- sub( "S4",  "T1R1", data_long$Name )
data_long$Name <- sub( "S5",  "T1R2", data_long$Name )
data_long$Name <- sub( "S6",  "T1R3", data_long$Name )
data_long$Name <- sub( "S7",  "T2R1", data_long$Name )
data_long$Name <- sub( "S8",  "T2R2", data_long$Name )
data_long$Name <- sub( "S9",  "T2R3", data_long$Name )

# =============================================================================
    print("Histogram")
# =============================================================================

if(TRUE) {

# Plot with full xlim range
ggplot(data_long, aes(x=FST)) +
    geom_histogram(bins=50) +
    xlab("FST") +
    xlim(0, 1) +
    facet_wrap(~ Pair, ncol=6)

# Need to specify white background here, for whatever reason...
ggsave("histogram-1.0.png", width=16, height=24, bg="white")


# Plot with better xlim range
ggplot(data_long, aes(x=FST)) +
    geom_histogram(bins=50) +
    xlab("FST") +
    xlim(0, 0.18) +
    facet_wrap(~ Pair, ncol=6)

# Need to specify white background here, for whatever reason...
ggsave("histogram-0.18.png", width=16, height=24, bg="white")

}

# =============================================================================
#     Boxplot function
# =============================================================================

# Make a plot. The first two arguments subset the columns by value:
# If col == "p", we use colnames to subset the Pair column;
# if col == "g", we use colnames to subset the Group column.
# Also, if groupcolor=TRUE, we use the group to color things.
make_boxplot <- function( col, colnames, title, filename, width=12, height=8, groupcolor=TRUE ) {

    if(col == "p") {
        data_sub <- data_long[data_long$Pair %in% colnames,]
    } else if(col == "g") {
        data_sub <- data_long[data_long$Group %in% colnames,]
    }

    if(TRUE) {

        # Get R to plot them in a nice order by times and replicates.

        #data_sub <- with(data_sub, data_sub[order( Group ),])
        #data_sub <- arrange(data_sub, Group)
        #data_sub$Pair <- factor(data_sub$Pair, levels = data_sub$Pair[order(data_sub$Group)])
        #fct_reorder(data_sub$Pair, data_sub$Group)

        # ... none of that works.
        # Fuck it. Manually ordering.

        group_order <- c(
            "S1.S4", "S1.S7", "S1.S10", "S4.S7", "S4.S10", "S7.S10", "S2.S5", "S2.S8", 
            "S2.S11", "S5.S8", "S5.S11", "S8.S11", "S3.S6", "S3.S9", "S3.S12", "S6.S9", 
            "S6.S12", "S9.S12", "S1.S2", "S1.S3", "S2.S3", "S4.S5", "S4.S6", "S5.S6", 
            "S7.S8", "S7.S9", "S8.S9", "S10.S11", "S10.S12", "S11.S12", "S1.S5", "S1.S6", 
            "S1.S8", "S1.S9", "S1.S11", "S1.S12", "S2.S4", "S2.S6", "S2.S7", "S2.S9", 
            "S2.S10", "S2.S12", "S3.S4", "S3.S5", "S3.S7", "S3.S8", "S3.S10", "S3.S11", 
            "S4.S8", "S4.S9", "S4.S11", "S4.S12", "S5.S7", "S5.S9", "S5.S10", "S5.S12", 
            "S6.S7", "S6.S8", "S6.S10", "S6.S11", "S7.S11", "S7.S12", "S8.S10", "S8.S12", 
            "S9.S10", "S9.S11"
        )
        data_sub$Pair <- factor(data_sub$Pair, levels = group_order)
        #data_sub$Name <- factor(data_sub$Name, levels = group_order)

        # Now, if we want to use that order, while also being able to use a different column
        # (our nice names!) for plotting, we need to to trickery to let R keep that order.
        # Let's turn the Pair into a numeric vector, that forcats and ggplot can
        # then use to order the rows...
        data_sub$PairOrder <- as.numeric(data_sub$Pair)

        # When using group to filter, we always want to color by group.
        groupcolor = TRUE
    }

    if(groupcolor) {
        # Do the plot, colored by Group
        #myplot <- ggplot(data_sub, aes(x=Pair, y=FST, color=Group))
        #myplot <- ggplot(data_sub, aes(x=reorder(Pair, Group), y=FST, color=Group))

        myplot <- ggplot(data_sub, aes(x=fct_reorder(Name, PairOrder), y=FST, color=Group))

    } else {
        myplot <- ggplot(data_sub, aes(x=Pair, y=FST))
    }   

    mean_sub <- mean(data_sub$FST, na.rm = TRUE)

    # Manual assignment of factors to colors, to keep it consistent across plots.
    colorscheme <- c(   
        #"R1" = "#64a1f4",
        #"R2" = "#4a91f2",
        #"R3" = "#3b7dd8",
        "R1" = "#A6CB45",
        "R2" = "#71B238",
        "R3" = "#6A8347",
        "T0" = "#f1c27d",
        "T1" = "#e0ac69",
        "T2" = "#c68642",
        "T3" = "#8d5524",
        "-"  = "#a7adba"
    )
    
    # Now if we just used the above list for plotting, we'd get ALL its colors all the time...
    # So let's subset to the ones that are actually in the plot.
    # AAAALSO, we need to go through them in the order of the factor level, but then filter
    # again by whether they actually appear in our selection, as otherwise, we mess up the order.
    # Wow, R, get your shit together. I guess, there is a deeper understanding to be gained here,
    # but at the moment, this is just annoying...
    mycolors = c()
    for( e in unique(levels(data_sub$Group)) ) {
        if( e %in% unique(data_sub$Group) ) {
            mycolors <- append(mycolors, colorscheme[e])
        }
    }

    # General plot settings for all.
    # We want to reverse the order of sample pair names, so that they are sorted
    # top to bottom. those are the x-axis in ggplot, even after applying the coordiate flip...
    # also, R wants special treatment of those...
    myplot <- myplot +
        geom_boxplot() +
        #scale_color_manual(values = wes_palette(grlen, "GrandBudapest1")) +
        scale_color_manual(values = mycolors) +
        coord_flip() +
        scale_x_discrete(limits=rev) +
        ylim(0.0, 0.3) +
        geom_hline(yintercept = mean_sub, color="gray") +
        annotate("text",  x=-Inf, y = Inf, label = paste("Mean FST:", format(mean_sub, digits=3)), vjust=-1, hjust=1, color="gray") +
        ggtitle(title) +
        xlab("Pair")

    # Facet the plot. Doesn't look good though.
    #if(col=="g"){
    #    myplot <- myplot + facet_wrap( ~ Names )
    #}

    ggsave(paste0(filename,".png"), width=width, height=height, bg="white")

}

# =============================================================================
    print("Boxplots")
# =============================================================================

if(TRUE) {

# Plot just the T0 sample pairs
smp_T0 <- c( "S1.S2", "S1.S3", "S2.S3" )
make_boxplot( "p", smp_T0, "FST between seed generation replicates", "box-T0" )

# Plot all the T1-3 sample pairs
smp_T1_3 <- c( 
    "S4.S5", "S4.S6", "S4.S7", "S4.S8", "S4.S9", "S4.S10", "S4.S11", "S4.S12", 
    "S5.S6", "S5.S7", "S5.S8", "S5.S9", "S5.S10", "S5.S11", "S5.S12", "S6.S7", 
    "S6.S8", "S6.S9", "S6.S10", "S6.S11", "S6.S12", "S7.S8", "S7.S9", "S7.S10", 
    "S7.S11", "S7.S12", "S8.S9", "S8.S10", "S8.S11", "S8.S12", "S9.S10", 
    "S9.S11", "S9.S12", "S10.S11", "S10.S12", "S11.S12"
)
make_boxplot( "p", smp_T1_3, "FST between pools of flowers at different times and in different replicates", "box-T1-3", height=12 )

# Plot all pairs of the same replicate
smp_R13 <- c( "R1", "R2", "R3" )
make_boxplot( "g", smp_R13, "FST between all pairs of pools in each replicate", "box-R13" )

# Plot all pairs of the same time
smp_T03 <- c( "T0", "T1", "T2", "T3" )
make_boxplot( "g", smp_T03, "FST between pairs of pools at each timepoint", "box-T03" )

# Plot all pairs of the same time, only 1-3
smp_T13 <- c( "T1", "T2", "T3" )
make_boxplot( "g", smp_T13, "FST between pairs of flower pools at each timepoint", "box-T13" )

# Plot all pairs of the same time, only 0.
# Should be identical to the above box-t0 plot, but let's check
#smp_T0 <- c( "T0" )
#make_boxplot( "g", smp_T0, "FST between pairs of pools at timepoint 0 (seeds)", "box-T0" )

# Plot all pairs of T0 to T1-3
smp_T0_13 <- c(
    "S1.S4", "S1.S5", "S1.S6", "S1.S7", "S1.S8", "S1.S9", "S1.S10", "S1.S11", 
    "S1.S12", "S2.S4", "S2.S5", "S2.S6", "S2.S7", "S2.S8", "S2.S9", "S2.S10", 
    "S2.S11", "S2.S12", "S3.S4", "S3.S5", "S3.S6", "S3.S7", "S3.S8", "S3.S9", 
    "S3.S10", "S3.S11", "S3.S12"
)
make_boxplot( "p", smp_T0_13, "FST between pairs of pools at timepoint 0 and later timepoints", "box-T0-13", height=12 )

# Plot all pairs of T0 to T1-3 that are of the same replicate
smp_T0_13_R <- c("S1.S4", "S1.S7", "S1.S10", "S2.S5", "S2.S8", "S2.S11", "S3.S6", "S3.S9", "S3.S12")
make_boxplot( "p", smp_T0_13_R, "FST between pairs of pools at timepoint 0 and later timepoints", "box-T0-13-R", height=12, groupcolor=TRUE )


}
