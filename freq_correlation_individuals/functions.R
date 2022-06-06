
make_corr_plot_quality <- function( df, suffix ) {
  
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
    base$ALTCNT1 <- data[, c(paste0(col1, ".ALT_CNT"))]
    base$ALTCNT2 <- data[, c(paste0(col2, ".ALT_CNT"))]
    #print(head(base))
    
    #print("baaaaaase")
    #print(head(base))
    
    # Drop rows with low cov or invariant freqs, based on both cols.
    # If one of them is seeds, those have already been filtered, but that's okay.
    base <- base[( base$COV1 >= 10 & base$FREQ1 > 0 & base$FREQ1 < 1),]
    base <- base[( base$COV2 >= 10 & base$FREQ2 > 0 & base$FREQ2 < 1 ),]
    #print(head(base))
    # add filter of MAC
    base <- base[( base$ALTCNT1 > 2 & base$ALTCNT1 > 2 ),]
    
    #print("baaaaaase")
    #print(head(base))
    
    # Get lower range subsets of coverage and compute correlations.
    for( minc in seq(10,200,5) ) {
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
  
  ggsave(paste0("correlations_bonafide2mac-", suffix, ".png"), bg="white", width=9, height=6)
  ggsave(paste0("correlations_bonafide2mac--", suffix, ".pdf"), bg="white", width=9, height=6)
  
}



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
  
  suppressMessages(ggsave(paste0(out, "-cov-hist_qualbonafide2mac.png"), bg="white"))
  #suppressMessages(ggsave(paste0(out, "-cov-hist.pdf"), bg="white"))
  #suppressMessages(ggsave(paste0(out, "-cov-hist.svg"), bg="white"))
  
}