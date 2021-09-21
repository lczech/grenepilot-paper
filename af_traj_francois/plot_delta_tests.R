library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyr)

# Read the table.
#francois <- read.csv( "frequency.csv", sep="\t" )
francois <- read.csv( "test.csv", sep="\t" )

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
francois$Repl1.Time1 <- francois$S4.FREQ  - francois$S1.FREQ
francois$Repl1.Time2 <- francois$S7.FREQ  - francois$S1.FREQ
francois$Repl1.Time3 <- francois$S10.FREQ - francois$S1.FREQ
francois$Repl2.Time1 <- francois$S5.FREQ  - francois$S2.FREQ
francois$Repl2.Time2 <- francois$S8.FREQ  - francois$S2.FREQ
francois$Repl2.Time3 <- francois$S11.FREQ - francois$S2.FREQ
francois$Repl3.Time1 <- francois$S6.FREQ  - francois$S3.FREQ
francois$Repl3.Time2 <- francois$S9.FREQ  - francois$S3.FREQ
francois$Repl3.Time3 <- francois$S12.FREQ - francois$S3.FREQ



# Turn it into the R long format, putting the previous column names into `Pair`,
# and the values into `FREQ.DIFF`, taking all columns with S..S as needed.
francois_long <- gather(francois, S, FREQ, c(`S1.FREQ`, `S2.FREQ`, `S3.FREQ`, `S4.FREQ`, `S5.FREQ`, `S6.FREQ`, `S7.FREQ`, `S8.FREQ`, `S9.FREQ`, `S10.FREQ`, `S11.FREQ`, `S12.FREQ`), factor_key=TRUE)
francois_long <- gather(francois_long, RT, DIFF, `Repl1.Time1`:`Repl3.Time3`, factor_key=TRUE)

#head(francois_long, n= 1000)
#quit()

# Turn it into the R long format, putting the previous column names into `Pair`,
# and the values into `FREQ.DIFF`, taking all columns with S..S as needed.
#francois_long <- gather(francois, Pair, FREQ.DIFF, `Repl1.Time1`:`Repl3.Time3`, factor_key=TRUE)

myplot <- ggplot(francois_long, aes(x=FREQ, y=DIFF)) + 
    geom_point(alpha = 1/100, size=0.05) + 
    facet_wrap( ~ RT ) +
    #facet_wrap( ~ S + RT ) +
    xlab("Frequency at T0") + 
    ylab("Frequency difference")  
    #ggtitle("R1T2")

# Plot with different sizes and resolutions, as those plots look very different,
# and we kind of have to explore what the "correct" visualization is...

ggsave("delta.png", width=12, height=8, dpi=600)

