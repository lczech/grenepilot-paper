#Basic plottingi
library(tidyr)
library(ggplot2)
library(cowplot)
#library(tidyverse)
theme_set(theme_cowplot())
suppressMessages(library(viridis))

# Input data
infile="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas/seed_fst_all2all/fst.csv"
#infile="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas/seed_fst_all2all/head.csv"
data = read.table( infile, sep="\t", header=TRUE )
#head(data)

# Turn it into the R long format, putting the previous column names into `Pair`,
# and the values into `FST`, taking all columns with S..S as needed.
data_long <- gather(data, Pair, FST, `S1.S2`:`S7.S8`, factor_key=TRUE)

# Plot with full xlim range
ggplot(data_long, aes(x=FST)) +
    geom_histogram(bins=50) +
    xlab("FST") +
    xlim(0, 1) +
    facet_wrap(~ Pair, ncol=7)

# Need to specify white background here, for whatever reason...
ggsave("histogram-1.0.png", width=16, height=8, bg="white")

# Plot with xlim set to a nice value for this dataset
ggplot(data_long, aes(x=FST)) + 
    geom_histogram(bins=50) +
    xlab("FST") +
    xlim(0, 0.075) +
    facet_wrap(~ Pair, ncol=7)

# Need to specify white background here, for whatever reason...
ggsave("histogram-0.075.png", width=16, height=8, bg="white")
