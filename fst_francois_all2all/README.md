We could not do fst on width 1, as the resulting table 
has too many rows for R to handle. How I love R!


awk '{ if( ($1 == 5) && ($2 >= 3167316) && ($2 <= 3185514) ){ print } }' fst-width-1.csv >> fst-width-1-FLC.csv
