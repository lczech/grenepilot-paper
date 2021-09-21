#!/bin/bash

head -n 1 frequency.csv > test.csv
egrep "^1" frequency.csv | head -n 10000 >> test.csv
egrep "^2" frequency.csv | head -n 10000 >> test.csv
egrep "^3" frequency.csv | head -n 10000 >> test.csv
egrep "^4" frequency.csv | head -n 10000 >> test.csv
egrep "^5" frequency.csv | head -n 10000 >> test.csv


