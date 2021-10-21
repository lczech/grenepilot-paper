#!/bin/bash

IN="frequency.csv"

head -n 1  $IN > test.csv
egrep "^1" $IN | head -n 10000 >> test.csv
egrep "^2" $IN | head -n 10000 >> test.csv
egrep "^3" $IN | head -n 10000 >> test.csv
egrep "^4" $IN | head -n 10000 >> test.csv
egrep "^5" $IN | head -n 10000 >> test.csv


