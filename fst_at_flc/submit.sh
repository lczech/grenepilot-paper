#!/bin/bash

for method in spence-nei spence-hudson kofler karlsson ; do

    cat run-template.sh | sed "s/#METHOD#/${method}/g" > run-${method}.sh
    chmod 755 run-${method}.sh
    sbatch run-${method}.sh
        
done
