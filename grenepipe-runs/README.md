command to just copy the config and sample files:

    rsync -rav -e ssh --exclude='/*/*/*/' --include '*/' --include='config.yaml' --include='samples.tsv' --exclude='*' lczech@memex.carnegiescience.edu:~/lscratch/grenepipe-runs .
