conda create env -n artic
conda activate artic
conda install -c conda-forge arctic
artic minion --min_depth 40  --medaka --medaka-model  r941_min_high_g360  --normalise 500 --threads 96 --scheme-directory ./primer_scheme/ --read-file ./nanopore/total.fq.gz  --scheme-version 1 MPXV total

