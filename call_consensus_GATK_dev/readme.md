## 1. creating conda env
```
conda create env -f SARS-CoV-2-Consensus_v1.yaml
```
## 2. activating conda env

```
conda activate  SARS-CoV-2-Consensus_v1
```

## 3. genome assembly

```
nohup ./call_consensus_from_GATK_quick.sh \
FT100000689_L01_UDB-397-400_1.fq.gz \
FT100000689_L01_UDB-397-400_2.fq.gz \
amplicon_G99_HD_recut 0 &
```

## 4. results

* amplicon_G99_HD_recut.cons.fasta #consensus genome
* amplicon_G99_HD_recut.multiallelic.txt #iSNV
* amplicon_G99_HD_recut.clean_snvs.txt #all SNV
* amplicon_G99_HD_recut.sorted.clean.bam #mapping file
* amplicon_G99_HD_recut.sorted.clean.depth #sequencing depth