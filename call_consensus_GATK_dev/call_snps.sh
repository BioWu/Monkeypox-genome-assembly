
#！/bin/bash
#mapping BAM to vcf
#软件, reference genome路径，根据实际情况修改
gatk=/gpfs1/wucc/anaconda3/envs/gatk4/bin/gatk
java_options="-Xmx40G"
samtools=/gpfs1/wucc/anaconda3/envs/biobase_py27/bin/samtools
bcftools=/gpfs1/wucc/anaconda3/envs/bcftools/bin/bcftools

ref_genome=/gpfs1/wucc/MonkeyPox-2022-05-21/02.sequencing/ref/MPXV.reference.fasta

#参数
BAM=$1 #mapping后的BAM文件
RGID=$2 #read group ID
library=$3 #测序文库编号
sample=$4 #sample名称
outdir=$5 #输出文件路径

#按样本设置目录
outdir=${outdir}/${sample}
#output diretory

if [ ! -d $outdir/variants ]
then mkdir -p $outdir/variants
fi
if [ ! -d $outdir/processBAM ]
then mkdir -p $outdir/processBAM
fi

#add RG & sort
time $gatk --java-options ${java_options} \
    AddOrReplaceReadGroups \
    --INPUT $BAM \
    --OUTPUT $outdir/processBAM/${sample}.RG.bam \
    --RGID $RGID \
    --RGLB $library \
    --RGPL illumina \
    --RGPU snpcall \
    --RGSM ${sample} && echo "** ADD RG done **"

#split N Trim
## do not markdup
# --I $outdir/processBAM/${sample}.RG.markdup.bam \
time $gatk  --java-options ${java_options} \
    SplitNCigarReads \
    --R $ref_genome \
    --I $outdir/processBAM/${sample}.RG.bam \
    --O $outdir/processBAM/${sample}.RG.markdup.split.bam && echo "** ${sample}.RG.markdup.bam split N done **"
#HaplotypeCaller

# time $gatk  --java-options ${java_options} \
#     HaplotypeCaller \
#     -R $ref_genome \
#     -I $outdir/processBAM/${sample}.RG.markdup.split.bam \
#     -O $outdir/variants/${sample}.raw.vcf.gz \
#     -dont-use-soft-clipped-bases \
#     -A AlleleFraction \
#     -A FisherStrand \
#     -A StrandOddsRatio \
#     -bamout $outdir/variants/${sample}.hap.bam \
#     --standard-min-confidence-threshold-for-calling 20 \
#     --force-active true   --adaptive-pruning true  \
#     --max-reads-per-alignment-start 500 \
#     --native-pair-hmm-threads 16 \
#     --max-assembly-region-size  200 \
#     --kmer-size 13 --kmer-size 15 \
#     --kmer-size 17 --kmer-size 19 \
#     --kmer-size 21  && echo "** ${sample}.vcf done **"

#--minimum-mapping-quality 0 Because of the ITR in the begining.

time $gatk  --java-options ${java_options} \
    Mutect2 \
    -R $ref_genome \
    -I $outdir/processBAM/${sample}.RG.markdup.split.bam \
    -O $outdir/variants/${sample}.raw.vcf.gz \
    -dont-use-soft-clipped-bases \
    -A AlleleFraction \
    -A FisherStrand \
    -A StrandOddsRatio \
    --max-reads-per-alignment-start 1000 \
    --max-assembly-region-size  200 \
    -mnp-dist 0 -min-AF 0.1 \
    --max-num-haplotypes-in-population 512 \
    --native-pair-hmm-threads 16 \
    --kmer-size 10 --kmer-size 15 \
	--minimum-mapping-quality 0 \
    --kmer-size 19 --kmer-size 25  && echo "** ${sample}.vcf done **"

    # -bamout $outdir/variants/${sample}.hap.bam \
    # --max-assembly-region-size  300 \
    # --force-active true  \
## split sites (multiple allelic)
time ${bcftools} norm -f ${ref_genome} -m - -O z $outdir/variants/${sample}.raw.vcf.gz -o $outdir/variants/${sample}.vcf.gz
time $gatk IndexFeatureFile -I $outdir/variants/${sample}.vcf.gz

#选择SNP，并filter
time $gatk SelectVariants \
    -select-type SNP \
    -V $outdir/variants/${sample}.vcf.gz \
    -O $outdir/variants/${sample}.snp.vcf.gz

#选择indel，并filter
time $gatk SelectVariants \
    -select-type INDEL \
    -V $outdir/variants/${sample}.vcf.gz \
    -O $outdir/variants/${sample}.indel.vcf.gz


time $gatk VariantFiltration \
    -V $outdir/variants/${sample}.snp.vcf.gz \
    -O $outdir/variants/${sample}.snp.filter.vcf.gz \
    -R $ref_genome \
    --filter-expression ' FS > 60.0  ' \
    --filter-name "my_filter" \
    -G-filter 'AF < 0.2' \
    -G-filter-name "rare_mutation" \
    -G-filter 'AF < 0.8 && AF > 0.2' \
    -G-filter-name "multiallelic" \
    -G-filter 'AF > 0.5' \
    -G-filter-name "major_allele" \
    -G-filter 'DP < 10' \
    -G-filter-name "low_depth"   && echo "** ${sample}.snp filter done **"


time $gatk VariantFiltration \
    -V $outdir/variants/${sample}.indel.vcf.gz \
    -O $outdir/variants/${sample}.indel.filter.vcf.gz \
    -R $ref_genome \
    --filter-expression ' FS > 60.0 ' \
    --filter-name "my_filter" \
    -G-filter 'AF < 0.2' \
    -G-filter-name "rare_mutation" \
    -G-filter 'AF < 0.8 && AF > 0.2' \
    -G-filter-name "multiallelic" \
    -G-filter 'AF > 0.5' \
    -G-filter-name "major_allele" \
    -G-filter 'DP < 10' \
    -G-filter-name "low_depth"   && echo "** ${sample}.indel filter done **"
