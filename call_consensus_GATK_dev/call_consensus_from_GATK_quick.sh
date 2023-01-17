#!/bin/sh
fq1=$1
fq2=$2
sample_name=$3
cut_primer=$4 #1: cut; 0: not-cut;

output_dir='./'
min_cov=10 #the min read depth >= 20

ref_fasta=MPXV.reference.fasta #genome file
primer_bed=MPXV.primer.bed #prime file in bed format
##01. maping file
#conda env  remove -n SARS-CoV-2-Consensus_v1
#conda env create -f SARS-CoV-2-Consensus_v1.yaml
#conda deactivate

bcftools=bcftools
#bwa wucc modified (multiple mapped reads MAPQ==10)
bwa=./call_consensus_GATK_dev/bwa/bwa

fastp -w 48 -i ${fq1} -I ${fq2} -o ${sample_name}.R1.fq.gz -O ${sample_name}.R2.fq.gz -l 60

${bwa} mem -t 96 ${ref_fasta} \
${sample_name}.R1.fq.gz ${sample_name}.R2.fq.gz 2>/dev/null | samtools view -hubS -F4 -@ 4  -F2048  - |samtools sort -@ 8 -o \
${sample_name}.sorted.bam

##02. Filtering low mapping quality reads
picard FilterSamReads I=${sample_name}.sorted.bam \
O=${sample_name}.sorted.before_cliped.bam \
JAVASCRIPT_FILE=./call_consensus_GATK_dev/filter_softcliping_reads.js \
FILTER=includeJavascript

##03. filtering amplicon PCR primers
if (( $cut_primer == 1 ))
then
    echo "Removing primers!"
	samtools ampliconclip -b ${primer_bed} --both-ends -@ 4 \
	-o ${sample_name}.sorted.after_cliped.bam  \
	${sample_name}.sorted.before_cliped.bam
else
    echo "Skipping primer removing!"
	mv ${sample_name}.sorted.before_cliped.bam ${sample_name}.sorted.after_cliped.bam
fi

samtools sort -@ 8 -o ${sample_name}.sorted.clean.bam ${sample_name}.sorted.after_cliped.bam 

## 04. Indexing
samtools index ${sample_name}.sorted.clean.bam

## 05. Calling Depth
samtools depth -a ${sample_name}.sorted.clean.bam > ${sample_name}.sorted.clean.depth
awk -v min_depth="${min_cov}" '{if($3 < min_depth){print $1"\t"$2-1"\t"$2}}' \
${sample_name}.sorted.clean.depth > ${sample_name}.sorted.clean.low_depth_mask.bed

## 06. call snps
mkdir GATK
cd GATK
./call_consensus_GATK_dev/call_snps_quick.sh ../${sample_name}.sorted.clean.bam \
aa aa ${sample_name} ./

cd ./${sample_name}/variants

## 06.1 used majority and high-quality snvs
zcat ${sample_name}.indel.filter.vcf.gz |grep -ivP '^#' | gzip | \
zcat  ${sample_name}.snp.filter.vcf.gz -  | \
grep -ivPw  'my_filter|low_depth|rare_mutation' | \
perl -ne 'if(/multiallelic/){if(/major_allele/){print $_}}else{print $_;}' | \
${bcftools} sort > ${sample_name}.snv.filter.clean.vcf

bgzip -f ${sample_name}.snv.filter.clean.vcf
tabix -f ${sample_name}.snv.filter.clean.vcf.gz
## 06.2 multiallelic sites
zcat ${sample_name}.indel.filter.vcf.gz |grep -ivP '^#' | gzip | \
zcat  ${sample_name}.snp.filter.vcf.gz -  \
|grep -ivPw  'my_filter|low_depth|rare_mutation' |grep -iPw "^#|multiallelic|CHROM" | \
${bcftools} sort > ${sample_name}.multiallelic.clean.vcf
bgzip -f ${sample_name}.multiallelic.clean.vcf
tabix -f ${sample_name}.multiallelic.clean.vcf.gz

## 07. call consensus from vcfs
cd ../../../

${bcftools} query -f '%CHROM\t%POS0\t%END\n' ./GATK/${sample_name}/variants/${sample_name}.snv.filter.clean.vcf.gz \
> ${sample_name}.clean_variants.bed
bedtools subtract -a ${sample_name}.sorted.clean.low_depth_mask.bed -b ${sample_name}.clean_variants.bed \
> ${sample_name}.sorted.clean.low_depth_mask.without_clean_variants.bed

${bcftools} consensus -f ${ref_fasta} \
-m ${sample_name}.sorted.clean.low_depth_mask.without_clean_variants.bed \
-o ${sample_name}.cons.fasta ./GATK/${sample_name}/variants/${sample_name}.snv.filter.clean.vcf.gz

##adding rename cmd
sed "s/>NC_045512.2/>${sample_name}-GATK/" ${sample_name}.cons.fasta > ${sample_name}.cons.fasta.rename
/bin/mv ${sample_name}.cons.fasta.rename ${sample_name}.cons.fasta

## 08. multiple alignment
cat ${ref_fasta}  \
${sample_name}.cons.fasta > commbined.fasta

mafft --auto commbined.fasta > commbined.mafft.fasta

## 09. summary of vaf and gaps

## 09.01 low coverage proportion
echo "Low-depth: `wc -l ${sample_name}.sorted.clean.low_depth_mask.without_clean_variants.bed`/29903"
echo "Variation-Number: `zcat ./GATK/${sample_name}/variants/${sample_name}.snv.filter.clean.vcf.gz |grep -v '^#'|wc -l`"
echo "Multiple allelic Variations-Number: `zcat ./GATK/${sample_name}/variants/${sample_name}.multiallelic.clean.vcf.gz |grep -v '^#'|wc -l`"

## 10 report multiple allelic sites
${bcftools} query -f '%CHROM\t%POS0\t%END\t%REF%END%ALT\t[%AF]\t[%AD]\n' ./GATK/${sample_name}/variants/${sample_name}.multiallelic.clean.vcf.gz \
> ${sample_name}.multiallelic.txt
${bcftools} query -f '%CHROM\t%POS0\t%END\t%REF%END%ALT\t[%AF]\t[%AD]\n' ./GATK/${sample_name}/variants/${sample_name}.snv.filter.clean.vcf.gz \
> ${sample_name}.clean_snvs.txt
