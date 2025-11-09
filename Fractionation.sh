#!/bin/bash
#set -euo pipefail

Exit the script if any command fails

###---Step1: Convert SRA to FASTQ-------------------------------###
:<<!
#sample.sra="SRR5048067 SRR5048068 SRR5048097 SRR5048098 SRR4421919 SRR4421920 SRR4421722 SRR4421723"

fasterq-dump is very slow, but consumes little memory, run in parallel

#MAX_JOBS=$(( $(sysctl -n hw.ncpu) - 2 ))
#find . -maxdepth 1 -name "*.sra" | parallel -j $MAX_JOBS "fasterq-dump --split-files {}"

####---Rename sample (add existence check)--------------------------####

Define the rename function

rename_file() {
  if [ -f "$1" ]; then
    mv "$1" "$2"
    echo "Renamed: $1 → $2"
  else
   echo "ERROR: File not found - $1" >&2
    exit 1
  fi
}

Use the function for safe renaming

rename_file "SRR4421722_1.fastq.gz" "K562.LabE.polyA.cyt.rep1.1.fastq.gz"
rename_file "SRR4421722_2.fastq.gz" "K562.LabE.polyA.cyt.rep1.2.fastq.gz"
rename_file "SRR4421723_1.fastq.gz" "K562.LabE.polyA.cyt.rep2.1.fastq.gz"
rename_file "SRR4421723_2.fastq.gz" "K562.LabE.polyA.cyt.rep2.2.fastq.gz"

rename_file "SRR4421919_1.fastq.gz" "K562.LabE.polyA.nuc.rep1.1.fastq.gz"
rename_file "SRR4421919_2.fastq.gz" "K562.LabE.polyA.nuc.rep1.2.fastq.gz"
rename_file "SRR4421920_1.fastq.gz" "K562.LabE.polyA.nuc.rep2.1.fastq.gz"
rename_file "SRR4421920_2.fastq.gz" "K562.LabE.polyA.nuc.rep2.2.fastq.gz"

###---Step2: Run FastQC----------------------------------------------------------###
Out_FastQC_dir="/Users/lixia/data/data/Fractionation_seq/FastQC_labE_trimmed"
mkdir -p "$Out_FastQC_dir"

Safe parallel settings (adjust parallel count based on M1 Max), dynamically calculate parallel jobs

MAX_JOBS=$(( $(sysctl -n hw.ncpu) - 2 ))
find . -maxdepth 1 -name "*.trimmed.fastq.gz" | parallel -j $MAX_JOBS "fastqc {} -o $Out_FastQC_dir"

Run MultiQC to aggregate FastQC results

multiqc "$Out_FastQC_dir" -o "$Out_FastQC_dir" -n "Fractionation_QC_Report"
echo " fastqc completed. Results in: $Out_FastQC_dir"

Loop through all samples###

sample="K562.LabE.polyA.cyt.rep1 K562.LabE.polyA.cyt.rep2 K562.LabE.polyA.nuc.rep1 K562.LabE.polyA.nuc.rep2"
for s in ${sample}
do

###---Step3: fastp removes adaptors-----------------------------------------###
Reads1="${s}.1.fastq.gz"
Reads2="${s}.2.fastq.gz"
fastp 

        -i "$Reads1" -I "$Reads2" \
        -o "${s}.1.trimmed.fastq.gz" \
        -O "${s}.2.trimmed.fastq.gz" 

        --detect_adapter_for_pe 

        --trim_poly_g 

        --correction 

        --length_required 30 

        --qualified_quality_phred 20 

        --unqualified_percent_limit 30 

        --thread 8 

        --html "${s}.fastp.html" \
        --json "${s}.fastp.json"
done

###---Step4: hisat2.alignment---------------------------------------------###

Genomedir provides the path and index filename prefix

Genomedir="/Users/lixia/Data/database/ref_genome/GRCh38.GENCODE/hisat2_index/grch38_tran/genome_tran"
hisat2 

     -x $Genomedir 

     -k 20 

     -p 8 

     --rna-strandness RF 

     --dta 

     -1 $s.1.trimmed.fastq.gz 

     -2 $s.2.trimmed.fastq.gz 

      --no-unal 

      --un-conc-gz $s.hisat2.un.gz 

      --summary-file $s.hisat2.log| 

      tee >(samtools flagstat - > $s.hisat2.flagstat) | 

      samtools sort -@ 2 -m 4G -O BAM | 

      tee $s.hisat2.bam
done

This next step is very slow, can it be processed in parallel?

####---Step5: samtools.picard.filter-----------------------------------------###

Merged steps: Exclude supplementary alignments + Extract unique alignments + Sort

samtools view -H $s.hisat2.bam > $s.head
samtools view -F 3332 $s.hisat2.bam -@ 8 | grep NH:i:1 |cat $s.head -|samtools view -Sb - |samtools sort - -o $s.uniq.bam -@ 8

picard runs single-threaded, consumes high memory; 100M reads need 10G RAM; thus the following steps run single-threaded, no parallel execution is used

Step: Add Read Group

java -jar ~/Data/Ubin/picard.jar AddOrReplaceReadGroups 

    I=$s.uniq.bam 

    O=$s.uniq.rg.bam 

    RGID=$s 

    RGLB=lib1 

    RGPL=ILLUMINA 

    RGPU=unit1 

    RGSM=$s 

    CREATE_INDEX=true

Step: Mark Duplicates

java -jar ~/Data/Ubin/picard.jar MarkDuplicates 

    INPUT=$s.uniq.rg.bam 

    OUTPUT=$s.uniq.nodup.bam 

    METRICS_FILE=$s.dup 

    ASSUME_SORTED=TRUE 

    REMOVE_DUPLICATES=true
 # Index the final file
samtools index -@ 8 $s.uniq.nodup.bam

Clean up intermediate files

rm -f "${s}.head" "${s}.uniq.bam" "${s}.uniq.rg.bam" "${s}.uniq.rg.bai"

###---Step6: Read quantification by HTSeq-count or featureCounts-------------------###
###---Step6: HTSeq-count###

HTSeq-count memory usage is typically 3-5 times the BAM file size.

4 sample BAM sizes: 1.2G, 1.2G, 2.1G, 2.5G; runtime completed in 1h.

Genome_dir="/Users/lixia/Data/database/ref_genome/GRCh38.ENSEMBL"
Genome_ref="Homo_sapiens.GRCh38.84.gtf.gz"
Out_dir_Readsquan="/Users/lixia/Data/data/Fractionation_seq/HTSeqcount_labE"
mkdir -p "${Out_dir_Readsquan}"

sample="K562.LabE.polyA.cyt.rep1 K562.LabE.polyA.cyt.rep2 K562.LabE.polyA.nuc.rep1 K562.LabE.polyA.nuc.rep2"
for s in ${sample}
do

htseq-count 

    --nprocesses 8 

    -f bam 

    -r pos 

    -m intersection-nonempty 

    -s reverse 

    -i gene_id 

    -t exon 

    "$s.uniq.nodup.bam" \
    "${Genome_dir}/${Genome_ref}" \
    1>"${Out_dir_Readsquan}/$s.uniq.exp" \
    2>"${Out_dir_Readsquan}/$s.uniq.err"

Add header --with-header \ # head: "$s.uniq.nodup.bam"

Remove counts for un-annotated genes#

grep __ -v ${Out_dir_Readsquan}/$s.uniq.exp > ${Out_dir_Readsquan}/$s.uniq.clean.exp
done
#rm -rf *.uniq.exp

<<!
###---Step6: featureCounts###

featureCounts memory requirement: loads the GTF file in each loop, consuming memory##

Genome_dir="/Users/lixia/Data/database/ref_genome/GRCh38.ENSEMBL"
Genome_ref="Homo_sapiens.GRCh38.84.gtf.gz"
Out_dir_Readsquan="/Users/lixia/Data/data/Fractionation_seq/featureCounts_labE"
mkdir -p "${Out_dir_Readsquan}"

Process all samples with a single command (no loop needed)##

featureCounts

    -p 

    --countReadPairs 

    -T 8 

    -s 2 

    -g gene_id 

    -t exon 

    -a "${Genome_dir}/${Genome_ref}" 

    -o "${Out_dir_Readsquan}/gene.counts.txt" 

    *.uniq.nodup.bam

Single sample loop processing###

featureCounts

    -p 

    --countReadPairs 

    -T 8 

    -s 2 

    -a "${Genome_dir}/${Genome_ref}" 

    -t exon 

    -g gene_id 

    -o "${Out_dir_Readsquan}/${s}.gene.counts.txt" 

    "${s}.uniq.nodup.bam"
!
###---Step7: Differential Gene Expression, DGE-------------------------------------- ###
###Step7: HTSeq-count + DESeq2###
Out_dir_Readsquan="/Users/lixia/Data/data/Fractionation_seq/HTSeqcount_labE"
cd "${Out_dir_Readsquan}"

echo -en "gene_id\tK562.LabE.polyA.cyt.rep1\tK562.LabE.polyA.cyt.rep2\tK562.LabE.polyA.nuc.rep1\tK562.LabE.polyA.nuc.rep2\n" > K562.LabE.polyA.CytvsNuc.HTseqcount.count.csv
paste K562.LabE.polyA.cyt.rep1.uniq.clean.exp K562.LabE.polyA.cyt.rep2.uniq.clean.exp K562.LabE.polyA.nuc.rep1.uniq.clean.exp K562.LabE.polyA.nuc.rep2.uniq.clean.exp |cut -f1,2,4,6,8 >> K562.LabE.polyA.CytvsNuc.HTseqcount.count.csv

Rscript HTSeq_count_DESeq2.R

#Rscript Fractionation.Deseq2.R
<<!
###Step7: featureCounts + DESeq2###

Counts for all samples

Out_dir_Readsquan="/Users/lixia/Data/data/Fractionation_seq/featureCounts_labE"
cd "${Out_dir_Readsquan}"
echo -en "gene_id\tgene_length\tK562.LabE.polyA.cyt.rep1\tK562.LabE.polyA.cyt.rep2\tK562.LabE.polyA.nuc.rep1\tK562.LabE.polyA.nuc.rep2\n" > K562.LabE.polyA.CytvsNuc.featureCounts.count.csv

Extract the 6 columns mentioned above

tail -n +3 gene.counts.txt | cut -f 1,6-10 >> K562.LabE.polyA.CytvsNuc.featureCounts.count.csv

Counts extraction can also be done using R

Next are the data organization steps

Confirm there are indeed differences between samples, and after data cleanup, start differential analysis

Please check the R script(s) starting with "Check_preparation" and "GeneExpressionDiff" in the "Other-R-code" repository.
