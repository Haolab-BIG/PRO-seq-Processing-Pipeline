#!/bin/bash

# --- Stricter Shell Safety ---
# Exit immediately if a command exits with a non-zero status.
set -e
set +H
# Fail a pipeline if any command fails, not just the last one.
set -o pipefail
# Set a sane default for the Internal Field Separator.
IFS=$'\n\t'

# Author: Liu Qian
#
# Description:
#   This script executes a complete TT-seq workflow with a single command. It
#   incorporates several improvements for robustness and flexibility, including
#   shell practices, and more robust file handling.

# --- 1. Default Parameters & Help Function ---
# MODIFICATION: Set a default for the SIF path, can be overridden by -c
SIF_PATH="$(dirname "$0")/proseq.sif"
HELP_MSG="Usage: $0 -o <out_dir> -r <ref_dir> [OPTIONS]

Required:
  -i  Sample name.
  -o  Main output directory.
  -d  Reference data directory.
  -r  Reference rDNA data directory.
  -f  Reference genome annotation file.
  -a  Reference genome file.
  -e  The exons (do not include the first exon of each transcript) position of the chromosome.
  -n  The introns position of the chromosome.
  -c  The singularity container.

Optional:
  -c  Path to the Singularity container (proseq.sif). (Default: same directory as script)
  -h  Display this help message.
"
# --- 2. Parse Command-line Arguments ---
T2G_MAP=""
# MODIFICATION: Added 'c:' to getopts string
while getopts "i:o:d:r:f:a:e:n:c:h" opt; do
  case ${opt} in
    i ) id="${OPTARG}" ;;
    o ) outputdir=$(realpath "${OPTARG}") ;;
    d ) BOWTIE2IDX="${OPTARG}" ;;
    r ) BOWTIE2_RDNAIDX="${OPTARG}" ;;
    f ) fn_gtf="${OPTARG}" ;;
    a ) fn_in="${OPTARG}" ;;
    e ) exon_name="${OPTARG}" ;;
    n ) intron_name="${OPTARG}" ;;
    # MODIFICATION: Handle the new -c argument for the SIF path
    c ) SIF_PATH=$(realpath "${OPTARG}") ;;
    h ) echo "${HELP_MSG}"; exit 0 ;;
    \? ) echo "Invalid option: -${OPTARG}" >&2; echo "${HELP_MSG}"; exit 1 ;;
  esac
done

# Check for mandatory arguments
# --- 3. Setup Environment and Directories ---
# MODIFICATION: The SIF_PATH is now validated here, whether it's the default or user-provided.
if [ ! -f "${SIF_PATH}" ]; then
    echo "Error: Singularity container not found at: ${SIF_PATH}" >&2
    exit 1
fi
SINGULARITY_BASE_CMD="singularity exec --cleanenv -B ${outputdir}:/root ${SIF_PATH}"
echo "================================================="
echo "====== PRO-seq Automated Pipeline Started ======="
echo "================================================="
echo "Output Directory: ${OUT_DIR}"
echo "Singularity Container: ${SIF_PATH}"
echo "================================================="
#id=$1

# BOWTIE2IDX="/mnt/liuq/genome/GRCh38.primary_assembly.genome.bowtie2_index/GRCh38.primary_assembly.genome"
# BOWTIE2_RDNAIDX="/mnt/liuq/genome/GRCh38.primary_assembly.rDNA.bowtie2_index/GRCh38.primary_assembly.rDNA"

#BOWTIE2IDX=$2
#BOWTIE2_RDNAIDX=$3
#outputdir=$4
if [ ! -d "qc/${id}" ]; then 
eval "$SINGULARITY_BASE_CMD mkdir -p qc/${id}"
fi
if [ ! -d "logs" ]; then
eval "$SINGULARITY_BASE_CMD mkdir -p logs"
fi
if [ ! -d "data" ]; then
eval "$SINGULARITY_BASE_CMD mkdir -p data"
fi
if [ ! -d "data" ]; then
eval "$SINGULARITY_BASE_CMD mkdir -p data"
fi
if [ ! -d "data" ]; then
eval "$SINGULARITY_BASE_CMD mkdir -p ${outputdir}"
fi
if [ ! -d "data/${id}" ]; then
eval "$SINGULARITY_BASE_CMD mkdir -p data/${id}"
fi
if [ ! -d "stat/" ]; then
eval "$SINGULARITY_BASE_CMD mkdir stat/"
fi
if [ ! -d "intermediate" ]; then
eval "$SINGULARITY_BASE_CMD mkdir intermediate"
fi
if [ ! -d "intermediate/stat" ]; then
eval "$SINGULARITY_BASE_CMD mkdir intermediate/stat"
fi
if [ ! -d "bed/" ]; then
eval "$SINGULARITY_BASE_CMD mkdir bed/"
fi
if [ ! -d "figure/" ]; then
eval "$SINGULARITY_BASE_CMD mkdir figure/"
fi
if [ ! -d "intermediate/${id}" ]; then
eval "$SINGULARITY_BASE_CMD mkdir intermediate/${id}"
fi
if [ ! -d "intermediate/${id}/HOMER_tag/" ]; then
eval "$SINGULARITY_BASE_CMD mkdir intermediate/${id}/HOMER_tag/"
fi
eval "$SINGULARITY_BASE_CMD fastqc -t 8 data/samples/${id}_R1.fastq.gz -o qc/${id}"
eval "$SINGULARITY_BASE_CMD trim_galore --fastqc --cores 8 -o qc/${id}/ data/samples/${id}_R1.fastq.gz"
eval "$SINGULARITY_BASE_CMD umi_tools extract --extract-method=string --bc-pattern=NNNNNNNN --stdin=qc/${id}/${id}_R1_trimmed.fq.gz --stdout=qc/${id}/${id}_R1_umi_extracted.fastq -L logs/${id}_umi_extracted.log"
eval "$SINGULARITY_BASE_CMD seqtk trimfq -b 1 qc/${id}/${id}_R1_umi_extracted.fastq > qc/${id}/${id}_R1_result_trimmed.fastq"
eval "$SINGULARITY_BASE_CMD seqtk seq -L 10 -r qc/${id}/${id}_R1_result_trimmed.fastq > qc/${id}/${id}_R1_processed.fastq"
eval "$SINGULARITY_BASE_CMD bowtie2 -p 16 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x ${BOWTIE2_RDNAIDX} --rg-id ${id} -U qc/${id}/${id}_R1_processed.fastq --un qc/${id}/${id}_R1_unmap.fastq | $SINGULARITY_BASE_CMD samtools view -bS - -@ 1 | $SINGULARITY_BASE_CMD samtools sort - -@ 1 -o data/prealignments_${id}.bam"
eval "$SINGULARITY_BASE_CMD fastqc -t 8 qc/${id}/${id}_R1_unmap.fastq -o data/${id}"
eval "$SINGULARITY_BASE_CMD bowtie2 -p 16 --very-sensitive -q --end-to-end --un-conc qc/${id}/${id}_hgRemovedTemp.fastq --rg-id ${id} -x ${BOWTIE2IDX} -U qc/${id}/${id}_R1_unmap.fastq -S data/${id}_temp.sam 2>&1 | $SINGULARITY_BASE_CMD tee stat/${id}_align_with_bt2.txt"
eval "$SINGULARITY_BASE_CMD umi_tools dedup --stdin data/${id}_temp.sam --stdout data/${id}_pair_deduped.bam --log logs/${id}_pair_umiTools_Log.log --umi-separator='_'"
eval "$SINGULARITY_BASE_CMD samtools view -h -o data/${id}_pair_deduped.sam data/${id}_pair_deduped.bam"
eval "$SINGULARITY_BASE_CMD awk '\$2 != 4 {print}' data/${id}_pair_deduped.sam | $SINGULARITY_BASE_CMD sed '/XS:/d' | $SINGULARITY_BASE_CMD samtools view -bS - -@ 1 | $SINGULARITY_BASE_CMD samtools sort - -@ 1 -o data/${id}_temp.bam"
eval "$SINGULARITY_BASE_CMD samtools index data/${id}_temp.bam"
eval "$SINGULARITY_BASE_CMD samtools view -q 10 -b -@ 16 -U data/${id}_fail_qc_dups.bam data/${id}_temp.bam > data/${id}_sort.bam"
eval "$SINGULARITY_BASE_CMD samtools index data/${id}_sort.bam"
eval "$SINGULARITY_BASE_CMD samtools idxstats data/${id}_sort.bam | $SINGULARITY_BASE_CMD cut -f 1-2 | $SINGULARITY_BASE_CMD awk '{print \$1, 0, \$2}' | $SINGULARITY_BASE_CMD grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > bed/${id}_chr_sizes.bed"
eval "$SINGULARITY_BASE_CMD samtools view -L bed/${id}_chr_sizes.bed -b -@ 16 data/${id}_sort.bam > data/${id}_noMT.bam"
eval "$SINGULARITY_BASE_CMD samtools index data/${id}_noMT.bam"
eval "$SINGULARITY_BASE_CMD samtools view -H data/${id}_noMT.bam | $SINGULARITY_BASE_CMD grep 'SN:' | $SINGULARITY_BASE_CMD awk -F':' '{print \$2,\$3}' | $SINGULARITY_BASE_CMD awk -F' ' -v OFS='\t' '{print \$1,\$3}' > intermediate/${id}_chr_order.txt"
eval "$SINGULARITY_BASE_CMD cut -f 1 intermediate/${id}_chr_order.txt > intermediate/${id}_chr_keep.txt"
eval "$SINGULARITY_BASE_CMD python scripts/pyTssEnrichment.py -a data/${id}_noMT.bam -b data/plus_TSS.tsv -p ends -c 16 -z -v -s 6 -o stat/${id}_minus_TssEnrichment.txt"
eval "$SINGULARITY_BASE_CMD python scripts/pyTssEnrichment.py -a data/${id}_noMT.bam -b data/minus_TSS.tsv -p ends -c 16 -z -v -s 6 -o stat/${id}_plus_TssEnrichment.txt"
eval "$SINGULARITY_BASE_CMD Rscript scripts/plot_tss_enrichment.R stat/${id}_plus_TssEnrichment.txt stat/${id}_minus_TssEnrichment.txt"
eval "$SINGULARITY_BASE_CMD mv intermediate/stat/${id}_TSSenrichment.pdf intermediate/stat/${id}_TSSenrichment.png figure/"
eval "$SINGULARITY_BASE_CMD bamCoverage -b data/${id}_noMT.bam -o data/${id}_plus.bw --filterRNAstrand forward --binSize 1 --ignoreDuplicates --minMappingQuality 30 --normalizeUsing CPM -p 20"
eval "$SINGULARITY_BASE_CMD bamCoverage -b data/${id}_noMT.bam -o data/${id}_minus.bw --filterRNAstrand reverse --binSize 1 --ignoreDuplicates --minMappingQuality 30 --normalizeUsing CPM -p 20"
eval "$SINGULARITY_BASE_CMD bedtools bamtobed -i data/${id}_noMT.bam > bed/${id}.converted.bed"
eval "$SINGULARITY_BASE_CMD bedtools sort -i bed/${id}.converted.bed > bed/${id}.sorted.bed"
eval "$SINGULARITY_BASE_CMD makeTagDirectory intermediate/${id}/HOMER_tag/ -format bed -forceBED  bed/${id}.sorted.bed"
eval "$SINGULARITY_BASE_CMD findPeaks intermediate/${id}/HOMER_tag/ -style groseq -o intermediate/${id}/transcript.txt -gtf intermediate/${id}/transcript.gtf"
eval "$SINGULARITY_BASE_CMD python scripts/pause_index.py . ${id} ${fn_gtf} ${fn_in}"
eval "$SINGULARITY_BASE_CMD awk '\$10>0 && \$10!=\"NA\"' data/${id}_raw.txt > intermediate/${id}_filtered_raw.txt"
eval "$SINGULARITY_BASE_CMD Rscript scripts/plot_pause_index.r intermediate/${id}_filtered_raw.txt"
eval "$SINGULARITY_BASE_CMD mv intermediate/${id}_pause_index.pdf intermediate/${id}_pause_index.png figure/"
eval "$SINGULARITY_BASE_CMD grep -wf intermediate/${id}_chr_keep.txt ${exon_name} | $SINGULARITY_BASE_CMD bedtools sort -i stdin -faidx intermediate/${id}_chr_order.txt > bed/${id}_exons_sort.bed"
eval "$SINGULARITY_BASE_CMD grep -wf intermediate/${id}_chr_keep.txt ${intron_name} | $SINGULARITY_BASE_CMD bedtools sort -i stdin -faidx intermediate/${id}_chr_order.txt | $SINGULARITY_BASE_CMD bedtools sort -i stdin -faidx intermediate/${id}_chr_order.txt > bed/${id}_introns_sort.bed"
eval "$SINGULARITY_BASE_CMD bedtools coverage -sorted -counts -s -a bed/${id}_exons_sort.bed -b data/${id}_noMT.bam -g intermediate/${id}_chr_order.txt > bed/${id}_exons_coverage.bed"
eval "$SINGULARITY_BASE_CMD bedtools coverage -sorted -counts -s -a bed/${id}_introns_sort.bed -b data/${id}_noMT.bam -g intermediate/${id}_chr_order.txt > bed/${id}_introns_coverage.bed"
ar=$(eval "$SINGULARITY_BASE_CMD samtools view -F 4 -c data/${id}_noMT.bam")
scaling_factor=$(echo "scale=6; $ar/1000000" | eval "$SINGULARITY_BASE_CMD bc")
eval "$SINGULARITY_BASE_CMD awk -v OFS='\t' -v scaling_factor=$scaling_factor '{chrom[\$4] = \$1; if(\$4!=prev4) {chromStart[\$4] = \$2} strand[\$4] = \$6; readCount[\$4] += \$7; exonCount[\$4] += 1; geneSizeKB[\$4] += (sqrt((\$3-\$2+0.00000001)^2)/1000); gene[\$4] = \$4; chromEnd[\$4]=\$3; prev4=\$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/scaling_factor)/geneSizeKB[a], strand[a] }}' bed/${id}_exons_coverage.bed | $SINGULARITY_BASE_CMD awk '\$5>0' | $SINGULARITY_BASE_CMD sort -k4 > bed/${id}_exons_rpkm.bed"
eval "$SINGULARITY_BASE_CMD awk -v OFS='\t' -v scaling_factor=$scaling_factor '{chrom[\$4] = \$1; if(\$4!=prev4) {chromStart[\$4] = \$2} strand[\$4] = \$6; readCount[\$4] += \$7; exonCount[\$4] += 1; geneSizeKB[\$4] += (sqrt((\$3-\$2+0.00000001)^2)/1000); gene[\$4] = \$4; chromEnd[\$4]=\$3; prev4=\$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/scaling_factor)/geneSizeKB[a], strand[a] }}' bed/${id}_introns_coverage.bed | $SINGULARITY_BASE_CMD awk '\$5>0' | $SINGULARITY_BASE_CMD sort -k4 > bed/${id}_introns_rpkm.bed"
eval "$SINGULARITY_BASE_CMD cat bed/${id}_exons_rpkm.bed bed/${id}_introns_rpkm.bed | $SINGULARITY_BASE_CMD awk '{print \$4}' > intermediate/${id}_all_genes.txt"
eval "$SINGULARITY_BASE_CMD awk -F';' '\$1 in first{print first[\$1] \$0; first[\$1]=\"\"; next} {first[\$1]=\$0 ORS}' intermediate/${id}_all_genes.txt | $SINGULARITY_BASE_CMD awk '!NF {print;next}; !(\$0 in a) {a[\$0];print}' > bed/${id}_shared_genes.txt"
# awk -F';' '$1 in first{print first[$1] $0; first[$1]=""; next} {first[$1]=$0 ORS}' intermediate/${id}_all_genes.txt | awk '!NF {print;next}; !($0 in a) {a[$0];print}' > bed/${id}_shared_genes.txt
eval "$SINGULARITY_BASE_CMD awk -F '\t' 'NR==FNR {id[\$1]; next} \$4 in id' bed/${id}_shared_genes.txt bed/${id}_exons_rpkm.bed > bed/${id}_shared_exons.bed"
eval "$SINGULARITY_BASE_CMD awk -F '\t' 'NR==FNR {id[\$1]; next} \$4 in id' bed/${id}_shared_genes.txt bed/${id}_introns_rpkm.bed > bed/${id}_shared_introns.bed"
eval "$SINGULARITY_BASE_CMD awk 'BEGIN{FS=OFS=\"\t\"} FNR>0 && FNR==NR{a[\$4]=\$4 OFS \$0; next} FNR>0{print \$0,a[\$4]?a[\$4]:\"\t\"}' bed/${id}_shared_exons.bed bed/${id}_shared_introns.bed | $SINGULARITY_BASE_CMD awk 'BEGIN { OFS=\"\t\" } {print \$8, \$9, \$10, \$11, (\$12/\$5), \$13}' | $SINGULARITY_BASE_CMD sort -k1,1 -k2,2n > bed/${id}_exon_intron_ratios.bed"
eval "$SINGULARITY_BASE_CMD Rscript scripts/PEPPRO.R mrna -i bed/${id}_exon_intron_ratios.bed --annotate"
eval "$SINGULARITY_BASE_CMD mv bed/${id}_mRNA_contamination.pdf bed/${id}_mRNA_contamination.png figure/"
eval "$SINGULARITY_BASE_CMD multiqc data/ -o data/multiqc_report/${id} --force"