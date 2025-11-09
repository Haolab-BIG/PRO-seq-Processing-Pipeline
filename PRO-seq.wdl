version 1.1

workflow main {
    input {
        String sample
        File FQ1
        File SIF
        String BOWTIE2_RDNAIDX
        String BOWTIE2IDX
        File plus_TSS
        File minus_TSS
        File fn_gtf
        File fa_in
        File exon_name
        File intron_name
        String PWD
    }

    call qc_before_filtered {
        input:
            sample=sample,
            FQ1=FQ1,
            PWD=PWD,
            SIF=SIF
    }

    call trim_fastq {
        input:
            html=qc_before_filtered.html,
            sample=sample,
            FQ1=FQ1,
            PWD=PWD,
            SIF=SIF
    }

    call umi_extracted {
        input:
            trimmed_FQ1=trim_fastq.trimmed_FQ1,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call remove_additional_nucleotide {
        input:
            umi_extracted_FQ1=umi_extracted.umi_extracted_FQ1,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call reverse_complement {
        input:
            remove_additional_nuc_FQ1=remove_additional_nucleotide.remove_additional_nuc_FQ1,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call prealignments {
        input:
            rev_compl_FQ1=reverse_complement.rev_compl_FQ1,
            BOWTIE2_RDNAIDX=BOWTIE2_RDNAIDX,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call qc_after_filtered {
        input:
            NON_RRNA_FQ1=prealignments.NON_RRNA_FQ1,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call primary_single_alignments {
        input:
            NON_RRNA_FQ1=prealignments.NON_RRNA_FQ1,
            sample=sample,
            BOWTIE2IDX=BOWTIE2IDX,
            PWD=PWD,
            SIF=SIF
    }

    call umi_dedup {
        input:
            temp_sam=primary_single_alignments.temp_sam,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call bam2sam {
        input:
            pair_deduped_bam=umi_dedup.pair_deduped_bam,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call Collapsing {
        input:
            pair_deduped_sam=bam2sam.pair_deduped_sam,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call index_temp_BAM {
        input:
            temp_BAM=Collapsing.temp_BAM,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call sorted_BAM {
        input:
            temp_BAM_index=index_temp_BAM.temp_BAM_index,
            sample=sample,
            temp_BAM=Collapsing.temp_BAM,
            PWD=PWD,
            SIF=SIF
    }

    call index_sorted_BAM {
        input:
            sort_BAM=sorted_BAM.sort_BAM,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call MT_Removed {
        input:
            sort_BAM_index=index_sorted_BAM.sort_BAM_index,
            sort_BAM=sorted_BAM.sort_BAM,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call GET_noMT_BAM {
        input:
            chr_sizes_bed=MT_Removed.chr_sizes_bed,
            sort_BAM=sorted_BAM.sort_BAM,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call index_noMT_BAM {
        input:
            noMT_BAM=GET_noMT_BAM.noMT_BAM,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call chr_order {
        input:
            noMT_BAM_index=index_noMT_BAM.noMT_BAM_index,
            noMT_BAM=GET_noMT_BAM.noMT_BAM,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call chr_keep {
        input:
            chr_order_txt=chr_order.chr_order_txt,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call TssEnrichment {
        input:
            noMT_BAM_index=index_noMT_BAM.noMT_BAM_index,
            noMT_BAM=GET_noMT_BAM.noMT_BAM,
            plus_TSS=plus_TSS,
            minus_TSS=minus_TSS,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call plot_tss_enrichment {
        input:
            plus_TssEnrichment=TssEnrichment.plus_TssEnrichment,
            minus_TssEnrichment=TssEnrichment.minus_TssEnrichment,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call tidy_tss_enrichment {
        input:
            pdf=plot_tss_enrichment.pdf,
            png=plot_tss_enrichment.png,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call Generate_BigWig {
        input:
            noMT_BAM_index=index_noMT_BAM.noMT_BAM_index,
            noMT_BAM=GET_noMT_BAM.noMT_BAM,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call converted_bed {
        input:
            noMT_BAM_index=index_noMT_BAM.noMT_BAM_index,
            noMT_BAM=GET_noMT_BAM.noMT_BAM,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call sorted_bed {
        input:
            convert_bed=converted_bed.convert_bed,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call HOMER_makeTagDirectory {
        input:
            sorted_bed=sorted_bed.sorted_bed,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call HOMER_findPeaks {
        input:
            HOMER_tag=HOMER_makeTagDirectory.HOMER_tag,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call pause_index {
        input:
            sorted_bed=sorted_bed.sorted_bed,
            convert_bed=converted_bed.convert_bed,
            fn_gtf=fn_gtf,
            fa_in=fa_in,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call filtered_raw_txt {
        input:
            raw_txt=pause_index.raw_txt,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call plot_pause_index {
        input:
            filtered_raw_txt=filtered_raw_txt.filtered_raw_txt,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call tidy_pause_index {
        input:
            pdf=plot_pause_index.pdf,
            png=plot_pause_index.png,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call exons_introns_sort {
        input:
            chr_keep_txt=chr_keep.chr_keep_txt,
            chr_order_txt=chr_order.chr_order_txt,
            exon_name=exon_name,
            intron_name=intron_name,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call exons_introns_coverage {
        input:
            exons_sort=exons_introns_sort.exons_sort,
            introns_sort=exons_introns_sort.introns_sort,
            noMT_BAM_index=index_noMT_BAM.noMT_BAM_index,
            chr_order_txt=chr_order.chr_order_txt,
            noMT_BAM=GET_noMT_BAM.noMT_BAM,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call exons_introns_rpkm {
        input:
            exons_coverage=exons_introns_coverage.exons_coverage,
            introns_coverage=exons_introns_coverage.introns_coverage,
            noMT_BAM_index=index_noMT_BAM.noMT_BAM_index,
            noMT_BAM=GET_noMT_BAM.noMT_BAM,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call all_genes {
        input:
            exons_rpkm=exons_introns_rpkm.exons_rpkm,
            introns_rpkm=exons_introns_rpkm.introns_rpkm,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call shared_genes {
        input:
            all_genes_txt=all_genes.all_genes_txt,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call shared_exons_introns {
        input:
            shared_genes_txt=shared_genes.shared_genes_txt,
            exons_rpkm_bed=exons_introns_rpkm.exons_rpkm,
            introns_rpkm_bed=exons_introns_rpkm.introns_rpkm,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call exon_intron_ratios {
        input:
            shared_exons=shared_exons_introns.shared_exons,
            shared_introns=shared_exons_introns.shared_introns,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call plot_mRNA_contamination {
        input:
            exon_intron_ratios=exon_intron_ratios.exon_intron_ratios,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call tidy_mRNA_contamination {
        input:
            pdf=plot_mRNA_contamination.pdf,
            png=plot_mRNA_contamination.png,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    call Final_MultiQC {
        input:
            noMT_BAM_index=index_noMT_BAM.noMT_BAM_index,
            sample=sample,
            PWD=PWD,
            SIF=SIF
    }

    output {
        File html=qc_after_filtered.html
        File tidy_mRNA_contamination_pdf=tidy_mRNA_contamination.mRNA_contamination_pdf
        File tidy_mRNA_contamination_png=tidy_mRNA_contamination.mRNA_contamination_png
        File tidy_tss_enrichment_pdf=tidy_tss_enrichment.TSSenrichment_pdf
        File tidy_tss_enrichment_png=tidy_tss_enrichment.TSSenrichment_png
        File plus_bw=Generate_BigWig.plus_bw
        File minus_bw=Generate_BigWig.minus_bw
        File HOMER_findPeaks_txt=HOMER_findPeaks.txt
        File HOMER_findPeaks_gtf=HOMER_findPeaks.gtf
        File tidy_pause_index_pdf=tidy_pause_index.pause_index_pdf
        File tidy_pause_index_png=tidy_pause_index.pause_index_png
        File multiqc_report_html=Final_MultiQC.multiqc_report_html
    }
}

task qc_before_filtered {
    input {
        String sample
        File FQ1
        File SIF
        String PWD
    }

    command <<<
        if [ ! -d "${PWD}/qc/${sample}" ]; then
        mkdir -p ${PWD}/qc/${sample}
        fi
        singularity exec --bind ${PWD}:/root ${SIF} fastqc -t 8 ${FQ1} -o qc/${sample}
    >>>

    output {
        File html="${PWD}/qc/${sample}/${sample}_R1_fastqc.html"
    }
}

task trim_fastq {
    input {
        String sample
        File FQ1
        File SIF
        File html
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} trim_galore --fastqc --cores 8 -o qc/${sample} ${FQ1}
    >>>

    output {
        File trimmed_FQ1="${PWD}/qc/${sample}/${sample}_R1_trimmed.fq.gz"
    }
}

task umi_extracted {
    input {
        String sample
        File trimmed_FQ1
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} umi_tools extract \
        --extract-method=string \
        --bc-pattern=NNNNNNNN \
        --stdin=${trimmed_FQ1} \
        --stdout=qc/${sample}/${sample}_R1_umi_extracted.fastq
    >>>

    output {
        File umi_extracted_FQ1="${PWD}/qc/${sample}/${sample}_R1_umi_extracted.fastq"
    }
}

task remove_additional_nucleotide {
    input {
        String sample
        File umi_extracted_FQ1
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} seqtk trimfq -b 1 ${umi_extracted_FQ1} > qc/${sample}/${sample}_R1_result_trimmed.fastq
    >>>

    output {
        File remove_additional_nuc_FQ1="${PWD}/qc/${sample}/${sample}_R1_result_trimmed.fastq"
    }
}

task reverse_complement {
    input {
        String sample
        File remove_additional_nuc_FQ1
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} seqtk seq -L 10 -r ${remove_additional_nuc_FQ1} > qc/${sample}/${sample}_R1_processed.fastq
    >>>

    output {
        File rev_compl_FQ1="${PWD}/qc/${sample}/${sample}_R1_processed.fastq"
    }
}

task prealignments {
    input {
        String sample
        File rev_compl_FQ1
        File SIF
        String BOWTIE2_RDNAIDX
        String PWD
    }

    command <<<
        if [ ! -d "${PWD}/data/" ]; then
        mkdir ${PWD}/data/
        fi
        singularity exec --bind ${PWD}:/root ${SIF} bowtie2 -p 16 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x ${BOWTIE2_RDNAIDX} --rg-id ${sample} -U ${rev_compl_FQ1} --un qc/${sample}/${sample}_R1_unmap.fastq | \
        singularity exec --bind ${PWD}:/root ${SIF} samtools view -bS - -@ 1 | \
        singularity exec --bind ${PWD}:/root ${SIF} samtools sort - -@ 1 -o data/prealignments_${sample}.bam
    >>>

    output {
        File NON_RRNA_FQ1="${PWD}/qc/${sample}/${sample}_R1_unmap.fastq"
    }
}

task qc_after_filtered {
    input {
        String sample
        File NON_RRNA_FQ1
        File SIF
        String PWD
    }

    command <<<
        if [ ! -d "${PWD}/data/${sample}" ]; then
        mkdir ${PWD}/data/${sample}
        fi
        singularity exec --bind ${PWD}:/root ${SIF} fastqc -t 8 ${NON_RRNA_FQ1} -o data/${sample}
    >>>

    output {
        File html="${PWD}/data/${sample}/${sample}_R1_unmap_fastqc.html"
    }
}

task primary_single_alignments {
    input {
        String sample
        File NON_RRNA_FQ1
        File SIF
        String BOWTIE2IDX
        String PWD
    }

    command <<<
        if [ ! -d "${PWD}/stat/" ]; then
        mkdir ${PWD}/stat/
        fi
        singularity exec --bind ${PWD}:/root ${SIF} bowtie2 -p 16 --very-sensitive -q --end-to-end --un-conc qc/${sample}/${sample}_hgRemovedTemp.fastq --rg-id ${sample} -x ${BOWTIE2IDX} -U ${NON_RRNA_FQ1} -S data/${sample}_temp.sam 2>&1 | \
        singularity exec --bind ${PWD}:/root ${SIF} tee stat/${sample}_align_with_bt2.txt
    >>>

    output {
        File temp_sam="${PWD}/data/${sample}_temp.sam"
    }
}

task umi_dedup {
    input {
        String sample
        File temp_sam
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} umi_tools dedup --stdin ${temp_sam} --stdout data/${sample}_pair_deduped.bam --umi-separator="_"
    >>>

    output {
        File pair_deduped_bam="${PWD}/data/${sample}_pair_deduped.bam"
    }
}

task bam2sam {
    input {
        String sample
        File pair_deduped_bam
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} samtools view -h -o ${pair_deduped_sam} data/${sample}_pair_deduped.sam
    >>>

    output {
        File pair_deduped_sam="${PWD}/data/${sample}_pair_deduped.sam"
    }
}

task Collapsing {
    input {
        String sample
        File pair_deduped_sam
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} awk '$2 != 4 {print}' ${pair_deduped_sam} | sed '/XS:/d' | \
        singularity exec --bind ${PWD}:/root ${SIF} samtools view -bS - -@ 1 | \
        singularity exec --bind ${PWD}:/root ${SIF} samtools sort - -@ 1 -o data/${sample}_temp.bam
    >>>

    output {
        File temp_BAM="${PWD}/data/${sample}_temp.bam"
    }
}

task index_temp_BAM {
    input {
        String sample
        File temp_BAM
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} samtools index ${temp_BAM}
    >>>

    output {
        File temp_BAM_index="${PWD}/data/${sample}_temp.bam.bai"
    }
}

task sorted_BAM {
    input {
        String sample
        File temp_BAM_index
        File temp_BAM
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} samtools view -q 10 -b -@ 16 -U data/${sample}_fail_qc_dups.bam ${temp_BAM} > data/${sample}_sort.bam
    >>>

    output {
        File sort_BAM="${PWD}/data/${sample}_sort.bam"
    }
}

task index_sorted_BAM {
    input {
        String sample
        File sort_BAM
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} samtools index ${sort_BAM}
    >>>

    output {
        File sort_BAM_index="${PWD}/data/${sample}_sort.bam.bai"
    }
}

task MT_Removed {
    input {
        String sample
        File sort_BAM_index
        File sort_BAM
        File SIF
        String PWD
    }

    command <<<
        if [ ! -d "${PWD}/bed" ]; then
        mkdir ${PWD}/bed
        fi
        singularity exec --bind ${PWD}:/root ${SIF} samtools idxstats ${sort_BAM} | cut -f 1-2 | \
        awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > bed/${sample}_chr_sizes.bed
    >>>

    output {
        File chr_sizes_bed="${PWD}/bed/${sample}_chr_sizes.bed"
    }
}

task GET_noMT_BAM {
    input {
        String sample
        File chr_sizes_bed
        File sort_BAM
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} samtools view -L ${chr_sizes_bed} -b -@ 16 ${sort_BAM} > data/${sample}_noMT.bam
    >>>

    output {
        File noMT_BAM="${PWD}/data/${sample}_noMT.bam"
    }
}

task index_noMT_BAM {
    input {
        String sample
        File noMT_BAM
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} samtools index ${noMT_BAM}
    >>>

    output {
        File noMT_BAM_index="${PWD}/data/${sample}_noMT.bam.bai"
    }
}

task chr_order {
    input {
        String sample
        File noMT_BAM_index
        File noMT_BAM
        File SIF
        String PWD
    }

    command <<<
        if [ ! -d "${PWD}/intermediate" ]; then
        mkdir ${PWD}/intermediate
        fi
        singularity exec --bind ${PWD}:/root ${SIF} samtools view -H ${noMT_BAM} | grep 'SN:' | awk -F':' '{print $2,$3}' | \
        awk -F' ' -v OFS='\t' '{print $1,$3}' > intermediate/${sample}_chr_order.txt
    >>>

    output {
        File chr_order_txt="${PWD}/intermediate/${sample}_chr_order.txt"
    }
}

task chr_keep {
    input {
        String sample
        File chr_order_txt
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} cut -f 1 ${chr_order_txt} > intermediate/${sample}_chr_keep.txt
    >>>

    output {
        File chr_keep_txt="${PWD}/intermediate/${sample}_chr_keep.txt"
    }
}

task TssEnrichment {
    input {
        String sample
        File noMT_BAM_index
        File noMT_BAM
        File plus_TSS
        File minus_TSS
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} python scripts/pyTssEnrichment.py -a ${noMT_BAM} -b ${plus_TSS} -p ends -c 16 -z -v -s 6 -o stat/${sample}_minus_TssEnrichment.txt
        singularity exec --bind ${PWD}:/root ${SIF} python scripts/pyTssEnrichment.py -a ${noMT_BAM} -b ${minus_TSS} -p ends -c 16 -z -v -s 6 -o stat/${sample}_plus_TssEnrichment.txt 
    >>>

    output {
        File plus_TssEnrichment="${PWD}/stat/${sample}_plus_TssEnrichment.txt"
        File minus_TssEnrichment="${PWD}/stat/${sample}_minus_TssEnrichment.txt"
    }
}

task plot_tss_enrichment {
    input {
        String sample
        File plus_TssEnrichment
        File minus_TssEnrichment
        File SIF
        String PWD
    }

    command <<<
        if [ ! -d "${PWD}/intermediate/stat/" ]; then
        mkdir ${PWD}/intermediate/stat/
        fi
        singularity exec --bind ${PWD}:/root ${SIF} Rscript scripts/plot_tss_enrichment.R stat/${sample}_plus_TssEnrichment.txt stat/${sample}_minus_TssEnrichment.txt
    >>>

    output {
        File pdf="${PWD}/intermediate/stat/${sample}_TSSenrichment.pdf"
        File png="${PWD}/intermediate/stat/${sample}_TSSenrichment.png"
    }
}

task tidy_tss_enrichment {
    input {
        String sample
        File pdf
        File png
        File SIF
        String PWD
    }

    command <<<
        if [ ! -d "${PWD}/figure" ]; then
        mkdir ${PWD}/figure
        fi
        singularity exec --bind ${PWD}:/root ${SIF} mv ${pdf} ${png} figure/
    >>>

    output {
        File TSSenrichment_pdf="${PWD}/figure/${sample}_TSSenrichment.pdf"
        File TSSenrichment_png="${PWD}/figure/${sample}_TSSenrichment.png"
    }
}

task Generate_BigWig {
    input {
        String sample
        File noMT_BAM_index
        File noMT_BAM
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} bamCoverage -b ${noMT_BAM} -o data/${sample}_plus.bw --filterRNAstrand forward --binSize 1 --ignoreDuplicates --minMappingQuality 30 --normalizeUsing CPM -p 20 
        singularity exec --bind ${PWD}:/root ${SIF} bamCoverage -b ${noMT_BAM} -o data/${sample}_minus.bw --filterRNAstrand reverse --binSize 1 --ignoreDuplicates --minMappingQuality 30 --normalizeUsing CPM -p 20
    >>>

    output {
        File plus_bw="${PWD}/data/${sample}_plus.bw"
        File minus_bw="${PWD}/data/${sample}_minus.bw"
    }
}

task converted_bed {
    input {
        String sample
        File noMT_BAM_index
        File noMT_BAM
        File SIF
        String PWD
    }

    command <<<
        if [ ! -d "${PWD}/bed" ]; then
        mkdir ${PWD}/bed
        fi
        singularity exec --bind ${PWD}:/root ${SIF} bedtools bamtobed -i ${noMT_BAM} > bed/${sample}.converted.bed
    >>>

    output {
        File convert_bed="${PWD}/bed/${sample}.converted.bed"
    }
}

task sorted_bed {
    input {
        String sample
        File convert_bed
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} bedtools sort -i ${convert_bed} > bed/${sample}.sorted.bed
    >>>

    output {
        File sorted_bed="${PWD}/bed/${sample}.sorted.bed"
    }
}

task HOMER_makeTagDirectory {
    input {
        String sample
        File sorted_bed
        File SIF
        String PWD
    }

    command <<<
        if [ ! -d "${PWD}/intermediate/${sample}/HOMER_tag/" ]; then
        mkdir -p ${PWD}/intermediate/${sample}/HOMER_tag/
        fi
        singularity exec --bind ${PWD}:/root ${SIF} makeTagDirectory intermediate/${sample}/HOMER_tag/ -format bed -forceBED  ${sorted_bed}
    >>>

    output {
        File HOMER_tag="${PWD}/intermediate/{sample}/HOMER_tag/tagInfo.txt"
    }
}

task HOMER_findPeaks {
    input {
        String sample
        File HOMER_tag
        File SIF
        String PWD
    }

    command <<<
        if [ ! -d "${PWD}/intermediate/${sample}/HOMER_tag/" ]; then
        mkdir -p ${PWD}/intermediate/${sample}/HOMER_tag/
        fi
        singularity exec --bind ${PWD}:/root ${SIF} findPeaks intermediate/${sample}/HOMER_tag/ -style groseq -o intermediate/${sample}/transcript.txt -gtf intermediate/${sample}/transcript.gtf
    >>>

    output {
        File txt="${PWD}/intermediate/${sample}/transcript.txt"
        File gtf="${PWD}/intermediate/${sample}/transcript.gtf"
    }
}

task pause_index {
    input {
        String sample
        File sorted_bed
        File convert_bed
        File fn_gtf
        File fa_in
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} python scripts/pause_index.py . ${sample} ${fn_gtf} ${fa_in}
    >>>

    output {
        File raw_txt="${PWD}/data/${sample}_raw.txt"
    }
}

task filtered_raw_txt {
    input {
        String sample
        File raw_txt
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} awk '$10>0 && $10!="NA"' ${raw_txt} > intermediate/${sample}_filtered_raw.txt
    >>>

    output {
        File filtered_raw_txt="${PWD}/intermediate/${sample}_filtered_raw.txt"
    }
}

task plot_pause_index {
    input {
        String sample
        File filtered_raw_txt
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} Rscript scripts/plot_pause_index.r intermediate/${sample}_filtered_raw.txt
    >>>

    output {
        File pdf="${PWD}/intermediate/${sample}_pause_index.pdf"
        File png="${PWD}/intermediate/${sample}_pause_index.png"
    }
}

task tidy_pause_index {
    input {
        String sample
        File pdf
        File png
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} mv ${pdf} ${png} figure/
    >>>

    output {
        File pause_index_pdf="${PWD}/figure/${sample}_pause_index.pdf"
        File pause_index_png="${PWD}/figure/${sample}_pause_index.png"
    }
}

task exons_introns_sort {
    input {
        String sample
        File chr_keep_txt
        File chr_order_txt
        File exon_name
        File intron_name
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} grep -wf ${chr_keep_txt} ${exon_name} | \
        singularity exec --bind ${PWD}:/root ${SIF} bedtools sort -i stdin -faidx ${chr_order_txt} > bed/${sample}_exons_sort.bed
        singularity exec --bind ${PWD}:/root ${SIF} grep -wf ${chr_keep_txt} ${intron_name} | \
        singularity exec --bind ${PWD}:/root ${SIF} bedtools sort -i stdin -faidx ${chr_order_txt} | \
        singularity exec --bind ${PWD}:/root ${SIF} bedtools sort -i stdin -faidx ${chr_order_txt} > bed/${sample}_introns_sort.bed

    >>>

    output {
        File exons_sort="${PWD}/bed/${sample}_exons_sort.bed"
        File introns_sort="${PWD}/bed/${sample}_introns_sort.bed"
    }
}

task exons_introns_coverage {
    input {
        String sample
        File exons_sort
        File introns_sort
        File noMT_BAM_index
        File chr_order_txt
        File noMT_BAM
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} bedtools coverage -sorted -counts -s -a bed/${sample}_exons_sort.bed -b ${noMT_BAM} -g ${chr_order_txt} > bed/${sample}_exons_coverage.bed
        singularity exec --bind ${PWD}:/root ${SIF} bedtools coverage -sorted -counts -s -a bed/${sample}_introns_sort.bed -b ${noMT_BAM} -g ${chr_order_txt} > bed/${sample}_introns_coverage.bed
    >>>

    output {
        File exons_coverage="${PWD}/bed/${sample}_exons_coverage.bed"
        File introns_coverage="${PWD}/bed/${sample}_introns_coverage.bed"
    }
}

task exons_introns_rpkm {
    input {
        String sample
        File exons_coverage
        File introns_coverage
        File noMT_BAM_index
        File noMT_BAM
        File SIF
        String PWD
    }

    command <<<
        ar=$(singularity exec --bind ${PWD}:/root ${SIF} samtools view -F 4 -c ${noMT_BAM})
        scaling_factor=$(echo "scale=6; $ar/1000000" | singularity exec --bind ${PWD}:/root ${SIF} bc)
        singularity exec --bind ${PWD}:/root ${SIF} awk -v OFS='\t' -v scaling_factor=$scaling_factor '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/scaling_factor)/geneSizeKB[a], strand[a] }}' bed/${sample}_exons_coverage.bed | awk '$5>0' | sort -k4 > bed/${sample}_exons_rpkm.bed
        singularity exec --bind ${PWD}:/root ${SIF} awk -v OFS='\t' -v scaling_factor=$scaling_factor '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/scaling_factor)/geneSizeKB[a], strand[a] }}' bed/${sample}_introns_coverage.bed | awk '$5>0' | sort -k4 > bed/${sample}_introns_rpkm.bed
    >>>

    output {
        File exons_rpkm="${PWD}/bed/${sample}_exons_rpkm.bed"
        File introns_rpkm="${PWD}/bed/${sample}_introns_rpkm.bed"
    }
}

task all_genes {
    input {
        String sample
        File exons_rpkm
        File introns_rpkm
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} cat bed/${sample}_exons_rpkm.bed bed/${sample}_introns_rpkm.bed | awk '{print $4}' > intermediate/${sample}_all_genes.txt
    >>>

    output {
        File all_genes_txt="${PWD}/intermediate/${sample}_all_genes.txt"
    }
}

task shared_genes {
    input {
        String sample
        File all_genes_txt
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} awk -F';' '$1 in first{print first[$1] $0; first[$1]=""; next} {first[$1]=$0 ORS}' ${all_genes_txt} | awk '!NF {print;next}; !($0 in a) {a[$0];print}' > bed/${sample}_shared_genes.txt
    >>>

    output {
        File shared_genes_txt="${PWD}/bed/${sample}_shared_genes.txt"
    }
}

task shared_exons_introns {
    input {
        String sample
        File shared_genes_txt
        File exons_rpkm_bed
        File introns_rpkm_bed
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} awk -F '\t' 'NR==FNR {id[$1]; next} $4 in id' ${shared_genes_txt} ${exons_rpkm_bed} > bed/${sample}_shared_exons.bed
        singularity exec --bind ${PWD}:/root ${SIF} awk -F '\t' 'NR==FNR {id[$1]; next} $4 in id' ${shared_genes_txt} ${introns_rpkm_bed} > bed/${sample}_shared_introns.bed
    >>>

    output {
        File shared_exons="${PWD}/bed/${sample}_shared_exons.bed"
        File shared_introns="${PWD}/bed/${sample}_shared_introns.bed"
    }
}

task exon_intron_ratios {
    input {
        String sample
        File shared_exons
        File shared_introns
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} awk 'BEGIN{FS=OFS="\t"} FNR>0 && FNR==NR{a[$4]=$4 OFS $0; next} FNR>0{print $0,a[$4]?a[$4]:"\t"}' bed/${sample}_shared_exons.bed bed/${sample}_shared_introns.bed | awk 'BEGIN { OFS="\t" } {print $8, $9, $10, $11, ($12/$5), $13}' | sort -k1,1 -k2,2n > bed/${sample}_exon_intron_ratios.bed
    >>>

    output {
        File exon_intron_ratios="${PWD}/bed/${sample}_exon_intron_ratios.bed"
    }
}

task plot_mRNA_contamination {
    input {
        String sample
        File exon_intron_ratios
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} Rscript scripts/PEPPRO.R mrna -i bed/${sample}_exon_intron_ratios.bed --annotate
    >>>

    output {
        File pdf="${PWD}/bed/${sample}_mRNA_contamination.pdf"
        File png="${PWD}/bed/${sample}_mRNA_contamination.png"
    }
}

task tidy_mRNA_contamination {
    input {
        String sample
        File pdf
        File png
        File SIF
        String PWD
    }

    command <<<
        singularity exec --bind ${PWD}:/root ${SIF} mv ${pdf} ${png} figure/
    >>>

    output {
        File mRNA_contamination_pdf="${PWD}/figure/${sample}_mRNA_contamination.pdf"
        File mRNA_contamination_png="${PWD}/figure/${sample}_mRNA_contamination.png"
    }
}

task Final_MultiQC {
    input {
        String sample
        File noMT_BAM_index
        File SIF
        String PWD
    }

    command <<<
        if [ ! -d "${PWD}/data/multiqc_report/" ]; then
        mkdir -p ${PWD}/data/multiqc_report/
        fi
        singularity exec --bind ${PWD}:/root ${SIF} multiqc data/ -o data/multiqc_report/${sample} --force
    >>>

    output {
        File multiqc_report_html="${PWD}/data/multiqc_report/${sample}/multiqc_report.html"
    }
}