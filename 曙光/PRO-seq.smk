import time
import os

import time

configfile : "config.yaml"
rule all:
    input:
        expand("data/{sample}/{sample}_R1_unmap_fastqc.html", sample=config["samples"]),
        expand("figure/{sample}_mRNA_contamination.{ext}", sample=config["samples"], ext=["pdf","png"]),
        expand("figure/{sample}_TSSenrichment.{ext}", sample=config["samples"], ext=["pdf","png"]),
        expand("data/{sample}_{strand}.bw", sample=config["samples"], strand=["plus","minus"]),
        expand("intermediate/{sample}/transcript.{ext}", sample=config["samples"], ext=["txt","gtf"]),
        expand("figure/{sample}_pause_index.{ext}", sample=config["samples"], ext=["pdf","png"]),
        expand("data/multiqc_report/{sample}/multiqc_report.html", sample=config["samples"]),

rule qc_before_filtered:
    output:
        html=temp("qc/{sample}/{sample}_R1_fastqc.html"),
    params:
        FQ1=lambda wildcards: config["samples"][wildcards.sample],
        threads=10,
    log:
        "logs/{sample}_qc_before_filtered.log"
    container:
        config["container"]
    shell:
        """
        set +e
        fastqc -t {params.threads} {params.FQ1} -o qc/{wildcards.sample} 2>{log}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tqc_before_filtered\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tqc_before_filtered\tdone\t$(date)" >> logs/logs

        """


rule trim_fastq:
    input:
        html="qc/{sample}/{sample}_R1_fastqc.html",
    output:
        trimmed_FQ1=temp("qc/{sample}/{sample}_R1_trimmed.fq.gz"),
    params:
        FQ1="data/samples/{sample}_R1.fastq.gz",
    threads: 
        10
    log:
        "logs/{sample}_trim_fastq.log"
    container:
        config["container"]
    shell:
        """
        set +e
        trim_galore --fastqc --cores {threads} -o qc/{wildcards.sample} {params.FQ1} 2>{log}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\ttrim_fastq\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\ttrim_fastq\tdone\t$(date)" >> logs/logs
        """
        

rule umi_extracted:
    input:
        trimmed_FQ1="qc/{sample}/{sample}_R1_trimmed.fq.gz",
    output:
        umi_extracted_FQ1=temp("qc/{sample}/{sample}_R1_umi_extracted.fastq"),
    log:
        "logs/{sample}_umi_extracted.log"
    container:
        config["container"]
    shell:
        """
        set +e
        umi_tools extract \
        --extract-method=string \
        --bc-pattern=NNNNNNNN \
        --stdin={input.trimmed_FQ1} \
        --stdout=qc/{wildcards.sample}/{wildcards.sample}_R1_umi_extracted.fastq \
        -L {log}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tumi_extracted\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tumi_extracted\tdone\t$(date)" >> logs/logs
        """


rule remove_additional_nucleotide:
    input:
        umi_extracted_FQ1="qc/{sample}/{sample}_R1_umi_extracted.fastq",
    output:
        remove_additional_nuc_FQ1=temp("qc/{sample}/{sample}_R1_result_trimmed.fastq"),
    log:
        "logs/{sample}_remove_additional_nucleotide.log"
    container:
        config["container"]
    shell:
        """
        set +e
        seqtk trimfq -b 1 {input.umi_extracted_FQ1} > {output.remove_additional_nuc_FQ1}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tremove_additional_nucleotide\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tremove_additional_nucleotide\tdone\t$(date)" >> logs/logs
        """


rule reverse_complement:
    input:
        remove_additional_nuc_FQ1="qc/{sample}/{sample}_R1_result_trimmed.fastq",
    output:
        rev_compl_FQ1=temp("qc/{sample}/{sample}_R1_processed.fastq"),
    log:
        "logs/{sample}_reverse_complement.log"
    container:
        config["container"]
    shell:
        """
        set +e
        seqtk seq -L 10 -r {input.remove_additional_nuc_FQ1} > {output.rev_compl_FQ1}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\treverse_complement\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\treverse_complement\tdone\t$(date)" >> logs/logs
        """

rule prealignments:
    input:
        rev_compl_FQ1="qc/{sample}/{sample}_R1_processed.fastq",
    output:
        NON_RRNA_FQ1=temp("qc/{sample}/{sample}_R1_unmap.fastq"),
    params:
        BOWTIE2_RDNAIDX = config["BOWTIE2_RDNAIDX"]
    log:
        "logs/{sample}_prealignments.log"
    container:
        config["container"]
    shell:
        """
        set +e
        bowtie2 -p 10 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x {params.BOWTIE2_RDNAIDX} --rg-id {wildcards.sample} -U {input.rev_compl_FQ1} --un qc/{wildcards.sample}/{wildcards.sample}_R1_unmap.fastq | \
        samtools view -bS - -@ 1 | \
        samtools sort - -@ 1 -o data/prealignments_{wildcards.sample}.bam
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tprealignments\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tprealignments\tdone\t$(date)" >> logs/logs
        """


rule qc_after_filtered:
    input:
        NON_RRNA_FQ1="qc/{sample}/{sample}_R1_unmap.fastq",
    output:
        html="data/{sample}/{sample}_R1_unmap_fastqc.html",
    threads: 
        10
    log:
        "logs/{sample}_qc_after_filtered.log"
    container:
        config["container"]
    shell:
        """
        set +e
        fastqc -t {threads} {input.NON_RRNA_FQ1} -o data/{wildcards.sample} 2>{log}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tqc_after_filtered\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tqc_after_filtered\tdone\t$(date)" >> logs/logs
        """


rule primary_single_alignments:
    input:
        NON_RRNA_FQ1="qc/{sample}/{sample}_R1_unmap.fastq",
    output:
        temp_sam=temp("data/{sample}_temp.sam"),
    threads: 
        10
    log:
        "logs/{sample}_primary_single_alignments.log"
    params:
        BOWTIE2IDX = config["BOWTIE2IDX"],
        home=config["home"],
    container:
        config["container"]
    shell:
        """ 
        set +e
        mkdir -p stat/
        bowtie2 -p {threads} --very-sensitive -q --end-to-end --un-conc qc/{wildcards.sample}/{wildcards.sample}_hgRemovedTemp.fastq --rg-id {wildcards.sample} -x {params.BOWTIE2IDX} -U {input.NON_RRNA_FQ1} -S {output.temp_sam} 2>&1 | \
        tee stat/{wildcards.sample}_align_with_bt2.txt
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tprimary_single_alignments\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tprimary_single_alignments\tdone\t$(date)" >> logs/logs
        """


rule umi_dedup:
    input:
        temp_sam="data/{sample}_temp.sam",
    output:
        pair_deduped_bam=temp("data/{sample}_pair_deduped.bam"),
    log:
        "logs/{sample}_pair_umiTools_Log.log"
    container:
        config["container"]
    shell:
        """
        set +e
        umi_tools dedup --stdin {input.temp_sam} --stdout data/{wildcards.sample}_pair_deduped.bam --log {log} --umi-separator="_"
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tumi_dedup\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tumi_dedup\tdone\t$(date)" >> logs/logs
        """


rule bam2sam:
    input:
        pair_deduped_bam="data/{sample}_pair_deduped.bam",
    output:
        pair_deduped_sam=temp("data/{sample}_pair_deduped.sam"),
    
    log:
        "logs/{sample}_bam2sam.log"
    container:
        config["container"]
    shell:
        """
        set +e
        samtools view -h -o {output.pair_deduped_sam} {input.pair_deduped_bam}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tbam2sam\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tbam2sam\tdone\t$(date)" >> logs/logs
        """


rule Collapsing:
    input:
        pair_deduped_sam="data/{sample}_pair_deduped.sam",
    output:
        temp_BAM="data/{sample}_temp.bam",
    log:
        "logs/{sample}_Collapsing.log"
    container:
        config["container"]
    shell:
        """
        set +e
        awk \'$2 != 4 {{print}}\' {input.pair_deduped_sam} | sed \'/XS:/d\' | samtools view -bS - -@ 1 | samtools sort - -@ 1 -o {output.temp_BAM}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tCollapsing\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tCollapsing\tdone\t$(date)" >> logs/logs
        """


rule index_temp_BAM:
    input:
        temp_BAM="data/{sample}_temp.bam",
    output:
        temp_BAM_index=temp("data/{sample}_temp.bam.bai"),
    log:
        "logs/{sample}_index_temp_BAM.log"
    container:
        config["container"]
    shell:
        """
        set +e
        samtools index {input.temp_BAM}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tindex_temp_BAM\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tindex_temp_BAM\tdone\t$(date)" >> logs/logs
        """


rule sorted_BAM:
    input:
        temp_BAM_index="data/{sample}_temp.bam.bai",
    output:
        sort_BAM="data/{sample}_sort.bam",
    params:
        temp_BAM="data/{sample}_temp.bam",
        threads=10,
    log:
        "logs/{sample}_sorted_BAM.log"
    container:
        config["container"]
    shell:
        """
        set +e
        samtools view -q 10 -b -@ {params.threads} -U data/{wildcards.sample}_fail_qc_dups.bam {params.temp_BAM} > {output.sort_BAM}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tsorted_BAM\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tsorted_BAM\tdone\t$(date)" >> logs/logs
        """


rule index_sorted_BAM:
    input:
        sort_BAM="data/{sample}_sort.bam",
    output:
        sort_BAM_index=temp("data/{sample}_sort.bam.bai"),
    container:
        config["container"]
    shell:
        """
        set +e
        samtools index {input.sort_BAM}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tindex_sorted_BAM\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tindex_sorted_BAM\tdone\t$(date)" >> logs/logs
        """

rule MT_Removed:
    input:
        sort_BAM_index="data/{sample}_sort.bam.bai",
    output:
        chr_sizes_bed=temp("bed/{sample}_chr_sizes.bed"),
    params:
        sort_BAM="data/{sample}_sort.bam",
    container:
        config["container"]
    shell:
        """
        set +e
        samtools idxstats {params.sort_BAM} | cut -f 1-2 | awk \'{{print $1, 0, $2}}\' | grep -vwe \'chrM\' -vwe \'chrMT\' -vwe \'M\' -vwe \'MT\' -vwe \'rCRSd\' -vwe \'rCRSd_3k\' > {output.chr_sizes_bed}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tMT_Removed\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tMT_Removed\tdone\t$(date)" >> logs/logs
        """

rule GET_noMT_BAM:
    input:
        chr_sizes_bed="bed/{sample}_chr_sizes.bed",
    output:
        noMT_BAM="data/{sample}_noMT.bam",
    params:
        sort_BAM="data/{sample}_sort.bam",
    threads: 
        10
    container:
        config["container"]
    shell:
        """
        set +e
        samtools view -L {input.chr_sizes_bed} -b -@ {threads} {params.sort_BAM} > {output.noMT_BAM}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tGET_noMT_BAM\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tGET_noMT_BAM\tdone\t$(date)" >> logs/logs
        """


rule index_noMT_BAM:
    input:
        noMT_BAM="data/{sample}_noMT.bam",
    output:
        noMT_BAM_index="data/{sample}_noMT.bam.bai",
    container:
        config["container"]
    shell:
        """
        set +e
        samtools index {input.noMT_BAM}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tindex_noMT_BAM\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tindex_noMT_BAM\tdone\t$(date)" >> logs/logs
        """


rule chr_order:
    input:
        noMT_BAM_index="data/{sample}_noMT.bam.bai",
    output:
        chr_order_txt="intermediate/{sample}_chr_order.txt",
    params:
        noMT_BAM="data/{sample}_noMT.bam",
    container:
        config["container"]
    shell:
        """
        set +e
        samtools view -H {params.noMT_BAM} | grep \'SN:\' | awk -F\':\' \'{{print $2,$3}}\' | awk -F\' \' -v OFS=\'\\t\' \'{{print $1,$3}}\' > {output.chr_order_txt}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tchr_order\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tchr_order\tdone\t$(date)" >> logs/logs
        """


rule chr_keep:
    input:
        chr_order_txt="intermediate/{sample}_chr_order.txt",
    output:
        chr_keep_txt=temp("intermediate/{sample}_chr_keep.txt"),
    container:
        config["container"]
    shell:
        """
        set +e
        cut -f 1 {input.chr_order_txt} > {output.chr_keep_txt}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tchr_keep\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tchr_keep\tdone\t$(date)" >> logs/logs
        """


rule TssEnrichment:
    input:
        noMT_BAM_index="data/{sample}_noMT.bam.bai",
    output:
        plus_TssEnrichment=temp("stat/{sample}_plus_TssEnrichment.txt"),
        minus_TssEnrichment=temp("stat/{sample}_minus_TssEnrichment.txt"),
    params:
        noMT_BAM="data/{sample}_noMT.bam",
        plus_TSS=config["plus_TSS"],
        minus_TSS=config["minus_TSS"],
        threads=10,
    log:
        "logs/{sample}_TssEnrichment.log"
    container:
        config["container"]
    shell:
        """
        set +e
        python /root/scripts/pyTssEnrichment.py -a {params.noMT_BAM} -b {params.plus_TSS} -p ends -c {params.threads} -z -v -s 6 -o stat/{wildcards.sample}_minus_TssEnrichment.txt
        python /root/scripts/pyTssEnrichment.py -a {params.noMT_BAM} -b {params.minus_TSS} -p ends -c {params.threads} -z -v -s 6 -o stat/{wildcards.sample}_plus_TssEnrichment.txt 
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tTssEnrichment\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tTssEnrichment\tdone\t$(date)" >> logs/logs
        """


rule plot_tss_enrichment:
    input:
        plus_TssEnrichment="stat/{sample}_plus_TssEnrichment.txt",
        minus_TssEnrichment="stat/{sample}_minus_TssEnrichment.txt",
    output:
        temp(multiext("intermediate/stat/{sample}_TSSenrichment", out1=".pdf", out2=".png")),
    log:
        "logs/{sample}_plot_tss_enrichment.log"
    container:
        config["container"]
    shell:
        """
        set +e
        Rscript /root/scripts/plot_tss_enrichment.R stat/{wildcards.sample}_plus_TssEnrichment.txt stat/{wildcards.sample}_minus_TssEnrichment.txt
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tplot_tss_enrichment\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tplot_tss_enrichment\tdone\t$(date)" >> logs/logs
        """


rule tidy_tss_enrichment:
    input:
        multiext("intermediate/stat/{sample}_TSSenrichment", in1=".pdf", in2=".png"),
    output:
        multiext("figure/{sample}_TSSenrichment", out1=".pdf", out2=".png"),
    container:
        config["container"]
    shell:
        """
        set +e
        mv {input.in1} {input.in2} figure/
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\ttidy_tss_enrichment\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\ttidy_tss_enrichment\tdone\t$(date)" >> logs/logs
        """



rule Generate_BigWig:
    input:
        noMT_BAM_index="data/{sample}_noMT.bam.bai",
    output:
        plus_bw="data/{sample}_plus.bw",
        minus_bw="data/{sample}_minus.bw",
    params:
        noMT_BAM="data/{sample}_noMT.bam",
        threads=16,
    log:
        "logs/{sample}_Generate_BigWig.log"
    container:
        config["container"]
    shell:
        """
        set +e
        bamCoverage -b {params.noMT_BAM} -o data/{wildcards.sample}_plus.bw --filterRNAstrand forward --binSize 1 --ignoreDuplicates --minMappingQuality 30 --normalizeUsing CPM -p {params.threads} 
        bamCoverage -b {params.noMT_BAM} -o data/{wildcards.sample}_minus.bw --filterRNAstrand reverse --binSize 1 --ignoreDuplicates --minMappingQuality 30 --normalizeUsing CPM -p {params.threads} 
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tGenerate_BigWig\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tGenerate_BigWig\tdone\t$(date)" >> logs/logs
        """  


rule converted_bed:
    input:
        noMT_BAM_index="data/{sample}_noMT.bam.bai",
    output:
        convert_bed="bed/{sample}.converted.bed",
    params:
        noMT_BAM="data/{sample}_noMT.bam",
    container:
        config["container"]
    shell:
        """
        set +e
        bedtools bamtobed -i {params.noMT_BAM} > {output.convert_bed}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tconverted_bed\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tconverted_bed\tdone\t$(date)" >> logs/logs
        """


rule sorted_bed:
    input:
        convert_bed="bed/{sample}.converted.bed",
    output:
        sorted_bed=temp("bed/{sample}.sorted.bed"),
    container:
        config["container"]
    shell:
        """
        set +e
        bedtools sort -i {input.convert_bed} > {output.sorted_bed}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tsorted_bed\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tsorted_bed\tdone\t$(date)" >> logs/logs
        """


rule HOMER_makeTagDirectory:
    input:
        sorted_bed="bed/{sample}.sorted.bed",
    output:
        HOMER_tag=temp("intermediate/{sample}/HOMER_tag/tagInfo.txt"),
    container:
        config["container"]
    shell:
        """
        set +e
        makeTagDirectory intermediate/{wildcards.sample}/HOMER_tag/ -format bed -forceBED  {input.sorted_bed}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tHOMER_makeTagDirectory\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tHOMER_makeTagDirectory\tdone\t$(date)" >> logs/logs
        """


rule HOMER_findPeaks:
    input:
        HOMER_tag="intermediate/{sample}/HOMER_tag/tagInfo.txt",
    output:
        multiext("intermediate/{sample}/transcript", out1=".txt", out2=".gtf"),
    container:
        config["container"]
    shell:
        """
        set +e
        findPeaks intermediate/{wildcards.sample}/HOMER_tag/ -style groseq -o {output.out1} -gtf intermediate/{wildcards.sample}/transcript.gtf
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tHOMER_findPeaks\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tHOMER_findPeaks\tdone\t$(date)" >> logs/logs
        """


rule pause_index:
    input:
        sorted_bed="bed/{sample}.sorted.bed",
    output:
        raw_txt=temp("data/{sample}_raw.txt"),
    params:
        convert_bed="bed/{sample}.converted.bed",
        fn_gtf= config["fn_gtf"],
        fa_in = config["fa_in"],
    log:
        "logs/{sample}_pause_index.log"
    container:
        config["container"]
    shell:
        """
        set +e
        python /root/scripts/pause_index.py . {wildcards.sample} {params.fn_gtf} {params.fa_in}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tpause_index\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tpause_index\tdone\t$(date)" >> logs/logs
        """


rule filtered_raw_txt:
    input:
        raw_txt="data/{sample}_raw.txt",
    output:
        filtered_raw_txt=temp("intermediate/{sample}_filtered_raw.txt"),

    container:
        config["container"]
    shell:
        """
        set +e
        awk \'$10>0 && $10!="NA"\' {input.raw_txt} > {output.filtered_raw_txt}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tfiltered_raw_txt\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tfiltered_raw_txt\tdone\t$(date)" >> logs/logs
        """
 

rule plot_pause_index:
    input:
        filtered_raw_txt="intermediate/{sample}_filtered_raw.txt",
    output:
        temp(multiext("intermediate/{sample}_pause_index", out1=".pdf", out2=".png")),
    log:
        "logs/{sample}_plot_pause_index.log"
    container:
        config["container"]
    shell:
        """
        set +e
        Rscript /root/scripts/plot_pause_index.r intermediate/{wildcards.sample}_filtered_raw.txt
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tplot_pause_index\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tplot_pause_index\tdone\t$(date)" >> logs/logs
        """


rule tidy_pause_index:
    input:
        multiext("intermediate/{sample}_pause_index", in1=".pdf", in2=".png"),
    output:
        multiext("figure/{sample}_pause_index", out1=".pdf", out2=".png"),
    container:
        config["container"]
    shell:
        """
        set +e
        mv {input.in1} {input.in2} figure/
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\ttidy_pause_index\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\ttidy_pause_index\tdone\t$(date)" >> logs/logs
        """


rule exons_introns_sort:
    input:
        chr_keep_txt="intermediate/{sample}_chr_keep.txt",
    output:
        exons_sort=temp("bed/{sample}_exons_sort.bed"),
        introns_sort=temp("bed/{sample}_introns_sort.bed"),
    params:
        chr_order_txt="intermediate/{sample}_chr_order.txt",
        exon_name=config["exon_name"],
        intron_name=config["intron_name"],
    container:
        config["container"]
    shell:
        """
        set +e
        grep -wf {input.chr_keep_txt} {params.exon_name} | bedtools sort -i stdin -faidx {params.chr_order_txt} > bed/{wildcards.sample}_exons_sort.bed
        grep -wf {input.chr_keep_txt} {params.intron_name} | bedtools sort -i stdin -faidx {params.chr_order_txt} | bedtools sort -i stdin -faidx {params.chr_order_txt} > bed/{wildcards.sample}_introns_sort.bed
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\texons_introns_sort\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\texons_introns_sort\tdone\t$(date)" >> logs/logs
        """


rule exons_introns_coverage:
    input:
        exons_sort="bed/{sample}_exons_sort.bed",
        introns_sort="bed/{sample}_introns_sort.bed",
    output:
        exons_coverage=temp("bed/{sample}_exons_coverage.bed"),
        introns_coverage=temp("bed/{sample}_introns_coverage.bed"),
    params:
        noMT_BAM_index="data/{sample}_noMT.bam.bai",
        chr_order_txt="intermediate/{sample}_chr_order.txt",
        noMT_BAM="data/{sample}_noMT.bam",
    log:
        "logs/{sample}_exons_introns_coverage.log"
    container:
        config["container"]
    shell:
        """
        set +e
        bedtools coverage -sorted -counts -s -a bed/{wildcards.sample}_exons_sort.bed -b {params.noMT_BAM} -g {params.chr_order_txt} > bed/{wildcards.sample}_exons_coverage.bed
        bedtools coverage -sorted -counts -s -a bed/{wildcards.sample}_introns_sort.bed -b {params.noMT_BAM} -g {params.chr_order_txt} > bed/{wildcards.sample}_introns_coverage.bed
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\texons_introns_coverage\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\texons_introns_coverage\tdone\t$(date)" >> logs/logs
        """


rule exons_introns_rpkm:
    input:
        exons_coverage="bed/{sample}_exons_coverage.bed",
        introns_coverage="bed/{sample}_introns_coverage.bed",
    output:
        exons_rpkm="bed/{sample}_exons_rpkm.bed",
        introns_rpkm="bed/{sample}_introns_rpkm.bed",
    params:
        noMT_BAM="data/{sample}_noMT.bam",
        noMT_BAM_index="data/{sample}_noMT.bam.bai",
    log:
        "logs/{sample}_exons_introns_rpkm.log"
    container:
        config["container"]
    shell:
        """
        set +e
        ar=$(samtools view -F 4 -c {params.noMT_BAM})
        scaling_factor=$(echo "scale=6; $ar/1000000" | bc)
        awk -v OFS=\'\\t\' -v scaling_factor=$scaling_factor \'{{chrom[$4] = $1; if($4!=prev4) {{chromStart[$4] = $2}} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4}} END {{ for (a in readCount) {{ print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/scaling_factor)/geneSizeKB[a], strand[a] }}}}\' bed/{wildcards.sample}_exons_coverage.bed | awk \'$5>0\' | sort -k4 > bed/{wildcards.sample}_exons_rpkm.bed
        awk -v OFS=\'\\t\' -v scaling_factor=$scaling_factor \'{{chrom[$4] = $1; if($4!=prev4) {{chromStart[$4] = $2}} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4}} END {{ for (a in readCount) {{ print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/scaling_factor)/geneSizeKB[a], strand[a] }}}}\' bed/{wildcards.sample}_introns_coverage.bed | awk \'$5>0\' | sort -k4 > bed/{wildcards.sample}_introns_rpkm.bed
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\texons_introns_rpkm\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\texons_introns_rpkm\tdone\t$(date)" >> logs/logs
        """
      

        
rule all_genes:
    input:
        exons_rpkm="bed/{sample}_exons_rpkm.bed",
        introns_rpkm="bed/{sample}_introns_rpkm.bed",
    output:
        all_genes_txt=temp("intermediate/{sample}_all_genes.txt"),
    container:
        config["container"]
    shell:
        """
        set +e
        cat bed/{wildcards.sample}_exons_rpkm.bed bed/{wildcards.sample}_introns_rpkm.bed | awk \'{{print $4}}\' > {output.all_genes_txt}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tall_genes\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tall_genes\tdone\t$(date)" >> logs/logs
        """


rule shared_genes:
    input:
        all_genes_txt="intermediate/{sample}_all_genes.txt",
    output:
        shared_genes_txt=temp("bed/{sample}_shared_genes.txt"),
    container:
        config["container"]
    shell:
        """
        set +e
        awk -F\';\' \'$1 in first{{print first[$1] $0; first[$1]=""; next}} {{first[$1]=$0 ORS}}\' {input.all_genes_txt} | awk \'!NF {{print;next}}; !($0 in a) {{a[$0];print}}\' > {output.shared_genes_txt}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tshared_genes\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tshared_genes\tdone\t$(date)" >> logs/logs
        """


rule shared_exons_introns:
    input:
        shared_genes_txt="bed/{sample}_shared_genes.txt",
    output:
        shared_exons=temp("bed/{sample}_shared_exons.bed"),
        shared_introns=temp("bed/{sample}_shared_introns.bed"),
    params:
        exons_rpkm_bed="bed/{sample}_exons_rpkm.bed",
        introns_rpkm_bed="bed/{sample}_introns_rpkm.bed",
    container:
        config["container"]
    shell:
        """
        set +e
        awk -F \'\\t\' \'NR==FNR {{id[$1]; next}} $4 in id\' {input.shared_genes_txt} {params.exons_rpkm_bed} > bed/{wildcards.sample}_shared_exons.bed
        awk -F \'\\t\' \'NR==FNR {{id[$1]; next}} $4 in id\' {input.shared_genes_txt} {params.introns_rpkm_bed} > bed/{wildcards.sample}_shared_introns.bed
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tshared_exons_introns\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tshared_exons_introns\tdone\t$(date)" >> logs/logs
        """


rule exon_intron_ratios:
    input:
        shared_exons="bed/{sample}_shared_exons.bed",
        shared_introns="bed/{sample}_shared_introns.bed",
    output:
        exon_intron_ratios=temp("bed/{sample}_exon_intron_ratios.bed"),
    container:
        config["container"]
    shell:
        """
        set +e
        awk \'BEGIN{{FS=OFS="\\t"}} FNR>0 && FNR==NR{{a[$4]=$4 OFS $0; next}} FNR>0{{print $0,a[$4]?a[$4]:"\\t"}}\' bed/{wildcards.sample}_shared_exons.bed bed/{wildcards.sample}_shared_introns.bed | awk \'BEGIN {{ OFS="\\t" }} {{print $8, $9, $10, $11, ($12/$5), $13}}\' | sort -k1,1 -k2,2n > {output.exon_intron_ratios}
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\texon_intron_ratios\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\texon_intron_ratios\tdone\t$(date)" >> logs/logs
        """


rule plot_mRNA_contamination:
    input:
        exon_intron_ratios="bed/{sample}_exon_intron_ratios.bed",
    output:
        temp(multiext("bed/{sample}_mRNA_contamination", out1=".pdf", out2=".png")),
    log:
        "logs/{sample}_plot_mRNA_contamination.log"
    container:
        config["container"]
    shell:
        """
        set +e
        Rscript /root/scripts/PEPPRO.R mrna -i bed/{wildcards.sample}_exon_intron_ratios.bed --annotate
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tplot_mRNA_contamination\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tplot_mRNA_contamination\tdone\t$(date)" >> logs/logs
        """


rule tidy_mRNA_contamination:
    input:
        multiext("bed/{sample}_mRNA_contamination", in1=".pdf", in2=".png"),
    output:
        multiext("figure/{sample}_mRNA_contamination", out1=".pdf", out2=".png"),
    container:
        config["container"]
    shell:
        """
        set +e
        mv {input.in1} {input.in2} figure/
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\ttidy_mRNA_contamination\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\ttidy_mRNA_contamination\tdone\t$(date)" >> logs/logs
        """


rule Final_MultiQC:
    input:
        noMT_BAM_index="data/{sample}_noMT.bam.bai",
    output:
        multiqc_report_html="data/multiqc_report/{sample}/multiqc_report.html",
    log:
        "logs/{sample}_Final_MultiQC.log"
    container:
        config["container"]
    shell:
        """
        set +e
        multiqc data/ -o data/multiqc_report/{wildcards.sample} --force
        if [ $? -ne 0 ]; then
            echo "{wildcards.sample}\tFinal_MultiQC\tfailed\t$(date)">> logs/logs
            exit 1
        fi
        echo "{wildcards.sample}\tFinal_MultiQC\tdone\t$(date)" >> logs/logs
        """

