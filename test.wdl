version 1.0

workflow main {
    input {
        String sample
        File FQ1
        File SIF
        String PWD
    }

    call qc_before_filtered {
        input:
            sample=sample,
            FQ1=FQ1,
            SIF=SIF,
            PWD=PWD
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
        if [ ! -d "qc/${sample}" ]; then
        mkdir -p qc/${sample}
        fi
        fastqc -t 8 ${FQ1} -o qc/${sample}
    >>>

    output {
        File html="${PWD}/qc/${sample}/${sample}_R1_fastqc.html"
    }

    runtime {
        docker: "crpi-6cdftxq0w4mkpmr1.cn-beijing.personal.cr.aliyuncs.com/wdl-docker/wdl-docker:latest"
    }
}