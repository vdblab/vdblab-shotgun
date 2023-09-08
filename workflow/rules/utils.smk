rule split_fastq:
    """
    These dummy inputs are intended to be overwritten when importing the rule
    """
    input:
        R1="",
        R2="",
    output:
        R1=[],
        R2=[],
    params:
        outdir="",
        nshards=1,
    container:
        "docker://pegi3s/seqkit:2.3.0"
    threads: 4  # see their docs
    resources:
        mem_mb=4000,
    log:
        e="logs/split_fastq.e",
        o="logs/split_fastq.o",
    shell:
        """
        seqkit split2 \
            --threads {threads} \
            --read1 {input.R1} \
            --read2 {input.R2} \
            --by-part {params.nshards} \
            --force \
            --out-dir {params.outdir}/ \
            > {log.o} 2>> {log.e}
        """


rule merge_shards:
    input:
        R1=[],
        R2=[],
    output:
        R1="",
        R2="",
    threads: 1
    resources:
        mem_mb=1024,
    shell:
        """
        cat {input.R1} > {output.R1} 2> {log.e}
        cat {input.R2} > {output.R2} 2>> {log.e}
        """


rule concat_R1_R2:
    """ Merges multilane fastqs and/or standardizes names
    There is some overhead with this rule but the alternatives are worse:
    - symlinking or linking doesn't work across worker nodes, so local jobs
      succeed but distributed jobs can't see input files
    - conditionally combining fastqs if more than one provided is possible via
      input functions, but that can break wildcard globbing due to the lack of
      any standard sample name /file name relationship.
    - this is also a "convenient" place to deal with non-gzipped files. default is frmt is "gz"
    """
    input:
        R1=[],
        R2=[],
    output:
        R1="",
        R2="",
    conda:
        "../envs/base.yaml"
    log:
        e="logs/concat_r1_r2_{sample}.e",
    shell:
        """
        case {input.R1[0]} in
        *gz )
            cat {input.R1} > {output.R1} 2>> {log.e}
            cat {input.R2} > {output.R2} 2>> {log.e}
        ;;
        *bz2 )
            bzcat {input.R1} | gzip -c > {output.R1} 2>> {log.e}
            bzcat {input.R2} | gzip -c > {output.R2} 2>> {log.e}
        ;;
        *fastq | *fq  )
            cat {input.R1} | gzip -c > {output.R1} 2>> {log.e}
            cat {input.R2} | gzip -c > {output.R2} 2>> {log.e}
        ;;
        * )
            echo "Supported formats are gz, bz2, or uncompressed"
            exit 1
        esac

        """
