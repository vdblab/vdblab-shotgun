rule split_fastq:
    """
    These dummy inputs are intended to be overwritten when importing the rule
    """
    input:
        R1=[],
    output:
        reads=[],
    params:
        outdir=lambda wc, output: os.path.dirname(output.reads[0]),
        inputstring=lambda wc, input: (
            f"--read1 {input['R1']} --read2 {input['R2']}"
            if is_paired()
            else f"--read1 {input['R1']}"
        ),
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
            {params.inputstring} \
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


rule concat_lanes_fix_names:
    """ Merges multilane fastqs and/or standardizes names
    There is some overhead with this rule but the alternatives are worse:
    - symlinking or linking doesn't work across worker nodes, so local jobs
      succeed but distributed jobs can't see input files
    - conditionally combining fastqs if more than one provided is possible via
      input functions, but that can break wildcard globbing due to the lack of
      any standard sample name /file name relationship.
    - this is also a "convenient" place to deal with non-gzipped files. default is frmt is "gz"
    """

    #
    input:
        fq=[],
    output:
        fq="out_{sample}.1.fq.gz",
    log:
        e="logs/concat_lanes_fix_names_{sample}.e",
    shell:
        """
        case {input.fq[0]} in
        *gz )
            cat {input.fq} > {output.fq} 2>> {log.e}
        ;;
        *bz2 )
            bzcat {input.fq} | gzip -c > {output.fq} 2>> {log.e}
        ;;
        *fastq | *fq  )
            cat {input.fq} | gzip -c > {output.fq} 2>> {log.e}
        ;;
        * )
            echo "Supported formats are gz, bz2, or uncompressed"
            exit 1
        esac

        """
