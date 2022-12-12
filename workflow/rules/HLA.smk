rule all:
    input:
        f"reports/optitype_result.tsv",


# rule sort_index:
#     input:
#         bam=config["bam"],
#     output:
#         results=temp(f"{os.path.basename(config['bam'])}.bai"),
#         bam=temp(f"{os.path.basename(config['bam'])}"),
#     container: config["docker_bowtie2"],
#     threads: 8
#     shell:
#         """
#         samtools view -bh -F 4 {input.bam} | samtools  sort -@ {threads} -o {output.bam} -O bam -
#         samtools index -@ {threads} {output.bam}
#         """

# rule xHLA:
#     """ Couldn't get this to work, got error where the tsv file was empty.  Perhaps too few reads?  Similar to https://github.com/humanlongevity/HLA/issues/62"""
#     input:
#         results=f"{os.path.basename(config['bam'])}.bai",
#         bam=f"{os.path.basename(config['bam'])}",
#     output:
#         results=f"xHLA/{os.path.basename(config['bam'])}.json",
#     container: config["docker_hla"]
#     threads: 4
#     resources:
#         mem_mb=1024*16,
#         runtime="12:00",
#     params:
#         sample=os.path.basename(config['bam']),
#     shell:
#         """ run.py \
#         --sample_id {params.sample} --input_bam_path {input.bam} \
#         --output_path xHLA
#         find .
#         """

# rule OptiType:
#     """ Couldn't get this to work from BAM, perhaps an issue with the reference or with BWA vs bowtie?  the Bam input is intended for rerunning optityper, not neccessaritly to go straight from an existing alignment. see allso https://github.com/FRED-2/OptiType/issues/64"""
#     input:
#         results=f"{os.path.basename(config['bam'])}.bai",
#         bam=f"{os.path.basename(config['bam'])}",
#     output:
#         results=f"optitype/{os.path.basename(config['bam'])}.json",
#     container:
#         "docker://fred2/optitype"
#     threads: 16
#     shell:"""
#     OptiTypePipeline.py -i {input.bam} --dna -o optitype -v
#     find optitype
#     """


rule OptiType_fastq:
    input:
        R1=config["R1"],
        R2=config["R2"],
    output:
        results=f"optitype_fastq/optitype_result.tsv",
    container:
        config["docker_optitype"]
    params:
        #        sample=os.path.splitext(os.path.basename(config['R1'][0]))[0]
        sample="optitype",
    threads: 16
    shell:
        """
        OptiTypePipeline.py -i {input.R1} {input.R2} --dna --prefix {params.sample} --outdir reports/ --verbose
        find reports
        """
