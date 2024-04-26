from pathlib import Path


def is_paired():
    if config["lib_layout"] not in ["paired", "single"]:
        raise ValueError("lib_layout must be specified as either paired or single")

    if config["lib_layout"] == "single":
        if "R2" in config.keys():
            raise ValueError("lib_layout specified as single, but R2 has been provided")
        else:
            pass
    else:
        if "R2" not in config.keys():
            raise ValueError(
                "lib_layout specified as paired, but R2 has not been provided"
            )
    return config["lib_layout"] == "paired"


def skip_dedup():
    return not config["dedup_reads"]


def get_pipeline_version():
    return (
        subprocess.check_output(["git", "describe", "--always"], cwd=workflow.basedir)
        .decode("utf-8")
        .strip()
    )


def make_shard_names(nshards):
    return [f"{x:03}" for x in range(1, config["nshards"] + 1)]


def make_assembly_split_names(nparts):
    split_names = []
    # deal with the 3 digit ones first
    for i in range(1, min(nparts, 100)):
        split_names.append(f"{i:03}")
    # values after 99 are not padded
    if nparts > 99:
        for i in range(100, nparts + 1):
            split_names.append(f"{i}")
    return split_names

def get_concat_input(wc):
    return config[f"R{wc.rd}"]

def files_to_split(wildcards):
    """full disclosure: I don't remember why this is returning a dict of tuples"""
    res = {"R1": None}
    if not skip_dedup():
        res["R1"] = ((f"dedup/{wildcards.sample}_R1.fastq.gz"),)
        if is_paired():
            res["R2"] = ((f"dedup/{wildcards.sample}_R2.fastq.gz"),)
    else:
        res["R1"] = ((f"concatenated/{wildcards.sample}_R1.fastq.gz"),)
        if is_paired():
            res["R2"] = ((f"concatenated/{wildcards.sample}_R2.fastq.gz"),)
    return res


def files_to_trim(wildcards, nshards=1):
    # if splitting into shards, use output from split_fastq
    # if not and you deduplicated reads, use that output;
    # otherwise, just trim the concatenated inputs
    res = {"R1": None}
    if config["nshards"] > 1:
        res["R1"] = (
            (f"split_fastq/{wildcards.sample}_R1.part_{wildcards.shard}.fastq.gz"),
        )
        if is_paired():
            res["R2"] = (
                (f"split_fastq/{wildcards.sample}_R2.part_{wildcards.shard}.fastq.gz"),
            )
    elif not skip_dedup():
        res["R1"] = ((f"dedup/{wildcards.sample}_R1.fastq.gz"),)
        if is_paired():
            res["R2"] = ((f"dedup/{wildcards.sample}_R2.fastq.gz"),)
    else:
        res["R1"] = ((f"concatenated/{wildcards.sample}_R1.fastq.gz"),)
        if is_paired():
            res["R2"] = ((f"concatenated/{wildcards.sample}_R2.fastq.gz"),)
    return res


# def get_input_fastqs(wildcards):
#     if config["lib_layout"] == "paired":
#         return({"R1": config["R1"],
#                 "R2": config["R2"]})
#     else:
#         return ({"R1": config["R1"]})

# def get_input_fastqs_wc(wildcards):
#     return(" if config["lib_layout"] == "paired":
#         return({"R1": config["R1"],
#                 "R2": config["R2"]})
#     else:
#         return ({"R1": config["R1"]})


def get_config_inputs(wc):
    if is_paired():
        return {"R1": config["R1"], "R2": config["R2"]}
    else:
        return {
            "R1": config["R1"],
        }


def bbmap_dedup_params_flags(wildcards):
    # normal optical duplicate removal (HiSeq, MiSeq, etc)
    # should use these flags
    flags = "dedupe optical"

    if "dedup_platform" in config:
        if config["dedup_platform"] == "NextSeq":
            # see their docs; this is the recommended command for NextSeq
            # https://www.biostars.org/p/225338/
            flags = "dedupe optical spany adjacent"
        # if SRA, they tossed the read names so we get errors in bbmap if they
        # try to parse the SRA names for optical deduplication
        if config["dedup_platform"] == "SRA":
            flags = "dedupe"
        else:
            dupedist = bbmap_dedup_params_dupedist(wildcards)
            flags += f" dupedist={dupedist}"
    return flags


def bbmap_dedup_params_dupedist(wildcards):
    # default for unspecified platform, HiSeq up to 2500, and NextSeq
    dupedist = 40

    if config["dedup_platform"] == "NovaSeq":
        dupedist = 12000

    elif config["dedup_platform"] == "HiSeq-3000":
        dupedist = 2500

    elif config["dedup_platform"] == "HiSeq-4000":
        dupedist = 2500

    return dupedist


def test_bbmap_dedup_params_flags():
    test_conditions = [
        [{"dedup_platform": "SRA"}, "dedupe"],
        [{"dedup_platform": "MiniSeq"}, "dedupe optical dupedist=40"],
        [{"dedup_platform": "NovaSeq"}, "dedupe optical dupedist=12000"],
        [{"dedup_platform": "NextSeq"}, "dedupe optical spany adjacent dupedist=40"],
        [{"dedup_platform": "HiSeq-3000"}, "dedupe optical dupedist=2500"],
        [{"dedup_platform": "HiSeq-4000"}, "dedupe optical dupedist=2500"],
    ]
    for i in test_conditions:
        global config
        config = i[0]
        outcome = i[1]
        assert bbmap_dedup_params_flags(None) == outcome
        print(f"passed dedup flags test for {config['dedup_platform']}")


# Kraken functions
def parse_read_lengths(wildcards):
    ck_output = checkpoints.get_read_len.get(**wildcards).output.readlen_dir
    (readlens,) = glob_wildcards(os.path.join(ck_output, "len{readlen,[0-9]+}"))
    if len(readlens) == 0:
        raise ValueError("No read Length files found!")
    # here is some unpleasant logic; we don't want to have to recompute kmers for every concievable read length.
    # Some of our data is 101bp reads, while the prebuilt indexes contain the kmer profiles for 100, 150, etc.
    # So here we fudge the numbers a bit; and this seems justifiable:
    # From the PeerJ paper:
    #  > Bracken can use these probabilities for any metagenomics data set,
    # including data with different read lengths, although the estimates
    # might be slightly improved by re-computing with a read length that
    # matches the experimental data.
    # see https://benlangmead.github.io/aws-indexes/k2
    prebuilt_lengths = [50, 75, 100, 150, 200, 250, 300]
    new_readlens = []
    outpaths = []
    # if we only glob 1, don't treat the string as a list

    if not isinstance(readlens, list):
        readlens = [readlens]

    final_basenames = []

    for readlen in readlens:
        # If we don't have a prebuilt db for this readlength, pick the closest
        if int(readlen) not in prebuilt_lengths:
            # from https://stackoverflow.com/questions/12141150
            closest_read_length = min(
                prebuilt_lengths, key=lambda x: abs(x - int(readlen))
            )
            logger.warning(
                f"Read length {readlen}bp detected, using {closest_read_length}bp for bracken instead."
            )
            # return the path to the file tagged with that
            # closest/actaual readlength, creating if needed
            thispath = os.path.join(ck_output, f"len{closest_read_length}approx")
            Path(thispath).touch()
            final_basenames.append(os.path.basename(thispath))
        else:
            final_basenames.append(f"len{readlen}")

    return expand(os.path.join(ck_output, "{basenames}"), basenames=final_basenames)


def get_dbs_needed(wildcards):
    ck_output = checkpoints.get_read_len.get(**wildcards).output.readlen_dir
    (readlens,) = glob_wildcards(os.path.join(ck_output, "len{readlen}"))
    path_pattern = os.path.join(
        config["kraken2_db"], "database{readlen}mers.kmer_distrib"
    )
    return expand(path_pattern, readlen=readlens)
