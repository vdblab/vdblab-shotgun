from pathlib import Path


def get_pipeline_version():
    return (
        subprocess.check_output(["git", "describe", "--always"], cwd=workflow.basedir)
        .decode("utf-8")
        .strip()
    )


def make_shard_names(nshards):
    return [f"{x:03}" for x in range(1, config["nshards"] + 1)]


def files_to_split(wildcards, dedup=False, read_dir=1):
    if dedup:
        return ((f"dedup/{wildcards.sample}_R{read_dir}.fastq.gz"),)
    else:
        return ((f"concatenated/{wildcards.sample}_R{read_dir}.fastq.gz"),)


def files_to_trim(wildcards, nshards=1, dedup=False, read_dir=1):
    # if splitting into shards, use output from split_fastq
    # if not and you deduplicated reads, use that output;
    # otherwise, just trim the concatenated inputs
    if nshards > 1:
        return (
            (
                f"split_fastq/{wildcards.sample}_R{read_dir}.part_{wildcards.shard}.fastq.gz"
            ),
        )
    elif dedup:
        return f"dedup/{wildcards.sample}_R{read_dir}.fastq.gz"
    else:
        return ((f"concatenated/{wildcards.sample}_R{read_dir}.fastq.gz"),)


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
    # So here we fudge the numbes a bit; and this seems justifiable:
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
            final_basenames.append(f"len{closest_read_length}")

    return expand(os.path.join(ck_output, "{basenames}"), basenames=final_basenames)


def get_dbs_needed(wildcards):
    ck_output = checkpoints.get_read_len.get(**wildcards).output.readlen_dir
    (readlens,) = glob_wildcards(os.path.join(ck_output, "len{readlen}"))
    path_pattern = os.path.join(
        config["kraken2_db"], "database{readlen}mers.kmer_distrib"
    )
    return expand(path_pattern, readlen=readlens)
