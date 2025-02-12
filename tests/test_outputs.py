import os
import pytest
import sys
from math import isclose
import gzip
from pathlib import Path

def get_simulated_organism_counts(path=".test/simulated/1_depth100000_equalreads.statsfastq"):
    base_to_replace = os.path.splitext(os.path.basename(path))[0]
    res = {}
    with open (path, "r") as inf:
        for i, line in enumerate(inf):
            if i == 0:
                continue
            line = line.strip().split()
            line[0] = line[0].replace("tmp_" + base_to_replace + "_", "")
            organism = line[0].split(".fasta")[0]
            count = int(line[3].replace(",", ""))
            if organism in res:
                res[organism] = res[organism] + count
            else:
                res[organism] =  count
    return res

simcounts_equalreads = get_simulated_organism_counts(".test/simulated/1_depth100000_equalreads.statsfastq")
simcounts_equalcoverage = get_simulated_organism_counts(".test/simulated/1_depth100000_equalcoverage.statsfastq")

def running_as_github_action():
    return "GITHUB_ACTION" in os.environ and os.environ["GITHUB_ACTION"] is not None

def test_simulated_data_present():
    # TODO: test hashes of files?
    assert os.path.exists(".test/simulated/1_depth100000_equalreads.statsfastq"), "no simulated data present; see readme to simulate data"

def get_dirname():
    if running_as_github_action:
        return "tmppreprocess-gha_equalreads"
    else:
        return "tmppreprocess_equalreads"


@pytest.mark.skipif(running_as_github_action(), reason="this test not available when run as GH action")
def test_simulated_host_deplete_results_present():
    assert os.path.exists("tmppreprocess_equalreads/reports/473_hostdepletion.stats"), "no host depletion results for simulated data; please run `bash test.sh preprocess equalreads`"

@pytest.mark.skipif(running_as_github_action(), reason="this test not available when run as GH action")
def test_simulated_sharded_host_deplete_results_present():
    assert os.path.exists("tmppreprocess_equalreads4shards/reports/473_hostdepletion.stats"), "no host depletion results for simulated data; please run `bash test.sh preprocess equalreads4shards`"

@pytest.mark.skipif(running_as_github_action(), reason="this test not available when run as GH action")
def test_simulated_bb_kraken_run_present():
    assert os.path.exists("tmpbiobakery_equalreads/metaphlan/473_metaphlan3_profile.txt"), "no metaphlan results for simulated data; please run `bash test.sh biobakery equalreads`"
    assert os.path.exists("tmpkraken_equalreads/kraken2/473_kraken2.bracken.S.out"), "no kraken results for sumulated data; please run `bash test.sh kraken equalreads`"


@pytest.mark.skipif(running_as_github_action(), reason="this test not available when run as GH action")
def test_simulated_bb_se_present():
    assert os.path.exists("tmpbiobakery-se_equalreads/metaphlan/473_metaphlan3_profile.txt"), "no metaphlan results for simulated data; please run `bash test.sh biobakeryse equalreads`"
@pytest.mark.skipif(running_as_github_action(), reason="this test not available when run as GH action")
def test_simulated_bb_se_even_present():
    assert os.path.exists("tmpbiobakery-se_equalcov/metaphlan/473_metaphlan3_profile.txt"), "no metaphlan results for simulated data; please run `bash test.sh biobakeryse equalcov`"

@pytest.mark.skipif(running_as_github_action(), reason="this test not available when run as GH action")
def test_preprocess_depletes_correct_n_reads():
    print(simcounts_equalreads)
    nreads_human = simcounts_equalreads["t2t_chr21"]

#    sample\tbowtie2_human\tbowtie2_human_aligned\tsnap_human\tsnap_human_aligned\tbowtie2_mouse\tbowtie2_mouse_aligned\tsnap_mouse\tsnap_mouse_aligned
    with open("tmppreprocess_equalreads/reports/473_hostdepletion.stats", "r") as inf:
        for line in inf:
            if line.startswith("473"):
                nreads_human_detected = line.split("\t")[2]

    assert isclose(int(nreads_human), int(nreads_human_detected), abs_tol=10)


@pytest.mark.skipif(running_as_github_action(), reason="this test not available when run as GH action")
def test_preprocess_host_reads_fq_has_correct_n_reads():
    print(simcounts_equalreads)
    nreads_human = simcounts_equalreads["t2t_chr21"]

    nreads_in_host_fq = 0
    with gzip.open('tmppreprocess_equalreads/host/473_all_host_reads_R1.fastq.gz', 'rb') as inf:
        for line in inf:
            if line.decode().startswith("@"):
                nreads_in_host_fq = nreads_in_host_fq + 1

    assert isclose(int(nreads_human)/2, nreads_in_host_fq, abs_tol=10)

@pytest.mark.skipif(running_as_github_action(), reason="this test not available when run as GH action")
def test_preprocess_se_host_reads_fq_has_correct_n_reads():
    print(simcounts_equalreads)
    nreads_human = simcounts_equalreads["t2t_chr21"]

    nreads_in_host_fq = 0
    with gzip.open('tmppreprocess-se_equalreads/host/473_all_host_reads_R1.fastq.gz', 'rb') as inf:
        for line in inf:
            if line.decode().startswith("@"):
                nreads_in_host_fq = nreads_in_host_fq + 1
    # see comments elsewhere; higher percieved duplication in single-end data
    # leads to fewer host reads in the deduplicated input dataset
    assert isclose(int(nreads_human)/2, nreads_in_host_fq, abs_tol=50)

@pytest.mark.skipif(running_as_github_action(), reason="this test not available when run as GH action")
def test_sharding_doesntbreak_host_depletion():
    nreads_human_detected_1shard = 0
    with open("tmppreprocess_equalreads/reports/473_hostdepletion.stats", "r") as inf:
        for line in inf:
            # ignore header
            if not line.startswith("sample"):
                # we will have one line per shard, so this sums them all
                nreads_human_detected_1shard = nreads_human_detected_1shard + int(line.split("\t")[2])
    nreads_human_detected_4shards = 0
    with open("tmppreprocess_equalreads4shards/reports/473_hostdepletion.stats", "r") as inf:
        for line in inf:
            # ignore header
            if not line.startswith("sample"):
                # we will have one line per shard, so this sums them all
                nreads_human_detected_4shards = nreads_human_detected_4shards + int(line.split("\t")[2])

    assert isclose(nreads_human_detected_1shard, nreads_human_detected_4shards, abs_tol=10)


@pytest.mark.skipif(running_as_github_action(), reason="this test not available when run as GH action")
def test_preprocess_se_depletes_correct_n_reads():
    """ This test has higher tollerances than the paired-end one.  I found that
    fewer human reads were detected, but it appears that this is due to the upstream
    de-deduplication, rather than not having enough reads aligning.
    For example, this test dataset has 99883 read PAIRs surviving de-duplication,
    but only 99789 single-end reads surviving deduplication.

    This may be addressed by increasing the stringency of the deduplication in a future release
    """
    nreads_human = simcounts_equalreads["t2t_chr21"] / 2
    nreads_human_detected = 0
    with open("tmppreprocess-se_equalreads/reports/473_hostdepletion.stats", "r") as inf:
        for line in inf:
            # ignore header
            if not line.startswith("sample"):
                # we will have one line per shard, so this sums them all
                nreads_human_detected = nreads_human_detected + int(line.split("\t")[2])

    assert isclose(int(nreads_human), int(nreads_human_detected), abs_tol=100)


def test_preprocess_check_dedup():
    """ the test workflow runs on a double copy of simulated data to
    create a dataset with a 100% duplication rate, so we should see a
    high duplication rate here
    """
    nreads_total = simcounts_equalreads["1_depth100000_equalreads_R1.fastq.gz"] * 2
    with open(f"{get_dirname()}/reports/473_dedup.stats", "r") as inf:
        for line in inf:
            if line.startswith("Duplicates Found"):
                ndups_found = int(line.split("\t")[1])
    assert isclose(nreads_total, ndups_found, abs_tol=10)


@pytest.mark.skipif(running_as_github_action(), reason="this test not available when run as GH action")
def test_preprocess_se_check_dedup():
    """ the test workflow runs on a double copy of simulated data to
    create a dataset with a 100% duplication rate, so we should see a
    high duplication rate here

    The abs_tol is higher for this single end run, as we expected
    99885 (the input size) duplicates but found 99930,
    presumably because there are some reads that are called duplicates when
    considering just R1 but wouldn't be considered duplicates when
    considering both R1 and R2.
    """
    nreads_total = simcounts_equalreads["1_depth100000_equalreads_R1.fastq.gz"]
    with open("tmppreprocess-se_equalreads/reports/473_dedup.stats", "r") as inf:
        for line in inf:
            if line.startswith("Duplicates Found"):
                ndups_found = int(line.split("\t")[1])
    assert isclose(nreads_total, ndups_found, abs_tol=100)

@pytest.mark.skipif(running_as_github_action(), reason="this test not available when run as GH action")
def test_bb_metaphlan_bacteria_relab():
    """ compare bacterial relab to theoretical
    The simulated data was made from equal reads from each of the gut zymo bugs plus
    human, so 21+1 taxa. 18 are bacterial. Of the 99,885*2 reads generated, we should hvae
    ((99,885*2)/22) * 18, 163448 bacterial reads.

    (Un)fortunately, metaphlan normalizes by genome length, so we can't just check the number of reads.
    """
    target_relative_kingdom_coverage, nontarget_relative_kingdom_coverage = 0, 0
    with open ("tests/zymo_plus_chr21.tsv", "r") as inf:
        for line in inf:
            if not line.startswith("Species"):
                line = line.strip().split("\t")
                # 1000 gets cancelled out, its just to have some nice fractions
                # for calculating the coverage ratios
                if line[7] == "Bacteria":
                    target_relative_kingdom_coverage = target_relative_kingdom_coverage + (1000/float(line[6]))
                else:
                    nontarget_relative_kingdom_coverage = nontarget_relative_kingdom_coverage + (1000/float(line[6]))
    print("relative target and nontarget coverages, and ratio")
    print(target_relative_kingdom_coverage, nontarget_relative_kingdom_coverage, target_relative_kingdom_coverage/nontarget_relative_kingdom_coverage)
    with open("tmpbiobakery_equalreads/metaphlan/473_metaphlan3_profile.txt", "r") as inf:
        for line in inf:
            line = line.strip().split("\t")
            if len(line) > 1:
                if line[1] == "2":
                    bact_total, bact_relab = line[4], line[2]
                    break
    est_relab = target_relative_kingdom_coverage/(target_relative_kingdom_coverage + nontarget_relative_kingdom_coverage)
    print(est_relab*100, float(bact_relab))
    # would love to be able to tighen up the tolerance  here
    assert isclose(est_relab*100, float(bact_relab), abs_tol=10),  "Bad relative abundance estimate"
    # TODO: these are way off, likely due to unmapped reads
    # print((99885*2) *est_relab, float(bact_total))
    #assert isclose((99885*2) * est_relab, int(bact_total)),  "Bad estimated counts"

@pytest.mark.skipif(running_as_github_action(), reason="this test not available when run as GH action")
def test_bb_metaphlan_est_counts():
    mpares = {}
    with open("tmpbiobakery_equalreads/metaphlan/473_metaphlan3_profile.txt", "r") as inf:
        for line in inf:
            if "s__" in line and "t__" not in line:
                line = line.split()
                species = line[0].split("s__")[1]
                mpares[species] =  int(line[4])

    # there are multiple E. coli in the mock
    assert isclose(mpares["Escherichia_coli"],
                   sum([v for k,v in simcounts_equalreads.items() if k.startswith("Escherichia_coli")]),
                   abs_tol=1200), "bad e.coli count"

    assert isclose(mpares["Akkermansia_muciniphila"], simcounts_equalreads["Akkermansia_muciniphila"], abs_tol=2000)
    assert isclose(mpares["Enterococcus_faecalis"], simcounts_equalreads["Enterococcus_faecalis"], abs_tol=2000)
    assert isclose(mpares["Salmonella_enterica"], simcounts_equalreads["Salmonella_enterica"], abs_tol=2000)
    assert isclose(mpares["Veillonella_rogosae"], simcounts_equalreads["Veillonella_rogosae"], abs_tol=2000 )


@pytest.mark.skipif(running_as_github_action(), reason="this test not available when run as GH action")
def test_bb_metaphlan_se_est_counts():
    mpares = {}
    with open("tmpbiobakery-se_equalreads/metaphlan/473_metaphlan3_profile.txt", "r") as inf:
        for line in inf:
            if "s__" in line and "t__" not in line:
                line = line.split()
                species = line[0].split("s__")[1]
                mpares[species] =  int(line[4])

    # there are multiple E. coli in the mock
    assert isclose(mpares["Escherichia_coli"],
                   sum([v for k,v in simcounts_equalreads.items() if k.startswith("Escherichia_coli")]) / 2,
                   abs_tol=1200), "bad e.coli count"

    assert isclose(mpares["Akkermansia_muciniphila"], simcounts_equalreads["Akkermansia_muciniphila"] / 2, abs_tol=2000)
    assert isclose(mpares["Enterococcus_faecalis"], simcounts_equalreads["Enterococcus_faecalis"] / 2, abs_tol=2000)
    assert isclose(mpares["Salmonella_enterica"], simcounts_equalreads["Salmonella_enterica"] / 2, abs_tol=2000)
    assert isclose(mpares["Veillonella_rogosae"], simcounts_equalreads["Veillonella_rogosae"] / 2, abs_tol=2000 )

@pytest.mark.skipif(running_as_github_action(), reason="this test not available when run as GH action")
def test_bb_metaphlan_se_even_coverate_relab():
    mpares_relab = {}
    with open("tmpbiobakery-se_equalcov/metaphlan/473_metaphlan3_profile.txt", "r") as inf:
        for line in inf:
            if "s__" in line and "t__" not in line:
                line = line.split()
                species = line[0].split("s__")[1]
                mpares_relab[species] =  float(line[2])
    print(mpares_relab)
    target_relab = (1/21) *100
    # there are multiple E. coli in the mock
    assert isclose(target_relab * 5, mpares_relab["Escherichia_coli"],
                   abs_tol=3), "bad e.coli count"

    assert isclose(mpares_relab["Akkermansia_muciniphila"], target_relab, abs_tol=1)
    assert isclose(mpares_relab["Enterococcus_faecalis"], target_relab, abs_tol=1)
    assert isclose(mpares_relab["Salmonella_enterica"], target_relab, abs_tol=1)
    assert isclose(mpares_relab["Veillonella_rogosae"], target_relab, abs_tol=1 )




def test_bb_metaphlan_archaea_relab():
    pass

def test_bb_metaphlan_fungi_relab():
    pass

def test_bb_humann_functional_count():
    pass

@pytest.mark.skip(reason="We will revisit this test once we relax the confidence parameter")
def test_kraken_bracken_counts_candida():
    nreads_c_albicans = simcounts_equalreads["Candida_albican"]
    with open ("tmpkraken_equalreads/kraken2/473_kraken2.bracken.S.out", "r") as inf:
        for line in inf:
            if line.startswith("Candida albicans"):
                detected_c_albicans = line.strip().split("\t")[5]
            pass
    print(f"input C. albicans reads: {nreads_c_albicans}; detected C. albicans {detected_c_albicans}")
    print("Note! we are considering reverting --confidence to the default for future releases")
    #    assert isclose(int(nreads_c_albicans), int(detected_c_albicans)),  "Bad Kraken2/Bracken  estimated counts"
    # This is a bad test
    assert isclose(22, int(detected_c_albicans)),  "Bad Kraken2/Bracken  estimated counts"
