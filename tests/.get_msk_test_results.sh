# run on cluster
# isabl get-results -fi application.name "Shotgun Biobakery 1m" -fi application.version "1.0.0" -fi status SUCCEEDED -fi projects 5 -r metaphlan_ab > ~/GitHub/vdblab-shotgun/.test/metaphlan_test_results.txt
# isabl get-results -fi application.name "Shotgun Biobakery 1m" -fi application.version "1.0.0" -fi status SUCCEEDED -fi projects 5 -r humann_gene_families  > ~/GitHub/vdblab-shotgun/.test/humann_gf_test_results.txt
# isabl get-results -fi application.name "Shotgun Biobakery 1m" -fi application.version "1.0.0" -fi status SUCCEEDED -fi projects 5 -r humann_ko  > ~/GitHub/vdblab-shotgun/.test/humann_ko_test_results.txt
# isabl get-results -fi application.name "Shotgun Kraken Mock 1m" -fi application.version "0.2.0"  -fi status SUCCEEDED -fi projects 5 -r bracken_species_report  > ~/GitHub/vdblab-shotgun/.test/bracken_test_results.txt

# sort | uniq | is because I forgot to filter on app version
mkdir -p .tests/results/
scp lilac:/home/watersn/GitHub/vdblab-shotgun/.test/metaphlan_test_results.txt ./
rsync -av --progress --partial-dir=/tmp --no-relative --files-from=metaphlan_test_results.txt  lilac:/ .tests/results/
cat metaphlan_test_results.txt | sed 's|.*\/|.tests/results/|'   > local_metaphlan_test_results.txt


scp lilac:/home/watersn/GitHub/vdblab-shotgun/.test/humann_gf_test_results.txt ./
rsync -av --progress --partial-dir=/tmp --no-relative --files-from=humann_gf_test_results.txt  lilac:/ .tests/results/
cat humann_gf_test_results.txt | sed 's|.*\/|.tests/results/|'  > local_humann_gf_test_results.txt

scp lilac:/home/watersn/GitHub/vdblab-shotgun/.test/humann_ko_test_results.txt ./
rsync -av --progress --partial-dir=/tmp --no-relative --files-from=humann_ko_test_results.txt  lilac:/ .tests/results/
cat humann_ko_test_results.txt | sed 's|.*\/|.tests/results/|'  > local_humann_ko_test_results.txt

scp lilac:/home/watersn/GitHub/vdblab-shotgun/.test/bracken_test_results.txt ./
rsync -av --progress --partial-dir=/tmp --no-relative --files-from=bracken_test_results.txt  lilac:/ .tests/results/
cat bracken_test_results.txt | sed 's|.*\/|.tests/results/|'  > local_bracken_test_results.txt


