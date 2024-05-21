set -o pipefail
echo "Running preprocessing on simulated data with equal reads per organism"
bash test.sh preprocess equalreads
rm -r  tmppreprocess_equalreads/sortmerna/tmp*
echo "Running preprocessing on simulated data with equal coverage per organism"
bash test.sh preprocess equalcov
rm -r  tmppreprocess_equalcov/sortmerna/tmp*
echo "Running preprocessing on simulated data with sharded data"
bash test.sh preprocess equalreads4shards
rm -r  tmppreprocess_equalreads4shards/sortmerna/tmp*
echo "Running biobakery on simulated data with equal reads per organism"
bash test.sh biobakery equalreads
echo "Running kraken on simulated data with equal reads per organism"
bash test.sh kraken equalreads

echo "Running biobakery on simulated single-end data with equal reads per organism"
bash test.sh preprocess-se equalreads
rm -r  tmppreprocess_equalreads/sortmerna/tmp*
echo "Running biobakery on simulated single-end data with equal reads per organism"
bash test.sh biobakery-se equalreads
echo "Running biobakery on simulated single-end data with equal coverage per organism"
bash test.sh biobakery-se equalcov


pytest
