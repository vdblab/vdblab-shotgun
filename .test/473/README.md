Shotgun test data generated as follows:

```
zcat /home/daia1/my_workdir/samples/Sample_473_IGO_12587_1/473_IGO_12587_1_S132_L003_R1_001.fastq.gz |head -n 4000 |gzip > /data/brinkvd/data/shotgun/test/473_IGO_12587_1_S132_L003_R1_001.fastq.gz
zcat /home/daia1/my_workdir/samples/Sample_473_IGO_12587_1/473_IGO_12587_1_S132_L003_R2_001.fastq.gz |head -n 4000 |gzip > /data/brinkvd/data/shotgun/test/473_IGO_12587_1_S132_L003_R2_001.fastq.gz
```

Ran through kneaddata to remove any possible human reads before committing.
This process mangles the fastq files which have been fixed like so:

```
zcat 473/473_IGO_12587_1_S132_L003_R1_001_FIXME.fastq.gz | sed 's|#0/1$||' | sed 's|:N:| 1:N:|' > 473/473_IGO_12587_1_S132_L003_R1_001.fastq
zcat 473_bak/473_IGO_12587_1_S132_L003_R2_001_FIXME.fastq.gz | sed 's|#0/2$||' | sed 's|:N:| 2:N:|' > 473/473_IGO_12587_1_S132_L003_R2_001.fastq
```
