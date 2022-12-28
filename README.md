# structural variation

## manta and AnnotSV

```shell
(base) hojin@ihojin-ui-MacBookPro scripts % sh manta_annotsv_Dec212022_hj.sh -h
manta_annotsv_Dec212022_hj.sh -i [bam or cram] -d [path to sample directory]
ex ) manta_annotsv_Dec212022_hj.sh -i ../../P1.cram -d ./P1
  -i : input file ex ) *.bam or *.cram
  -d : directory with sampleID as name
  -h : show this message
```

![image](https://user-images.githubusercontent.com/121307215/209815820-a832d3e2-83e2-4406-9390-97e4c23e48d4.png)


this script makes directory with sampleID as name. finally, the following files are obtained for each sample.  
```shell
(base) hojin@ihojin-ui-MacBookPro P1 % tree
.
├── manta_annotsv
│   ├── P1_manta.pass.mantadel.dom.filtered.txt
│   ├── P1_manta.pass.mantadel.rec.filtered.txt
│   ├── P1_manta.pass.mantadup.dom.filtered.txt
│   ├── P1_manta.pass.mantadup.rec.filtered.txt
│   ├── P1_manta.pass.mantains.dom.filtered.txt
│   ├── P1_manta.pass.mantains.rec.filtered.txt
│   └── diploidSV_edit.annotated.tsv
├── results
│   ├── evidence
│   ├── stats
│   │   ├── alignmentStatsSummary.txt
│   │   ├── svCandidateGenerationStats.tsv
│   │   ├── svCandidateGenerationStats.xml
│   │   └── svLocusGraphStats.tsv
│   └── variants
│       ├── candidateSV.vcf.gz
│       ├── candidateSV.vcf.gz.tbi
│       ├── candidateSmallIndels.vcf.gz
│       ├── candidateSmallIndels.vcf.gz.tbi
│       ├── diploidSV.vcf.gz
│       └── diploidSV.vcf.gz.tbi
```

and then, to get the *geneCounts.txt, the results of all samples are concatenated using the R script. 
```shell
(base) hojin@ihojin-ui-MacBookPro manta_geneCounts % ls /Volumes/hjdrive/thyroiditis/thyroiditis_sv/manta_geneCounts/data/*/manta_annotsv/*manta.pass.mantadel.dom.filtered.txt > 20221228_manta_annotsv_dom.del.list.txt
(base) hojin@ihojin-ui-MacBookPro manta_geneCounts % Rscript 20221228_manta_annotsv_geneCount_deletion_Dec282022_hj.R
```

