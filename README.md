SeqMinerCmd Tutorial
========================================================
A command line tool for [SeqMiner][SeqMiner_link]

(Updated October 1, 2013)

[SeqMinerCmd][SeqMinerCmd_link] is a convenient data extraction tools for [VCF][VCF_Format_link] and [BCF][BCF_Format_link] files. The software, `SeqMinerCmd`, can be downloaded in [this link][SeqMinerCmd_zip_link] and its source code can be obtained from [github][SeqMinerCmd_link].

This tutorial demonstrates by-region and by-gene based methods, and the extraction results can be stored as a R matrix or a list. 

Extract gentoype matrix from VCF file
--------------------------------------

1. extract by region

```bash
./seqMiner -r "1:196621007-196634467" vcf/all.anno.filtered.extract.bcf.gz
```

```
## 1 region to be extracted.
## ----- 1:196621007-196634467 -----
## Position	NA12286	NA12341	NA12342
## 1:196623337	1	2	0
## 1:196632129	0	2	0
## 1:196632470	1	2	0
## 1:196633606	2	2	2
```

2. extract by gene

```bash
./seqMiner --geneFile vcf/refFlat_hg19_6col.txt.gz -n CFH vcf/all.anno.filtered.extract.vcf.gz
```

```
## 1 region to be extracted.
## ----- CFH -----
## Position	NA12286	NA12341	NA12342
## 1:196623337	1	2	0
## 1:196632129	0	2	0
## 1:196632470	1	2	0
## 1:196633606	2	2	2
## 1:196623337	1	2	0
## 1:196632129	0	2	0
## 1:196632470	1	2	0
## 1:196633606	2	2	2
```


Extract arbitrary fields from VCF file
--------------------------------------

1. extract by region

```bash
./seqMiner -r "1:196621007-196634467" -e CHROM,POS:DP:GT,GD vcf/all.anno.filtered.extract.vcf.gz
```

```
## CHROM	POS	DP	NA12286:GD	NA12341:GD	NA12342:GD	NA12286:GT	NA12341:GT	NA12342:GT
## 1	196623337	23	0	2	21	0/1	1/1	0/0
## 1	196632129	25	1	7	17	0/0	1/1	0/0
## 1	196632470	28	3	7	18	0/1	1/1	0/0
## 1	196633606	26	2	5	19	1/1	1/1	1/1
```

2. extract by gene

```bash
./seqMiner --geneFile vcf/refFlat_hg19_6col.txt.gz -n CFH -e CHROM,POS:DP vcf/all.anno.filtered.extract.vcf.gz
```

```
## CHROM	POS	DP
## 1	196623337	23
## 1	196632129	25
## 1	196632470	28
## 1	196633606	26
## 1	196623337	23
## 1	196632129	25
## 1	196632470	28
## 1	196633606	26
```


Contact
-------

SeqMiner is developed by [Xiaowei Zhan][zhanxw_link] and [Dajiang Liu][dajiang_link].
We welcome your questions and feedbacks.

[SeqMiner_link]: http://cran.r-project.org/web/packages/seqminer/index.html
[Vcf2geno_link]: http://cran.r-project.org/web/packages/vcf2geno/index.html
[VCF_Format_link]: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
[BCF_Format_link]: http://www.1000genomes.org/wiki/analysis/variant-call-format/bcf-binary-vcf-version-2
[Anno_link]: https://github.com/zhanxw/anno
[TabAnno_link]: https://github.com/zhanxw/anno
[TASER_link]: http://zhanxw.com/taser/
[Tabix_link]: http://sourceforge.net/projects/samtools/files/tabix/
[zhanxw_link]: mailto:zhanxw@gmail.com
[dajiang_link]: mailto:dajiang.liu@gmail.com
[1kg_link]: http://www.1000genomes.org/
[1kg_population_link]: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20111108_samples_pedigree/20111108_1000genomes_samples.xls
[SKAT_link]: http://www.hsph.harvard.edu/skat/
[SeqMinerCmd_link]: https://github.com/zhanxw/SeqMinerCmd
[SeqMinerCmd_zip_link]: https://github.com/zhanxw/SeqMinerCmd/archive/master.zip


[![Bitdeli Badge](https://d2weczhvl823v0.cloudfront.net/zhanxw/seqminercmd/trend.png)](https://bitdeli.com/free "Bitdeli Badge")

