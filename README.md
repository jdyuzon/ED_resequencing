# Genome Integrity of Exophiala dermatitidis adapted to ionizing radiation 
E. dermatitidis WT and pks1 were subjected to LD99 doses of acute gamma radiation for 15 consecutive weeks to produce several independently adapted strains. These samples were sequenced at Yale Center for Genome Analysis (YCGA). Evolved isolates and unevolved isolates were exposed to a final dose of 8000 Gy and sequenced at the Joint Genome Institute (JGI).

## Filtering and mapping FASTQ reads
FASTQ files from the JGI were processed using BBTools. Read files from YCGA were also processed in a similar manner.
```
Yale_sequences.06302022.pbs
```

## Single Nucleotide Polymorphism Calling
The following files show the samples and BAM used in SNP calling.
```
sample_names.jgi_resequencing.Jul052022.txt
bam_list.jgi_resequencing.Jul052022.txt
```

Samtools/BCFtools was used to call and filter SNPs. 
```
Samtools_JGI_resequencing.Jul052022.pbs
```

## Structural Variant calling
sv-callers was used to identify Structural Variants. The following files were used to configure the workflow:
```
WT_PKS.paired.samples.Jul042022.csv
analysis.yaml
```

## Identification of Parallel Mutation
Mutations that arose multiple times independently in evolved lineages were identified. We focused on mutations that persisted in evolved isolates and were associated with differential gene expression.
```
JGI_resequencing.R
```
