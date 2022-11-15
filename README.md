# smsk_snps: A Snakemake pipeline for SNP calling


## 1. Description

This is a repo that contains installers and snakemake scripts to execute the pipelines described by Kofler et al. in popoolation 1 and 2:

- QC with Trimmomatic

- Mapping with bwa

- BAM wrangling and SNP calling with samtools

- SNP calling with bcftools


## 2. First steps


1. Install ([ana](https://www.continuum.io/downloads)|[mini](http://conda.pydata.org/miniconda.html))conda

2. Clone and install the software

    ```sh
    git clone https://github.com/jlanga/smsk_snps.git smsk_snps
    cd smsk_snps
    snakemake --use-conda --conda-frontend mamba --create-envs-only
    ```

3. Run the test dataset:

    ```sh
    snakemake --use-conda -
    ```

4. Modify the following files:

    - `features.yaml`: the path to the genome reference, and the names of every chromosome to be processed.

    - `samples.tsv`: paths and library information of each of the samples.

    - `params.yml`: execution parameters of different tools.

    - `cluster.json`: memory limits for some of the steps-

5. Execute the pipeline:

    ```sh
    snakemake --use-conda --conda-frontend mamba -j 40
    ```

## Representation of the pipeline

![smsk_popoolation pipeline](https://cdn.rawgit.com/jlanga/smsk_popoolation/master/rulegraph.svg)

## Bibliography

- [A Quick Guide to Organizing Computational Biology Projects](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424)

- [Snakemake—a scalable bioinformatics workflow engine](http://bioinformatics.oxfordjournals.org/content/28/19/2520)

- [popoolation](https://sourceforge.net/p/popoolation/wiki/Main/)

- [popoolation2](https://sourceforge.net/p/popoolation2/wiki/Main/)

- [samtools](http://www.htslib.org/)
