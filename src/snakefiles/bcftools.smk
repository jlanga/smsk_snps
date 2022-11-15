def get_library_files_from_sample(wildcards):
    """TODO: needs improvement/simplification
    Return the list of libraries corresponding to a population and chromosome.
    """
    population = wildcards.population
    libraries = samples[samples["population"] == population]["library"].values.tolist()
    files = [MAP + population + "." + library + ".cram" for library in libraries]
    return files


rule bcftools_mpileup:
    """Compute the mpileup and compress it"""
    input:
        cram=get_library_files_from_sample,
        fa=RAW + "genome.fa",
        fai=RAW + "genome.fa.fai",
    output:
        mpileup_gz=BCFTOOLS + "{population}.vcf.gz",
    log:
        BCFTOOLS + "{population}.log",
    benchmark:
        BCFTOOLS + "{population}.bmk"
    conda:
        "bcftools.yml"
    shell:
        """
        (samtools merge \
            -u \
            --reference {input.fa} \
            - \
            {input.cram} \
        | bcftools mpileup \
            --output-type u \
            --fasta-ref {input.fa} \
            - \
        | bcftools call \
            --variants-only \
            --multiallelic-caller \
            --output-type z \
            --output {output.mpileup_gz} \
        ) 2> {log}
        """


rule bcftools:
    input:
        [BCFTOOLS + f"{population}.vcf.gz" for population in POPULATIONS],
