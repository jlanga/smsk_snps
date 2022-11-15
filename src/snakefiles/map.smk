rule map_bwa_index:
    """Index with bwa"""
    input:
        fa=RAW + "genome.fa",
    output:
        mock=touch(MAP_INDEX + "genome"),
        buckets=expand(
            MAP_INDEX + "genome.{suffix}", suffix="amb ann bwt pac sa".split()
        ),
    log:
        MAP_INDEX + "bwa_index.log",
    benchmark:
        MAP_INDEX + "bwa_index.bmk"
    conda:
        "map.yml"
    shell:
        "bwa index -p {output.mock} {input.fa} > {log} 2>&1"


def compose_rg_tag(wildcards):
    identifier = "ID:" + wildcards.population + "_" + wildcards.library
    library = "LB:truseq_" + wildcards.library
    platform = "PL:Illumina"
    sample = "SM:" + wildcards.population
    return f"@RG\\t{identifier}\\t{library}\\t{platform}\\t{sample}"


rule map_bwa_map:
    """Map population with bowtie2, sort with samtools, compress to cram"""
    input:
        fwd=QC + "{population}.{library}_1.fq.gz",
        rev=QC + "{population}.{library}_2.fq.gz",
        fwd_unp=QC + "{population}.{library}_3.fq.gz",
        rev_unp=QC + "{population}.{library}_4.fq.gz",
        index=MAP_INDEX + "genome",
        reference=RAW + "genome.fa",
    output:
        cram=protected(MAP + "{population}.{library}.cram"),
    params:
        extra=params["bwa"]["extra"],
        rg_tag=compose_rg_tag,
    threads: params["bwa"]["threads"]
    log:
        MAP + "{population}.{library}.bwa_mem.log",
    benchmark:
        MAP + "{population}.{library}.bwa_mem.bmk"
    conda:
        "map.yml"
    shell:
        """
        (bwa mem \
            -M \
            -R '{params.rg_tag}' \
            -t {threads} \
            {params.extra} \
            {input.index} \
            {input.fwd} \
            {input.rev} \
        | samtools fixmate - - \
        | samtools sort \
            -l 9 \
            -o {output.cram} \
            --reference {input.reference} \
            --output-fmt CRAM \
            -@ {threads} \
            /dev/stdin \
        ) 2> {log}
        """


rule map:
    input:
        [
            MAP + f"{population}.{library}.cram"
            for population, library in (
                samples[["population", "library"]].values.tolist()
            )
        ],
