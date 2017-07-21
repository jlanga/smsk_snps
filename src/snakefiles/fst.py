rule fst_sliding_chromosome:
    input:
        sync_gz = SYNC_SUB + "{chromosome}.sync.gz"
    output:
        tsv_gz  = TABLE_FST + "{chromosome}.tsv.gz"
    params:
        sync = SYNC_SUB + "{chromosome}.sync",
        tsv  = TABLE_FST + "{chromosome}.tsv",
        window_size = config["popoolation2_params"]["fst_sliding"]["window_size"],
        step_size = config["popoolation2_params"]["fst_sliding"]["step_size"],
        min_covered_fraction = config["popoolation2_params"]["fst_sliding"]["min_covered_fraction"],
        min_coverage = config["popoolation2_params"]["fst_sliding"]["min_coverage"],
        max_coverage = config["popoolation2_params"]["fst_sliding"]["max_coverage"],
        pool_size = config["popoolation2_params"]["fst_sliding"]["pool_size"],
        min_count = config["popoolation2_params"]["fst_sliding"]["min_count"]
    threads:
        1
    log:
        TABLE_FST + "{chromosome}.log"
    benchmark:
        TABLE_FST + "{chromosome}.json"
    shell:
        "pigz --decompress --stdout {input.sync_gz} > {params.sync} 2> {log} ; "
        "perl src/popoolation2_1201/fst-sliding.pl "
            "--window-size {params.window_size} "
            "--step-size {params.step_size} "
            "--suppress-noninformative "
            "--input {params.sync} "
            "--min-covered-fraction {params.min_covered_fraction} "
            "--min-coverage {params.min_coverage} "
            "--max-coverage {params.max_coverage} "
            "--min-count {params.min_count} "
            "--output {params.tsv} "
            "--pool-size {params.pool_size} "
        "2>> {log} ; "
        "pigz --best {params.tsv} 2>> {log}"
 



rule fst_plot:
    input:
        tsvs_gz = expand(
            TABLE_FST + "{chromosome}.tsv.gz",
            chromosome = CHROMOSOMES
        )
    output:
        merged_tsv_gz = PLOT_FST + "all.tsv.gz",
        z_pdfs = [
            PLOT_FST + str(i) + "_" + str(j) +"_z.pdf" 
            for i in range(1, len(POPULATIONS)) 
            for j in range(i+1, len(POPULATIONS)+1)
        ],
        pdfs = [ 
            PLOT_FST + str(i) + "_" + str(j) +".pdf" 
            for i in range(1, len(POPULATIONS))
            for j in range(i+1, len(POPULATIONS)+1)
        ]
    params:
        plot_fst = PLOT_FST,
        merged_tsv = PLOT_FST + "all.tsv",
        n_pop = len(POPULATIONS) 
    threads:
        1
    log:
        PLOT_FST + "plot.log"
    benchmark:
        PLOT_FST + "plot.json"
    shell:
        "pigz --decompress --stdout {input.tsvs_gz} > {params.merged_tsv} 2> {log} ; "
        "for i in `seq 1 {params.n_pop}`; do "
            "for j in `seq $(($i + 1)) {params.n_pop}`; do "
                "cat {params.merged_tsv} "
                    "| python src/fst_to_genomic_score.py "
                        "$(( $i - 1 )) "
                        "$(( $j - 1 )) "
                "> {params.plot_fst}/${{i}}_${{j}}.fst "
                "2>> {log} ; "
                "Rscript src/plot_score.R "
                    "none "
                    "{params.plot_fst}/${{i}}_${{j}}.fst "
                    "{params.plot_fst}/${{i}}_${{j}}.pdf "
                "2>> {log} ; "
                "Rscript src/plot_score.R "
                    "z "
                    "{params.plot_fst}/${{i}}_${{j}}.fst "
                    "{params.plot_fst}/${{i}}_${{j}}_z.pdf "
                "2>> {log} ; "
            "done ; "
        "done ; "
        "pigz --best {params.merged_tsv} "