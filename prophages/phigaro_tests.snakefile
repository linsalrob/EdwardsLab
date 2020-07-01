phispydir = "/home3/redwards/GitHubs/PhiSpy/test_genbank_files"
GENOMES, = glob_wildcards(os.path.join(phispydir, '{genome}.gb.gz'))

phigarodir = "phigaro_tests"


rule all:
    input:
        expand(os.path.join(phigarodir, "{genome}_phigaro_tptn.tsv"), genome=GENOMES)

rule convert_gb_to_fna:
    input:
        gen = os.path.join(phispydir, "{genome}.gb.gz")
    output:
        fna = os.path.join(phigarodir, "{genome}.fna")
    shell:
        """
        python3 ~/GitHubs/EdwardsLab/bin/genbank2sequences.py -g {input.gen} -n {output.fna}
        """

rule run_phigaro:
    input:
        fna = os.path.join(phigarodir, "{genome}.fna")
    output:
        tsv = "{genome}_phigaro.tsv"
    params:
        tsv = "{genome}_phigaro" # phigaro adds the .tsv extension
    benchmark:
        os.path.join(phigarodir, "benchmarks", "{genome}_phigaro.txt")
    conda:
        "/home3/redwards/anaconda3/envs/phigaro"
    shell:
        """
        phigaro -f {input.fna} -e tsv -o {params.tsv} --delete-shorts
        """

rule phigaro_to_tbl:
    input:
        tsv = "{genome}_phigaro.tsv"
    output:
        os.path.join(phigarodir, "{genome}_phigaro_locs.tsv")
    shell:
        """
        if [ $(stat -c %s {input}) == 41 ]; then
            touch {output}
        else
            grep -v scaffold {input.tsv} | cut -f 1,2,3 > {output}
        fi
        """

rule count_tp_tn:
    input:
        gen = os.path.join(phispydir, "{genome}.gb.gz"),
        tbl = os.path.join(phigarodir, "{genome}_phigaro_locs.tsv")
    output:
        tp = os.path.join(phigarodir, "{genome}_phigaro_tptn.tsv")
    shell:
        """
        python3 ~/GitHubs/PhiSpy/scripts/compare_predictions_to_phages.py -t {input.gen} -r {input.tbl} > {output.tp}
        """
