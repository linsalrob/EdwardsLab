########################################################################
#                                                                      #
# Run and test PhageBoost                                              #
#                                                                      #
# PhageBoost is available from                                         #
# https://github.com/ku-cbd/PhageBoost                                 #
# https://www.biorxiv.org/content/10.1101/2020.08.09.243022v1.full.pdf #
#                                                                      #
#                                                                      #
# For this test, we need phageboost and phispy installed in conda      #
#                                                                      #
########################################################################



phispydir = "/home3/redwards/GitHubs/PhiSpy/test_genbank_files"
GENOMES, = glob_wildcards(os.path.join(phispydir, '{genome}.gb.gz'))

phageboostdir = "phageboost_tests"


rule all:
    input:
        expand(os.path.join(phageboostdir, "{genome}_phageboost_tptn.tsv"), genome=GENOMES)


rule run_phageboost:
    input:
        gen = os.path.join(phispydir, "{genome}.gb.gz")
    output:
        tsv = os.path.join(phageboostdir, "{genome}_phageboost.tsv")
    benchmark:
        os.path.join(phageboostdir, "benchmarks", "{genome}_phageboost.txt")
    shell:
        """
        python3 ~/GitHubs/EdwardsLab/prophages/phageboost_genbank.py -g {input.gen} -o {output.tsv} -m /home3/redwards/phage/prophage/model_delta_std_hacked.pickled.silent.gz
        """

rule phageboost_to_tbl:
    input:
        tsv = os.path.join(phageboostdir, "{genome}_phageboost.tsv")
    output:
        os.path.join(phageboostdir, "{genome}_phageboost_locs.tsv")
    shell:
        """
        if [ $(stat -c %s {input}) -lt 50 ]; then
            touch {output}
        else
            grep -v probability {input.tsv} | cut -f 3,4,5 > {output}
        fi
        """

rule count_tp_tn:
    input:
        gen = os.path.join(phispydir, "{genome}.gb.gz"),
        tbl = os.path.join(phageboostdir, "{genome}_phageboost_locs.tsv")
    output:
        tp = os.path.join(phageboostdir, "{genome}_phageboost_tptn.tsv")
    shell:
        """
        python3 ~/GitHubs/PhiSpy/scripts/compare_predictions_to_phages.py -t {input.gen} -r {input.tbl} > {output.tp}
        """
