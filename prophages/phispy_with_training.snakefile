

import sys



def get_trainingset(wildcards):
    ts = {    
        "Achromobacter_denitrificans_strain_PR1": "data/trainSet_32002.17.txt",
        "Bacillus_halodurans_C-125": "data/trainSet_272558.23.txt",
        "Bacillus_subtilis_subsp._subtilis_str._168": "data/trainSet_224308.360.txt",
        "Bacteroides_uniformis_ATCC_8492_strain_81A2": "data/trainSet_411479.31.txt",
        "Bifidobacterium_longum_NCC2705": "data/trainSet_206672.37.txt",
        "Brucella_melitensis_16M": "data/trainSet_224914.79.txt",
        "Caulobacter_crescentus_CB15": "data/trainSet_190650.21.txt",
        "Clostridium_perfringens_str._13": "data/trainSet_195102.53.txt",
        "Clostridium_tetani_E88": "data/trainSet_212717.31.txt",
        "Deinococcus_radiodurans_R1": "data/trainSet_243230.96.txt",
        "Enterococcus_faecalis_strain_V583": "data/trainSet_1351.557.txt",
        "Enterococcus_faecalis_V583": "data/trainSet_226185.9.txt",
        "Escherichia_coli_CFT073": "data/trainSet_199310.168.txt",
        "Escherichia_coli_K12": "data/trainSet_83333.998.txt",
        "Escherichia_coli_O157-H7_EDL933": "data/trainSet_155864.289.txt",
        "Escherichia_coli_O157-H7": "data/trainSet_83334.295.txt",
        "Haemophilus_influenzae_Rd_KW20": "data/trainSet_71421.45.txt",
        "Lactococcus_lactis_subsp._lactis_Il1403": "data/trainSet_272623.42.txt",
        "Listeria_innocua_Clip11262": "data/trainSet_272626.22.txt",
        "Listeria_monocytogenes_EGD-e": "data/trainSet_169963.176.txt",
        "Mesorhizobium_loti_MAFF303099": "data/trainSet_266835.41.txt",
        "Mycobacterium_tuberculosis_CDC1551": "data/trainSet_83331.121.txt",
        "Mycobacterium_tuberculosis_H37Rv": "data/trainSet_83332.460.txt",
        "Neisseria_meningitidis_MC58": "data/trainSet_122586.26.txt",
        "Neisseria_meningitidis_Z2491": "data/trainSet_122587.18.txt",
        "Paracoccus_aminophilus_JCM_7686": "data/trainSet_1367847.3.txt",
        "Paracoccus_denitrificans_PD1222": "data/trainSet_318586.5.txt",
        "Paracoccus_sanguinis_5503": "data/trainSet_1525717.3.txt",
        "Paracoccus_sp._SCN_68-21": "data/trainSet_1660154.3.txt",
        "Paracoccus_yeei_TT13": "data/trainSet_147645.106.txt",
        "Pasteurella_multocida_subsp._multocida_str._Pm70": "data/trainSet_272843.53.txt",
        "Pseudomonas_aeruginosa_PAO1": "data/trainSet_208964.452.txt",
        "Pseudomonas_putida_KT2440": "data/trainSet_160488.79.txt",
        "Ralstonia_solanacearum_GMI1000": "data/trainSet_267608.42.txt",
        "Salmonella_enterica_subsp._enterica_serovar_Typhi_str._CT18": "data/trainSet_220341.87.txt",
        "Shewanella_oneidensis_MR-1": "data/trainSet_211586.69.txt",
        "Shigella_flexneri_2a_str._301": "data/trainSet_198214.txt",
        "Staphylococcus_aureus_strain_Sa_Newman_UoM": "data/trainSet_1280.10152.txt",
        "Staphylococcus_aureus_subsp._aureus_Mu50": "data/trainSet_158878.38.txt",
        "Staphylococcus_aureus_subsp._aureus_MW2": "data/trainSet_196620.15.txt",
        "Streptococcus_pyogenes_M1_GAS": "data/trainSet_160490.61.txt",
        "Streptococcus_pyogenes_MGAS315": "data/trainSet_198466.10.txt",
        "Streptococcus_pyogenes_MGAS8232": "data/trainSet_186103.26.txt",
        "Vibrio_cholerae_O1_biovar_eltor_str._N16961": "data/trainSet_243277.252.txt",
        "Xanthomonas_axonopodis_pv._citri_str._306": "data/trainSet_190486.46.txt",
        "Xylella_fastidiosa_9a5c": "data/trainSet_160492.65.txt",
        "Xylella_fastidiosa_Temecula1": "data/trainSet_183190.38.txt",
        "Yersinia_pestis_CO92": "data/trainSet_214092.200.txt",
        "Yersinia_pestis_KIM": "data/trainSet_187410.24.txt",
    }
    if wildcards.genome in ts:
        return ts[wildcards.genome]
    else:
        sys.stderr.write(f"FATAL: No training set for {wildcards.genome}\n")
        sys.exit(-2)






outputdir = "phispy_with_training"

phispydir = "/home3/redwards/GitHubs/PhiSpy/test_genbank_files"
phispydata = "/home3/redwards/GitHubs/PhiSpy/PhiSpyModules"
GENOMES, = glob_wildcards(os.path.join(phispydir, '{genome}.gb.gz'))

rule all:
    input:
        expand(os.path.join(outputdir, "{genome}_phispy_tptn.tsv"), genome=GENOMES)

rule run_phispy:
    input:
        g = os.path.join(phispydir, "{genome}.gb.gz")
    params:
        t = get_trainingset,
        o = os.path.join(outputdir, "{genome}.phispy")
    benchmark:
        os.path.join(outputdir, "benchmarks", "{genome}.benchmarks.txt")
    output:
        temporary(os.path.join(outputdir, "{genome}.phispy", "bacteria.fasta")),
        temporary(os.path.join(outputdir, "{genome}.phispy", "bacteria.gbk")),
        temporary(os.path.join(outputdir, "{genome}.phispy", "phage.fasta")),
        os.path.join(outputdir, "{genome}.phispy", "phage.gbk"),
        os.path.join(outputdir, "{genome}.phispy", "phispy.log"),
    shell:
        """
        PhiSpy.py -t {params.t} -o {params.o} --output_choice 4 {input.g}
        """

rule count_tp_tn:
    input:
        gen = os.path.join(phispydir, "{genome}.gb.gz"),
        phg = os.path.join(outputdir, "{genome}.phispy", "phage.gbk")
    output:
        tp = os.path.join(outputdir, "{genome}_phispy_tptn.tsv")
    shell:
        """
        python3 ~/GitHubs/PhiSpy/scripts/compare_predictions_to_phages.py -t {input.gen} -p {input.phg} > {output.tp}
        """
