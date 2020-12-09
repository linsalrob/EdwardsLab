"""
Snakefile to process all of Liz's metagenomes
"""

# Name of our shark samples. These should be the same as the names of the folders
SHARKS = ['Angelshark', 'Batray', 'Butterflyray', 'Duskyshark',
          'Galapagosshark', 'Leopardshark', 'Roundray', 'Thornbackray',
          'Tigershark', 'Whaleshark']

# READDIR = config['Paths']['Reads']
READDIR = "fastq"
PSEQDIR = "prinseq"
SUPFDIR = "superfocus"
FOCSDIR = "focus"


SAMPLES,EXTENSIONS = glob_wildcards(os.path.join(READDIR, '{sample}_R1{extentions}'))
# we just get the generic extension. This is changed in Step 1
file_extension = EXTENSIONS[0]

# a convenience so we don't need to use '{sample}_R1' all the time
PATTERN_R1 = '{sample}_R1'
PATTERN_R2 = '{sample}_R2'



rule all:
    input:
        expand(os.path.join(PSEQDIR, "{sample}_good_out_R1.fastq"), sample=SAMPLES),
        os.path.join(FOCSDIR, "output_All_levels.csv"),
        os.path.join(SUPFDIR, "all_taxonomy.tsv"),
        os.path.join(SUPFDIR, "output_all_levels_and_function.xls"),
        "statistics.tsv"


rule prinseq:
    input:
        r1 = os.path.join(READDIR, "{sample}_R1" + file_extension),
        r2 = os.path.join(READDIR, "{sample}_R2" + file_extension)
    output:
        r1 = os.path.join(PSEQDIR, "{sample}_good_out_R1.fastq"),
        r2 = temporary(os.path.join(PSEQDIR, "{sample}_good_out_R2.fastq")),
        s1 = temporary(os.path.join(PSEQDIR, "{sample}_single_out_R1.fastq")),
        s2 = temporary(os.path.join(PSEQDIR, "{sample}_single_out_R2.fastq")),
        b1 = temporary(os.path.join(PSEQDIR, "{sample}_bad_out_R1.fastq")),
        b2 = temporary(os.path.join(PSEQDIR, "{sample}_bad_out_R2.fastq"))
    conda: "liz_metagenomes_env.yaml"
    params:
        o = os.path.join(PSEQDIR, "{sample}")
    shell:
        """
            prinseq++ -min_len 60 -min_qual_mean 25 -ns_max_n 1 -derep 1 \
                    -out_format 0 -trim_tail_left 5 -trim_tail_right 5 \
                    -ns_max_n 5  -trim_qual_type min -trim_qual_left 30 \
                    -trim_qual_right 30 -trim_qual_window 10 \
                    -threads {threads} \
                    -out_name {params.o} \
                    -fastq {input.r1} \
                    -fastq2 {input.r2};
        """

rule stats:
    input:
        expand(os.path.join(PSEQDIR, "{smps}_good_out_R1.fastq"), smps=SAMPLES)
    output:
        "statistics.tsv"
    conda: "liz_metagenomes_env.yaml"
    params:
        i = PSEQDIR,
    shell:
        """
        echo -e "name\tnum sequences\ttotal len\tshortest\tlongest\tn50\tn75" > {output} &&
        countfastq.py -t -d {PSEQDIR} >> {output}
        """


rule focus:
    input:
        expand(os.path.join(PSEQDIR, "{smps}_good_out_R1.fastq"), smps=SAMPLES)
    output:
        os.path.join(FOCSDIR, "output_All_levels.csv"),
        temporary(os.path.join(FOCSDIR, "output_Family_tabular.csv")),
        temporary(os.path.join(FOCSDIR, "output_Kingdom_tabular.csv")),
        temporary(os.path.join(FOCSDIR, "output_Phylum_tabular.csv")),
        temporary(os.path.join(FOCSDIR, "output_Strain_tabular.csv")),
        temporary(os.path.join(FOCSDIR, "output_Class_tabular.csv")),
        temporary(os.path.join(FOCSDIR, "output_Genus_tabular.csv")),
        temporary(os.path.join(FOCSDIR, "output_Order_tabular.csv")),
        temporary(os.path.join(FOCSDIR, "output_Species_tabular.csv"))
    conda: "liz_metagenomes_env.yaml"
    params:
        i = PSEQDIR,
        o = FOCSDIR
    shell:
        """
        focus -q {params.i} -o {params.o} -t {threads}
        """



rule superfocus:
    input:
        expand(os.path.join(PSEQDIR, "{smps}_good_out_R1.fastq"), smps=SAMPLES)
    output:
        os.path.join(SUPFDIR, "output_all_levels_and_function.xls"),
        expand(os.path.join(SUPFDIR, "{smps}_good_out_R1.fastq_alignments.m8"), smps=SAMPLES),
        temporary(os.path.join(SUPFDIR, "output_subsystem_level_1.xls")),
        temporary(os.path.join(SUPFDIR, "output_subsystem_level_2.xls")),
        temporary(os.path.join(SUPFDIR, "output_subsystem_level_3.xls")),
    conda: "liz_metagenomes_env.yaml"
    params:
        i = PSEQDIR,
        o = SUPFDIR
    shell:
        """
        superfocus -q {params.i} -dir {params.o} -a diamond -t 12
        """

# NOTE: This uses SQL, so we need to move the database to a local database
#
rule superfocus_taxonomy:
    input:
        os.path.join(SUPFDIR, "{smps}_good_out_R1.fastq_alignments.m8"),
    output:
        os.path.join(SUPFDIR, "{smps}_good_out_R1.taxonomy")
    shell:
        """
        python3 ~/GitHubs/EdwardsLab/superfocus_all/superfocus_to_taxonomy.py -f {input} --tophit > {output}
        """

rule count_sf_taxonomy:
    input:
        os.path.join(SUPFDIR, "{smps}_good_out_R1.taxonomy")
    output:
        os.path.join(SUPFDIR, "{smps}_tax_counts.tsv")
    params:
        s = "{smps}"
    shell:
        """
        cut -d$'\t' -f 2- {input} | sed -e 's/\t/|/g' | \
        awk -v F={params.s} 'BEGIN {{print "Superkingdom|Phylum|Class|Order|Family|Genus|Species|Taxid\t"F}} s[$0]++ {{}} END {{ for (i in s) print i"\t"s[i] }}' \
        > {output}
        """

rule join_superfocus_taxonomy:
    input:
        expand(os.path.join(SUPFDIR, "{smps}_tax_counts.tsv"), smps=SAMPLES)
    output:
        os.path.join(SUPFDIR, "all_taxonomy.tsv")
    shell:
        """
        joinlists.pl -h -z {input} > {output}
        """

