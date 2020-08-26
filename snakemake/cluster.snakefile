import os

# set this to the directory that has the fastq files
readdir = "step_5"

# set this to the path to the indexed host genome
# that you want to filter out. This should be indexed
# with `bowtie2-build`, and should be just the base 
# filename of the .1.bt2, .2.bt2, .3.bt2 indexes etc.
host_bt_index = "/home3/redwards/IBD/CERVAID/databases_20200604/databases/human_masked/human_virus_masked"

# set this to the host name that we will use in the output names
hostname = "human"


# set this to the location where you would like the results
# written. We will make the directory if it doesn't exist
outdir = "step_6"



SAMPLES,EXTENSIONS = glob_wildcards(os.path.join(readdir, '{sample}_R1{extentions}'))

if not EXTENSIONS:
    sys.stderr.write("""
        FATAL: We could not parse the sequence file names.
        We are expecting {sample}_R1{extension}, and so your files
        should contain the characters '_R1' in the fwd reads
        and '_R2' in the rev reads
        """)
    sys.exit()
# we just get the generic extension. This is changed in Step 1

file_extension = EXTENSIONS[0]
# a convenience so we don't need to use '{sample}_R1' all the time
PATTERN_R1 = '{sample}_R1'
PATTERN_R2 = '{sample}_R2'



clustdir = "clusters"


rule all:
    input:
        expand([
        os.path.join(clustdir, '{sample}_R1.cdhit'),
        os.path.join(clustdir, "{sample}_R1.best.fasta")
        ], sample=SAMPLES)



rule fq2fa:
    input:
        r1 = os.path.join(readdir, '{sample}_R1' + file_extension)
    output:
        r1 = os.path.join(clustdir, '{sample}_R1.fasta')
    shell:
        """
        seqtk seq -A {input.r1} > {output.r1}
        """

rule cdhit:
    input:
        r1 = os.path.join(clustdir, '{sample}_R1.fasta')
    output:
        i = os.path.join(clustdir, '{sample}_R1.cdhit'),
        c = os.path.join(clustdir, '{sample}_R1.cdhit.clstr')
    benchmark:
        "benchmark/cdhit.{sample}.txt"
    shell:
        """
        ~/anaconda3/bin/cd-hit -d 0 -T 0 -M 0 -i {input.r1} -o {output.i}
        """


rule deduplicate:
    """
    Step 10: Dereplicate
    """
    input:
        r1 = os.path.join(readdir, '{sample}_R1' + file_extension)
    output:
        fa = os.path.join(clustdir, "{sample}_R1.best.fasta"),
        stats = os.path.join(clustdir, "{sample}_R1.stats.txt")
    benchmark:
        "benchmark/dedupe.{sample}.txt"
    shell:
        """
        dedupe.sh in={input} \
            csf={output.stats} out={output.fa} \
            ow=t s=4 rnc=t pbr=t 
        """

        
