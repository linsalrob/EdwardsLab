"""
Run phispy with VOGs on a set of IDs listed in a file.

We expect a file with the genome assembly ID to process, 
one per line

e.g.

snakemake -s phispy_vogs_download.snakefile --config filelist=20220130/needed.10 gbk=20220130/gbk output=20220130/phispy assembly=assembly_summary_20220130.txt  --profile sge

See phispy_vogs_download_submit.sh for a script to do this

We will download the sequence if required using either rsync or curl,
and then process it using phispy

NOTES:

    1. unalias rsync before proceeding. 
    2. best way to run this is on the cluster with ~100 genomes per file
       otherwise the DAG gets confusing for snakemake.
    3. Provide all the configs on the command line because then you can 
       code them in the slurm/sge script!

Configs we need:

    1. filelist: the list of accessions to process (aim for 100-1000 accessions)
    2. output: directory base where to write the files (e.g. phispy_20220130)
    3. gbk: directory where to store the genomes and/or look for them
    4. assembly: the latest assembly summary from genbank. We use this to get the URLs
    5. vogs: full path to the VOGdb hmms. e.g  /home3/redwards/VOGs/VOGs.hmm (we use version 99 because that was what was available in 2020)


"""


import os
import sys
import gzip


if 'filelist' not in config:
    sys.stderr.write("FATAL: Please provide a file with a list of genbank accessions, one per line\n")
    sys.stderr.write("as the config filelist. e.g. --config filelist=files.1\n")
    sys.exit(0)

if 'output' not in config:
    sys.stderr.write("FATAL: Please provide an output directory to write the data to\n")
    sys.exit(0)

if 'gbk' not in config:
    sys.stderr.write("FATAL: Please provide the location of the genbank files\n")
    sys.exit(0)

if 'assembly' not in config:
    sys.stderr.write("FATAL: Please provide the genbank assembly summary file\n")
    sys.exit(0)

if 'vogs' not in config:
    sys.stderr.write("FATAL: Please provide the location of the vogs. E.g. /home3/redwards/VOGs/VOGs.hmm\n")
    sys.exit(0)


def get_url(wildcards):
    return URLS[wildcards.sample]


def get_directories(samps):
    d = []
    for s in samps:
        d.append(os.path.join(s[0:9], s[0:13]))
    return d

def get_samples():
    samples = {}
    with open(config['filelist'], 'r') as f:
        for l in f:
            samples[l.strip()] = ""

    if config['assembly'].endswith('.gz'):
        f = gzip.open(config['assembly'], 'rt')
    else:
        f = open(config['assembly'], 'r')

    for l in f:
        if l.startswith("#"):
            continue
        p = l.strip().split("\t")
        # we strip off the protocol because we prefer rsync but fall
        # fall back to curl if that doesn't work
        if p[0] not in samples:
            continue
        if p[19] == "na":
            samples.pop(p[0])
            continue
        # remove the protocol
        p[19] = p[19].replace("https://", "", 1)
        p[19] = p[19].replace("http://", "", 1)
        p[19] = p[19].replace("ftp://", "", 1)
        ass_id = p[19].split("/")[-1]
        samples[p[0]] = f"{p[19]}/{ass_id}_genomic.gbff.gz"

    f.close()

    return list(samples.keys()), samples

SAMPLES, URLS = get_samples()
DIRS = get_directories(SAMPLES)

rule all:
    input:
        expand(
            [
                os.path.join(config['output'], "{directory}", "VOGS", "{sample}_VOGS_phispy.log.gz"),
                os.path.join(config['output'], "{directory}", "VOGS", "{sample}_VOGS_prophage_coordinates.tsv.gz"),
                os.path.join(config['output'], "{directory}", "VOGS", "{sample}_VOGS_phage.gbk.gz"),
            ], zip, sample=SAMPLES, directory=DIRS)
                
# removed this because we don't need              
# os.path.join(config['output'], "{directory}", "VOGS", "{sample}_VOGS_protein_functions.txt.gz")


rule download_genbank:
    """
    We start with rsync (it is the preferred protocol per 
    https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#protocols
    but that seems to fail quite frequently (perhaps gets swamped by
    multiple requests?) and so we fall back to curl/https.

    Note that Rob has an alias to rsync to use ssh by default, 
    so unaliasing that unsets it
    """
    output:
        gbf = os.path.join(config['gbk'], "{directory}", "{sample}_genomic.gbff.gz")
    resources:
        load_onehundred=10
    params:
        url = get_url,
        gbd = os.path.join(config['gbk'], "{directory}")
    shell:
        """
        set +e
        if alias rsync 2>/dev/null; then unalias rsync; fi
        mkdir --parents {params.gbd}
        rsync --copy-links --recursive --times --verbose rsync://{params.url} {output.gbf}
        exitcode=$?
        if [ $exitcode == 10 ]
        then
            curl -Lo {output} https://{params.url}
        fi
        """





rule run_phispy:
    input:
        os.path.join(config['gbk'], "{directory}", "{sample}_genomic.gbff.gz")
    output:
        temporary(os.path.join(config['output'], '{directory}', "VOGS", "{sample}_VOGS_bacteria.fasta")),
        temporary(os.path.join(config['output'], '{directory}', "VOGS", "{sample}_VOGS_phage.fasta")),
        os.path.join(config['output'], '{directory}',"VOGS", "{sample}_VOGS_phispy.log"),
        temporary(os.path.join(config['output'], '{directory}', "VOGS", "{sample}_VOGS_bacteria.gbk")),
        os.path.join(config['output'], '{directory}', "VOGS", "{sample}_VOGS_phage.gbk"),
        os.path.join(config['output'], '{directory}', "VOGS", "{sample}_VOGS_prophage_coordinates.tsv"),
    params:
        phispydir = os.path.join(config['output'], '{directory}', "VOGS"),
        sample = "{sample}_VOGS",
	vogs=config['vogs']
    shell:
        """
        set +e
        mkdir --parents phispydir;
        PhiSpy.py -o {params.phispydir} --quiet --output_choice 5 -p {params.sample} --phmms {params.vogs} {input}
        exitcode=$?
        if [ $exitcode -gt 19 ]
        then
            for o in {output}; do touch $o; done
        fi
        """


rule gzip_log:
    input:
        os.path.join(config['output'], '{directory}', "VOGS", "{sample}_VOGS_phispy.log")
    output:
        os.path.join(config['output'], '{directory}', "VOGS", "{sample}_VOGS_phispy.log.gz")
    shell:
        "gzip {input}"

rule gzip_gbk:
    input:
        gbk = os.path.join(config['output'], '{directory}', "VOGS", "{sample}_VOGS_phage.gbk"),
        pc = os.path.join(config['output'], '{directory}', "VOGS", "{sample}_VOGS_prophage_coordinates.tsv.gz"),
    output:
        os.path.join(config['output'], '{directory}', "VOGS", "{sample}_VOGS_phage.gbk.gz")
    shell:
        "gzip {input.gbk}"

rule gzip_tsv:
    input:
        os.path.join(config['output'], '{directory}', "VOGS", "{sample}_VOGS_prophage_coordinates.tsv")
    output:
        os.path.join(config['output'], '{directory}', "VOGS", "{sample}_VOGS_prophage_coordinates.tsv.gz")
    shell:
        "gzip {input}"
        
