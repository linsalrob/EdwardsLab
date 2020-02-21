"""

A simple test

"""

VAR = "hello world"

rule all:
    input:
        "ls.txt"

rule lsd:
    input:
        "bac_giant_unique_species"
    output:
        "ls.txt"
    shell:
        #"/bin/ls {input} > {output}"
        "echo {VAR} > {output}"
