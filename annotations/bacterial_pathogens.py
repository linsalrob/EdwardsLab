"""

"""

import os
import sys
import argparse

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='input file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    pathogs = {
        "Genera": [
            {
                "genus": "Streptococcus",
                "family": "Streptococcaceae",
                "phylum": "Firmicutes"
            },
            {
                "genus": "Staphylococcus",
                "family": "Staphylococcaceae",
                "phylum": "Firmicutes"
            },
            {
                "genus": "Haemophilus",
                "family": "Pasteurellaceae",
                "phylum": "Proteobacteria"
            },
            {
                "genus": "Mycobacterium",
                "family": "Mycobacteriaceae",
                "phylum": "Actinobacteria"
            },
            {
                "genus": "Pseudomonas",
                "family": "Pseudomonadaceae",
                "phylum": "Proteobacteria"
            },
            {
                "genus": "Klebsiella",
                "family": "Enterobacteriaceae",
                "phylum": "Proteobacteria"
            },
            {
                "genus": "",
                "family": "Moraxellaceae",
                "phylum": "Proteobacteria"
            },
            {
                "genus": "Bordetella",
                "family": "Alcaligenaceae",
                "phylum": "Proteobacteria"
            },
            {
                "genus": "Legionella",
                "family": "Legionellaceae",
                "phylum": "Proteobacteria"
            },
            {
                "genus": "Corynebacterium",
                "family": "Corynebacteriaceae",
                "phylum": "Actinobacteria"
            },
            {
                "genus": "Chlamydia",
                "family": "Chlamydiaceae",
                "phylum": "Chlamydiota"
            },
            {
                "genus": "Mycoplasma",
                "family": "Mycoplasmataceae",
                "phylum": "Tenericutes (Mollicutes)"
            },
            {
                "genus": "Neisseria",
                "family": "Neisseriaceae",
                "phylum": "Proteobacteria"
            },
            {
                "genus": "Burkholderia",
                "family": "Burkholderiaceae",
                "phylum": "Proteobacteria"
            },
            {
                "genus": "Acinetobacter",
                "family": "Moraxellaceae",
                "phylum": "Proteobacteria"
            },
            {
                "genus": "Francisella",
                "family": "Francisellaceae",
                "phylum": "Proteobacteria"
            },
            {
                "genus": "Escherichia",
                "family": "Enterobacteriaceae",
                "phylum": "Proteobacteria"
            },
            {
                "genus": "Pasteurella",
                "family": "Pasteurellaceae",
                "phylum": "Proteobacteria"
            },
            {
                "genus": "Nocardia",
                "family": "Nocardiaceae",
                "phylum": "Actinobacteria"
            },
            {
                "genus": "Actinomyces",
                "family": "Actinomycetaceae",
                "phylum": "Actinobacteria"
            },
            {
                "family": "Bacillaceae",
                "notes": "Includes Bacillus anthracis, which can cause inhalational anthrax."
            },
            {
                "family": "Enterobacteriaceae",
                "notes": "Beyond Klebsiella and Escherichia, genera like Enterobacter and Serratia can cause respiratory infections."
            },
            {
                "family": "Rickettsiaceae",
                "notes": "Primarily vector-borne, but can have pulmonary involvement in immunocompromised hosts."
            },
            {
                "family": "Brucellaceae",
                "notes": "Brucella species can sometimes involve the respiratory system through inhalation."
            }
        ]
    }
