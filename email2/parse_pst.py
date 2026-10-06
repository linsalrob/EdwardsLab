"""
Parse a .pst file downloaded from outlook
"""

import os
import sys
import argparse
from roblib import bcolors
__author__ = 'Rob Edwards'

import pypff

def parse_message(message):
    print(f"{message.sender_name}\t{message.subject}\t{message.client_submit_time}")
    
# Recursively parse the PST
def parse_folder(folder):
    print(f"Parsing folder: {folder.name}", file=sys.stderr)
    for i in range(folder.number_of_sub_items):
        sub_item = folder.get_sub_item(i)
        if isinstance(sub_item, pypff.message):
            parse_message(sub_item)
        elif isinstance(sub_item, pypff.folder):
            parse_folder(sub_item)

def parse_pst(file_path):
    pst_file = pypff.file()
    pst_file.open(file_path)
    
    # Access the root folder
    root_folder = pst_file.get_root_folder()
    
    
    parse_folder(root_folder)
    pst_file.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-p', help='pst outlook output file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if args.v:
        print(f"{bcolors.GREEN}Parsing {args.p}{bcolors.ENDC}")

    if not os.path.exists(args.p):
        print(f"{bcolors.RED}ERROR: {args.p} does not exist{bcolors.ENDC}")
        sys.exit(1)

    parse_pst(args.p)


