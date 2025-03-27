"""
If you open the outlook app you can export the emails to a .pst file. This script will extract the emails from that file

File -> Open & Export -> Import/Export -> Export to a file -> Outlook Data File (.pst) -> Select the folder -> save the file
"""

import os
import sys
import argparse

import pypff


def parse_message(message):
    print(f"{message.sender_name}\t{message.subject}\t{message.client_submit_time}")

def parse_folder(folder):
    print(f"Parsing folder: {folder.name}")
    for i in range(folder.number_of_sub_items):
        sub_item = folder.get_sub_item(i)
        if sub_item.is_message():
            parse_message(sub_item)
        elif sub_item.is_folder():
            parse_folder(sub_item)

def parse_pst(file_path):
    pst_file = pypff.file()
    pst_file.open(file_path)

    # Access the root folder
    root_folder = pst_file.get_root_folder()

    # Recursively parse the PST
    parse_folder(root_folder)
    pst_file.close()




__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='pst file exported from outlook', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    parse_pst(args.f)