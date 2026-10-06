"""
Convert date in column 2 to day of the week.
"""

import os
import sys
import argparse
from roblib import bcolors
from datetime import datetime
import pytz
__author__ = 'Rob Edwards'



def process_date(date_str):
    # Parse the date string into a datetime object
    try:
        utc_time = datetime.strptime(date_str, "%Y-%m-%d %H:%M:%S")
    except ValueError as e:
        print(f"Error parsing date {date_str}: {e}", file=sys.stderr)
        return
    
    # Ensure the date is in 2024
    if utc_time.year != 2024:
        print(f"The date {date_str} is not in 2024.", file=sys.stderr)
        return
    
    # Define UTC and Adelaide timezones
    utc_tz = pytz.utc
    adelaide_tz = pytz.timezone("Australia/Adelaide")
    
    # Localize the UTC time and convert to Adelaide time
    utc_time = utc_tz.localize(utc_time)
    adelaide_time = utc_time.astimezone(adelaide_tz)
    
    # Print the day of the week
    print(adelaide_time.strftime("%A"))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='input file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    with open(args.f, 'r') as f:
        for line in f:
            process_date(line.strip().split("\t")[2])

