"""
Get some data information about a repo
"""

import os
import sys
import argparse
import re
import requests

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-r', '--repo', help='repository url')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not args.repo:
        args.repo = 'https://github.com/Vini2/phables'
        print(f"WARNING: Using {args.repo} as the repository", file=sys.stderr)

    # https://api.github.com/repos/Vini2/phables/branches/develop

    # api = re.sub(r'^.*github.com', 'https://api.github.com', args.repo)
    if args.repo.endswith('/'):
        args.repo = args.repo[:-1]
    pieces = args.repo.split("/")
    username = pieces[-2]
    reponame = pieces[-1]
    api = f"https://api.github.com/repos/{username}/{reponame}/branches/master"

    print(f"Getting: {api}", file=sys.stderr)
    r = requests.get(api)
    if r.status_code != 200:
        print(f"WARNING: Status code: {r.status_code}", file=sys.stderr)
        sys.exit(r.status_code)

    js = r.json()
    # print(js)
    author = js['commit']['author']['login']
    branch = js['name']
    dt = js['commit']['commit']['author']['date']
    lk =  js['_links']['html']

    print("\t".join([username, reponame, author, branch, dt, lk]))




