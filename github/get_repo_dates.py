"""
Get some data information about a repo
"""

import os
import sys
import argparse
import re
import requests

__author__ = 'Rob Edwards'

github_user = None
github_token = None

if 'GHTOKEN' in os.environ:
    github_token = os.environ['GHTOKEN']
if 'GHUSER' in os.environ:
    github_user = os.environ['GHUSER']



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-r', '--repo', help='repository url')
    parser.add_argument('-f', '--file', help='file with one repo per line')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not github_user:
        print("Please specify a github user using the GHUSER environment variable", file=sys.stderr)
        sys.exit(1)
    if not github_token:
        print("Please specify a github token using the GHTOKEN environment variable", file=sys.stderr)
        sys.exit(1)

    headers = {
        'Authorization': f'token {github_token}'
    }

    if not args.repo and not args.file:
        args.repo = 'https://github.com/Vini2/phables'
        print(f"WARNING: Using {args.repo} as the repository", file=sys.stderr)

    repos = []

    if args.repo:
        repos.append(args.repo)
    if args.file:
        with open(args.file, 'r') as f:
            for l in f:
                repos.append(l.strip())

    for repo in repos:
        # https://api.github.com/repos/Vini2/phables/branches/develop
        # api = re.sub(r'^.*github.com', 'https://api.github.com', args.repo)
        if repo.endswith('/'):
            repo = repo[:-1]
        pieces = repo.split("/")
        username = pieces[-2]
        reponame = pieces[-1]
        api = f"https://api.github.com/repos/{username}/{reponame}"

        if args.verbose:
            print(f"Getting: {api}", file=sys.stderr)
        r = requests.get(api, headers=headers)
        if r.status_code != 200:
            print(f"WARNING: Status code: {r.status_code} for {api}", file=sys.stderr)
        else:
            js = r.json()
            if js:
                author = js['owner']['login']
                rname = js['name']
                cat = js['created_at']
                uat = js['updated_at']
                lk =  js['html_url']

                print("\t".join([repo, username, reponame, rname, author, cat, uat, lk]))
            else:
                print(f"ERROR: No json for {api}", file=sys.stderr)




