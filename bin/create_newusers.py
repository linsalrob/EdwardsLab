"""
Create a new users file.

To use this, run the code and then use the `newusers` command to add these to an AWS instance.

"""
import argparse
from random import shuffle

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-b', '--base', help='base of the username (default="user")', default="user")
    parser.add_argument('-n', '--number', help='number of accounts to create. (default=100)', type=int, default=100)
    parser.add_argument('-s', '--servers', help='number of servers to share users among', default=6)
    parser.add_argument('-u', '--users', help='users file for bash to write (default="users.tsv")')
    parser.add_argument('-a', '--accounts', help='accounts file to give to the users (default="accounts.tsv")')
    args = parser.parse_args()


    abt = [x for x in 'ABCDEFGHKMNPQRSTWXYZabcdefghkmnpqrstwxyz23456789@#$%^*']

    with open(args.users, 'w') as users, open(args.accounts, 'w') as accounts:
        for i in range(1, args.number):
            server = ( i % args.servers) + 1
            shuffle(abt)
            pwd = "".join(abt[0:8])
            print(f"{server}\t{args.base}{i:0>4}\t{pwd}", file=accounts)
            print(f"{args.base}{i:0>4}:{pwd}:::Rob Added User:/home/{args.base}{i:0>4}:/bin/bash", file=users)
