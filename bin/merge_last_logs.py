"""
Merge output from multiple lastlog outputs to give a single last log
"""

import os
import sys

import argparse
import re
import dateutil.parser

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-l', help='lastlog file(s). You can provide multiple -l at once', action='append', required=True)
    args = parser.parse_args()

    # the string is white space separated and looks like
    # redwards         pts/0    anthill          Thu Feb 11 16:39:56 -0800 2016
    # redwards         pts/42   rohan.sdsu.edu   Wed Apr 13 17:14:03 -0700 2016


    special_users = ['named', 'webalizer', 'nut', 'squid', 'sync', 'netdump', 'rpc', 'smmsp', 'shutdown', 'mysql',
                     'operator', 'rpm', 'gdm', 'postgres', 'ftp', 'uucp', 'mailman', 'mailnull', 'distcache', 'pvm',
                     'hsqldb', 'beaglidx', 'lp', 'ldap', 'mail', 'ganglia', 'adm', 'bin', 'prime', 'fax', 'nslcd',
                     'halt', 'abrt', 'ntp', 'nobody', 'sshd', 'gopher', 'privoxy', 'pegasus', 'dovecot', 'apache',
                     'news', 'xfs', 'pcap', 'haldaemon', 'biobike3', 'biobike2', 'daemon', 'ident', 'tomcat', 'rpcuser',
                     'postfix', 'dbus', 'biobike4', 'ais', 'nfsnobody', 'games', 'vcsa', 'nscd', 'avahi', 'saslauth',
                     'fig', 'jottotoo', 'root', 'jotto', 'rtkit', 'biobike', 'pulse', 'tcpdump', 'avahi-autoipd',
                     'dockerroot']
    special_users.sort()
    special_users_set = set(special_users)


    access = {}
    log = {}
    never = {}

    for f in args.l:
        with open(f, 'r') as fin:
            for l in fin:
                p=l.strip().split()
                if p[0] in special_users_set:
                    continue
                if 'Never logged in' in l:
                    never[p[0]] = l
                    continue
                m = re.match('(\S+)\s+(\S+)\s+(\S+)\s+(.*?)$', l)
                if not m:
                    sys.stderr.write("Can't parse: '{}'\n".format(l))
                    continue
                if m.group(1) not in access:
                    try:
                        access[m.group(1)] = dateutil.parser.parse(m.group(4))
                        log[m.group(1)] = l
                    except:
                        sys.stderr.write("Can't parse date from '{}' in '{}'\n".format(m.group(4), l.strip()))
                elif dateutil.parser.parse(m.group(4)) > access[m.group(1)]:
                    access[m.group(1)] = dateutil.parser.parse(m.group(4))
                    log[m.group(1)] = l

# delete the never logged in if they did
sys.stderr.write("{} never logged in\n".format(len(never.keys())))
for k in log:
    if k in never:
        never.pop(k)
sys.stderr.write("{} never logged in\n".format(len(never.keys())))

for s in special_users:
    sys.stdout.write("{}\tspecial user. Do not deactivate\n".format(s))
sys.stdout.write("\n\nNEVER LOGGED IN\n")

for n in never.values():
    sys.stdout.write(n)

sys.stdout.write("\n\nLOGGED IN USERS\n")

ke = sorted(log, key=access.__getitem__)
for l in ke:
    sys.stdout.write(log[l])
