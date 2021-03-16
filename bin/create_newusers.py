"""
Create a new users file.

To use this, run the code and then use the `newusers` command to add these to an AWS instance.

"""

from random import shuffle

abt = [x for x in 'ABCDEFGHKMNPQRSTWXYZabcdefghkmnpqrstwxyz23456789@#$%^*']

with open("users.tsv", 'w') as out:
    for i in range(1, 76):
        server = ( i % 6) + 1
        shuffle(abt)
        pwd = "".join(abt[0:8])
        if i < 10:
            out.write(f"{server}\tsagc0{i}:{pwd}:::Rob Added User:/home/sagc0{i}:/bin/bash\n")
        else:
            out.write(f"{server}\tsagc{i}:{pwd}:::Rob Added User:/home/sagc{i}:/bin/bash\n")
