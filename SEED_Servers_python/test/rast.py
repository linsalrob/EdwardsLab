import sys
import os
from servers.RAST import RASTserver
import getpass

username = getpass.getuser()
sys.stdout.write("Please enter your username [" + username + "]: ")
testun = sys.stdin.readline()
testun = testun.strip()
if (len(testun) > 0):
    username = testun

password = getpass.getpass()




rast = RASTserver(username, password)

status = rast.status_of_RAST_job({'-job' : 246902})
