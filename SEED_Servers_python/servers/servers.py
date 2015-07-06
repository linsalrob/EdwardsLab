'''
A python implementation of the Annotation server, Anno server. 

Written by Rob Edwards, Daniel Cuevas, and Genivaldo Gueiros


'''

#Note: this requires requests. See http://www.python-requests.org
# requests is now often installed directly

import requests
import sys
import json

class server:
    def __init__(self):
        #self.urlbase = 'http://servers.nmpdr.org/'
        self.urlbase = 'http://pubseed.theseed.org/'
        self.urlservice = None
        self.urlcgi  = '/server.cgi'
        self.params(email = 'redwards@mcs.anl.gov', source = 'Robs python implementation', encoding='json')

    def service(self, name):
        self.urlservice = name

    def params(self, function=None, email=None,  source=None, username=None, password=None, encoding=None):
        self.data = {}
        self.data['encoding'] = encoding

        # construct the query
        if email:
            self.data['email']=email

        if source:
            self.data['source']=source

        if username:
            self.data['username']=username

        if password:
            self.data['password']=password

    def function(self):
        '''Get the current function'''
        return self.data['function']

    def function(self, function=None):
        '''Define the function to call'''
        self.data['function']=function

    def getData(self):
        return self.data

    def addData(self, data):
        '''add more data to the call. Note that this does not delete the existing data (in case you define things separately'''
        for x in data:
            self.data[x]=data[x] 

    def clearData(self):
        '''delete all the data'''
        self.data = {}


    def retrieve(self, data):
        '''make the call and get the response'''

        if not self.urlservice:
            sys.stderr.write("You did not define a service to call. Please use one of the servers like SAPserver, ANNOserver, FBAMODELserver, etc\n")
            sys.exit(-1)

        url = self.urlbase + self.urlservice + self.urlcgi
        sys.stderr.write("Connecting to " + url + " for function " + self.data['function'] + "\n")

        dataToProcess = self.data.copy()
        ## convert the data to unicode
        if data == None:
            data = {"":""}

        dataToProcess['args'] = json.dumps(data)
        # print("TRYING: " + str(data))

        #print("ARGS: ||>>", dataToProcess['args'], "<<||\n")
        # sys.stderr.write("URL: " + url + "\nData: " + str(dataToProcess) + "\n");
        resp = requests.post(url, data=dataToProcess)

        if resp is None:
            sys.stderr.write("A connection to " + url + " could not be established. Please check your internets and try again\n")
            sys.exit(-1)

        jsonData = None
        try:
            jsonData = resp.json()
        except(ValueError):
            sys.stderr.write("\nERROR: No json object could be detected\nContent:\n")
            sys.stderr.write(resp.text + "\n\n")

        return jsonData

