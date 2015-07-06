'''
A python implementation of the Annotation server, Anno server. 

Written by Rob Edwards, Daniel Cuevas, and Genivaldo Gueiros

'''

import urllib
import json

class server:
    def __init__(self, url = 'http://servers.nmpdr.org/anno/server.cgi', email = 'redwards@mcs.anl.gov', source = 'Robs python implementation', username=None, password=None, encoding='json', method=None, data={}):
        self.url = url
        self.method = method
        self.data = data
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
        


    def method(self):
        '''Get the current method'''
        return self.data['method']

    def method(self, method=None):
        '''Define the method to call'''
        self.data['method']=method

    def addData(self, data):
        '''add more data to the call. Note that this does not delete the existing data (in case you define things separately'''
        for x in data:
            self.data[x]=data[x] 

    def clearData(self):
        '''delete all the data'''
        self.data = {}


    def query():
        '''make the call and get the response'''
        query = urllib.parse.urlencode(data).encode('utf-8')
        resp = None
        try:
            resp = urllib.request.urlopen(self.url, query)
        except urllib.error.URLError as e:
            print(e.reason)
        except urllib.error.HTTPError as e:
            print(e.code)
            print(e.read())

        data = json.load(resp.read())
        return data








