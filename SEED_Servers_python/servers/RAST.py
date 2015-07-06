from .servers import server

class RASTserver(server):
    '''
    The RASTserver connects to RAST and will allow you to submit or retrieve
    RAST jobs.

    To make a connection, instantiate the SAPserver:
        rast_server - servers.RASTserver()

    All functions take a dictionary object, that contains the
    appropriate query information you are looking for. 

    For example, to retrieve a specific job, you can use:
        data = [ '-job': ['1234', '5678'] ]
        tarball = rast_server.retrieve_RAST_job(data)

    or

        tarball = rast_server.retrieve_RAST_job( [ '-data' : [''1234', '5678'] ]

    The objects returned are also almost always dict objects, so you can
    retrieve the appropriate information.

    By using dict objects you can send and retrieve multuple queries at 
    once.

    For more detailed discussion of both the arguments and the returned
    data, see the sapling documention: http://servers.nmpdr.org/sapling/server.cgi?pod=RAST

    You can add arbitrary function calls as they are added to the RAST
    servers. For examplle if a function is listed on the documentation
    website and is not available directly, then you can call it using
    the function:
        rast_server.function('newfunction')
        rast_server.retrieveData(dataObj)


    '''


    def __init__(self, username, password):
        server.__init__(self)
        self.urlbase = 'http://servers.nmpdr.org/'
        self.urlservice = "rast"
        self.urlcgi = '/server.cgi'
        self.params(username = username, password = password)



    def submit_RAST_job(self,data):
        '''
        This method submits a job to the RAST server. It returns a hash of:
        '''

        data['username']=self.username
        data['password']=self.password

        self.function('submit_RAST_job')
        return self.retrieve(data)

    def status_of_RAST_job(self,data):
        '''
        Where data is a list of jobs
        The return value is a hash keyed by Jobid of   
        '''

        self.function('status_of_RAST_job')
        return self.retrieve(data)

    def retrieve_RAST_job(self,data):
        '''
        where $jobid is the RAST id of the job and
        $format is  one of:
        The return is a hash of 
        '''

        self.function('retrieve_RAST_job')
        return self.retrieve(data)

    def kill_RAST_job(self,data):
        '''
        where @jobids is an array of RAST job ids to kill.
        Return is a hash keyed by Job ID of 
        '''

        self.function('kill_RAST_job')
        return self.retrieve(data)

    def delete_RAST_job(self,data):
        '''
        where @jobids is an array of RAST job ids to kill.
        Return is a hash keyed by Job ID of 
        '''

        self.function('delete_RAST_job')
        return self.retrieve(data)

    def get_job_metadata(self,data):
        '''
        where $jobid is the RAST id  of a RAST job
        Return is a hash of 
        '''

        self.function('get_job_metadata')
        return self.retrieve(data)

    def copy_to_RAST_dir(self,data):
        '''
        where $jobid is the RAST id of the job,
        The return is a hash of 
        '''

        self.function('copy_to_RAST_dir')
        return self.retrieve(data)

