from servers.servers import server

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


    def __init__(self):
        server.__init__(self)
        self.urlservice = "sapling"


