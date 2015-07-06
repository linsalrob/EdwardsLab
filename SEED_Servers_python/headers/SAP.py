

from servers.servers import server

class SAPserver(server):
    '''
    The SAPserver connects to the SAPLING database that contains all genome information. You can use this class to retrieve data from the SEED database. 

    To make a connection, instantiate the SAPserver:
        sapServer = servers.SAPserver()

    All functions take a dictionary object, that contains the appropriate query information you are looking for. 

    For example, to find the sequence of a specific protein you can call the idsToSequences function with the id of the protein:
        data = { '-id':'fig|83333.1.peg.1234' }
        seqObject = sapServer.idsToSequences(data)

    or
        seqObject = sapServer.idsToSequences({ '-id':'fig|83333.1.peg.1234' })

    The objects returned are also almost always dict objects, so you can retrieve the appropriate information.

    By using dict objects you can send and retrieve multuple queries at once, such as retrieving all the protein sequences for all the ids in a genome.

    For more detailed discussion of both the arguments and the returned data, see the sapling documention: http://servers.nmpdr.org/sapling/server.cgi?pod=SAP

    You can add arbitrary function calls as they are added to the SAP servers. For examplle if a function is listed on the documentation website and is not available directly, then you can call it using the function:
        sapServer.function('newfunction')
        sapServer.retrieveData(dataObj)


    '''


    def __init__(self):
        server.__init__(self)
        self.urlservice = "sapling"


