from .servers import server

class ANNOserver(server):
    '''
    The ANNOserver connects to the annotation server and will allow you
    to make annotations.

    To make a connection, instantiate the SAPserver:
        anno_server - servers.ANNOserver()

    All functions take a dictionary object, that contains the
    appropriate query information you are looking for. 

    The objects returned are also almost always dict objects, so you can
    retrieve the appropriate information.

    By using dict objects you can send and retrieve multuple queries at 
    once.

    For more detailed discussion of both the arguments and the returned
    data, see the sapling documention: http://servers.nmpdr.org/sapling/server.cgi?pod=RAST

    You can add arbitrary function calls as they are added to the ANNO
    servers. For examplle if a function is listed on the documentation
    website and is not available directly, then you can call it using
    the function:
        anno_server.function('newfunction')
        anno_server.retrieveData(dataObj)


    '''


    def __init__(self):
        server.__init__(self)
        self.urlservice = "sapling"


    def metabolic_reconstruction(self,data):
        '''
        This method will find for each subsystem, the subsystem variant that contains a
        maximal subset of the roles in an incoming list, and output the ID of the
        variant and a list of the roles in it.
        '''

        self.function('metabolic_reconstruction')
        return self.retrieve(data)

    def find_special_proteins(self,data):
        '''
        This method searches for special proteins in a list of contigs. The method is
        specifically designed to find selenoproteins and pyrrolysoproteins, but custom
        protein templates can be specified to allow searching for any type of protein
        family.
        '''

        self.function('find_special_proteins')
        return self.retrieve(data)

    def assign_function_to_prot(self,data):
        '''
        For each incoming protein sequence, attempt to assign a function.
        There are two ways functions can get assigned.  The first is based on
        kmers, and these are normally viewed as the most reliable (at least
        they give a consistent vocabulary!).  If no kmer match is made,
        you can optionally try to make an assignment based on similarity
        computations.
        The attempt is made using kmer-technology.  A pass through the sequence
        will locate "signature kmers", and scores will be computed.
        The scores are based on the number of nonoverlapping hits, the number of
        overlapping hits, and the difference in counts between hits against the
        most probable function's kmer-set and the next most probable function's 
        kmer set.  Basically, we compute all matching kmers.  Then, we split them
        into sets based on the predictions each would make (each kmer, in effect,
        predicts a single function).  One threshhold (the B<scoreThreshold>) is
        the difference between total number of overlapping hits for the "best function"
        versus the total number for the "next best".  B<hitTheshold> is the number
        of overlapping hits required for the "best function".  Similarly,
        B<seqHitThreshold> is the minimum number of non-overlapping hits.
        Now, to add complexity, these thresholds are based on counting "1" for
        each matched Kmer.  That is, the B<scoreThreshold> is normally thought of
        as a difference in the number of occurrences.  However, you may wish to
        "normalize" this threshold by dividing the count by the length of the sequence.
        This then gives scores between 0 and 1, rather than between 0 and the length
        of the sequence (-K if you wish to be pedantic).
        '''

        self.function('assign_function_to_prot')
        return self.retrieve(data)

    def call_genes(self,data):
        '''
        Call the protein-encoding genes for the specified DNA sequences. 
        '''

        self.function('call_genes')
        return self.retrieve(data)

    def find_rnas(self,data):
        '''
        Call the RNAs for the specified DNA sequences.
        '''

        self.function('find_rnas')
        return self.retrieve(data)

    def assign_functions_to_dna(self,data):
        '''
        Analyze DNA sequences and output regions that probably belong to FIGfams.
        The selected regions will be high-probability candidates for protein
        encoding sequences.
        '''

        self.function('assign_functions_to_dna')
        return self.retrieve(data)

