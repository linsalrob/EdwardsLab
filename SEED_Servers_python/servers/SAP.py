
from .servers import server

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



    def methods(self,data):
        '''
        Return a reference to a list of the methods allowed on this object.
        '''

        self.function('methods')
        return self.retrieve(data)

    def exists(self,data):
        '''
        Return a hash indicating which of the specified objects of the given type exist
        in the database. This method is used as a general mechanism for finding what
        exists and what doesn't exist when you know the ID. In particular, you can use
        it to check for the presence or absence of subsystems, genomes, features,
        or FIGfams.
        '''

        self.function('exists')
        return self.retrieve(data)

    def get(self,data):
        '''
        Query the Sapling database. This is a variant of the L</query> method in
        which a certain amount of power is sacrificed for ease of use. Instead of
        a full-blown filter clause, the caller specifies a filter hash that maps
        field identifiers to values.
        '''

        self.function('get')
        return self.retrieve(data)

    def query(self,data):
        '''
        This method queries the Sapling database and returns a reference to a list of
        lists. The query is specified in the form of an object name string, a filter
        string, an optional list of parameter values, and a list of desired output
        fields. The result document can be thought of as a two-dimensional array, with
        each row being a record returned by the query and each column representing an
        output field.
        This function buys a great deal of flexibility as the cost of ease of use.
        Before attempting to formulate a query, you will need to look at the
        L<ERDB> documentation.
        '''

        self.function('query')
        return self.retrieve(data)

    def select(self,data):
        '''
        Query the Sapling database. This is a variant of the L</get> method in
        which a further amount of power is sacrificed for ease of use. The
        return is a list of lists, and the criteria are always in the form of
        lists of possible values.
        '''

        self.function('select')
        return self.retrieve(data)

    def equiv_precise_assertions(self,data):
        '''
        Return the assertions for all genes in the database that match the
        identified gene. The gene can be specified by any prefixed gene
        identifier (e.g. C<uni|AYQ44>, C<gi|85841784>, or
        C<fig|360108.3.peg.1041>).
        '''

        self.function('equiv_precise_assertions')
        return self.retrieve(data)

    def equiv_sequence_assertions(self,data):
        '''
        Return the assertions for all genes in the database that match the
        identified protein sequences. A protein sequence can be identified by a
        protein MD5 code or any prefixed gene identifier (e.g. C<uni|AYQ44>,
        C<gi|85841784>, or C<fig|360108.3.peg.1041>).
        '''

        self.function('equiv_sequence_assertions')
        return self.retrieve(data)

    def feature_assignments(self,data):
        '''
        Return all features of the specified type for the specified genome along
        with their assignments.
        '''

        self.function('feature_assignments')
        return self.retrieve(data)

    def ids_to_assertions(self,data):
        '''
        Return the assertions associated with each prefixed ID.
        '''

        self.function('ids_to_assertions')
        return self.retrieve(data)

    def ids_to_annotations(self,data):
        '''
        Return the annotations associated with each prefixed ID. Annotations are
        comments attached to each feature (gene), and include past functional
        assignments as well as more general information.
        '''

        self.function('ids_to_annotations')
        return self.retrieve(data)

    def ids_to_functions(self,data):
        '''
        Return the functional assignment for each feature in the incoming list.
        '''

        self.function('ids_to_functions')
        return self.retrieve(data)

    def occ_of_role(self,data):
        '''
        Search for features in a specified genome with the indicated roles or
        functions.
        '''

        self.function('occ_of_role')
        return self.retrieve(data)

    def all_complexes(self,data):
        '''
        Return a list of all the complexes in the database.
        '''

        self.function('all_complexes')
        return self.retrieve(data)

    def all_models(self,data):
        '''
        Return a hash of all the models in the database, mapping each one to the relevant
        genome.
        '''

        self.function('all_models')
        return self.retrieve(data)

    def all_reactions(self,data):
        '''
        Return a list of all the reactions in the database.
        '''

        self.function('all_reactions')
        return self.retrieve(data)

    def all_roles_used_in_models(self,data):
        '''
        Return a list of all the roles used in models.
        '''

        self.function('all_roles_used_in_models')
        return self.retrieve(data)

    def complex_data(self,data):
        '''
        Return the specified data items for each incoming reaction complex.
        '''

        self.function('complex_data')
        return self.retrieve(data)

    def coupled_reactions(self,data):
        '''
        For each of a set of reactions, get the adjacent reactions in the metabolic network.
        Two reactions are considered I<adjacent> if they share at least one compound that
        is neither a cofactor or a ubiquitous compound (like water or oxygen). The compounds
        that relate the adjacent reactions are called the I<connecting compounds>. In most cases,
        each pair of adjacent reactions will have only one connecting compound, but this is
        not guaranteed to be true.
        '''

        self.function('coupled_reactions')
        return self.retrieve(data)

    def models_to_reactions(self,data):
        '''
        Return the list of reactions in each specified model.
        '''

        self.function('models_to_reactions')
        return self.retrieve(data)

    def reaction_neighbors(self,data):
        '''
        Return a list of the reactions in the immediate neighborhood of the specified reactions.
        A separate neighborhood list will be generated for each incoming reaction; the neighborhood
        will consist of reactions connected to the incoming reaction and reactions connected to those
        reactions up to the specified depth. (Two reactions are I<connected> if they have a compound
        in common that is not a cofactor or a ubiquitous chemical like water or ATP).
        '''

        self.function('reaction_neighbors')
        return self.retrieve(data)

    def reaction_path(self,data):
        '''
        Find the shortest reaction path that represents as many of the specified roles
        as possible. Note that since the a reaction may be associated with multiple
        roles, it is possible for a single role to be represented more than once in
        the path.
        The search is artificially limited to paths under a maximum length that can
        be specified in the parameters.
        '''

        self.function('reaction_path')
        return self.retrieve(data)

    def reaction_strings(self,data):
        '''
        Return the display string for each reaction. The display string contains the compound IDs
        (as opposed to the atomic formulas) and the associated stoichiometries, with the
        substrates on the left of the arrow and the products on the right.
        '''

        self.function('reaction_strings')
        return self.retrieve(data)

    def reactions_to_complexes(self,data):
        '''
        Return the complexes containing each reaction. Note that most reactions
        are in more than one complex, so the complexes for each reaction are returned
        as a list.
        '''

        self.function('reactions_to_complexes')
        return self.retrieve(data)

    def reactions_to_roles(self,data):
        '''
        Return the roles associated with each reaction.
        '''

        self.function('reactions_to_roles')
        return self.retrieve(data)

    def role_neighbors(self,data):
        '''
        For each role, return a list of roles in the immediate chemical neighborhood. A role is
        in the immediate chemical neighborhood of another role if the two roles are associated with
        reactions that share a compound that is not ubiquitous or a cofactor.
        '''

        self.function('role_neighbors')
        return self.retrieve(data)

    def role_reactions(self,data):
        '''
        Return a list of all the reactions associated with each incoming role.
        '''

        self.function('role_reactions')
        return self.retrieve(data)

    def roles_to_complexes(self,data):
        '''
        Return the complexes (sets of related reactions) associated with each role in
        the incoming list. Roles trigger many complexes, and a complex may be triggered
        by many roles. A given role is considered either I<optional> or I<necessary>
        to the complex, and an indication of this will be included in the output.
        '''

        self.function('roles_to_complexes')
        return self.retrieve(data)

    def dlits_for_ids(self,data):
        '''
        Find the PUBMED literature references for a list of proteins. The
        proteins can be specified either was FIG feature IDs or protein
        sequence MD5s.
        '''

        self.function('dlits_for_ids')
        return self.retrieve(data)

    def equiv_ids_for_sequences(self,data):
        '''
        Find all the identifiers in the database that produce the specified proteins.
        '''

        self.function('equiv_ids_for_sequences')
        return self.retrieve(data)

    def find_closest_genes(self,data):
        '''
        Find the closest genes to the specified sequences in the specified genome.
        Each indicated sequence will be converted to a DNA sequence and then the contigs
        of the specified genome will be searched for the sequence. The genes in closest
        proximity to the sequence will be returned. The sequences are named; in the return
        hash, the genes found will be associated with the appropriate sequence name.
        '''

        self.function('find_closest_genes')
        return self.retrieve(data)

    def ids_to_sequences(self,data):
        '''
        Compute a DNA or protein string for each incoming feature ID.
        '''

        self.function('ids_to_sequences')
        return self.retrieve(data)

    def locs_to_dna(self,data):
        '''
        Return the DNA sequences for the specified locations.
        '''

        self.function('locs_to_dna')
        return self.retrieve(data)

    def roles_to_proteins(self,data):
        '''
        Return a list of the proteins associated with each of the incoming functional
        roles.
        '''

        self.function('roles_to_proteins')
        return self.retrieve(data)

    def upstream(self,data):
        '''
        Return the DNA sequences for the upstream regions of the specified
        features. The nucleotides inside coding regions are displayed in upper
        case; others are displayed in lower case.
        '''

        self.function('upstream')
        return self.retrieve(data)

    def all_experiments(self,data):
        '''
        Return a list of all the experiment names.
        '''

        self.function('all_experiments')
        return self.retrieve(data)

    def atomic_regulon_vectors(self,data):
        '''
        Return a map of the expression levels for each specified atomic regulon. The
        expression levels will be returned in the form of vectors with values C<-1>
        (suppressed), C<1> (expressed), or C<0> (unknown) in each position. The positions
        will correspond to the experiments in the order returned by L</genome_experiments>.
        '''

        self.function('atomic_regulon_vectors')
        return self.retrieve(data)

    def atomic_regulons(self,data):
        '''
        Return a map of the atomic regulons for the specified genome. Each atomic
        regulon is a set of genes that are always regulated together. The map will
        connect each regulon ID to a list of those genes. A given gene can only be
        in one atomic regulon.
        '''

        self.function('atomic_regulons')
        return self.retrieve(data)

    def coregulated_correspondence(self,data):
        '''
        Given a gene, return genes that may be coregulated because they correspond to
        coregulated genes in genomes for which we have expression data (an
        I<expression-analyzed genome>). For each incoming gene, a corresponding
        gene will be found in each expression-analyzed genome. The coregulated
        genes for the corresponding gene will be determined, and then these will be
        mapped back to the original genome. The resulting genes can be considered
        likely candidates for coregulation in the original genome.
        '''

        self.function('coregulated_correspondence')
        return self.retrieve(data)

    def coregulated_fids(self,data):
        '''
        Given a gene, return the coregulated genes and their pearson coefficients.
        Two genes are considered coregulated if there is some experimental evidence
        that their expression levels are related: the pearson coefficient indicates
        the strength of the relationship.
        '''

        self.function('coregulated_fids')
        return self.retrieve(data)

    def experiment_fid_levels(self,data):
        '''
        Given an experiment, return the on/off levels for all genes in that
        experiment. An on/off level is either C<1> (expressed), C<-1> (inhibited),
        or C<0> (unknown).
        '''

        self.function('experiment_fid_levels')
        return self.retrieve(data)

    def experiment_regulon_levels(self,data):
        '''
        Given an experiment, return the on/off levels for all atomic regulons
        affected by that experiment. An on/off level is either C<1> (expressed), C<-1>
        (inhibited), or C<0> (unknown).
        '''

        self.function('experiment_regulon_levels')
        return self.retrieve(data)

    def expressed_genomes(self,data):
        '''
        List the IDs of genomes for which expression data exists in the database.
        '''

        self.function('expressed_genomes')
        return self.retrieve(data)

    def fid_experiments(self,data):
        '''
        Return the expression levels for the specified features in all experiments for which they
        have results.
        '''

        self.function('fid_experiments')
        return self.retrieve(data)

    def fid_vectors(self,data):
        '''
        Return a map of the expression levels for each specified feature (gene). The
        expression levels will be returned in the form of vectors with values C<-1>
        (suppressed), C<1> (expressed), or C<0> (unknown) in each position. The positions
        will correspond to the experiments in the order returned by L</genome_experiments>.
        '''

        self.function('fid_vectors')
        return self.retrieve(data)

    def fids_expressed_in_range(self,data):
        '''
        Return for each genome the genes that are expressed in a given fraction of the experiments
        for that ganome.
        '''

        self.function('fids_expressed_in_range')
        return self.retrieve(data)

    def fids_to_regulons(self,data):
        '''
        Return the atomic regulons associated with each incoming gene.
        '''

        self.function('fids_to_regulons')
        return self.retrieve(data)

    def genome_experiments(self,data):
        '''
        Return a list of the experiments for each indicated genome.
        '''

        self.function('genome_experiments')
        return self.retrieve(data)

    def genome_experiment_levels(self,data):
        '''
        Return the expression levels for the specified features in all experiments for which they
        have results.
        '''

        self.function('genome_experiment_levels')
        return self.retrieve(data)

    def regulons_to_fids(self,data):
        '''
        Return the list of genes in each specified atomic regulon.
        '''

        self.function('regulons_to_fids')
        return self.retrieve(data)

    def compared_regions(self,data):
        '''
        Return information about the context of a focus gene and the corresponding genes in
        other genomes (known as I<pinned genes>). The information returned can be used to
        create a compare-regions display.
        The return information will be in the form of a reference to a list of contexts,
        each context containing genes in a region surrounding the pinned gene on a particular
        genome. The genome containing the focus gene will always be the first in the list.
        '''

        self.function('compared_regions')
        return self.retrieve(data)

    def equiv_sequence_ids(self,data):
        '''
        Return all identifiers for genes in the database that are
        protein-sequence-equivalent to the specified identifiers. In this case, the
        identifiers are assumed to be in their natural form (without prefixes). For
        each identifier, the identified protein sequences will be found and then
        for each protein sequence, all identifiers for that protein sequence or for
        genes that produce that protein sequence will be returned.
        Alternatively, you can ask for identifiers that are precisely equivalent, that is,
        that identify the same location on the same genome.
        '''

        self.function('equiv_sequence_ids')
        return self.retrieve(data)

    def fid_correspondences(self,data):
        '''
        Return the corresponding genes for the specified features in the specified genomes.
        The correspondences are determined in the same way as used by L</gene_correspondence_map>,
        but this method returns substantially less data.
        '''

        self.function('fid_correspondences')
        return self.retrieve(data)

    def fid_locations(self,data):
        '''
        Return the DNA locations for the specified features.
        '''

        self.function('fid_locations')
        return self.retrieve(data)

    def fid_map_for_genome(self,data):
        '''
        Find FIG IDs corresponding to caller-provided genes in a specific genome.
        In some situations you may have multiple external identifiers for
        various genes in a genome without knowing which ones are present in the Sapling
        database and which are not. The external identifiers present in the Sapling
        database are culled from numerous sources, but different genomes will tend to
        have coverage from different identifier types: some genomes are represented
        heavily by CMR identifiers and have no Locus Tags, others have lots of Locus
        Tags but no CMR identifiers, and so forth. This method allows you to throw everything
        you have at the database in hopes of finding a match.
        '''

        self.function('fid_map_for_genome')
        return self.retrieve(data)

    def fid_possibly_truncated(self,data):
        '''
        For each specified gene, return C<stop> if its end is possibly truncated,
        C<start> if its beginning is possibly truncated, and an empty string
        otherwise. Truncation occurs if the gene is located near either edge of a
        contig.
        '''

        self.function('fid_possibly_truncated')
        return self.retrieve(data)

    def fids_to_ids(self,data):
        '''
        Find all aliases and/or synonyms for the specified FIG IDs. For each FIG
        ID, a hash will be returned that maps each ID type to a list of the IDs
        of that type.
        '''

        self.function('fids_to_ids')
        return self.retrieve(data)

    def fids_to_proteins(self,data):
        '''
        Return the ID or amino acid sequence associated with each specified gene's protein. If the gene
        does not produce a protein, it will not be included in the output.
        '''

        self.function('fids_to_proteins')
        return self.retrieve(data)

    def fids_with_evidence_codes(self,data):
        '''
        Return the ID, assignment, and evidence for all features having an
        evidence code of one of the specified types. The output can be restricted
        to one or more specified genomes.
        '''

        self.function('fids_with_evidence_codes')
        return self.retrieve(data)

    def genes_in_region(self,data):
        '''
        Return a list of the IDs for the features that overlap the specified
        regions on a contig.
        '''

        self.function('genes_in_region')
        return self.retrieve(data)

    def ids_to_data(self,data):
        '''
        Return the specified data items for the specified features.
        '''

        self.function('ids_to_data')
        return self.retrieve(data)

    def ids_to_fids(self,data):
        '''
        Return a list of the FIG IDs corresponding to each of the specified
        identifiers. The correspondence can either be gene-based (same feature)
        or sequence-based (same protein).
        '''

        self.function('ids_to_fids')
        return self.retrieve(data)

    def ids_to_genomes(self,data):
        '''
        Return the genome information for each incoming gene ID.
        '''

        self.function('ids_to_genomes')
        return self.retrieve(data)

    def ids_to_lengths(self,data):
        '''
        Return the DNA or protein length of each specified gene.
        '''

        self.function('ids_to_lengths')
        return self.retrieve(data)

    def make_runs(self,data):
        '''
        Look at sequences of feature IDs and separate them into operons. An
        operon contains features that are close together on the same contig going
        in the same direction.
        '''

        self.function('make_runs')
        return self.retrieve(data)

    def proteins_to_fids(self,data):
        '''
        Return the FIG feature IDs associated with each incoming protein. The protein can be
        specified as an amino acid sequence or MD5 protein ID.
        '''

        self.function('proteins_to_fids')
        return self.retrieve(data)

    def all_figfams(self,data):
        '''
        Return a list of all the FIGfams along with their functions. Optionally, you
        can specify a role or a function, and only FIGfams with that role or function
        will be returned.
        '''

        self.function('all_figfams')
        return self.retrieve(data)

    def discriminating_figfams(self,data):
        '''
        Determine the FIGfams that discriminate between two groups of genomes.
        A FIGfam discriminates between genome groups if it is common in one group and
        uncommon in the other. The degree of discrimination is assigned a score based
        on statistical significance, with 0 being insignificant and 2 being extremely
        significant. FIGfams with a score greater than 1 are returned by this method.
        '''

        self.function('discriminating_figfams')
        return self.retrieve(data)

    def figfam_fids(self,data):
        '''
        Return a list of all the protein encoding genes in a FIGfam. The genes
        can be returned as IDs or as FASTA strings.
        '''

        self.function('figfam_fids')
        return self.retrieve(data)

    def figfam_fids_batch(self,data):
        '''
        Return a list of all the protein encoding genes in one or more FIGfams. This
        method is an alternative to L</figfam_fids> that is faster when you need the
        feature IDs but not the protein sequences.
        '''

        self.function('figfam_fids_batch')
        return self.retrieve(data)

    def figfam_function(self,data):
        '''
        For each incoming FIGfam ID, return its function, that is, the common
        functional assignment of all its members.
        '''

        self.function('figfam_function')
        return self.retrieve(data)

    def genome_figfams(self,data):
        '''
        Compute the list of FIGfams represented in each specific genome.
        '''

        self.function('genome_figfams')
        return self.retrieve(data)

    def ids_to_figfams(self,data):
        '''
        This method returns a hash mapping each incoming feature to its FIGfam.
        '''

        self.function('ids_to_figfams')
        return self.retrieve(data)

    def related_figfams(self,data):
        '''
        This method takes a list of FIGfam IDs. For each FIGfam, it returns a
        list of FIGfams related to it by functional coupling.
        '''

        self.function('related_figfams')
        return self.retrieve(data)

    def roles_to_figfams(self,data):
        '''
        For each incoming role, return a list of the FIGfams that implement
        the role, that is, whose functional assignments include the role.
        '''

        self.function('roles_to_figfams')
        return self.retrieve(data)

    def clusters_containing(self,data):
        '''
        This method takes as input a list of FIG feature IDs. For each feature, it
        returns the IDs and functions of other features in the same cluster of
        functionally-coupled features.
        '''

        self.function('clusters_containing')
        return self.retrieve(data)

    def co_occurrence_evidence(self,data):
        '''
        For each specified pair of genes, this method returns the evidence that
        the genes are functionally coupled (if any); that is, it returns a list
        of the physically close homologs for the pair.
        '''

        self.function('co_occurrence_evidence')
        return self.retrieve(data)

    def conserved_in_neighborhood(self,data):
        '''
        This method takes a list of feature IDs. For each feature ID, it will
        return the set of other features to which it is functionally coupled,
        along with the appropriate score.
        '''

        self.function('conserved_in_neighborhood')
        return self.retrieve(data)

    def pairsets(self,data):
        '''
        This method takes as input a list of functional-coupling pair set IDs
        (such as those returned in the output of L</conserved_in_neighborhood>). For
        each pair set, it returns the set's score (number of significant couplings) and
        a list of the coupled pairs in the set.
        '''

        self.function('pairsets')
        return self.retrieve(data)

    def related_clusters(self,data):
        '''
        This method returns the functional-coupling clusters related to the specified
        input features. Each cluster contains features on a single genome that are
        related by functional coupling.
        '''

        self.function('related_clusters')
        return self.retrieve(data)

    def all_features(self,data):
        '''
        Return a list of the IDs for all features of a specified type in a specified
        genome.
        '''

        self.function('all_features')
        return self.retrieve(data)

    def all_genomes(self,data):
        '''
        Return a list of the IDs for all the genomes in the system.
        '''

        self.function('all_genomes')
        return self.retrieve(data)

    def all_proteins(self,data):
        '''
        Return the protein sequences for all protein-encoding genes in the specified
        genome.
        '''

        self.function('all_proteins')
        return self.retrieve(data)

    def close_genomes(self,data):
        '''
        Find the genomes functionally close to the input genomes.
        Functional closeness is determined by the number of FIGfams in common. As a result,
        this method will not produce good results for genomes that do not have good FIGfam
        coverage.
        '''

        self.function('close_genomes')
        return self.retrieve(data)

    def contig_sequences(self,data):
        '''
        Return the DNA sequences for the specified contigs.
        '''

        self.function('contig_sequences')
        return self.retrieve(data)

    def contig_lengths(self,data):
        '''
        Return the lengths for the specified contigs.
        '''

        self.function('contig_lengths')
        return self.retrieve(data)

    def gene_correspondence_map(self,data):
        '''
        Return a map of genes in the specified second genome that correspond to genes in
        the specified first genome.
        '''

        self.function('gene_correspondence_map')
        return self.retrieve(data)

    def genome_contig_md5s(self,data):
        '''
        For each incoming genome, return a hash mapping its contigs to their MD5 identifiers.
        '''

        self.function('genome_contig_md5s')
        return self.retrieve(data)

    def genome_contigs(self,data):
        '''
        For each incoming genome, return a list of its contigs.
        '''

        self.function('genome_contigs')
        return self.retrieve(data)

    def genome_data(self,data):
        '''
        Return the specified data items for the specified genomes.
        '''

        self.function('genome_data')
        return self.retrieve(data)

    def genome_domain(self,data):
        '''
        Return the domain for each specified genome (e.g. C<Archaea>, C<Bacteria>, C<Plasmid>).
        '''

        self.function('genome_domain')
        return self.retrieve(data)

    def genome_fid_md5s(self,data):
        '''
        For each incoming genome, return a hash mapping its genes to their MD5 identifiers.
        '''

        self.function('genome_fid_md5s')
        return self.retrieve(data)

    def genome_ids(self,data):
        '''
        Find the specific genome ID for each specified genome name or taxonomic number.
        This method helps to find the correct version of a given genome when only the
        species and strain are known.
        '''

        self.function('genome_ids')
        return self.retrieve(data)

    def genome_metrics(self,data):
        '''
        For each incoming genome ID, returns the number of contigs, the total
        number of base pairs in the genome's DNA, and the genome's default genetic
        code.
        '''

        self.function('genome_metrics')
        return self.retrieve(data)

    def genome_names(self,data):
        '''
        Return the name of the genome containing each specified feature or genome.
        '''

        self.function('genome_names')
        return self.retrieve(data)

    def genomes_by_md5(self,data):
        '''
        Find the genomes associated with each specified MD5 genome identifier. The MD5
        genome identifier is computed from the DNA sequences of the genome's contigs; as
        a result, two genomes with identical sequences arranged in identical contigs
        will have the same MD5 identifier even if they have different genome IDs.
        '''

        self.function('genomes_by_md5')
        return self.retrieve(data)

    def intergenic_regions(self,data):
        '''
        Return a list of L</Location Strings> for the regions in the specified genome that are
        not occupied by genes of the specified types. All of these will be construed to be on
        the forward strand, and sorted by contig ID and start location within contig.
        '''

        self.function('intergenic_regions')
        return self.retrieve(data)

    def is_prokaryotic(self,data):
        '''
        For each incoming genome ID, returns 1 if it is prokaryotic and 0
        otherwise.
        '''

        self.function('is_prokaryotic')
        return self.retrieve(data)

    def mapped_genomes(self,data):
        '''
        For each incoming genome, return a list of the genomes that have an existing
        gene correspondence map (see L<ServerThing/Gene Correspondence List>). Gene
        correspondence maps indicate which genes in the target genome are the best hit
        of each gene in the source genome. If a correspondence map does not yet exist,
        it will be created when you ask for it, but this is an expensive process and it
        is sometimes useful to find an alternate genome that will give you a faster
        result.
        '''

        self.function('mapped_genomes')
        return self.retrieve(data)

    def otu_members(self,data):
        '''
        For each incoming genome, return the name and ID of each other genome in the same
        OTU.
        '''

        self.function('otu_members')
        return self.retrieve(data)

    def representative(self,data):
        '''
        Return the representative genome for each specified incoming genome ID.
        Genomes with the same representative are considered closely related, while
        genomes with a different representative would be considered different
        enough that similarities between them have evolutionary significance.
        '''

        self.function('representative')
        return self.retrieve(data)

    def representative_genomes(self,data):
        '''
        Compute mappings for the genome sets (OTUs) in the database. This method will
        return a mapping from each genome to its genome set ID and from each
        genome set ID to a list of the genomes in the set. For the second
        mapping, the first genome in the set will be the representative.
        This method does not require any parameters.
        '''

        self.function('representative_genomes')
        return self.retrieve(data)

    def submit_gene_correspondence(self,data):
        '''
        Submit a set of gene correspondences to be stored on the server.
        '''

        self.function('submit_gene_correspondence')
        return self.retrieve(data)

    def taxonomy_of(self,data):
        '''
        Return the taxonomy of each specified genome. The taxonomy will start at
        the domain level and moving down to the node where the genome is
        attached.
        '''

        self.function('taxonomy_of')
        return self.retrieve(data)

    def scenario_names(self,data):
        '''
        Return the names of all the scenarios for the specified subsystem. Each scenario
        has an internal ID number and a common name. This method returns both.
        '''

        self.function('scenario_names')
        return self.retrieve(data)

    def all_subsystems(self,data):
        '''
        Return a list of all subsystems in the system. For each subsystem, this
        method will return the ID, curator, the classifications, and roles.
        '''

        self.function('all_subsystems')
        return self.retrieve(data)

    def classification_of(self,data):
        '''
        Return the classification for each specified subsystem.
        '''

        self.function('classification_of')
        ss = self.retrieve(data)
        for k in ss:
            if (len(ss[k]) == 1):
                print(k)
                ss[k].append(None)
        return ss

    def genomes_to_subsystems(self,data):
        '''
        Return a list of the subsystems participated in by each of the specified
        genomes.
        '''

        self.function('genomes_to_subsystems')
        return self.retrieve(data)

    def get_subsystems(self,data):
        '''
        Get a complete description of each specified subsystem. This will include
        the basic subsystem properties, the list of roles, and the spreadsheet.
        '''

        self.function('get_subsystems')
        return self.retrieve(data)

    def ids_in_subsystems(self,data):
        '''
        Return the features of each specified subsystems in the specified genome, or
        alternatively, return all features of each specified subsystem.
        '''

        self.function('ids_in_subsystems')
        return self.retrieve(data)

    def ids_to_publications(self,data):
        '''
        Return the PUBMED ID and title of each publication relevant to the specified
        feature IDs.
        '''

        self.function('ids_to_publications')
        return self.retrieve(data)

    def ids_to_subsystems(self,data):
        '''
        Return the subsystem and role for each feature in the incoming list. A
        feature may have multiple roles in a subsystem and may belong to multiple
        subsystems, so the role/subsystem information is returned in the form of
        a list of ordered pairs for each feature.
        '''

        self.function('ids_to_subsystems')
        return self.retrieve(data)

    def is_in_subsystem(self,data):
        '''
        Return the subsystem and role for each specified feature.
        '''

        self.function('is_in_subsystem')
        return self.retrieve(data)

    def is_in_subsystem_with(self,data):
        '''
        For each incoming feature, returns a list of the features in the same genome that
        are part of the same subsystem. For each other feature returned, its role,
        functional assignment, subsystem variant, and subsystem ID will be returned as
        well.
        '''

        self.function('is_in_subsystem_with')
        return self.retrieve(data)

    def pegs_implementing_roles(self,data):
        '''
        Given a subsystem and a list of roles, return a list of the subsystem's
        features for each role.
        '''

        self.function('pegs_implementing_roles')
        return self.retrieve(data)

    def pegs_in_subsystems(self,data):
        '''
        This method takes a list of genomes and a list of subsystems and returns
        a list of the roles represented in each genome/subsystem pair.
        '''

        self.function('pegs_in_subsystems')
        return self.retrieve(data)

    def pegs_in_variants(self,data):
        '''
        This method takes a list of genomes and a list of subsystems and returns
        a list of the pegs represented in each genome/subsystem pair.
        The main difference between this method and L</pegs_in_subsystems> is in
        the organization of the output, which is more like a subsystem spreadsheet.
        '''

        self.function('pegs_in_variants')
        return self.retrieve(data)

    def roles_exist_in_subsystem(self,data):
        '''
        Indicate which roles in a given list belong to a specified subsystem.
        '''

        self.function('roles_exist_in_subsystem')
        return self.retrieve(data)

    def roles_to_subsystems(self,data):
        '''
        Return the subsystems containing each specified role.
        '''

        self.function('roles_to_subsystems')
        return self.retrieve(data)

    def rows_of_subsystem(self,data):
        '''
        Return the subsystem row for each subsystem/genome pair. A row in this
        case consists of a reference to a hash mapping role names to a list of
        the FIG feature IDs for the features in the genome performing that
        role.
        In the Sapling database, a subsystem row is represented by the
        B<MolecularMachine> entity. The strategy of this method is therefore
        to find the molecular machine for each subsystem/genome pair, and
        then use its ID to get the roles and features.
        '''

        self.function('rows_of_subsystem')
        return self.retrieve(data)

    def subsystem_data(self,data):
        '''
        For each incoming subsystem ID, return the specified data field. This
        method can be used to find the curator, description, or version of the
        specified subsystems.
        '''

        self.function('subsystem_data')
        return self.retrieve(data)

    def subsystem_genomes(self,data):
        '''
        For each subsystem, return the genomes that participate in it and their
        associated variant codes.
        '''

        self.function('subsystem_genomes')
        return self.retrieve(data)

    def subsystem_names(self,data):
        '''
        Return a list of all subsystems in the database.
        '''

        self.function('subsystem_names')
        return self.retrieve(data)

    def subsystem_roles(self,data):
        '''
        Return the list of roles for each subsystem, in order.
        '''

        self.function('subsystem_roles')
        return self.retrieve(data)

    def subsystem_spreadsheet(self,data):
        '''
        This method takes a list of subsystem IDs, and for each one returns a
        list of the features in the subsystem. For each feature, it will include
        the feature's functional assignment, the subsystem name and variant
        (spreadsheet row), and its role (spreadsheet column).
        '''

        self.function('subsystem_spreadsheet')
        return self.retrieve(data)

    def subsystem_type(self,data):
        '''
        For each incoming subsystem, return TRUE if it has the specified
        characteristic, else FALSE.
        '''

        self.function('subsystem_type')
        return self.retrieve(data)

    def subsystems_for_role(self,data):
        '''
        For each role, return a list of the subsystems containing that role. The
        results can be filtered to include unusable subsystems or exclude subsystems
        of certain exotic types.
        '''

        self.function('subsystems_for_role')
        return self.retrieve(data)

