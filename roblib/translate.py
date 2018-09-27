


@aa_1_letter_order = qw( A C D E F G H I K L M N P Q R S T V W Y );  # Alpha by 1 letter
@aa_3_letter_order = qw( A R N D C Q E G H I L K M F P S T W Y V );  # PAM matrix order
@aa_n_codon_order  = qw( L R S A G P T V I C D E F H K N Q Y M W );

%genetic_code = (

    # DNA version

    TTT => 'F',  TCT => 'S',  TAT => 'Y',  TGT => 'C',
    TTC => 'F',  TCC => 'S',  TAC => 'Y',  TGC => 'C',
    TTA => 'L',  TCA => 'S',  TAA => '*',  TGA => '*',
    TTG => 'L',  TCG => 'S',  TAG => '*',  TGG => 'W',
    CTT => 'L',  CCT => 'P',  CAT => 'H',  CGT => 'R',
    CTC => 'L',  CCC => 'P',  CAC => 'H',  CGC => 'R',
    CTA => 'L',  CCA => 'P',  CAA => 'Q',  CGA => 'R',
    CTG => 'L',  CCG => 'P',  CAG => 'Q',  CGG => 'R',
    ATT => 'I',  ACT => 'T',  AAT => 'N',  AGT => 'S',
    ATC => 'I',  ACC => 'T',  AAC => 'N',  AGC => 'S',
    ATA => 'I',  ACA => 'T',  AAA => 'K',  AGA => 'R',
    ATG => 'M',  ACG => 'T',  AAG => 'K',  AGG => 'R',
    GTT => 'V',  GCT => 'A',  GAT => 'D',  GGT => 'G',
    GTC => 'V',  GCC => 'A',  GAC => 'D',  GGC => 'G',
    GTA => 'V',  GCA => 'A',  GAA => 'E',  GGA => 'G',
    GTG => 'V',  GCG => 'A',  GAG => 'E',  GGG => 'G',

    #  The following ambiguous encodings are not necessary,  but
    #  speed up the processing of some ambiguous triplets:

    TTY => 'F',  TCY => 'S',  TAY => 'Y',  TGY => 'C',
    TTR => 'L',  TCR => 'S',  TAR => '*',
                 TCN => 'S',
    CTY => 'L',  CCY => 'P',  CAY => 'H',  CGY => 'R',
    CTR => 'L',  CCR => 'P',  CAR => 'Q',  CGR => 'R',
    CTN => 'L',  CCN => 'P',               CGN => 'R',
    ATY => 'I',  ACY => 'T',  AAY => 'N',  AGY => 'S',
                 ACR => 'T',  AAR => 'K',  AGR => 'R',
                 ACN => 'T',
    GTY => 'V',  GCY => 'A',  GAY => 'D',  GGY => 'G',
    GTR => 'V',  GCR => 'A',  GAR => 'E',  GGR => 'G',
    GTN => 'V',  GCN => 'A',               GGN => 'G'
);




%amino_acid_codons_DNA = (
         L  => [ qw( TTA TTG CTA CTG CTT CTC ) ],
         R  => [ qw( AGA AGG CGA CGG CGT CGC ) ],
         S  => [ qw( AGT AGC TCA TCG TCT TCC ) ],
         A  => [ qw( GCA GCG GCT GCC ) ],
         G  => [ qw( GGA GGG GGT GGC ) ],
         P  => [ qw( CCA CCG CCT CCC ) ],
         T  => [ qw( ACA ACG ACT ACC ) ],
         V  => [ qw( GTA GTG GTT GTC ) ],
         I  => [ qw( ATA ATT ATC ) ],
         C  => [ qw( TGT TGC ) ],
         D  => [ qw( GAT GAC ) ],
         E  => [ qw( GAA GAG ) ],
         F  => [ qw( TTT TTC ) ],
         H  => [ qw( CAT CAC ) ],
         K  => [ qw( AAA AAG ) ],
         N  => [ qw( AAT AAC ) ],
         Q  => [ qw( CAA CAG ) ],
         Y  => [ qw( TAT TAC ) ],
         M  => [ qw( ATG ) ],
         U  => [ qw( TGA ) ],
         W  => [ qw( TGG ) ],

         l  => [ qw( tta ttg cta ctg ctt ctc ) ],
         r  => [ qw( aga agg cga cgg cgt cgc ) ],
         s  => [ qw( agt agc tca tcg tct tcc ) ],
         a  => [ qw( gca gcg gct gcc ) ],
         g  => [ qw( gga ggg ggt ggc ) ],
         p  => [ qw( cca ccg cct ccc ) ],
         t  => [ qw( aca acg act acc ) ],
         v  => [ qw( gta gtg gtt gtc ) ],
         i  => [ qw( ata att atc ) ],
         c  => [ qw( tgt tgc ) ],
         d  => [ qw( gat gac ) ],
         e  => [ qw( gaa gag ) ],
         f  => [ qw( ttt ttc ) ],
         h  => [ qw( cat cac ) ],
         k  => [ qw( aaa aag ) ],
         n  => [ qw( aat aac ) ],
         q  => [ qw( caa cag ) ],
         y  => [ qw( tat tac ) ],
         m  => [ qw( atg ) ],
         u  => [ qw( tga ) ],
         w  => [ qw( tgg ) ],

        '*' => [ qw( TAA TAG TGA ) ]
);

