Scavager - a proteomics post-search validation tool
---------------------------------------------------------------

The .pep.xml or .mzid files are required for basic operation of the script. Currently supported search engines:
Identipy, X!Tandem, Comet, MSFragger, msgf+, Morpheus.

.fasta file is required for calculation NSAF (label-free quantitation index), protein sequence coverage and amino acid statistics.

For msgf+ and morpheus search engines it is desirable to provide cleavage rules used in search (These search engines do not report number of missed cleavages for peptides).

The output of Scavager contains:


    tab-separated table with unfiltered peptide-spectrum matches (ends with _PSMs_full.tsv)

    tab-separated table with identified peptide-spectrum matches at 1% PSM FDR (ends with _PSMs.tsv)

    tab-separated table with identified peptides at 1% peptide FDR (ends with _peptides.tsv)

    tab-separated table with identified proteins without grouping at 1% protein FDR (ends with _proteins.tsv)

    tab-separated table with identified protein groups at 1% protein FDR (ends with _protein_groups.tsv)

    png figure with PSM, peptide and protein features distributions

Installation
-----
Using the pip:

    pip install Scavager


Usage
-----
Algorithm can be run with following command (works with Python2.7/Python3+):

    scavager path_to_pepXML/MZID

    OR

    scavager -h

Links
-----

- BitBucket repo & issue tracker: https://bitbucket.org/markmipt/scavager
- Mailing list: pyteomics@googlegroups.com, markmipt@gmail.com