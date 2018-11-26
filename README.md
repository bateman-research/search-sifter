# Search-Sifter
Search-Sifter uses hashes to rapidly compare protein families, in order to
analyse their relationship to each other.

## Installation

Search-Sifter requires Python 3.3 or greater. It's recommended that Search-Sifter
is installed into a virtual environment.

To install Search-Sifter, clone this repository:

    git clone URLHERE

Install the Search-Sifter package:

    pip install [path to Search-Sifter]

## Usage

### Generating Pfam hashes
First download a copy of the Pfam database in Stockholm file format. For example:

    wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.full.uniprot.gz

Run Search-Sifter with appropriate arguments depending on desired hash length
and window size:

    python -m searchsifter.scripts.generate_residue_hashes -n [hash length]
    -w [window size(s)] -o [output directory] -p [path to Pfam file]
    -t stockholm

For each window size specified, a file `rhashes_[w].json` will be created in the
output directory. This will contain a hash for each of the families in the
input file. The hashes are stored in a JSON dictonary keyed by family
accession. Each hash is a JSON list. The elements of the list are in the
following format:
    [hash, [protein accession, chunk number]]

### Running analysis

Two scripts are provided.

#### Accuracy

To analyse the accuracy of generated hashes:

    python -m searchsifter.scripts.performance -s [test families]
    -p [path to Pfam file] -t stockholm -a [hash files] [-n [hash lengths]]

Where `hash files` should be paths to one or more of the hashes generated by
`generate_residue_hashes`.

In `test families`, supply a file containing a JSON list of accessions.

In `hash lengths`, specify one or more lengths of hashes to test. Note that the
script can shorten longer hashes. So if you have generated hashes of length 800,
you could use `-n 200 400 800` to test hashes of length 200, 400 and 800.

If the n flag is not used, the script will not use hashes to estimate the
index and containment, but will instead compute the true index and containment.

The script will write to standard output a TSV file with the following columns:

    family_A                Accession of the family being compared to
    family_B                Accession from test families
    jaccard_index           Estimated (or true if -n is not used) Jaccard Index
    jaccard_containment	    Estimated (or true if -n is not used) Jaccard Containment
    n                       Length of hash (or 1 if -n is not used)
    w                       Window size (or 1 if -n is not used)
    type                    "estimated" if -n is used, "exact" otherwise

#### Time

To analyse performance:

    python -m searchsifter.scripts.time -s [test families]
    -p [path to Pfam file] -t stockholm -a [hash files] -n [hash lengths]

Arguments are as for searchsifter.scripts.performance.

The script will write to standard output a TSV with the following columns:

    test_acc    Accesssion from the test families
    ji_time     Time to estimate Jaccard index
    jc_time     Time to estimate Jaccard containment
    n           Length of hash
    w           Window size
    size        Family size in number of proteins

### Further usage

The file `searchsifter/Family.py` provides functions for creating objects to
represent protein families, and comparing them to each other.

The file `searchsifter/relationships/hmmer.py` provides the ability to create
Family objects from HMMER searches, given a Stockholm format output file.

The file `searchsifter/relationships/jaccard.py` provides functions for
calculating Jaccard index and containment.
`searchsifter/relationships/minhash.py` provides functions for estimating these
using MinHash.

Further documentation and usage is given as docstrings.