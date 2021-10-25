User guide
==========

Installation
------------

To install, clone the repository or `download and unzip
<https://github.com/Singh-Lab/DeMaSk/archive/master.zip>`_.  To get
Python dependencies and be able to run or import the modules from any
directory, install with pip::

  pip install -e DeMaSk/

Unless you're supplying your own aligned homologs, you'll need to have
the ``blastp`` program.  It can be downloaded as part of `BLAST+
<https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`_.

The ``blastp`` step requires a formatted sequence database.  To use
UniRef90, which demask.princeton.edu uses, download the `zipped fasta
file <https://www.uniprot.org/downloads>`_, unzip, and use
``makeblastdb`` from BLAST+ to format the database.

To avoid having to specify the location of the ``blastp`` binary and
the database for every DeMaSk run, put them in a config file, e.g.::

  blastp=/usr/local/ncbi/blast/bin/blastp
  db=/path/to/database/uniref90.fasta

By default, DeMaSk will look for the config file
``DeMaSk/config.ini``, which can also contain any other command line
arguments, such as ``nseqs``, ``threads``, and ``matrix``.

Running from the command line
-----------------------------

Once the demask package is installed, you can run run it from
anywhere.  If you don't have aligned homologs for your query yet, run
the demask.homologs module::
 
 python3 -m demask.homologs -s myquery.fa -o myquery_homologs.a2m

The above command will also produce a file myquery.blast.json
containing the intermediate blastp output.

Then, generate fitness impact predictions for all single-residue
variants of the query sequence::
 
 python3 -m demask.predict -i myquery_homologs.a2m -o myquery_predictions.txt

Full options for finding homologs:

.. code-block:: none

   usage: homologs.py [-h] -s SEQFILE [-c FILE] --blastp BLASTP --db PATH
                     [-t THREADS] [-e EVALUE_CUTOFF] [-b BITSCORE_CUTOFF]
                     [-n NSEQS] (-o OUTFILE | -d OUTDIR)

   -h, --help            show this help message and exit
   -s SEQFILE, --seqfile SEQFILE
                           Name of a fasta file with one or more query sequences.
   -c FILE, --config FILE
                           Configuration file (e.g. the same configuration file
                           used for running DeMaSk). Defaults to 'config.ini' in
                           the demask directory.
   --blastp BLASTP       The full path of your blastp executable.
   --db PATH             Path of the BLAST database in which supporting
                           sequences will be searched. Must be formatted using
                           the makeblastdb program.
   -t THREADS, --threads THREADS
                           Number of threads (CPUs) blastp will use. Defaults to
                           1.
   -e EVALUE_CUTOFF, --evalue_cutoff EVALUE_CUTOFF
                           E-value threshold for blastp hits. Defaults to 1e-5.
   -b BITSCORE_CUTOFF, --bitscore_cutoff BITSCORE_CUTOFF
                           Bits per query residue threshold for blastp hits.
                           Default is no threshold
   -n NSEQS, --nseqs NSEQS
                           Maximum number of top hits to include per query
                           sequence. Defaults to 500.
   -o OUTFILE, --outfile OUTFILE
                           Name of output file (only used for single-query runs).
   -d OUTDIR, --outdir OUTDIR
                           Name of directory to write output files in a2m (FASTA)
                           format. The sequence titles (after the '>' and before
                           any whitespace) will be used for file names.

   For each query sequence, an output file in a2m (FASTA) format will be produced
   containing the query sequence and the most similar homologs up to a specified
   number (default 500).

Full options for getting DeMaSk predictions:

.. code-block:: none

   usage: predict.py [-h] -i INFILE [-o OUTFILE] [-c FILE] [-m FILE]
                     [--coefs FILE] [-n NSEQS] [-w WEIGHT_THRESHOLD]

   -h, --help            show this help message and exit
   -i INFILE, --infile INFILE
                           Name of the file containing a sequence alignment in
                           A2M (FASTA) format, with the query protein as the
                           first sequence.
   -o OUTFILE, --outfile OUTFILE
                           Name of new file to write scores to.
   -c FILE, --config FILE
                           Configuration file. Defaults to 'config.ini' in the
                           demask directory.
   -m FILE, --matrix FILE
                           File containing the directional substitution matrix.
                           Defaults to the file at 'DeMaSk/data/matrix.txt'.
   --coefs FILE          File containing the intercept, entropy, log2f_var, and
                           matrix coefficients. Defaults to the file at
                           'DeMaSk/data/coefficients.txt'.
   -n NSEQS, --nseqs NSEQS
                           Maximum number of supporting sequences per query
                           sequence. Defaults to 500.
   -w WEIGHT_THRESHOLD, --weight_threshold WEIGHT_THRESHOLD
                           Sequence identity threshold used for sequence
                           weighting, e.g. 0.8. Sequences are weighted by the
                           inverse of the number of sequences within this percent
                           identity. If None (default), sequence weighting is not
                           used.

Running in Python
^^^^^^^^^^^^^^^^^

Corresponding functions can be run in Python code:

.. autofunction:: demask.homologs.find_homologs

.. autofunction:: demask.predict.run_demask

User-generated matrix and coefficients
--------------------------------------

DeMaSk comes with a directional substitution matrix computed from a
collection of deep mutational scanning datasets, as well as
corresponding linear model coefficients.  Additional commands are
included in case you want to fit the model to a custom matrix, or even
calculate a matrix from a custom data collection and then fit the
linear model to it.

For example, the default matrix was generated like this::

 python3 -m demask.matrix -d DeMaSk/data/datasets -o DeMaSk/data/matrix.txt

Then, linear model coefficients were calculated::

 python3 -m demask.fit -d DeMaSk/data/datasets -a DeMaSk/data/alignments -m DeMaSk/data/matrix.txt -o DeMaSk/data/coefficients.txt

Full options for computing a matrix:

.. code-block:: none

   usage: matrix.py [-h] -d DATADIR -o OUTFILE [--columns COLUMNS] [-c FILE]

   -h, --help            show this help message and exit
   -d DATADIR, --datadir DATADIR
                           Name of directory containing (only) DMS data files,
                           one per protein. If a file name is supplied, only that
                           file will be used.
   -o OUTFILE, --outfile OUTFILE
                           Name of file to write matrix to.
   --columns COLUMNS     An optional comma-separated list of column names to
                           read in place of 'pos', 'WT', 'var', and 'score',
                           respectively. For example,
                           'Position,wt,mutant,Fitness'. Must apply to all files
                           provided.
   -c FILE, --config FILE
                           Configuration file. Defaults to 'config.ini' in the
                           DeMaSk directory.

   Loads a collection of DMS datasets and averages scores by WT + variant
   identity. Data should already be normalized as desired. For the default DeMaSk
   matrix, variant scores were rank-normalized per protein, and then the wild-
   type fitness levels subtracted, so that scores represent delta fitness.
   Dataset files must be tab-delimited with column names. Columns containing
   residue position, WT residue, variant residue, and score must be present in
   any order, and other columns will be ignored. Entries for which the WT and
   variant residues are not among the 20 canonical amino acids or are the same
   are ignored.

Full options for calculating coefficients:

.. code-block:: none

   usage: fit.py [-h] -d DATADIR -a ALNDIR [-m FILE] -o OUTFILE
               [--columns COLUMNS] [-c FILE] [-n NSEQS] [-w WEIGHT_THRESHOLD]

   -h, --help            show this help message and exit
   -d DATADIR, --datadir DATADIR
                           Name of directory containing (only) DMS data files,
                           one per protein.
   -a ALNDIR, --alndir ALNDIR
                           Name of the directory containing alignment files,
                           eachof which must be named as the basename of a DMS
                           data filefollowed by the '.a2m' extension.
   -m FILE, --matrix FILE
                           File containing the directional substitution matrix.
                           Defaults to the file at 'demask/matrix/matrix.txt'.
   -o OUTFILE, --outfile OUTFILE
                           Name of the file in which to save the coefficients.
   --columns COLUMNS     An optional comma-separated list of DMS data column
                           names to read in place of 'pos', 'WT', 'var', and
                           'score', respectively. For example,
                           'Position,wt,mutant,Fitness'. Must apply to all files
                           provided.
   -c FILE, --config FILE
                           Configuration file. Defaults to 'config.ini' in the
                           demask directory.
   -n NSEQS, --nseqs NSEQS
                           Maximum number of supporting sequences per query
                           sequence. Defaults to 500.
   -w WEIGHT_THRESHOLD, --weight_threshold WEIGHT_THRESHOLD
                           Sequence identity threshold used for sequence
                           weighting, e.g. 0.8. Sequences are weighted by the
                           inverse of the number of sequences within this percent
                           identity. If None (default), sequence weighting is not
                           used.

   For a set of DMS datasets, load the DMS scores, aligned homologs, and pre-
   computed substitution matrix, and fit a linear model, saving the coefficients
   to a file for later prediction on new proteins. Dataset files must be tab-
   delimited with column names. Columns containing residue position, WT residue,
   variant residue, and score must be present in any order, and other columns
   will be ignored. Entries for which the WT and variant residues are not among
   the 20 canonical amino acids or are the same are ignored. DMS scores should
   already be normalized as desired. For the default DeMaSk matrix, variant
   scores were rank-normalized per protein, and then the wild-type fitness levels
   subtracted, so that scores represent delta fitness.

Running in Python
^^^^^^^^^^^^^^^^^

Corresponding functions can be run in Python code:

.. autofunction:: demask.matrix.prepare_matrix

.. autofunction:: demask.fit.fit_model

Full documentation
==================

You won't need to use most of these functions directly, but they are
described here for custom use and curiosity.

.. automodule:: demask.matrix
   :members:

.. automodule:: demask.homologs
   :members:

.. automodule:: demask.profiles
   :members:

.. automodule:: demask.fit
   :members:

.. automodule:: demask.predict
   :members:
      
.. automodule:: demask.utils
   :members:

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
