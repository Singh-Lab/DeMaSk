DeMaSk
******

Prediction of amino acid substitution impact using homologs and a
directional substitution matrix.

DeMaSk predictions can be easily obtained for any protein sequence
using the `web tool <https://demask.princeton.edu>`_.  This package
can be downloaded for customized usage.  See the `full documentation
<https://demask.readthedocs.io>`_ for more detailed instructions.

Installation
============

Install from PyPI with pip:

  pip install demask

To install from GitHub, clone the repository or `download and unzip
<https://github.com/Singh-Lab/DeMaSk/archive/master.zip>`_.  To get
Python dependencies and be able to run or import the modules from any
directory, install with pip (or pip3 if pip = pip2)::

  pip install -e DeMaSk/

Unless you're supplying your own aligned homologs, you'll need to have
the ``blastp`` program.  It can be downloaded as part of `BLAST+
<https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`_.

The ``blastp`` step requires a formatted sequence database.  To use
UniRef90, which demask.princeton.edu uses, download the zipped fasta
file from `here <https://www.uniprot.org/downloads>`_, unzip, and use
``makeblastdb`` from BLAST+ to format the database.

To avoid having to specify the location of the ``blastp`` binary and
the database for every DeMaSk run, put them in a config file, e.g.::

  blastp=/usr/local/ncbi/blast/bin/blastp
  db=/path/to/database/uniref90.fasta

By default, DeMaSk will look for the config file
``DeMaSk/config.ini``, which can also contain any other command line
arguments, such as ``nseqs``, ``threads``, and ``matrix``.

Usage
=====

Run any of the command modules with ``-h`` to see all options, e.g.::

 python3 -m demask.homologs -h

Get predictions for a query sequence
------------------------------------

Once the demask package is installed, you can run run it from
anywhere.  If you don't have aligned homologs for your query yet, run
the demask.homologs module::
 
 python3 -m demask.homologs -s myquery.fa -o myquery_homologs.a2m

The above command will also produce a file myquery.blast.json
containing the intermediate blastp output.

Then, generate fitness impact predictions for all single-residue
variants of the query sequence::
 
 python3 -m demask.predict -i myquery_homologs.a2m -o myquery_predictions.txt

The output looks like this::

  pos   WT      var     score   entropy log2f_var       matrix
  1     M       A       -0.3019 1.0666  -21.9207        -0.2641
  1     M       C       -0.3074 1.0666  -21.9207        -0.2713
  1     M       D       -0.4134 1.0666  -21.9207        -0.4100
  1     M       E       -0.4036 1.0666  -21.9207        -0.3972
  1     M       F       -0.2183 1.0666  -4.5455         -0.2707
  1     M       G       -0.3828 1.0666  -21.9207        -0.3700
  ...

Corresponding functions can be run in Python code by importing
``demask.homologs.find_homologs`` and ``demask.predict.run_demask``.

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

Corresponding functions can be run in Python code by importing
``demask.matrix.prepare_matrix`` and ``demask.fit.fit_model``.
