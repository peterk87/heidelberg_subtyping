*****************************************************************
``heidelberg_subtyping``: Subtype *Salmonella* Heidelberg genomes
*****************************************************************

Subtype *Salmonella enterica* subsp. enterica serovar Heidelberg genomes using an *in-silico* 33 bp k-mer subtyping method developed by Genevieve Labbe et al.

Subtype *Salmonella* Heidelberg genome assemblies (FASTA files) and/or whole-genome sequencing reads (FASTQ files)!

Citation
========

If you find this tool useful, please cite as:

.. epigraph::

    A robust genotyping scheme for *Salmonella enterica* serovar Heidelberg clones circulating in North America.
    Geneviève Labbé, James Robertson, Peter Kruczkiewicz, Chad R. Laing, Kim Ziebell, Aleisha R. Reimer, Lorelee Tschetter, Gary Van Domselaar, Sadjia Bekal, Kimberley A. MacDonald, Linda Hoang, Linda Chui, Danielle Daignault, Durda Slavic, Frank Pollari, E. Jane Parmley, Philip Mabon, Elissa Giang, Lok Kan Lee, Jonathan Moffat, Marisa Rankin, Joanne MacKinnon, Roger Johnson, John H.E. Nash.
    [Manuscript in preparation]


Requirements and Dependencies
=============================

This tool has only been tested on Linux (specifically Arch Linux). It may or may not work on OSX.

These are the external dependencies required for ``heidelberg_subtyping``:

- Python (>=v3.5)
- BLAST+ (>=v2.6)
    + for subtyping of FASTA input (assemblies)
- JELLYFISH (v2.0) (http://www.cbcb.umd.edu/software/jellyfish/)
    + for subtyping of FASTQ input (raw reads)


Installation
============

Ensure that BLAST+ and/or JELLYFISH are installed and accessible in your ``$PATH``.

Install ``heidelberg_subtyping`` from PyPI:

.. code-block:: bash

    pip install heidelberg_subtyping

Or install the latest master branch version directly from Github:

.. code-block:: bash

    pip install git+https://github.com/peterk87/heidelberg_subtyping.git@master


Usage
=====

If you run ``heidelberg_subtyping -h``, you should see the following usage statement:

.. code-block:: none

    usage: heidelberg_subtyping [-h] [-p forward_reads reverse_reads]
                                [-i fasta_path genome_name] [-D INPUT_DIRECTORY]
                                [-o OUTPUT_SUMMARY] [-O OUTPUT_TILE_RESULTS]
                                [--min-kmer-freq MIN_KMER_FREQ]
                                [--max-kmer-freq MAX_KMER_FREQ] [-t THREADS]
                                [-T TMP_DIR] [-K] [-v] [-V]
                                [F [F ...]]

    Subtype Salmonella Heidelberg genomes using a 33bp k-mer typing scheme
    Developed by Genevieve Labbe, Roger Johnson, PHAC-NML Guelph

    positional arguments:
      F                     Input genome FASTA/FASTQ files

    optional arguments:
      -h, --help            show this help message and exit
      -p forward_reads reverse_reads, --paired-reads forward_reads reverse_reads
                            FASTQ paired-end reads
      -i fasta_path genome_name, --input-fasta-genome-name fasta_path genome_name
                            fasta file path to genome name pair
      -D INPUT_DIRECTORY, --input-directory INPUT_DIRECTORY
                            directory of input fasta files (.fasta|.fa|.fna) or
                            FASTQ files (paired FASTQ should have same basename
                            with "_\d\.(fastq|fq)" postfix to be automatically
                            paired)
      -o OUTPUT_SUMMARY, --output-summary OUTPUT_SUMMARY
                            Subtyping summary output path (tab-delimited)
      -O OUTPUT_TILE_RESULTS, --output-tile-results OUTPUT_TILE_RESULTS
                            Subtyping tile matching output path (tab-delimited)
      --min-kmer-freq MIN_KMER_FREQ
                            Min k-mer freq/coverage
      --max-kmer-freq MAX_KMER_FREQ
                            Max k-mer freq/coverage
      -t THREADS, --threads THREADS
                            Number of parallel threads to run analysis (default=1)
      -T TMP_DIR, --tmp-dir TMP_DIR
                            Base temporary working directory for intermediate
                            analysis files
      -K, --keep-tmp        Keep temporary analysis files
      -v, --verbose         Logging verbosity level (-v == show warnings; -vvv ==
                            show debug info)
      -V, --version         show program's version number and exit



Example Usage
=============

Analysis of a single FASTA file
-------------------------------

.. code-block:: bash

    heidelberg_subtyping -vv -o results.tab -O match_results.tab /path/to/SRR1002850.fasta


Contents of ``results.tab``:

.. code-block:: none

    sample      subtype      all_subtypes                                    tiles_matching_subtype                                         are_subtypes_consistent  inconsistent_subtypes  n_tiles_matching_all  n_tiles_matching_positive  n_tiles_matching_subtype  file_path
    SRR1002850  2.2.2.2.1.4  2; 2.2; 2.2.2; 2.2.2.2; 2.2.2.2.1; 2.2.2.2.1.4  1037658-2.2.2.2.1.4; 3785187-2.2.2.2.1.4; 2154958-2.2.2.2.1.4  True                                            212                   17                         3                         SRR1002850.fasta


Contents of ``match_results.tab``:

.. code-block:: none

    tilename                     stitle                                 pident  length  mismatch  gapopen  qstart  qend  sstart  send    evalue   bitscore  qlen  slen    seq                                coverage  is_trunc  refposition      subtype      is_pos_tile  sample      file_path
    775920-2.2.2.2               NODE_3_length_511571_cov_26.9963_ID_5  100.0   33      0         0        1       33    475240  475272  1.5e-11  62.1      33    511571  GTTCAGGTGCTACCGAGGATCGTTTTTGGTGCG  1.0       False     775920           2.2.2.2      True         SRR1002850  SRR1002850.fasta
    negative3113857-1.2          NODE_4_length_474326_cov_28.1591_ID_7  100.0   33      0         0        1       33    84804   84836   1.5e-11  62.1      33    474326  TTCATGACGTCATCCCAGTCTTTTTCCGTGAAA  1.0       False     negative3113857  1.2          False        SRR1002850  SRR1002850.fasta
    negative3159204-2.2.1.1.3    NODE_4_length_474326_cov_28.1591_ID_7  100.0   33      0         0        1       33    130145  130177  1.5e-11  62.1      33    474326  CCGCCTCGCCAACCTGCGGCGGAGTCGCGAGCT  1.0       False     negative3159204  2.2.1.1.3    False        SRR1002850  SRR1002850.fasta
    negative3187428-2.2.3.1.1    NODE_4_length_474326_cov_28.1591_ID_7  100.0   33      0         0        1       33    158369  158401  1.5e-11  62.1      33    474326  CTTTATCAGCGCGCAGTGTCCCATTCCATCATC  1.0       False     negative3187428  2.2.3.1.1    False        SRR1002850  SRR1002850.fasta
    negative3200083-2.1          NODE_4_length_474326_cov_28.1591_ID_7  100.0   33      0         0        1       33    171024  171056  1.5e-11  62.1      33    474326  ACCCGGTCTACCGCAAAATGGAAAGCGATATGC  1.0       False     negative3200083  2.1          False        SRR1002850  SRR1002850.fasta
    negative3204925-2.2.3.1.5    NODE_4_length_474326_cov_28.1591_ID_7  100.0   33      0         0        1       33    175866  175898  1.5e-11  62.1      33    474326  CTCGCTGGCAAGCAGTGCGGGTACTATCGGCGG  1.0       False     negative3204925  2.2.3.1.5    False        SRR1002850  SRR1002850.fasta
    negative3230678-2.2.2.1.1.1  NODE_4_length_474326_cov_28.1591_ID_7  100.0   33      0         0        1       33    201619  201651  1.5e-11  62.1      33    474326  AGCGGTGCGCCAAACCACCCGGAATGATGAGTG  1.0       False     negative3230678  2.2.2.1.1.1  False        SRR1002850  SRR1002850.fasta
    negative3233869-2.1.1.1.1    NODE_4_length_474326_cov_28.1591_ID_7  100.0   33      0         0        1       33    204810  204842  1.5e-11  62.1      33    474326  CAGCGCTGGTATGTGGCTGCACCATCGTCATTA  1.0       False     negative3233869  2.1.1.1.1    False        SRR1002850  SRR1002850.fasta
    negative3254229-2.2.3.1.3    NODE_4_length_474326_cov_28.1591_ID_7  100.0   33      0         0        1       33    225170  225202  1.5e-11  62.1      33    474326  CGCCACCACGCGGTTAGCGTCACGCTGACATTC  1.0       False     negative3254229  2.2.3.1.3    False        SRR1002850  SRR1002850.fasta


Analysis of a single FASTQ readset
----------------------------------

.. code-block:: bash

    heidelberg_subtyping -vv -t 4 -o results.tab -O match_results.tab -p SRR5646583_1.fastq SRR5646583_2.fastq


Contents of ``results.tab``:

.. code-block:: none

    sample      subtype      all_subtypes                                    tiles_matching_subtype                                         are_subtypes_consistent  inconsistent_subtypes  n_tiles_matching_all  n_tiles_matching_positive  n_tiles_matching_subtype  file_path
    SRR5646583  2.2.1.1.1.1  2; 2.2; 2.2.1; 2.2.1.1; 2.2.1.1.1; 2.2.1.1.1.1  1983064-2.2.1.1.1.1; 4211912-2.2.1.1.1.1; 4568600-2.2.1.1.1.1  True                                            212                   21                         3                         SRR5646583_1.fastq; SRR5646583_2.fastq


Contents of ``match_results.tab``:

.. code-block:: none

    seq                                freq  sample      file_path                                tilename         is_pos_tile  subtype      refposition        is_kmer_freq_okay
    ACGGTAAAAGAGGACTTGACTGGCGCGATTTGC  68    SRR5646583  SRR5646583_1.fastq;  SRR5646583_2.fastq  21097-2.2.1.1.1      True     2.2.1.1.1    21097              True
    AACCGGCGGTATTGGCTGCGGTAAAAGTACCGT  77    SRR5646583  SRR5646583_1.fastq;  SRR5646583_2.fastq  157792-2.2.1.1.1     True     2.2.1.1.1    157792             True
    CCGCTGCTTTCTGAAATCGCGCGTCGTTTCAAC  67    SRR5646583  SRR5646583_1.fastq;  SRR5646583_2.fastq  293728-2.2.1.1       True     2.2.1.1      293728             True
    GAATAACAGCAAAGTGATCATGATGCCGCTGGA  91    SRR5646583  SRR5646583_1.fastq;  SRR5646583_2.fastq  607438-2.2.1         True     2.2.1        607438             True
    CAGTTTTACATCCTGCGAAATGCGCAGCGTCAA  87    SRR5646583  SRR5646583_1.fastq;  SRR5646583_2.fastq  691203-2.2.1.1       True     2.2.1.1      691203             True
    CAGGAGAAAGGATGCCAGGGTCAACACGTAAAC  33    SRR5646583  SRR5646583_1.fastq;  SRR5646583_2.fastq  944885-2.2.1.1.1     True     2.2.1.1.1    944885             True
    GCGAACTGGCGAAACGCCTTGGCGTGGAACAAC  77    SRR5646583  SRR5646583_1.fastq;  SRR5646583_2.fastq  1047714-2.2.1.1.1    True     2.2.1.1.1    1047714            True
    ACAACACCGGGGTGGAGGCGCTGATTGTGCAGG  1     SRR5646583  SRR5646583_1.fastq;  SRR5646583_2.fastq  1697637-2.2.2.2.2.1  True     2.2.2.2.2.1  1697637            False
    GCCTGCGTTCAGTCGCTTGGGCGATATGCTGGA  65    SRR5646583  SRR5646583_1.fastq;  SRR5646583_2.fastq  1983064-2.2.1.1.1.1  True     2.2.1.1.1.1  1983064            True


Analysis of all FASTA/FASTQ files in a directory
------------------------------------------------

.. code-block:: bash

    heidelberg_subtyping -vv --threads <n_cpu> -o results.tab -O match_results.tab -D /path/to/fastas_or_fastqs/


``heidelberg_subtyping`` will only attempt to analyze the FASTA/FASTQ files within the specified directory and will not descend into any subdirectories!


License
=======

Copyright 2017 Public Health Agency of Canada

Distributed under the GNU Public License version 3.0
