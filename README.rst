*****************************************************************
``heidelberg_subtyping``: Subtype *Salmonella* Heidelberg genomes
*****************************************************************

Subtype *Salmonella enterica* subsp. enterica serovar Heidelberg genomes using an *in-silico* 33 bp k-mer subtyping method developed by Genevieve Labbe, Roger Johnson, et al from the Public Health Agency of Canada.

Citation
========

If you find this tool useful, please cite as:

.. epigraph::

    [article citation goes here]


Requirements and Dependencies
=============================

This tool has only been tested on Linux (specifically Arch Linux). It may or may not work on OSX.

These are the external dependencies required for ``heidelberg_subtyping``:

- Python (>=v3.5)
- BLAST+ (>=v2.6)


Installation
============

Install BLAST+ using your package manager

.. code-block:: bash

    pip install heidelberg_subtyping


Usage
=====

If you run ``heidelberg_subtyping -h``, you should see the following usage statement:

.. code-block:: none

    usage: heidelberg_subtyping [-h] [-i fasta_path genome_name]                                                           
                                [-D INPUT_DIRECTORY] [-o OUTPUT_SUMMARY]                                                   
                                [-O OUTPUT_TILE_RESULTS] [-t THREADS] [-T TMP_DIR]                                         
                                [-K] [-v] [-V]                 
                                [F [F ...]]                    

    Subtype Salmonella Heidelberg genomes using a 33bp k-mer typing scheme                                                 
    Developed by Genevieve Labbe, Roger Johnson, PHAC-NML Guelph                                                           

    positional arguments:                                      
      F                     Input genome FASTA file            

    optional arguments:                                        
      -h, --help            show this help message and exit    
      -i fasta_path genome_name, --input-fasta-genome-name fasta_path genome_name                                          
                            fasta file path to genome name pair                                                            
      -D INPUT_DIRECTORY, --input-directory INPUT_DIRECTORY    
                            directory of input fasta files (.fasta|.fa|.fna)                                               
      -o OUTPUT_SUMMARY, --output-summary OUTPUT_SUMMARY       
                            Subtyping summary output path (tab-delimited)                                                  
      -O OUTPUT_TILE_RESULTS, --output-tile-results OUTPUT_TILE_RESULTS                                                    
                            Subtyping tile matching output path (tab-delimited)                                            
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

Analysis of a single FASTA file:

.. code-block:: bash

    heidelberg_subtyping -o results.tab -O match_results.tab /path/to/SRR1002850.fasta


Contents of ``results.tab``:

.. code-block:: none

    sample	subtype	all_subtypes	tiles_matching_subtype	are_subtypes_consistent	inconsistent_subtypes	n_tiles_matching_all	n_tiles_matching_positive	n_tiles_matching_subtype	file_path
    SRR1002850	2.2.2.2.1.4	2; 2.2; 2.2.2; 2.2.2.2; 2.2.2.2.1; 2.2.2.2.1.4	1037658-2.2.2.2.1.4; 3785187-2.2.2.2.1.4; 2154958-2.2.2.2.1.4	True		212	17	3	/path/to/SRR1002850.fasta


Contents of ``match_results.tab``:

.. code-block:: none

    tilename	stitle	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore	qlen	slen	sseq	coverage	is_trunc	refposition	subtype	is_pos_tile	sample	file_path
    775920-2.2.2.2	NODE_3_length_511571_cov_26.9963_ID_5	100.0	33	0	0	1	33	475240	475272	1.5e-11	62.1	33	511571	GTTCAGGTGCTACCGAGGATCGTTTTTGGTGCG	1.0	False	775920	2.2.2.2	True	SRR1002850	/path/to/SRR1002850.fasta
    negative3113857-1.2	NODE_4_length_474326_cov_28.1591_ID_7	100.0	33	0	0	1	33	84804	84836	1.5e-11	62.1	33	474326	TTCATGACGTCATCCCAGTCTTTTTCCGTGAAA	1.0	False	negative3113857	1.2	False	SRR1002850	/path/to/SRR1002850.fasta
    negative3159204-2.2.1.1.3	NODE_4_length_474326_cov_28.1591_ID_7	100.0	33	0	0	1	33	130145	130177	1.5e-11	62.1	33	474326	CCGCCTCGCCAACCTGCGGCGGAGTCGCGAGCT	1.0	False	negative3159204	2.2.1.1.3	False	SRR1002850	/path/to/SRR1002850.fasta



License
=======

Copyright 2017 Public Health Agency of Canada

Distributed under the GNU Public License version 3.0
