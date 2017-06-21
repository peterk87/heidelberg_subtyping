#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import attr
import logging
import re

from . import program_name, program_desc, __version__
from .subtyper import subtype_fasta, SUBTYPE_SUMMARY_COLS
from .utils import genome_name_from_fasta_path

LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'


def init_console_logger(logging_verbosity=3):
    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    if logging_verbosity > (len(logging_levels) - 1):
        logging_verbosity = 3
    lvl = logging_levels[logging_verbosity]

    logging.basicConfig(format=LOG_FORMAT, level=lvl)


def init_parser():
    parser = argparse.ArgumentParser(prog=program_name,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=program_desc)
    parser.add_argument('fastas',
                        metavar='F',
                        nargs='*',
                        help='Input genome FASTA file')
    parser.add_argument('-i', '--input-fasta-genome-name',
                        nargs=2,
                        metavar=('fasta_path', 'genome_name'),
                        action='append',
                        help='fasta file path to genome name pair')
    parser.add_argument('-D', '--input-directory',
                        help='directory of input fasta files (.fasta|.fa|.fna)')
    parser.add_argument('-o', '--output-summary',
                        help='Subtyping summary output path (tab-delimited)')
    parser.add_argument('-O', '--output-tile-results',
                        help='Subtyping tile matching output path (tab-delimited)')
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=1,
                        help='Number of parallel threads to run analysis (default=1)')
    parser.add_argument('-T', '--tmp-dir',
                        default='/tmp',
                        help='Base temporary working directory for intermediate analysis files')
    parser.add_argument('-K', '--keep-tmp',
                        default=False,
                        action='store_true',
                        help='Keep temporary analysis files')
    parser.add_argument('-v', '--verbose',
                        action='count',
                        default=0,
                        help='Logging verbosity level (-v == show warnings; -vvv == show debug info)')
    parser.add_argument('-V', '--version',
                        action='version',
                        version='%(prog)s {}'.format(__version__))
    return parser


def main():
    parser = init_parser()
    args = parser.parse_args()
    init_console_logger(args.verbose)
    output_summary_path = args.output_summary
    output_tile_results = args.output_tile_results

    input_genomes = []
    logging.debug(args)
    if args.fastas:
        logging.debug('# of input fastas %s', len(args.fastas))
        for fasta_path in args.fastas:
            fasta_path = os.path.abspath(fasta_path)
            if os.path.exists(fasta_path):
                genome_name = genome_name_from_fasta_path(fasta_path)
                input_genomes.append((fasta_path, genome_name))
            else:
                logging.error('Input fasta "%s" does not exist!', fasta_path)

    if args.input_fasta_genome_name:
        for fasta_path, genome_name in args.input_fasta_genome_name:
            input_genomes.append((os.path.abspath(fasta_path), genome_name))

    if args.input_directory:
        logging.info('Searching dir "%s" for fasta files', args.input_directory)
        for x in os.listdir(args.input_directory):
            full_file_path = os.path.abspath(os.path.join(args.input_directory, x))
            if os.path.isfile(full_file_path) and re.match(r'^.+\.(fasta|fa|fna)$', x):
                genome_name = genome_name_from_fasta_path(full_file_path)
                input_genomes.append((full_file_path, genome_name))

    if len(input_genomes) == 0:
        raise Exception('No input files specified!')

    n_threads = args.threads
    tmp_dir = args.tmp_dir
    if n_threads == 1:
        logging.info('Serial single threaded run mode on %s input genomes', len(input_genomes))
        outputs = [subtype_fasta(input_fasta, genome_name, tmp_dir=tmp_dir) for input_fasta, genome_name in
                   input_genomes]
    else:
        from multiprocessing import Pool
        logging.info('Initializing thread pool with %s threads', n_threads)
        pool = Pool(processes=n_threads)
        logging.info('Running analysis asynchronously on %s input genomes', len(input_genomes))
        res = [pool.apply_async(subtype_fasta, (input_fasta, genome_name, tmp_dir)) for
               input_fasta, genome_name in input_genomes]

        logging.info('Parallel analysis complete! Retrieving analysis results')
        outputs = [x.get() for x in res]

    import pandas as pd

    subtype_results = []
    dfs = []
    for subtype, df in outputs:
        if df is not None:
            dfs.append(df)
        subtype_results.append(attr.asdict(subtype))
    dfall = pd.concat(dfs)
    dfsummary = pd.DataFrame(subtype_results)
    dfsummary = dfsummary[SUBTYPE_SUMMARY_COLS]

    if output_summary_path:
        dfsummary.to_csv(output_summary_path, sep='\t', index=None)
        logging.info('Wrote subtyping output summary to %s', output_summary_path)
    else:
        print(dfsummary.to_csv(sep='\t'))

    if output_tile_results:
        dfall.to_csv(output_tile_results, sep='\t', index=None)


if __name__ == '__main__':
    main()
