# -*- coding: utf-8 -*-

import os
import re
import attr
import logging
from datetime import datetime
from pkg_resources import resource_filename

from . import program_name
from .subtype import Subtype
from .utils import find_inconsistent_subtypes
from .kmer_count import Jellyfisher
from .blast_wrapper import BlastRunner, BlastReader

TILES_FASTA = resource_filename('heidelberg_subtyping', 'data/tiles.fasta')

SUBTYPE_SUMMARY_COLS ="""
sample
subtype
all_subtypes
tiles_matching_subtype
are_subtypes_consistent
inconsistent_subtypes
n_tiles_matching_all
n_tiles_matching_positive
n_tiles_matching_subtype
file_path""".strip().split('\n')


def subtype_fasta(fasta_path, genome_name, tmp_dir='/tmp'):
    dtnow = datetime.now()
    genome_name_no_spaces = re.sub(r'\W', '_', genome_name)
    genome_tmp_dir = os.path.join(tmp_dir, dtnow.strftime("%Y%m%d%H%M%S") + '-' + program_name + '-' + genome_name_no_spaces)
    with BlastRunner(fasta_path=fasta_path, tmp_work_dir=genome_tmp_dir) as brunner:
        blast_outfile = brunner.blast_against_query(TILES_FASTA, word_size=33)
        with BlastReader(blast_outfile=blast_outfile) as breader:
            df = breader.parse()

    st = Subtype(sample=genome_name, file_path=fasta_path)
    if df is None or df.shape[0] == 0:
        logging.warning('No Heidelberg subtyping tile matches for "%s"', fasta_path)
        st.are_subtypes_consistent = False
        return st, None

    df.rename(columns={'qseqid':'tilename',
                       'sseq':'seq'},
              inplace=True)

    refpositions = [x for x, y in df.tilename.str.split('-')]
    subtypes = [y for x, y in df.tilename.str.split('-')]
    df['refposition'] = refpositions
    df['subtype'] = subtypes
    df['is_pos_tile'] = ~df.tilename.str.contains('negative')
    logging.debug('df: %s', df)
    dfpos = df[df.is_pos_tile]
    logging.debug('dfpos: %s', dfpos)
    subtype_lens = dfpos.subtype.apply(len)
    max_subtype_strlen = subtype_lens.max()
    logging.debug('max substype str len: %s', max_subtype_strlen)
    dfpos_highest_res = dfpos[subtype_lens == max_subtype_strlen]
    logging.debug('dfpos_highest_res: %s', dfpos_highest_res)
    pos_subtypes = [[int(y) for y in x.split('.')] for x in dfpos.subtype.unique()]
    pos_subtypes.sort(key=lambda a: len(a))
    logging.debug('pos_subtypes: %s', pos_subtypes)
    inconsistent_subtypes = find_inconsistent_subtypes(pos_subtypes)
    logging.debug('inconsistent_subtypes: %s', inconsistent_subtypes)
    st.n_tiles_matching_all = df.tilename.unique().size
    st.n_tiles_matching_positive = dfpos.tilename.unique().size
    st.n_tiles_matching_subtype = dfpos_highest_res.tilename.unique().size
    pos_subtypes_str = [x for x in dfpos.subtype.unique()]
    pos_subtypes_str.sort(key=lambda x: len(x))
    st.all_subtypes = '; '.join(pos_subtypes_str)
    st.subtype = '; '.join([x for x in dfpos_highest_res.subtype.unique()])
    st.tiles_matching_subtype = '; '.join([x for x in dfpos_highest_res.tilename.unique()])

    if len(inconsistent_subtypes) > 0:
        st.are_subtypes_consistent = False
        st.inconsistent_subtypes = inconsistent_subtypes

    logging.info(st)

    df['sample'] = genome_name
    df['file_path'] = fasta_path
    return st, df


def subtype_reads(reads, genome_name, tmp_dir='/tmp', threads=1, min_kmer_freq=10, max_kmer_freq=200):
    dtnow = datetime.now()
    genome_name_no_spaces = re.sub(r'\W', '_', genome_name)
    genome_tmp_dir = os.path.join(tmp_dir,
                                  dtnow.strftime("%Y%m%d%H%M%S") + '-' + program_name + '-' + genome_name_no_spaces)
    with Jellyfisher(genome_name=genome_name,
                       reads=reads,
                       min_kmer_freq=min_kmer_freq,
                       max_kmer_freq=max_kmer_freq,
                       tmp_dir=genome_tmp_dir,
                       threads=threads) as jfer:
        st, df = jfer.summary()
        return st, df