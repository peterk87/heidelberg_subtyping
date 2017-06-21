# -*- coding: utf-8 -*-
import os

import attr
import logging

from _datetime import datetime

import re
from pkg_resources import resource_filename

from . import program_name
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

@attr.s
class Subtype(object):
    sample = attr.ib(validator=attr.validators.instance_of(str))
    file_path = attr.ib(validator=attr.validators.instance_of(str))
    subtype = attr.ib(default=None, validator=attr.validators.optional(attr.validators.instance_of(str)))
    all_subtypes = attr.ib(default=None, validator=attr.validators.optional(attr.validators.instance_of(str)))
    inconsistent_subtypes = attr.ib(default=None, validator=attr.validators.optional(attr.validators.instance_of(str)))
    tiles_matching_subtype = attr.ib(default=None, validator=attr.validators.optional(attr.validators.instance_of(str)))
    are_subtypes_consistent = attr.ib(default=True, validator=attr.validators.instance_of(bool))
    n_tiles_matching_all = attr.ib(default=0, validator=attr.validators.instance_of(int))
    n_tiles_matching_positive = attr.ib(default=0, validator=attr.validators.instance_of(int))
    n_tiles_matching_subtype = attr.ib(default=0, validator=attr.validators.instance_of(int))


def compare_subtypes(a, b):
    for x, y in zip(a, b):
        if x != y:
            return False
    return True

def find_inconsistent_subtypes(subtypes):
    from collections import Counter
    incon = []
    for i in range(len(subtypes) - 1):
        a = subtypes[i]
        for j in range(i + 1, len(subtypes)):
            b = subtypes[j]
            is_consistent = compare_subtypes(a, b)
            if not is_consistent:
                incon.append((a, b))
    l = []
    for a,b in incon:
        astr = '.'.join([str(x) for x in a])
        bstr = '.'.join([str(x) for x in b])
        l += [astr, bstr]
    c = Counter(l)
    incon_subtypes = []
    for subtype,freq in c.most_common():
        if freq > 1:
            incon_subtypes.append(subtype)
        else:
            break
    return incon_subtypes


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

    df.rename(columns={'qseqid':'tilename'}, inplace=True)

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
    st.n_tiles_matching_all = df.shape[0]
    st.n_tiles_matching_positive = dfpos.shape[0]
    st.n_tiles_matching_subtype = dfpos_highest_res.shape[0]
    pos_subtypes_str = [x for x in dfpos.subtype.unique()]
    pos_subtypes_str.sort(key=lambda x: len(x))
    st.all_subtypes = '; '.join(pos_subtypes_str)
    st.subtype = '; '.join([x for x in dfpos_highest_res.subtype.unique()])
    st.tiles_matching_subtype = '; '.join([x for x in dfpos_highest_res.tilename])

    if len(inconsistent_subtypes) > 0:
        st.are_subtypes_consistent = False
        st.inconsistent_subtypes = inconsistent_subtypes

    logging.info(st)

    df['sample'] = genome_name
    df['file_path'] = fasta_path
    return st, df
