# -*- coding: utf-8 -*-
import shutil

import attr
import os
from datetime import datetime
import logging
import pandas as pd
from pkg_resources import resource_filename
from ..utils import exc_exists, run_command, find_inconsistent_subtypes
from ..blast_wrapper.helpers import parse_fasta, revcomp
from ..subtype import Subtype


@attr.s
class Jellyfisher(object):
    genome_name = attr.ib()
    reads = attr.ib()
    min_kmer_freq = attr.ib(default=10, validator=attr.validators.instance_of(int))
    max_kmer_freq = attr.ib(default=200, validator=attr.validators.instance_of(int))
    tmp_dir = attr.ib(default='/tmp', validator=attr.validators.instance_of(str))
    threads = attr.ib(default=1, validator=attr.validators.instance_of(int))
    kmer_fasta = attr.ib(default=resource_filename('heidelberg_subtyping', 'data/tiles.fasta'),
                         validator=attr.validators.instance_of(str))
    jellyfish_exc = attr.ib(default='jellyfish', validator=attr.validators.instance_of(str))
    jf_file = attr.ib(default=None)
    jf_query_tiles_file = attr.ib(default=None)
    df_results = attr.ib(default=None, validator=attr.validators.optional(attr.validators.instance_of(pd.DataFrame)))

    @min_kmer_freq.validator
    def _check_min_kmer_freq(self, attribute, value):
        if value <= 0:
            raise ValueError('Cannot have zero or negative k-mer frequency threshold!')
        if value > self.max_kmer_freq:
            raise ValueError('Min k-mer freq threshold must be less than max k-mer freq threshold!')

    @max_kmer_freq.validator
    def _check_max_kmer_freq(self, attribute, value):
        if value <= 0:
            raise ValueError('Cannot have zero or negative k-mer frequency threshold!')
        if value < self.min_kmer_freq:
            raise ValueError('Max k-mer freq threshold must be greater than min k-mer freq threshold!')

    @jellyfish_exc.validator
    def _check_jellyfish_exists(self, attribute, value):
        if not exc_exists(value):
            raise OSError('{} does not exist in the user $PATH'.format(value))

    @reads.validator
    def _reads_exist(self, attribute, value):
        if isinstance(value, str):
            if not os.path.exists(value):
                raise FileNotFoundError('Reads file {} does not exist!'.format(value))
            if not os.path.isfile(value):
                raise OSError('{} is not a valid reads file'.format(value))
        elif isinstance(value, list):
            for x in value:
                if not isinstance(x, str):
                    raise Exception(
                        'Reads file not specified as string or list of string: type={} "{}"'.format(type(x), x))
                if not os.path.exists(x):
                    raise FileNotFoundError('Reads file {} does not exist!'.format(x))
                if not os.path.isfile(x):
                    raise OSError('{} is not a valid reads file'.format(x))
        else:
            raise Exception(
                'Reads file(s) not specified as string or list of string: type={} "{}"'.format(type(value), value))

    def _create_tmp_folder(self):
        count = 1
        tmp_dir = self.tmp_dir
        while True:
            try:
                logging.info('Trying to create analysis directory at: %s', tmp_dir)
                os.makedirs(tmp_dir)
                break
            except OSError as e:
                logging.warning('Error on creation of tmp analysis directory "{}"! {}'.format(
                    tmp_dir,
                    e
                ))
                tmp_dir = '{}_{}'.format(self.tmp_dir, count)
                count += 1
        self.tmp_dir = tmp_dir
        return self.tmp_dir

    def cleanup(self):
        shutil.rmtree(self.tmp_dir)

    def kmer_count(self, k=33):
        """Jellyfish k-mer counting of FASTQ reads

        Args:
            k (int): k-mer size

        Returns:
            str: Jellyfish k-mer counts outfile path
        """
        if not os.path.exists(self.tmp_dir):
            self._create_tmp_folder()
        timestamp = '{:%Y%b%d_%H_%M_%S}'.format(datetime.now())
        self.jf_file = os.path.join(self.tmp_dir, '{}-{}.jf'.format(timestamp, self.genome_name))
        cmd_list = [self.jellyfish_exc,
                    'count',
                    '-m', '{}'.format(k),  # kmer size
                    '-s', '100M',  # initial hash size
                    '-C',  # Count both strand, canonical representation
                    '-t', '{}'.format(self.threads),
                    '-o', self.jf_file,
                    ]
        if isinstance(self.reads, list):
            cmd_list += self.reads
        elif isinstance(self.reads, str):
            cmd_list.append(self.reads)

        logging.info('Running Jellyfish k-mer counting with command: %s', ' '.join(cmd_list))
        exit_code, stdout, stderr = run_command(cmd_list)
        if exit_code == 0:
            if os.path.exists(self.jf_file):
                logging.info('Jellyfish k-mer counts written to %s', self.jf_file)
                return self.jf_file
            else:
                raise Exception('Jellyfish count file does not exist at %s', self.jf_file)
        else:
            raise Exception('Could not run "jellyfish count"! Exit code {}; stderr: {}'.format(
                exit_code,
                stderr))

    def kmer_query(self):
        """Jellyfish query k-mers against Jellyfish k-mer counts database

        Returns:
            str: Jellyfish k-mer counts outfile path
        """

        if self.jf_file is None:
            self.kmer_count()

        cmd_list = [self.jellyfish_exc,
                    'query',
                    '-s', self.kmer_fasta,
                    self.jf_file,
                    ]

        logging.info('Querying tiles.fasta against Jellyfish k-mer count DB. Running command: %s', ' '.join(cmd_list))
        exit_code, stdout, stderr = run_command(cmd_list)
        if exit_code == 0:
            self.ran_query = True
            if stdout is None or stdout == '':
                return None
            self.jf_query_tiles_file = self.jf_file + '-vs-tiles.fasta.txt'
            with open(self.jf_query_tiles_file, 'w') as fout:
                fout.write(stdout)
            return self.jf_query_tiles_file
        else:
            raise Exception('Could not run "jellyfish query"! Exit code {}; stderr: {}'.format(
                exit_code,
                stderr))

    def _reads_to_str(self):
        reads_str = ''
        if isinstance(self.reads, str):
            reads_str = self.reads
        elif isinstance(self.reads, (list, tuple)):
            reads_str = '; '.join(self.reads)
        else:
            reads_str = '{}'.format(self.reads)
        return reads_str

    def parse_query(self):
        jf_query_file = self.jf_query_tiles_file
        if jf_query_file is None:
            jf_query_file = self.kmer_query()

        if jf_query_file is None:
            return None

        df = pd.read_table(self.jf_query_tiles_file, sep=' ')
        df.columns = ['seq', 'freq']
        df = df[df.freq > 0]
        if df.shape[0] == 0:
            return None
        seq_tiles = {s: h for h, s in parse_fasta(self.kmer_fasta)}
        tiles = []
        for x in df.seq:
            if x in seq_tiles:
                tiles.append(seq_tiles[x])
                continue
            x_rc = revcomp(x)
            if x_rc in seq_tiles:
                tiles.append(seq_tiles[x_rc])
                continue
            logging.error('%s/%s not found in tiles.fasta!', x, x_rc)
        df['sample'] = self.genome_name

        df['file_path'] = self._reads_to_str()
        df['tilename'] = tiles
        df['is_pos_tile'] = [not x.startswith('negative') for x in tiles]
        df['subtype'] = [y for x, y in df.tilename.str.split('-')]
        df['refposition'] = [x for x, y in df.tilename.str.split('-')]
        df['is_kmer_freq_okay'] = (df.freq >= self.min_kmer_freq) & (df.freq <= self.max_kmer_freq)
        logging.info('n=%s k-mers with freq within thresholds of %s and %s',
                     df.is_kmer_freq_okay.sum(),
                     self.min_kmer_freq,
                     self.max_kmer_freq)
        logging.info('n=%s k-mers with freq NOT within thresholds of %s and %s',
                     (~df.is_kmer_freq_okay).sum(),
                     self.min_kmer_freq,
                     self.max_kmer_freq)
        self.df_results = df
        return df

    def summary(self):
        if self.df_results is None:
            self.parse_query()
        df = self.df_results
        st = Subtype(sample=self.genome_name, file_path=self._reads_to_str())
        self.subtype = st
        if df is None or df.shape[0] == 0:
            logging.warning('No Heidelberg subtyping tile matches for "%s"', self.reads)
            st.are_subtypes_consistent = False
            return st, None
        dfgood = df[df.is_kmer_freq_okay]
        dfpos = dfgood[dfgood.is_pos_tile]
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
        st.n_tiles_matching_all = dfgood.shape[0]
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
        return st, df

    def __enter__(self):
        if self.genome_name is None \
                or self.reads is None:
            return self
        self._create_tmp_folder()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()
