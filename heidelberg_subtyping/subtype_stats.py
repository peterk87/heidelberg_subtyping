# -*- coding: utf-8 -*-

import attr
from collections import defaultdict

from pkg_resources import resource_filename

from .blast_wrapper.helpers import parse_fasta

TILES_FASTA = resource_filename('heidelberg_subtyping', 'data/tiles.fasta')


@attr.s
class SubtypeCounts:
    subtype = attr.ib()
    subtype_tile_count = attr.ib(default=0, validator=attr.validators.instance_of(int))
    positive_tile_count = attr.ib(default=0, validator=attr.validators.instance_of(int))
    all_tile_count = attr.ib(default=0, validator=attr.validators.instance_of(int))

    @subtype.validator
    def _check_subtype(self, attribute, value):
        if value is None or value == '':
            raise ValueError('Subtype cannot be None or empty string')
        if len(value) > 1:
            if '.' not in value:
                raise ValueError(
                    'Invalid subtype specified! "{}" does not numbers delimited by "." (periods)'.format(value))

    @subtype_tile_count.validator
    def _check_subtype_tile_count(self, attribute, value):
        if value == 0:
            raise ValueError('Subtype tile count cannot be zero!')

    @positive_tile_count.validator
    def _check_positive_tile_count(self, attribute, value):
        if value < self.subtype_tile_count:
            raise ValueError(
                'Subtype {}: Number of all subtype positive tiles (n={}) cannot be less than the number of subtype specific tiles (n={})'.format(
                    self.subtype,
                    value,
                    self.subtype_tile_count))
        if value > self.all_tile_count:
            raise ValueError(
                'Subtype {}: Number of all subtype positive tiles (n={}) cannot exceed number of all matching tiles for subtype (n={})'.format(
                    self.subtype,
                    value,
                    self.all_tile_count))
        if value == 0:
            raise ValueError('Subtype {}: Number of all subtype positive tiles cannot be zero!'.format(self.subtype))


def _tiles():
    tiles = defaultdict(list)
    neg_tiles = defaultdict(list)
    for h, s in parse_fasta(TILES_FASTA):
        _, st = h.split('-')
        if 'negative' not in h:
            tiles[st].append(h)
        else:
            neg_tiles[st].append(h)
    return tiles, neg_tiles


def subtype_counts():
    subtype_counts = {}
    tiles, neg_tiles = _tiles()
    n_tiles_total = sum([len(vs) for vs in tiles.values()])
    ks = [x for x in tiles.keys()]
    ks.sort()
    for k in ks:
        st_pos_count_rest = 0
        subtypes_set = {k}
        if len(k) > 1:
            *r, _ = k.split('.')
            for i in range(len(r)):
                sub_k = '.'.join(r[:i + 1])
                subtypes_set.add(sub_k)
                pos_tiles = tiles[sub_k]
                st_pos_count_rest += len(pos_tiles)
        st_count = len(tiles[k])
        st_count_pos = st_pos_count_rest + st_count
        st_neg_count = sum([len(v) for k,v in neg_tiles.items() if k not in subtypes_set])
        subtype_count = SubtypeCounts(subtype=k,
                                      subtype_tile_count=st_count,
                                      positive_tile_count= st_count_pos,
                                      all_tile_count=st_neg_count + st_count_pos)
        subtype_counts[k] = subtype_count
    return subtype_counts


SUBTYPE_COUNTS = subtype_counts()