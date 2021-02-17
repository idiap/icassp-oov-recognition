# Copyright (c) 2021 Idiap Research Institute, http://www.idiap.ch/
# Written by Rudolf A. Braun <rbraun@idiap.ch>
#
# This file is part of icassp-oov-recognition
#
# icassp-oov-recognition is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 3 as
# published by the Free Software Foundation.
#
# icassp-oov-recognition is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with icassp-oov-recognition. If not, see <http://www.gnu.org/licenses/>.

import wrappedfst
import os


def expand_fst(fst, to_expand, isyms):
    for state in fst.states():
        arc_iterator = wrappedfst.ArcIterator(fst, state)
        arcs_to_keep = []
        arcs_to_expand = []
        while not arc_iterator.Done():
            arc = arc_iterator.Value()
            olabel = arc.olabel
            if isyms[olabel] in to_expand:
                cs = isyms[olabel].split('|')
                arcs_to_expand.append((arc.ilabel, arc.weight, arc.nextstate, cs,))
            else:
                arcs_to_keep.append((arc.ilabel, olabel, arc.weight, arc.nextstate,))
            arc_iterator.Next()
        fst.delete_arcs(state)
        for arc in arcs_to_keep:
            fst.add_arc(state, arc[3], arc[0], arc[1], arc[2])
        for ilabel, weight, nextstate, cs in arcs_to_expand:
            newstate = fst.add_state()
            assert len(cs) == 2
            for i, c in enumerate(cs):
                if i == 0:
                    fst.add_arc(state, newstate, ilabel, isyms[c], weight)
                else:
                    fst.add_arc(newstate, nextstate, 0, isyms[c], 0.)
    return fst

def add_start_end(fst, isyms):
    finals = []
    for state in fst.states():
        if fst.final(state) != float('inf'):
            finals.append(state)
    l = isyms['<b>']
    new_final = fst.add_state()
    for state in finals:
        fst.add_arc(state, new_final, l, l, 0.)
        fst.set_final(state, float('inf'))
    fst.set_final(new_final)
    return fst
        

def main(infsts, to_expand_f, isym_f, outfsts, noexpand: ('', 'flag', None) = False):

    isyms = {}
    with open(isym_f) as fh:
        for line in fh:
            w, i = line.split()
            i = int(i)
            isyms[w] = i
            isyms[i] = w

    to_expand = set()
    if not noexpand:
        with open(to_expand_f) as fh:
            to_expand = set(fh.read().splitlines())

    lst = wrappedfst.WrappedFst.read_ark_entries(infsts)
    if os.path.isfile(outfsts): os.remove(outfsts)
    for key, fst in lst:
        if not noexpand:
            fst = expand_fst(fst, to_expand, isyms)        
        fst = add_start_end(fst, isyms)
        fst.arc_sort("olabel")
        fst.write_ark_entry(key, outfsts)

import plac; plac.call(main)
