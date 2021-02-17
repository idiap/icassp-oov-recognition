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

import plac
from collections import defaultdict


def main(vocab_f, idir: ('Should be a kaldi style data directory'), test_f, train_f):
    """ Will aim to create a test set with at least one OOV per utterance, at least 1000 utterances big """
    vocab = set()
    with open(vocab_f) as fh:
        for line in fh:
            vocab.add(line.strip())

    utts_oov = []
    utts_other = []
    with open(f'{idir}/text') as fh:
        for line in fh:
            utt, *words = line.split()
            add = False
            for word in words:
                if word not in vocab:
                    #print(word)
                    add = True
            if add:
                utts_oov.append(utt)
            else:
                utts_other.append(utt)
    print(f'Num utts oov {len(utts_oov)} out of {len(utts_other)}')
    cid_to_utts = defaultdict(list)
    utt_to_cid = {}
    with open(f'{idir}/utt2spk') as fh:
        for line in fh:
            uttid, cid = line.split()
            cid_to_utts[cid].append(uttid)
            utt_to_cid[uttid] = cid

    train_cids = defaultdict(int) 
    for utt in utts_other:
        cid = utt_to_cid[utt]
        train_cids[cid] += 1

    test_utts = []
    for utt in utts_oov:
        if utt_to_cid[utt] not in train_cids:
            test_utts.append(utt)

    if len(test_utts) < 1000:
        oov_cids = set()
        for utt in utts_oov:
            oov_cids.add(utt_to_cid[utt])
        oov_cids_cnt = {oov_cid: train_cids.get(oov_cid, 10000) for oov_cid in oov_cids}
        oov_cids_sorted = sorted(oov_cids_cnt.items(), key=lambda x: x[1])
        train_utts_to_remove = []
        to_add = 1000 - len(test_utts)
        if to_add > len(oov_cids_sorted) // 5:
            raise RuntimeError("Looks like there's not enough OOVs / train set too small :(")
        for i in range(to_add):
            cid = oov_cids_sorted[i][0]
            del train_cids[cid]
            utts = cid_to_utts[cid]
            train_utts_to_remove.extend(utts)
        utts_other = set(utts_other) - set(train_utts_to_remove)
        for utt in utts_oov:
            if utt_to_cid[utt] not in train_cids and utt not in test_utts:
                test_utts.append(utt)
    print(f'Num test utts final {len(test_utts)}')
    with open(test_f, 'w') as fh:
        for utt in test_utts:
            fh.write(f'{utt}\n')

    with open(train_f, 'w') as fh:
        for utt in utts_other:
            fh.write(f'{utt}\n')

plac.call(main)
