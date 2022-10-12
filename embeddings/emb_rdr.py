#!/usr/bin/env python3

import argparse
import sys
import time
from pathlib import Path

import h5py
import numpy as np

from Evaluation import split_dataset
from Evaluation.lib import FastaThreeLineIterator


class EmbReader:

    def __init__(self, _dir=Path(__file__).resolve().parents[1], as_script=False, **kwargs):
        if as_script:
            # parse the commandline
            parser = argparse.ArgumentParser()
            parser.add_argument('-e', '--embeddings',
                                action='store', dest='h5',
                                help='h5 datei', default=_dir / 'data' / 'baseline_embeddings_transmembrane.h5')
            parser.add_argument('-f', '--fasta',
                                action='store', dest='fasta',
                                help='fasta datei', default=_dir / 'data' / 'transmembrane.fasta')
            self.options = parser.parse_args()
        else:
            # make the keyword arguments that were passed if this is a function call look like an options object
            self.options = argparse.Namespace(**kwargs)

        print(f'now using: {self.options.fasta}')
        self.fasta_dict = {}
        with open(self.options.fasta, 'r') as stream:
            for r in FastaThreeLineIterator(stream):
                self.fasta_dict[r.id.split('|')[0]] = r

        # start parsing the embeddings
        f = h5py.File(self.options.h5, 'r')
        # print(list(f.keys()))
        # print(f['1afo_A'])

        if 'splits' in self.options:
            # fetch splits passed as argument if this right now is a function call
            self.splits = self.options.splits
        else:
            try:
                # or create the splits now
                print('fetching new splits')
                self.splits = [list(d.keys()) for d in split_dataset.find_optimal_split(seed=4000)]
            except:
                # let's just pretend
                self.splits = [[], [], [], [], [], []]
                [self.splits[i % 6].append(_id) for i, _id in enumerate(list(f.keys()))]

        num_residues_overall = 0
        proteins_clustered = [[], [], [], [], [], []]
        self.all_protein_ids_clustered = [[], [], [], [], [], []]
        self.all_labels_clustered = [[], [], [], [], [], []]

        print(f'proteins: {len(f)}')
        start = time.time()
        # iterate over the H5 file in the outer loop, assuming access there is slower than in self.splits
        for counter, key in enumerate(f.keys()):

            for i, split in enumerate(self.splits):
                if key not in split:
                    # skip to the next split
                    continue
                try:
                    labels = self.fasta_dict[key].labels
                    # only pick the rows where the label is not unknown
                    ar = f[key][[j for j, l in enumerate(labels) if l != 'U'], :]
                    # drop the U's now, not earlier!
                    labels = labels.replace('U', '')
                    # get a numpy column representing the labels
                    col = np.array(['H12SU'.index(l) for l in labels]).reshape(-1, 1)
                    # glue the label column to the left side of the array and save it
                    proteins_clustered[i].append(np.hstack((col, ar)))
                    # append to this loooong string of labels for this split
                    self.all_labels_clustered[i].append(labels)
                    self.all_protein_ids_clustered[i].append(key)
                    num_residues_overall += len(f[key])

                except Exception as ex:
                    print(ex)
                    print(f'The protein {key} is too short. Let\'s ignore it!', file=sys.stderr)
                # we're done for this protein, so don't bother checking the other splits:
                break

            # if not counter % 20:  # for every 20th protein, print the current progress
            #     print(f'{counter}/{len(f.keys())}')
        f.close()  # close the file!
        print(f'time taken: {(time.time() - start):.2f}sec')
        print(f'saw {num_residues_overall} residues overall')

        for i, cluster in enumerate(proteins_clustered):
            np.save(file=_dir / 'data' / f'ppvector_pp1_{i + 1}.npy',
                    arr=np.vstack(cluster), allow_pickle=False)  # takes 90% less hard drive space!
            # with open(_dir / 'data' / f'ppvector_pp1_{i + 1}.pkl', 'wb') as per_protein:
            #     pickle.dump(cluster, per_protein)

        with open(_dir / 'data' / 'split_labels.txt', 'w') as fh:
            for i, (split_protein_ids, split_labels) in enumerate(
                    zip(self.all_protein_ids_clustered, self.all_labels_clustered)):
                fh.write(f'>split {i + 1}\n')
                for _id, _labels in zip(split_protein_ids, split_labels):
                    fh.write(f'{_id}: {_labels}\n')
        print('embedding reader is done')


if __name__ == '__main__':
    EmbReader(as_script=True)
