#!/usr/bin/env python

import os
import sys

class Qvality():
    '''A class to handle qvality output (and in the future run it)'''
    
    def __init__(self, target_path=None, decoy_path='', output_path=None):
        self.target_path = target_path
        self.decoy_path = decoy_path
        self.output_path = output_path

    def run(self, reverse=False, include_negatives=True):
        '''Run qvality'''
        if not self.target_path:
            sys.exit('No target input file were set while initiating Qvality object')
        command = 'qvality '
        if reverse:
            command = command + '-r '  # Lower score is better
        if include_negatives:
            command = command + '-d '  # Decoy probabilities are added
        if self.output_path:
            command = command + '-o %s ' % (self.output_path)
        command = command + '%s %s' % (self.target_path, self.decoy_path)
        print command
        # Run command
        os.system(command)

    def read_output(self):
        '''Given a path to a qvality output, return a dict of scores and PEPs'''
        qvality_dict = {}
        for line in open(self.output_path):
            words = line.split()
            if words[0] == 'Score':
                # Skip first line
                continue
            score = float(words[0])
            pep = float(words[1])
            q = float(words[2])
            qvality_dict[score] = {'q':q, 'PEP':pep}
        return qvality_dict

    def count_confident(self, threshold=0.01, metric='qvalue'):
        '''Count the number of occurrences with q-value (or PEP) lower than the input threshold'''
        # Choose column
        if metric == 'qvalue':
            column = 2
        elif metric.lower() == 'pep':
            column = 1
        else:
            raise IOError('Did not recognize metric, should be qvalue or pep')
        # Count occurrences
        count = 0
        for line in open(self.output_path):
            if line.startswith('Score'):
                continue  # Skip header line
            words = line.split()
            if float(words[column]) < threshold:
                count += 1
            else:
                return count


    def remap_scores(self, peptides, qvality):
        ''' (Made by Jorrit)
        Remaps qvality scores using two dictionaries (where scores are float)
        peptides has form peptides[seq]['score'] = score
        qvality has form qvality[score]['q'] and qvality[score]['PEP']
        '''
        def check_score_and_update(score, seq):
            try:
                peptides[seq]['q'] = qvality[score]['q']
                peptides[seq]['PEP'] = qvality[score]['PEP']
                return True
            except KeyError:
                return False
        print 'Remapping scores'
        for seq in peptides:
            score = peptides[seq]['score']
            # TODO: Are all of these checks useful?
            if check_score_and_update(score, seq):
                continue
            elif check_score_and_update(round(score,6), seq):
                continue
            elif check_score_and_update(float('{:.5e}'.format(score)), seq):
                continue
            elif check_score_and_update(round(score+0.0000001,6), seq):
                continue
            elif check_score_and_update(round(score-0.0000001,6), seq):
                continue
            else:
                print 'Score: ', score
                raise Exception, 'error in remapping qvality scores'
        return peptides


def main():
    print 'A module to parse qvality output'

if __name__ == "__main__":
    main()
