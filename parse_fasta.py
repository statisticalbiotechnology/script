#!/usr/bin/env python

from collections import defaultdict

class Fasta():
    '''Imports a fasta file, and stores proteins and genes in dict'''

    def __init__(self, pathname, unique_proteins=False, ensembl=False):
        '''Reads and stores the proteins'''
        '''If unique_proteins is True, identical proteins are stored as one'''
        self.proteins = {}  # Protein to sequence (or protein to list if trypsinize has been run)
        self.genes = defaultdict(list)     # Genes to list of proteins
        self.protein_to_gene_dict = {}
        self.previous_sequences = {}  # To check whether a given protein is identical to previous protein
        # Read fasta
        sequence = ''
        for line in open(pathname):
            # Header line, first time (or when sequence is empty in fasta-file)
            if line[0] == '>' and sequence == '':
                # Obtain protein name
                protein = line.split()[0][1:]
                if ensembl:
                    gene = [word[5:] for word in line.split() if word.find('gene:') == 0][0]
                else:
                    gene = 'N/A'
                self.protein_to_gene_dict[protein] = gene
            # Header line, other times
            elif line[0] == '>' and sequence != '':
                self.store_previous_protein(protein, gene, sequence, unique_proteins)
                # Update values
                sequence = ''
                protein = line.split()[0][1:]
                if ensembl:
                    gene = [word[5:] for word in line.split() if word.find('gene:') == 0][0]
                else:
                    gene = 'N/A'
                self.protein_to_gene_dict[protein] = gene
            # Sequence line
            else:
                sequence = ''.join([sequence, line.strip()])  # Should be fast... 
        # Store the last protein of the fasta-file
        self.store_previous_protein(protein, gene, sequence, unique_proteins)


    def store_previous_protein(self, protein, gene, sequence, unique_proteins):
        '''Takes protein and gene names and aa-sequence, and boolean about how to treat repeated proteins'''
        # Store previous protein, (if unique, store only new proteins)
        if unique_proteins and sequence not in self.previous_sequences:
            self.proteins[protein] = sequence
            self.genes[gene].append(protein)
            self.previous_sequences[sequence] = 1
        # Identical protein sequence has already been stored, so pass
        elif unique_proteins:
            pass
        # Or... store all proteins
        else:
            self.proteins[protein] = sequence
            self.genes[gene].append(protein)  


    def trypsinize(self, min_len=6, max_len=40):
        '''
        Trypsinizes the imported protein sequences.
        Hence, it makes the string in self.proteins to a list of strings, with peptides fulfilling the length criteria
        '''
        self.peptide_to_protein_dict = defaultdict(list)
        self.peptide_to_gene_dict = defaultdict(list)
        for protein in self.proteins:
            gene = self.protein_to_gene_dict[protein]
            sequence = self.proteins[protein]
            peptide_list = []
            peptide = ''
            for i, aa in enumerate(sequence):
                # Get next amino acid
                try:
                    next_aa = sequence[i+1]
                except IndexError:  # Last aa
                    next_aa = ''
                # Update peptide
                peptide = ''.join([peptide, aa])
                # Cleavage?
                if aa in 'RK' and next_aa != 'P':
                    if min_len <= len(peptide) <= max_len:
                        peptide_list.append(peptide)
                        self.peptide_to_protein_dict[peptide].append(protein)
                        self.peptide_to_gene_dict[peptide].append(gene)
                    peptide = ''
            # Last peptide in protein
            if min_len <= len(peptide) <= max_len:
                peptide_list.append(peptide)
                self.peptide_to_protein_dict[peptide].append(protein)
                self.peptide_to_gene_dict[peptide].append(gene)
            # Replace the sequence string with a list of all valid tryptic peptides
            self.proteins[protein] = peptide_list
        # Uniqify the lists of unique peptide mappings
        for peptide in self.peptide_to_protein_dict:
            self.peptide_to_protein_dict[peptide] = list(set(self.peptide_to_protein_dict[peptide]))
            self.peptide_to_gene_dict[peptide] = list(set(self.peptide_to_gene_dict[peptide]))




def main():
    print 'A module with a class to parse fasta-files'




if __name__ == "__main__":
    main()
