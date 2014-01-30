# @Created by L. Moruz 
# May 12th, 2011
# Edited by V. Granholm, 2014

'''
Script that provides different utilities for peptides
'''

class PeptideUtilities:
  def __init__(self):
    pass 

  def ComputeMonoisotopicMassWithoutCharge(self, peptide_sequence):
    '''
    Compute the monoisotopic mass of a peptide using the masses given below 
    '''

    monoisotopicMasses = [71.0371103, 0.0, 160.0306443, 115.0269385, 129.0425877, \
                         147.0684087, 57.0214611, 137.0589059, 113.0840579, 0.0, \
                         128.0949557, 113.0840579, 131.0404787, 114.0429222, 0.0, \
                         97.0527595, 128.0585714, 156.1011021, 87.0320244, 101.0476736, \
                         0.0, 99.0684087, 186.0793065, 0.0, 163.0633228, 0.0]
    monoisotopicmassWater = 18.0105633
 
    # check if the peptide is given in the format X.Y.Z
    if peptide_sequence.find(".") != -1:
      peptide = peptide_sequence.split(".")[1]
    else: 
      peptide = peptide_sequence
    # compute the mass 
    mm = 0.0
    for i in range(len(peptide)):
      mm += monoisotopicMasses[ord(peptide[i]) - ord('A')]

    return (mm + monoisotopicmassWater)


  def ComputeMonoisotopicMassWithCharge(self, peptide_sequence):
    '''
    Compute the monoisotopic mass of a peptide including the charge (one proton)
    '''

    monoisotopicmassProton = 1.00727646677
    mm = self.ComputeMonoisotopicMassWithoutCharge(peptide_sequence)
 
    return (mm + monoisotopicmassProton)

def estimate_pI(peptide):
    iso_alphabet = "DECYHKR"
    pKiso = [-3.86, -4.25, -8.33, -10.0, 6.0, 10.5, 12.4 ]  # Lehninger
    pKN = 9.69
    pKC = 2.34
    epsilon = 0.01;
    numAA = 7*[0]
    for i in range(len(iso_alphabet)):
        numAA[i] = peptide.count(iso_alphabet[i])
    pH = 6.5
    pHlow = 2.0
    pHhigh = 13.0
    while (((pH - pHlow) > epsilon) or ((pHhigh - pH) > epsilon)):
        NQ = 1 / (1 + pow(10, (pH - pKN))) - 1 / (1 + pow(10,(pKC - pH)))
        for i in range(len(numAA)):
            if numAA[i] == 0:
                pass
            elif pKiso[i] > 0:
                NQ += numAA[i] / (1 + pow(10, (pH - pKiso[i])))
            else:
                NQ -= numAA[i] / (1 + pow(10, (-pKiso[i] - pH)))
        if NQ < 0:
            pHhigh = pH
            pH -= ((pH - pHlow) / 2)
        else:
            pHlow = pH;
            pH += ((pHhigh - pH) / 2);
    return pH
