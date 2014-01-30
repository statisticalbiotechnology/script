#! /usr/bin/env python

'''
Module with a lot of useful functions for parsing pin and pout xml files.
If percolator_xml.py is in the python path, example usage:

>> from percolator_xml import Pin
>> pin = Pin('path/to/pin.xml')
>> target_feature_value = pin.get_feature_values(2, isDecoy=False)  # Gets third feature

/Viktor
'''

import sys
from lxml import etree

class Pout():
    '''Class that reads pout file, takes pathname as input'''

    def __init__(self, pathname, ns_num=14, is_huge_tree=False):
        '''Initialize with pathname to pout-file'''        
        self.pathname = pathname
        self.ns = 'http://per-colator.com/percolator_out/%s' % (ns_num)
        parser = etree.XMLParser(ns_clean=False, huge_tree=is_huge_tree)
        self.tree = etree.parse(self.pathname, parser)

    def get_qvalues(self, level, isDecoy=False):
        '''Output a list of a qvalues of a pout file, for peptides or psms'''
        qvalues = []
        for element in self.tree.findall('{%(n)s}%(l)ss/{%(n)s}%(l)s' % dict(n=self.ns, l=level)):
            try:
                if (element.attrib['{%s}decoy' % (self.ns)] == 'true') == isDecoy:
                    qvalue = float(element.find('{%s}q_value' % (self.ns)).text)
                    qvalues.append(qvalue)
            except KeyError:
                # Seems Percolator was run without the -Z flag, so there's no p:decoy attribute
                if isDecoy == True:
                    sys.exit('KeyError: There might be no p:decoy attribute')
                qvalue = float(element.find('{%s}q_value' % (self.ns)).text)
                qvalues.append(qvalue)
        return qvalues

    def get_psm_ids(self, isDecoy=False):
        '''Output a list of a qvalues of a pout file, for peptides or psms'''
        psm_ids = []
        for element in self.tree.findall('{%s}psms/{%s}psm' % (self.ns, self.ns)):
            try:
                if (element.attrib['{%s}decoy' % (self.ns)] == 'true') == isDecoy:
                    psm_id = element.attrib['{%s}psm_id' % (self.ns)]
                    psm_ids.append(psm_id)
            except KeyError:
                # Seems Percolator was run without the -Z flag, so there's no p:decoy attribute
                if isDecoy == True:
                    sys.exit('KeyError: There might be no p:decoy attribute')
                psm_id = element.attrib['{%s}psm_id' % (self.ns)]
                psm_ids.append(psm_id)
        return psm_ids

    def get_scores(self, level, isDecoy=False):
        '''Output a list of a scores of a pout file, for peptides or psms'''
        scores = []
        for element in self.tree.findall('{%s}%ss/{%s}%s' % (self.ns, level, self.ns, level)):
            try:
                if (element.attrib['{%s}decoy' % (self.ns)] == 'true') == isDecoy:
                    score = float(element.find('{%s}svm_score' % (self.ns)).text)
                    scores.append(score)
            except KeyError:
                # Seems Percolator was run without the -Z flag, so there's no p:decoy attribute
                score = float(element.find('{%s}svm_score' % (self.ns)).text)
                scores.append(score)
        return scores

    def get_peptides(self, threshold=1, value_name='svm_score', ptms=True, isDecoy=False):
        '''Output dictionary of peptides below the q-value threshold, value_name decides the type of value in dict'''
        peptides = {}
        for element in self.tree.findall('{%s}peptides/{%s}peptide' % (self.ns, self.ns)):
            try:
                if (element.attrib['{%s}decoy' % (self.ns)] == 'true') == isDecoy:
                    peptide = element.attrib['{%s}peptide_id' % (self.ns)]
                    # Remove UNIMOD bits, if not considering PTMs
                    if not ptms:
                        peptide = self.strip_mods(peptide)
                    qvalue = float(element.find('{%s}q_value' % (self.ns)).text)
                    value = float(element.find('{%s}%s' % (self.ns, value_name)).text)
                    # Add peptide, unless it's already added
                    if qvalue <= threshold and peptide not in peptides:
                        peptides[peptide] = value
            except KeyError:
                # Seems Percolator was run without the -Z flag, so there's no p:decoy attribute
                if isDecoy == True:
                    sys.exit('KeyError: There might be no p:decoy attribute')
                peptide = element.attrib['{%s}peptide_id' % (self.ns)]
                # Remove UNIMOD bits, if not considering PTMs
                if not ptms:
                    peptide = self.strip_mods(peptide)
                qvalue = float(element.find('{%s}q_value' % (self.ns)).text)
                value = float(element.find('{%s}%s' % (self.ns, value_name)).text)
                # Add peptide, unless it's already added
                if qvalue <= threshold and peptide not in peptides:
                    peptides[peptide] = value
        return peptides

    def get_proteins(self, threshold=1, value_name='q_value', isDecoy=False):
        '''Output dictionary of proteins below the q-value threshold, value_name decides the type of value in dict'''
        proteins = {}
        for element in self.tree.findall('{%s}proteins/{%s}protein' % (self.ns, self.ns)):
            try:
                if (element.attrib['{%s}decoy' % (self.ns)] == 'true') == isDecoy:
                    protein = element.attrib['{%s}protein_id' % (self.ns)]
                    qvalue = float(element.find('{%s}q_value' % (self.ns)).text)
                    value = float(element.find('{%s}%s' % (self.ns, value_name)).text)
                    if qvalue <= threshold:
                        proteins[protein] = value
            except KeyError:
                # Seems Percolator was run without the -Z flag, so there's no p:decoy attribute
                if isDecoy == True:
                    sys.exit('KeyError: There might be no p:decoy attribute')
                protein = element.attrib['{%s}protein_id' % (self.ns)]
                qvalue = float(element.find('{%s}q_value' % (self.ns)).text)
                value = float(element.find('{%s}%s' % (self.ns, value_name)).text)
                if qvalue <= threshold:
                    proteins[protein] = value
        return proteins

    def get_peptide_to_protein_dict(self, threshold=1, ptms=True, isDecoy=False):
        '''Output dictionary of peptides to protein dictionary, for those below the q-value threshold'''
        peptides = {}
        for element in self.tree.findall('{%s}peptides/{%s}peptide' % (self.ns, self.ns)):
            proteins = []
            if (element.attrib['{%s}decoy' % (self.ns)] == 'true') == isDecoy:
                peptide = element.attrib['{%s}peptide_id' % (self.ns)]
                # Remove UNIMOD bits, if not considering PTMs
                if not ptms:
                    peptide = self.strip_mods(peptide)
                qvalue = float(element.find('{%s}q_value' % (self.ns)).text)
                # Get protein_ids
                for protein_id in element.findall('{%s}protein_id' % (self.ns)):
                    proteins.append(protein_id.text)
                # Store peptide, unless it's already stored
                if qvalue <= threshold and peptide not in peptides:
                    peptides[peptide] = proteins
        return peptides

    def strip_mods(self, seq):
        '''Take a peptide sequence, strip all occurences of [UNIMOD:XX]'''
        seq = seq.replace('[UNIMOD:', '')
        seq = seq.replace(']', '')
        seq = ''.join([letter for letter in seq if letter.isalpha()])
        return seq



class Pin():
    '''Class that reads pin file, takes a pathname as input'''

    def __init__(self, pathname, ns_num=13, is_huge_tree=False):
        '''Initialize with pathname to pin-file'''
        self.pathname = pathname
        self.ns = 'http://per-colator.com/percolator_in/%s' % (ns_num)
        parser = etree.XMLParser(ns_clean=False, huge_tree=is_huge_tree, remove_blank_text=True)  # Remove blank enables pretty print
        self.tree = etree.parse(self.pathname, parser)

    def get_feature_names(self):
        '''Output a list of all the feature names in the pin file'''
        feature_names = []
        for element in self.tree.findall('{%s}featureDescriptions/{%s}featureDescription' % (self.ns, self.ns)):
            name = element.attrib['name']
            feature_names.append(name)
        return feature_names

    def get_feature_values(self, feature_index, isDecoy=False, type_func=float):
        '''Extract all features values for a given feature, and output a list (first feature has index 0)'''
        values = []
        for element in self.tree.findall('{%s}fragSpectrumScan/{%s}peptideSpectrumMatch' % (self.ns, self.ns)):
            element_isDecoy = element.attrib['isDecoy'] == 'true'  # Convert element string boolean to boolean
            if isDecoy == element_isDecoy:
                feature_element = element.findall('{%s}features/{%s}feature' % (self.ns, self.ns))[feature_index]
                feature_value = feature_element.text
                values.append(type_func(feature_value))  # type_func: float(), int(), etc...
        return values

    def get_sequences(self, isDecoy=False):
        '''Extract peptide sequences for each psm, and output these in a list'''
        sequences = []
        for element in self.tree.findall('{%s}fragSpectrumScan/{%s}peptideSpectrumMatch' % (self.ns, self.ns)):
            element_isDecoy = element.attrib['isDecoy'] == 'true'  # Convert element string boolean to boolean
            if isDecoy == element_isDecoy:
                sequence = element.find('{%s}peptide/{%s}peptideSequence' % (self.ns, self.ns)).text
                sequences.append(sequence)
        return sequences

    def get_psm_attributes(self, attrib_name, isDecoy=False, type_func=str):
        '''Extract values of the attrib_name attributes of the peptideSpectrumMatch elements'''
        attributes = []
        for element in self.tree.findall('{%s}fragSpectrumScan/{%s}peptideSpectrumMatch' % (self.ns, self.ns)):
            element_isDecoy = element.attrib['isDecoy'] == 'true'  # Convert element string boolean to boolean
            if isDecoy == element_isDecoy:
                attribute = type_func(element.attrib[attrib_name])  # type_func is the function (float, int, etc.)
                attributes.append(attribute)
        return attributes

    def get_peptide_to_protein_dictionary(self, isDecoy=False):
        '''Store the peptides in a directory as keys, with their proteins in lists as values'''
        peptide_to_protein_dict = {}
        for element in self.tree.findall('{%s}fragSpectrumScan/{%s}peptideSpectrumMatch' % (self.ns, self.ns)):
            element_isDecoy = element.attrib['isDecoy'] == 'true'  # Convert element string boolean to boolean
            if isDecoy == element_isDecoy:
                peptide = element.find('{%s}peptide/{%s}peptideSequence' % (self.ns, self.ns)).text
                # Get proteins
                protein_list = []
                for occ in element.findall('{%s}occurence' % (self.ns)):
                    protein_list.append(occ.attrib['proteinId'])
                # Store proteins in dictionary
                try:
                    peptide_to_protein_dict[peptide].extend(protein_list)
                except KeyError:
                    peptide_to_protein_dict[peptide] = protein_list
                # Uniqify
                peptide_to_protein_dict[peptide] = list(set(peptide_to_protein_dict[peptide]))
        return peptide_to_protein_dict

    def remove_features(self, feature_names, add_value_to_ind=None, type_func=float):
        '''Remove all instances of certain feature(s) from the pin-file'''
        # Make a list of the feature if it is only one
        if type(feature_names) == str:
            feature_names = [str]
        # Find indeces
        indeces = []
        descriptions = self.tree.find('{%s}featureDescriptions' % (self.ns))
        for description in descriptions.findall('{%s}featureDescription' % (self.ns)):
            if description.attrib['name'] in feature_names:
                indeces.append(descriptions.index(description))
        indeces.sort(reverse=True)
        # Remove feature descriptions
        for i in indeces:
            descriptions.remove(descriptions[i])
        # Remove features from PSM elements
        for scan in self.tree.findall('{%s}fragSpectrumScan' % (self.ns)):
            for psm in scan.findall('{%s}peptideSpectrumMatch' % (self.ns)):
                value = type_func(0)
                features = psm.find('{%s}features' % (self.ns))
                for i in indeces:
                    value += type_func(features[i].text)
                    features.remove(features[i])
                if add_value_to_ind is not None:
                    features[add_value_to_ind].text = str(value)

    def add_feature(self, feature_name, target_dict, decoy_dict):
        '''Add a feature to a pin-file, requires dictionaries between psm-id and feature values'''
        # Add feature name to featureDescriptions
        feat_descrip = self.tree.find('{%s}featureDescriptions' % (self.ns))
        new_element = feat_descrip.makeelement('{%s}featureDescription' % (self.ns))
        new_element.set('name', feature_name)
        feat_descrip.append(new_element)
        # Add feature values to PSMs
        for isDecoy, id_to_value in [(False, target_dict), (True, decoy_dict)]:
            for psm in self.tree.findall('{%s}fragSpectrumScan/{%s}peptideSpectrumMatch' % (self.ns, self.ns)):
                element_isDecoy = psm.attrib['isDecoy'] == 'true'  # Convert element string boolean to boolean
                if isDecoy == element_isDecoy:
                    psm_id = psm.attrib['id']
                    try:
                        feature_value = id_to_value[psm_id]
                        features = psm.find('{%s}features' % (self.ns))
                        new_feature = features.makeelement('{%s}feature' % (self.ns))
                        new_feature.text = str(feature_value)
                        features.append(new_feature)
                    except KeyError:
                        sys.exit('KeyError: %s was not present in id dictionary' % (psm_id))

    def get_ids(self, isDecoy=False):
        '''Extract psm id's output these to list'''
        ids = []
        for element in self.tree.findall('{%s}fragSpectrumScan/{%s}peptideSpectrumMatch' % (self.ns, self.ns)):
            element_isDecoy = element.attrib['isDecoy'] == 'true'  # Convert element string boolean to boolean
            if isDecoy == element_isDecoy:
                psm_id = element.attrib['id']
                ids.append(psm_id)
        return ids        

    def write(self, outpath):
        '''Write the current tree to outpath'''
        self.tree.write(outpath, pretty_print=True)


def main():
    print "This is a module with a few functions for parsing percolator xml files"

if __name__ == "__main__":
    main()


