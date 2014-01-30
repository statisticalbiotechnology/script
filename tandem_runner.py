#!/usr/bin/env python

import os
import glob

class Tandem():
    
    def __init__(self, parameters, parameter_file_path):
        '''Initializes class
        Arguments:
        parameters - a dictionary with key: xtandem bioml labels, value: the value of the labels
        '''
        self.TANDEM_PATH = '/afs/pdc.kth.se/home/v/vgr/vol01/other/xtandem/tandem-linux-13-09-01-1/bin/32-bit/all_static/tandem.exe'
        default_path = '/afs/pdc.kth.se/home/v/vgr/vol01/other/xtandem/tandem-linux-13-09-01-1/bin/default_input.xml'
        taxonomy_path = '/afs/pdc.kth.se/home/v/vgr/vol01/other/xtandem/tandem-linux-13-09-01-1/bin/taxonomy.xml'
        default_element = '<note type="input" label="list path, default parameters">%s</note>' % (default_path)
        taxonomy_element = '<note type="input" label="list path, taxonomy information">%s</note>' % (taxonomy_path)
        self.BIOML_HEADER = '<?xml version="1.0"?>\n<?xml-stylesheet type="text/xsl" href="tandem-input-style.xsl"?>\n\
<bioml>\n%s\n%s\n\n' % (default_element, taxonomy_element)
        self.BIOML_FOOTER = '</bioml>'
        self.parameters = parameters
        self.parameter_file_path = parameter_file_path

    def run(self):
        '''Write a tandem input bioml file, and execite tandem'''
        # Make input file
        input_file = open(self.parameter_file_path, 'w')
        input_file.write(self.BIOML_HEADER)
        for label, value in self.parameters.items():
            input_file.write('<note type="input" label="%s">%s</note>\n' % (label, value))
        input_file.write(self.BIOML_FOOTER)
        input_file.close()    
        # Run X!Tandem
        os.system('%s %s' % (self.TANDEM_PATH, self.parameter_file_path))


class TandemJob():

    def __init__(self):
        self.spectra_path = ''
        self.parameters = {}
        self.database = ''
        self.output_dir = ''
        self.job_dir = ''
        self.input_dir = ''
        self.time = 120
        self.log_file_dir = ''
        # Private
        self._default_path = '/afs/pdc.kth.se/home/v/vgr/vol01/other/xtandem/tandem-linux-13-09-01-1/bin/default_input.xml'
        self._taxonomy_path = '/afs/pdc.kth.se/home/v/vgr/vol01/other/xtandem/tandem-linux-13-09-01-1/bin/taxonomy.xml'
        default_element = '<note type="input" label="list path, default parameters">%s</note>' % (self._default_path)
        taxonomy_element = '<note type="input" label="list path, taxonomy information">%s</note>' % (self._taxonomy_path)
        self._BIOML_HEADER = '<?xml version="1.0"?>\n<?xml-stylesheet type="text/xsl" href="tandem-input-style.xsl"?>\n\
<bioml>\n%s\n%s\n\n' % (default_element, taxonomy_element)
        self._BIOML_FOOTER = '</bioml>'
        self._TANDEM_PATH = '/afs/pdc.kth.se/home/v/vgr/vol01/other/xtandem/tandem-linux-13-09-01-1/bin/32-bit/all_static/tandem.exe'

    def run(self):
        # Make input files
        for spectra_filepath in glob.glob(self.spectra_path):
            input_xml_path = self._make_input_xml(spectra_filepath)
            job_path = self._make_jobscript(input_xml_path)
            self._submit_job(job_path)

    def _make_input_xml(self, spectra_filepath):
        # Make path
        basename = self._get_basename(spectra_filepath)
        input_path = '%s/%s.xml' % (self.input_dir, basename)
        # Write beginning of file
        input_file = open(input_path, 'w')
        input_file.write(self._BIOML_HEADER)
        input_file.write(self._make_bioml_note('spectrum, path', spectra_filepath))
        input_file.write(self._make_bioml_note('output, path', '%s/%s.xml' % (self.output_dir, basename)))
        input_file.write(self._make_bioml_note('protein, taxon', self.database))
        # Write parameters
        input_file.write('\n')
        for label, value in self.parameters.items():
            input_file.write(self._make_bioml_note(label, value))
        # Write end
        input_file.write(self._BIOML_FOOTER)
        input_file.close()
        return input_path     

    def _make_bioml_note(self, label, value):
        return '<note type="input" label="%s">%s</note>\n' % (label, value)

    def _get_basename(self, spectra_filepath):
        filename = spectra_filepath.split('/')[-1]
        basename_list = filename.split('.')[:-1]
        return '.'.join(basename_list)

    def _make_jobscript(self, input_xml_path):
        # Make jobscript path
        basename = self._get_basename(input_xml_path)
        job_path = '%s/%s.sh' % (self.job_dir, basename)
        # Write jobscript
        jobfile = open(job_path, 'w')
        jobfile.write('#!/bin/bash\n\n')
        if not self.log_file_dir == '':
            log_path = '%s/%s.log' % (self.log_file_dir, basename)
            jobfile.write('logfile=%s\nexec > $logfile 2>&1\n\n' % (log_path))
        jobfile.write('%s %s\n' % (self._TANDEM_PATH, input_xml_path))
        jobfile.close()
        os.system('chmod +x %s' % (job_path))
        return job_path

    def _submit_job(self, job_path):
        print 'Submit job %s' % (job_path)
        os.system('esubmit -n 1 -t %s %s' % (self.time, job_path))


def main():
    '''Doc:'''
    print 'This is a module for running X!Tandem'

if __name__ == "__main__":
    main()
