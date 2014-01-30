#!/usr/bin/env

import sys
import os
import glob
import numpy as np
import random


class Local():
    '''Class for the local computer'''
    def __init__(self):
        self.name = 'local'
        self.core_number = None  # A search engine job will be submitted to two cores... (target and decoy)
        self.just_test = False

    def submit(self, job_path):
        '''Just submits the job'''
        os.system('%s' % (job_path))


class Ferlin():
    '''Class for the ferlin cluster'''
    def __init__(self):
        # Name
        self.name = 'ferlin'
        # Parameters
        self.time = 300
        self.core_number = None  # If None, no parallelization (only target & decoy)
        self.just_test = False
        os.system('module add easy')

    def submit(self, job_path):
        '''Submits the job, time is needed'''
        os.system('esubmit -n 1 -t %s %s' % (self.time, job_path))


class Kalkyl():
    '''Class for Kalkyl on Uppmax'''
    def __init__(self):
        # Name
        self.name = 'kalkyl'
        # Parameters
        self.time = 300
        self.project = 's00111-298'
        self.partition = 'node'
        self.core_number = None    # If None, no parallelization (only target & decoy)
        self.share_index_script = '/bubo/home/h14/viktorg/files/src/uppmax/share_crux-index_on_node.py'
        self.priority = False
        self.fat = False
        self.just_test = False

    def submit(self, job_path):
        os.system('sbatch %s' % (job_path))


class Search_engine():
    '''Parent class of all search_engine procedures'''
    def __init__(self):
        '''Initiate variables'''
        self.job_dir = None
        self.spectra_path = None
        self.target_database = None
        self.decoy_database = None
        self.output_dir = None
        self.job_label = 'search'  # This partly decides the label of the job-file

    def check_required_parameters(self):
        '''Check whether some crucial parameters were not set'''
        # Some parameters required by all search engines
        required_for_all = [self.job_dir, self.spectra_path, self.target_database, self.decoy_database, self.output_dir]
        names_for_all = ['job_dir', 'spectra_path', 'target_database', 'decoy_database', 'output_dir']
        # This function adds parameters required by a specific search engine
        required, names = self.specific_required_params(required_for_all, names_for_all)
        # Check whether a required parameter is None, if so, error
        for i, value in enumerate(required):
            if value == None:
                sys.exit('Error: Parameter %s not set, required by search engine' % names[i])

    def prepare_jobs(self):
        '''Counts how many cores, and prepare that many jobs'''
        jobs = []
        random_tag = ''.join(random.sample('123456789ABCDEFGHJKLMPQRSTUVWXYZabcdefghijkmoprqrstuvwxyz', 4))
        self.databases = dict(target=self.target_database, decoy=self.decoy_database)
        self.check_required_parameters()
        # For each node/core, write the commands to write
        filenames_per_partition = self.split_list_of_files()
        for i, file_set in enumerate(filenames_per_partition):
            # Separate jobs for targets and decoys
            for db in ('target', 'decoy'):
                job = Job(self.computer, self.name)
                job.path = '%s/%s_%s_%s_%.4d.sh' % (self.job_dir, self.job_label, db, random_tag, i+1)
                job.commands.append('outdir=$(readlink -f %s)' % (self.output_dir))
                job.commands.append('mkdir $outdir/%s' % (db))
                job.commands.append('workdir=$(mktemp -dp /scratch/)')
                job.commands.append('cd $workdir')
                job.commands.extend(self.installation_commands)
                job = self.search_engine_commands(job, file_set, db)  # Different code for MSGF+ and Crux, etc.
                job.finish_job()
                jobs.append(job)
        return jobs

    def split_list_of_files(self):
        '''The input spectra are separated into smaller lists of spectra to run on each node/core'''
        all_files = glob.glob(self.spectra_path)
        filenames_per_partition = []
        if all_files == []:
            sys.exit('ERROR: No spectra files found for path %s' % self.spectra_path)
        if self.computer.core_number == None:
            # Minimum parallelization, target & decoy are separated though
            filenames_per_partition.append(all_files)
            return filenames_per_partition
        if self.computer.name == 'kalkyl':  # FIXME: Add possibility to reserve whole node
            # Extra option if on Kalkyl, node or core
            if self.computer.partition == 'node':
                partition_number = int(self.computer.core_number) / 8
            else:
                partition_number = self.computer.core_number
        else:
            # If not kalkyl, use core
            partition_number = self.computer.core_number
        # Split files into sublists in a list
        num_files_per_partition = int(np.ceil(len(all_files) / float(partition_number)))
        filenames_per_partition = []
        for i in xrange(0, len(all_files), num_files_per_partition):
            subset_files = all_files[i:i+num_files_per_partition]
            if subset_files != []:
                filenames_per_partition.append(subset_files)
        return filenames_per_partition


class MSGF(Search_engine):
    '''Contains methods etc. for running MSGF'''
    def __init__(self, computer):
        self.computer = computer
        # Inherit from Search_engine
        Search_engine.__init__(self)
        # Name
        self.name = 'msgf'
        # Set parameters
        self.memory = 3500  # Mb
        self.inst = 0       # Instrument, 0: Low-res LTQ
        self.t = '10ppm'    # Tolerance
        self.protocol = '0' # 0: No modifications
        self.ti = '0,2'     # Isotope range
        self.m = '1'        # Activation method, 1: CID, 2: ETD
        self.e = '1'        # Enzyme, 0: No enzyme, 1: Trypsin
        self.ntt = '2'      # Number of tolerable termini
        self.n = '1'        # Number of matches per spectrum
        self.mod = None     # Modification filename
        # For local machine Feynman
        if self.computer.name == 'local':
            self.version = '9540'
            self.installation_commands = []
            self.path = '/afs/pdc.kth.se/home/v/vgr/vol01/other/MSGFDB/v%s/MSGFPlus.jar' % (self.version)
        # For Ferlin
        elif self.computer.name == 'ferlin':
            self.version = '9540'
            self.installation_commands = [
                'module add jdk']
            self.path = '/afs/pdc.kth.se/home/v/vgr/vol01/other/MSGFDB/v%s/MSGFPlus.jar' % (self.version)
        # For kalkyl
        elif self.computer.name == 'kalkyl':
            self.version = '9012'
            self.installation_commands = ['module add java']
            self.path = '/bubo/home/h14/viktorg/files/msgfplus/v%s/MSGFPlus.jar' % (self.version)
        print 'MSGF+ v%s, for %s computer' % (self.version, self.computer.name)

    def search_engine_commands(self, job, file_set, db):
        '''Run the commands of the specific search engine'''
        job.commands.append('')
        job.commands.append('cp %s %s.fasta' % (self.databases[db], db))
        for f in file_set:
            spectra_filename = f.split('/')[-1]
            mzid_filename = ''.join(spectra_filename.split('.')[:-1]) + '.mzid'
            job.commands.append('cp %s %s' % (f, spectra_filename))
            cmd = 'java -Xmx%sM -jar %s -s %s -d %s.fasta' % (self.memory, self.path, spectra_filename, db)
            cmd = cmd + ' -inst %s -t %s -protocol %s -ti %s' % (self.inst, self.t, self.protocol, self.ti)
            cmd = cmd + ' -m %s -e %s -n %s -addFeatures 1 -o %s' % (self.m, self.e, self.n, mzid_filename)
            if self.mod is not None:
                cmd = cmd + ' -mod %s' % (self.mod)
            if self.ntt != 2:
                cmd = cmd + ' -ntt %s' % (self.ntt)
            job.commands.append(cmd)
            job.commands.append('cp %s $outdir/%s/%s' % (mzid_filename, db, mzid_filename))
        return job  

    def specific_required_params(self, required, names):
        '''Adds the parameters specifically for MSGF+'''  # No parameters, in this case...
        required.extend([])
        names.extend([])
        return required, names


class Crux(Search_engine):
    '''Contains methods etc. for running Crux'''
    def __init__(self, computer):
        # Inherit from Search_engine
        Search_engine.__init__(self)
        # Name
        self.name = 'crux'
        # Set parameters
        self.computer = computer
        self.parameter_file = None  # Absolute path
        self.share_index = False
        # For Feynman (local)
        if self.computer.name == 'local':
#            crux_version = '1.39'
#            self.installation_commands = []
#            self.path = '/afs/pdc.kth.se/home/v/vgr/vol01/other/Crux/crux_%(v)s/crux-%(v)s.1-Linux/bin/crux' % dict(v=crux_version)
            crux_version = '1.40'
            self.installation_commands = []
            self.path = '/afs/pdc.kth.se/home/v/vgr/vol01/other/Crux/crux_%(v)s/crux-%(v)s.Linux.x86_64/bin/crux' % dict(v=crux_version)
        # For ferlin cluster
        elif self.computer.name == 'ferlin':
            version_number = 40
            crux_version = '1.%s' % (version_number)  # Does not handle mgf-files
            tar_path = '/afs/pdc.kth.se/home/v/vgr/vol01/other/Crux/crux_%(v)s/crux_%(v)s.tar.gz' % dict(v=crux_version)
            if version_number < 39:
                self.installation_commands = [
                    'cp %s .' % (tar_path),
                    'tar -zxvf crux_%s.tar.gz' % crux_version,
                    'cd crux_%s' % crux_version,
                    'package_dir=$(pwd)',
                    './configure --prefix $package_dir/install-dir',
                    'make',
                    'make install',
                    'cd ..']
            elif version_number >= 40:
                # These installations commands are valid for crux that uses cmake
                self.installation_commands = [
                    'module add cmake/2.8.8',
                    'module add subversion',
                    'cp %s .' % (tar_path),
                    'tar -zxvf crux_%s.tar.gz' % crux_version,
                    'cd crux-%s.Source' % crux_version,
                    'package_dir=$(pwd)',
                    'mkdir install_dir',
                    'cmake -DCMAKE_INSTALL_PREFIX:PATH=$package_dir/install_dir .',
                    'make',
                    'make install',
                    'cd ..']
            self.path = '$package_dir/install_dir/bin/crux'
        # For kalkyl cluster
        elif self.computer.name == 'kalkyl':
            crux_version = '1.39'
            self.installation_commands = []
            self.path = '/bubo/home/h14/viktorg/files/crux/crux_%s_patch/bin/crux' % (crux_version)  # TODO: PATCH!
        print 'Crux version %s for %s computer' % (crux_version, self.computer.name)

    def search_engine_commands(self, job, file_set, db):
        '''Run the commands of the specific search engine'''
        if self.parameter_file == None:
            sys.exit('ERROR: parameter_file not set')
        job.commands.append('params=%s' % (self.parameter_file))
        job.commands.append('')
        # Smart crux search, a single database per node, but multiple jobs
        if self.computer.name == 'kalkyl' and self.share_index:
            job.commands.append('cp %s db.fasta' % self.databases[db])
            job.commands.append('%s create-index --parameter-file $params db.fasta db.index' % (self.path))
            spectra_files = []
            for f in file_set:
                spectra_filename = f.split('/')[-1]
                job.commands.append('cp %s %s' % (f, spectra_filename))
                spectra_files.append(spectra_filename)
                spectra_file_string = ' '.join(spectra_files)
            job.commands.append('du -smh $workdir/*')
            job.commands.append('ls')
            job.commands.append('%s %s $params db.index %s' % (self.computer.share_index_script, self.path, spectra_file_string))
            job.commands.append('cp *.sqt $outdir/%s/.' % (db))
        # Normal crux search
        else:
            # index database
            job.commands.append('cp %s db.fasta' % self.databases[db])
            job.commands.append('%s create-index --parameter-file $params db.fasta db.index' % (self.path))
            # Make searches
            for f in file_set:
                spectra_filename = f.split('/')[-1]
                spectra_base = spectra_filename.split('.')[0]
                job.commands.append('cp %s %s' % (f, spectra_filename))
                job.commands.append('%s sequest-search --overwrite T --parameter-file $params %s db.index' % (self.path, spectra_filename))
                job.commands.append('cp crux-output/sequest.target.sqt $outdir/%s/%s.sqt' % (db, spectra_base))
        return job

    def specific_required_params(self, required, names):
        '''Adds the parameters specifically for Crux'''
        required.extend([self.parameter_file])
        names.extend(['parameter_file'])
        return required, names
  

class Tandem(Search_engine):
    '''Contains methods etc. for running X!Tandem'''
    def __init__(self, computer):
        self.computer = computer
        # Inherit from Search_engine
        Search_engine.__init__(self)
        # Name
        self.name = 'tandem'
        # Set parameters
        self.parameters = {}
        self.parameters
        self.memory = 3500  # Mb
        self.inst = 0       # Instrument, 0: Low-res LTQ
        self.t = '10ppm'    # Tolerance
        self.protocol = '0' # 0: No modifications
        self.ti = '0,2'     # Isotope range
        self.m = '1'        # Activation method, 1: CID, 2: ETD
        self.e = '1'        # Enzyme, 0: No enzyme, 1: Trypsin
        self.ntt = '2'      # Number of tolerable termini
        self.n = '1'        # Number of matches per spectrum
        self.mod = None     # Modification filename
        # For local machine Feynman
        if self.computer.name == 'local':
            self.installation_commands = []
            self.path = '/afs/pdc.kth.se/home/v/vgr/vol01/other/xtandem/tandem-linux-13-09-01-1/bin/32-bit/all_static/tandem.exe'
        # For Ferlin
        elif self.computer.name == 'ferlin':
            self.installation_commands = []
            self.path = '/afs/pdc.kth.se/home/v/vgr/vol01/other/xtandem/tandem-linux-13-09-01-1/bin/32-bit/all_static/tandem.exe'
        # For kalkyl
        elif self.computer.name == 'kalkyl':
            sys.exit('Tandem is not installed for Kalkyl yet')
        print 'Tandem, for %s computer' % (self.version, self.computer.name)

    def search_engine_commands(self, job, file_set, db):
        '''Run the commands of the specific search engine'''
        job.commands.append('')
        job.commands.append('cp %s %s.fasta' % (self.databases[db], db))
        for f in file_set:
            spectra_filename = f.split('/')[-1]
            mzid_filename = ''.join(spectra_filename.split('.')[:-1]) + '.mzid'
            job.commands.append('cp %s %s' % (f, spectra_filename))
            cmd = 'java -Xmx%sM -jar %s -s %s -d %s.fasta' % (self.memory, self.path, spectra_filename, db)
            cmd = cmd + ' -inst %s -t %s -protocol %s -ti %s' % (self.inst, self.t, self.protocol, self.ti)
            cmd = cmd + ' -m %s -e %s -n %s -addFeatures 1 -o %s' % (self.m, self.e, self.n, mzid_filename)
            if self.mod is not None:
                cmd = cmd + ' -mod %s' % (self.mod)
            if self.ntt != 2:
                cmd = cmd + ' -ntt %s' % (self.ntt)
            job.commands.append(cmd)
            job.commands.append('cp %s $outdir/%s/%s' % (mzid_filename, db, mzid_filename))
        return job  

    def specific_required_params(self, required, names):
        '''Adds the parameters specifically for MSGF+'''  # No parameters, in this case...
        required.extend([])
        names.extend([])
        return required, names



class Job():
    '''Class that describes a job script to be written'''

    def __init__(self, computer, name):
        '''Initialize a list of commands, and path to the script'''
        self.path = 'default_job_path.sh'
        self.commands = []
        self.commands.extend(['#!/bin/bash', ''])
        # Computer specific addition to beginning of job-file
        if computer.name == 'local':
            pass
        elif computer.name == 'ferlin':
            self.commands.append('hostname')
            self.commands.append('module add python')
        elif computer.name == 'kalkyl':
            if computer.partition == 'core':
                cores = 1 # computer.core_number
            elif computer.partition == 'node':
                cores = 8           
            self.commands.append('#SBATCH -A %s' % (computer.project))
            self.commands.append('#SBATCH -p %s' % (computer.partition))
            self.commands.append('#SBATCH -n %s' % (cores))  # If -p node, cores is always 8
            self.commands.append('#SBATCH -t %s' % (computer.time))
            self.commands.append('#SBATCH -J %s' % (name))
            if computer.priority:
                self.commands.append('#SBATCH --qos=short')
            if computer.fat:
                self.commands.append('#SBATCH -C fat')
            self.commands.append('')
            self.commands.append('hostname')
        # First job-line
        self.commands.append('')
        self.commands.append('date')

    def finish_job(self):
        '''Add last command and write to file'''
        self.commands.append('')
        self.commands.append('date')
        job_string = '\n'.join(self.commands)
        f = open(self.path, 'w')
        f.write(job_string)
        f.close()
        # change permisisons
        os.system('chmod +x %s' % (self.path))
        

class Submitter():
    '''
    Class with everything needed to submit
    Usage:
      search = Submitter('crux', 'kalkyl')
      search.program.target_database('target.fasta') ...
      search.submit()
    '''

    def __init__(self, program='crux', computer='local'):
        '''Sets the type of search to be submitted'''
        # Set computer
        if computer.lower() == 'local':
            self.computer = Local()
        elif computer.lower() == 'ferlin':
            self.computer = Ferlin()
        elif computer.lower() == 'kalkyl':
            self.computer = Kalkyl()
        else:
            sys.exit('ERROR: Computer %s not known' % computer)
        # Set program
        if program.lower() == 'crux':
            self.program = Crux(self.computer)
        elif program.lower() == 'msgf':
            self.program = MSGF(self.computer)
        elif program.lower() == 'tandem':
            self.program = Tandem(self.computer)
        else:
            sys.exit('ERROR: Program %s not known' % program)
        
    def submit(self):
        '''Writes the actual job scripts, to be submitted'''
        jobs = self.program.prepare_jobs()
        for job in jobs:
            if self.computer.just_test:
                print 'Did not submit job %s, because \'just_test = True\'' % (job.path)
            else:
                self.computer.submit(job.path)
                print 'Submitted job %s' % job.path
        if jobs == []:
            sys.exit('ERROR: No jobs were submitted')
                                         

if __name__ == '__main__':
    main()
