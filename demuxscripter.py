#!/usr/bin/env python

import multiprocessing
import os
import pandas as pd
import stat
import xml.dom.minidom
import re
import pipelineparams as params

class DemuxScripter:
    def __init__(self, run_id, run_folder_path, out_path, demux_type,
                 nextera, args_bcl2fastq=None):
        self.run_id = run_id
        self.in_path = run_folder_path
        self.out_path = out_path
        self.out_file = os.path.join(out_path, 'demuxer.sh')
        self.type = demux_type
        check = demux_type in ['standard', '10x', '10x-atac', 'patchpcr',
                               'nano']
        if not check:
            raise NameError('Type of demultiplexing not available.')
        self.args = args_bcl2fastq
        # boolean True if nextera adapters used for miseq nano
        self.nextera = nextera
        self.nprocs = str(multiprocessing.cpu_count() - 1)
        transfer_path = os.path.join(params.alademux_path, 'templates')
        self.transfer = {'standard': os.path.join(transfer_path,
                                                  'standard_transfer.sh'),
                         '10x': os.path.join(transfer_path,
                                             '10x_transfer.sh'),
                         '10x-atac': os.path.join(transfer_path,
                                                  '10x_transfer.sh'),
                         '10x-long': os.path.join(transfer_path,
                                                  '10x_transfer.sh'),
                         'patchpcr': os.path.join(transfer_path,
                                                  'standard_transfer.sh'),
                         'nano': os.path.join(transfer_path,
                                                  'nano_report.sh')
                         }
        self.script = '#!/bin/bash\n\nset -e\n'
        self.script += 'IN_BCL2FASTQ=' + self.in_path + '\n'
        self.script += 'OUT_BCL2FASTQ=' + self.out_path + '\n'
        self.script += 'RUN=' + self.run_id + '\n\n'
        self.script += 'cd $OUT_BCL2FASTQ\n\n'
        if self.type == 'nano':
            if self.nextera:
                print("Trimming Nextera adapters...")
                self.script += 'ADAPTER=CTGTCTCTTATACACATCT\n'
            else:
                print("Trimming Illumina adapters...")
                self.script += 'ADAPTER=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC\n'
    def count_index_reads(self):
        configfile = os.path.join(self.in_path, 'RunInfo.xml')
        if not os.path.exists(configfile):
          raise Exception('Flowcell File RunInfo.xml does not exist.')
        doc = xml.dom.minidom.parse(configfile)
        reads=doc.getElementsByTagName("Read")
        key = 'IsIndexedRead'
        n_index = 0
        for read in reads:
            if read.getAttribute(key) == 'Y':
                n_index = n_index + 1
        return(n_index)
    def tag_swapper(self):
        fname_samples = os.path.join(self.out_path, 'SampleSheet.csv')
        samples = pd.read_csv(fname_samples)
        # same barcodes for expression/VDJ as longranger
        if self.type == '10x' or self.type == '10x-long':
            fname_tags = os.path.join(params.alademux_path,'10x_tags.tsv')
        elif self.type == '10x-atac':
            fname_tags = os.path.join(params.alademux_path,'10x_atac_tags.tsv')
        else:
            raise Exception('Only 10x tags should be swapped for indices.')
        tags = pd.read_table(fname_tags)
        # left join samples with tags
        updated = pd.merge(samples,tags, on='index', how='left')
        # report samples missing tags
        missing = updated['tag'].isna()
        if any(missing):
            miss_samples=updated[missing]
            print('Samples missing 10x Genomics Barcodes are listed below: ')
            for ms in miss_samples['Sample_ID']:
                print(ms)
        # the nucleotide information is no longer needed
        updated = updated.drop(['index'], axis=1)
        # rename the tag to index
        updated.rename(columns={'tag': 'index'}, inplace=True)
        # set to original order
        cols=list(samples.columns)
        updated = updated[cols]
        # re-name original IEM sample sheet
        os.path.basename
        fname_previous = os.path.join(self.out_path, 'GNomEx_SampleSheet.csv')
        os.rename(fname_samples, fname_previous)
        # over-write filename
        updated.to_csv(path_or_buf=fname_samples, index=False)
    def check_use_bases(self):
        """Make sure the use-bases-mask argument has been specified."""
        if self.args:
            patt = '--use-bases-mask='
            is_matched = [re.search(patt, arg) for arg in self.args]
            if not any(is_matched):
                raise Exception("""User must specify --use-bases-mask=
                                for demultiplexing 10x data.""")
        else:
            raise Exception("""User must specify --use-bases-mask=
                            for demultiplexing 10x data.""")
        return(None)
    def standard(self):
        cmd = [params.bcl2fastq_path,
               '--runfolder-dir', '$IN_BCL2FASTQ',
               '--sample-sheet', os.path.join('$OUT_BCL2FASTQ',
                                              'SampleSheet.csv'),
               '--output-dir', '$OUT_BCL2FASTQ',
               '--processing-threads', self.nprocs,
               ]
        return(cmd)
    def tenx(self):
        """Demultiplexing 10x Genomics libraries."""
        if self.type == '10x':
            tool = params.cellranger_path
        elif self.type == '10x-atac':
            tool = params.cellranger_atac_path
        elif self.type == '10x-long':
            tool = params.longranger_path
        else:
            raise Exception("Demux type incorrectly specified.")
        cmd = [tool,
               'mkfastq',
               '--id=10x_mkfastq_log',
               '--qc',
               '--run=$IN_BCL2FASTQ',
               '--output-dir=$OUT_BCL2FASTQ',
               '--samplesheet=', os.path.join('$OUT_BCL2FASTQ',
                                              'SampleSheet.csv'),
               ]
        # ignore the dual index if it exists for 5',3',VDJ libraries
        # or long ranger
        n_index = self.count_index_reads()
        flag = self.type == '10x' or self.type == '10x-long'
        if flag and n_index == 2:
            cmd.append('--ignore-dual-index')
        # check for required argument
        self.check_use_bases()
        # swap the 1st nucleotide index with the 10x tag
        self.tag_swapper()
        return(cmd)
    def keep_short_reads(self):
        """Add bcl2fastq options to keep short read lengths."""
        cmd = self.standard()
        short = ['--minimum-trimmed-read-length 0',
               '--mask-short-adapter-reads 0']
        cmd.extend(short)
        return(cmd)
    def patchpcr(self):
        # check for required argument
        self.check_use_bases()
        # remove NNNNs from SampleSheet index2
        fname_samples = os.path.join(self.out_path,'SampleSheet.csv')
        f = open(fname_samples,'r')
        filedata = f.read()
        f.close()
        newdata = re.sub("NNNN+", "", filedata)
        f = open(fname_samples,'w')
        f.write(newdata)
        f.close()
        return(self.keep_short_reads())
    def nano(self):
        return(self.keep_short_reads())
    def get_cmd(self):
        """Construct the demultiplexing command based on the type."""
        if self.type == 'standard':
            cmd = self.standard()
        elif self.type == '10x' or self.type == '10x-atac':
            cmd = self.tenx()
        elif self.type == 'patchpcr':
            cmd = self.patchpcr()
        elif self.type == 'nano':
            cmd = self.nano()
        # add the additional arguments to the command.
        if self.args:
            cmd.extend(self.args)
        return(cmd)
    def write_demux_script(self):
        """Combine demultiplexing command with the appropriate transfer script."""
        self.script += ' '.join(self.get_cmd())
        with open(self.out_file, 'w+') as outfile:
            outfile.write(self.script)
            outfile.write('\n\n# Transfer data to GNomEx\n')
            with open(self.transfer[self.type]) as infile:
                for line in infile:
                    outfile.write(line)
        # change to executible file
        st = os.stat(self.out_file)
        os.chmod(self.out_file, st.st_mode | stat.S_IEXEC)
