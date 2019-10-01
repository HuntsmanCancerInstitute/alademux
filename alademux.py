#!/usr/bin/env python

import argparse
import os
from datetime import datetime
from iemwriter import IEMWriter,MissingBarcodeError
from demuxscripter import DemuxScripter
import pipelineparams as params

p = argparse.ArgumentParser()
p.add_argument('-r', '--run_id', required=True,
               type=str,
               help = '''Name of Illumina Instrument Run
                    (e.g. 190624_A00421_0081_AHC7G3DRXX)''',
                metavar='run')
p.add_argument('-l', '--lanes', nargs='*', type = int,
                    help = '''List of lanes to demultiplex.
                    e.g. -l 1 4 8 or -l 2''', metavar='lane')
p.add_argument('-t', '--type', nargs='?', type=str,
                    help = '''Library type, valid options include:
                    standard, 10x, 10x-atac, 10x-long, patchpcr, nano''',
                    choices = ['standard', '10x', '10x-atac',
                               '10x-long', 'patchpcr', 'nano'],
                    default = 'standard')
p.add_argument('-i', '--run_path', nargs='?',
                    help = 'Path to directory with Illumina Instrument Run',
                    default = params.illumina_run_path)
p.add_argument('-o', '--out_path', nargs='?',
                    help = 'Path to demultiplexing results',
                    default = params.demux_result_path)
p.add_argument('-b','--bcl2fastq', nargs=argparse.REMAINDER,
               help='''Any argument to pass to bcl2fastq. Note, these args
               are not sanitized.''')
p.add_argument('-n','--nextera', action = 'store_true')
               help='''Add flag to indicate nextera adapters for MiSeq Nano Runs.
               Otherwise, Illumina adapters are assumed to be the case.''')
args = p.parse_args()

# sanitize lanes
if args.lanes:
    valid_lanes = list(range(1,9))
    args.lanes = [lane for lane in args.lanes if lane in valid_lanes]
    # drop duplicate lanes
    args.lanes = list(dict.fromkeys(args.lanes))
else:
    args.lanes = None

# check the run folder exists
run_folder_path = os.path.join(args.run_path, args.run_id)
if not os.path.exists(run_folder_path):
    raise OSError('Run path does not exist: %s' % run_folder_path)
# create unique output directory
date_folder = datetime.now().strftime("%Y%m%d-%H%M%S")
out_demux_path = os.path.join(args.out_path, args.run_id, date_folder)
os.makedirs(out_demux_path, exist_ok=True)

if args.type == 'nano':
    print("Nano: No SampleSheet.csv being written. User must specify SampleSheet.csv")
else:
    # write SampleSheet.CSV
    iem = IEMWriter(args.run_id, out_demux_path, args.lanes)
    write_header = args.type == 'standard' or args.type == 'patchpcr'
    success = iem.write_sample_sheet(iem_header=write_header)
    if not success:
        raise Exception('Do not execute run until GNomEx records are available.')

# write script for demultiplexing demux.sh
demuxer = DemuxScripter(args.run_id, run_folder_path, out_demux_path,
                        args.type, args.bcl2fastq, args.nextera)
demuxer.write_demux_script()

print("Start demultiplexing with the following commands: ")
print("cd " + out_demux_path)
print("nohup ./demuxer.sh & ")
