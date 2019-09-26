
# Demultiplexing with alademux

The tool `alademux` is meant to allow you to customize your demultiplexing
needs in a la carte manner, hence the name. What follows are examples of how to
demultiplex the 5 types of commonly encountered scenarios in the bioinformatics
core.

Log in to hci-zion. The software `alademux` software was built and tested
 on python3.

```bash
module load python3
```

## Overview

In this tutorial assume the git repository has been cloned to your home directory.
Review the usage of the main `alademux` command:

```bash
python3 alademux/alademux.py -h
```

You will see the basis of demultiplexing starts with a required argument of the
Illumina instrument run id.  Every other argument is optional; for example,
you may specify a subset of lanes to demultiple. By default the type of
demultiplexing defaults to the standard scenario.  One can override the
input and output file paths. Importantly, the last argument allows you to
pass any valid `bcl2fastq` argument to the demultiplexing run.

```
usage: alademux.py [-h] -r run [-l [lane [lane ...]]]
                   [-t [{standard,10x,10x-atac,patchpcr}]] [-i [RUN_PATH]]
                   [-o [OUT_PATH]] [-b ...]

optional arguments:
  -h, --help            show this help message and exit
  -r run, --run_id run  Name of Illumina Instrument Run (e.g.
                        190624_A00421_0081_AHC7G3DRXX)
  -l [lane [lane ...]], --lanes [lane [lane ...]]
                        List of lanes to demultiplex. e.g. -l 1 4 8 or -l 2
  -t [{standard,10x,10x-atac,patchpcr}], --type [{standard,10x,10x-atac,patchpcr}]
                        Library type, valid options include: standard, 10x,
                        10x-atac, patchpcr
  -i [RUN_PATH], --run_path [RUN_PATH]
                        Path to directory with Illumina Instrument Run
  -o [OUT_PATH], --out_path [OUT_PATH]
                        Path to demultiplexing results
  -b ..., --bcl2fastq ...
                        Any argument to pass to bcl2fastq. Note, these args
                        are not sanitized.
```

##  Preview run contents

Suppose we want to demultiplex a recent NovaSeq run. Preview the library types
on a (library,project) basis.

```bash
python3 alademux/preview_demux.py 190920_A00421_0112_BHKTM3DMXX
```

The preview shows

```
2019-09-25 17:12 Flow cell configuration:
2019-09-25 17:12 Number | NumCycles | IndexRead
2019-09-25 17:12 1      | 151       | N
2019-09-25 17:12 2      | 10       | Y
2019-09-25 17:12 3      | 10       | Y
2019-09-25 17:12 4      | 151       | N


>Lane,Application,Sample_Project
1,10X Genomics Single Cell 3' Gene Expression Library Prep v2,17385R,
1,10X Genomics Single Cell 3' Gene Expression Library Prep v3,17334R,
1,10X Genomics Single Cell Mouse T Cell V(D)J Enrichment,17350R,
2,Illumina TruSeq DNA Nano LIbrary Prep (450 bp mean insert size) with UDI,16274R,
2,Illumina TruSeq DNA Nano LIbrary Prep (450 bp mean insert size) with UDI,17353R,
2,Illumina TruSeq Stranded mRNA Library Prep (RIN 8-10) with UDI,17354R,
2,Illumina TruSeq Stranded Total RNA Library Prep Ribo-Zero Gold (human, mouse, rat etc) ,16180R,
2,NEBNext Ultra II Directional RNA Library Prep with rRNA Depletion Kit (human,mouse,rat),16093R,
2,NEBNext Ultra II Directional RNA Library Prep with rRNA Depletion Kit (human,mouse,rat),16298R,
```
In this first example, we encounter a run with 10x Genomics libraries on lane 1
as well as other Illumina/NEB libraries fitting the standard scenario on lane 2. Let's explore how to demultiplex lane 2 first.

## Standard demultiplex

Set up the standard demultiplexing scripts for lane 2 with the following command.

```bash
python3 alademux/alademux.py -r 190920_A00421_0112_BHKTM3DMXX -l 2
```

The output will give you the commands to start demultiplexing from the
target output directory.

```
Start demultiplexing :
cd /path/to/demux/190920_A00421_0112_BHKTM3DMXX/20190925-173229
nohup ./demux.sh &
```

The reason for not starting the demultiplexing in a single command is if
the user wants to inspection the script or sample sheet before launching the demultiplexing job.


## 10x Genomics Libraries

### [5'/3' expr. and V(D)J]

The 10x 5',3' and V(D)J libraries are handled a little differently than the
standard scenario. Trim the reads based on the 10x Genomics docs [here]((https://support.10xgenomics.com/permalink/3IQFKIvEuskMoEWkWUis2s)) and [here](https://support.10xgenomics.com/single-cell-vdj/sequencing/doc/specifications-sequencing-requirements-for-single-cell-vdj):

+ **R1** Keep the first 28 bases. The 3' gene expression libraries (v3 chemistry) uses 28 bases for R1, 5' gene expression and V(D)J uses  26 bases.
+ **R2** Keep 8 bp for the sample index.
+ **R3** Do not use read; the `alademux` software will add the `--ignore-dual-index`
option automatically to the script.
+ **R4** Use the full length except for the last cycle which is generally lower
quality.

The `--use-bases-mask` argument is how the trimming instructions get passed
along. Recall these instructions are specific to the flow cell, which is why
we first preview the flow cell's contents to get the number of cycles correct.


```bash
python3 alademux/alademux.py -r 190920_A00421_0112_BHKTM3DMXX -l 1 \
  -t 10x \
  -b "--use-bases-mask='Y28n*,i8nn,n*,Y150n'"
```

Behind the scenes, this command generates the sample sheet from GNomEx (see file GNomEx_SampleSheet.csv) and swaps the nucleotide based index with the 10x Genomics tag (see file SampleSheet.csv).

```
Start demultiplexing :
cd /path/to/demux/190920_A00421_0112_BHKTM3DMXX/20190925-174609
nohup ./demux.sh &
```

Listing the files in the demultiplexing path you can see the demultiplexing
script as well as the 2 sample sheets.
```bash
ls /path/to/demux/190920_A00421_0112_BHKTM3DMXX/20190925-174609
```

```
demuxer.sh  GNomEx_SampleSheet.csv  SampleSheet.csv
```

### ATAC

While there are no recent 10x ATAC libraries to demonstrate the software,
we can still show what commands to run in this demultiplexing scenario.

```bash
python3 alademux/preview_demux.py 190710_A00421_0084_AHC7K3DRXX
```
Because the run has been archived to tape, the flow cell configuration is
hidden from view.
å
```
2019-09-25 18:05 Directory /path/to/illumina_runs/190710_A00421_0084_AHC7K3DRXX does not exist.


Lane,Application,Sample_Project
1,10X Genomics Chromium Single Cell ATAC Library Prep,16164R1,
1,10X Genomics Chromium Single Cell ATAC Library Prep,16237R1,
2,10X Genomics Chromium Single Cell ATAC Library Prep,16164R1,
2,10X Genomics Chromium Single Cell ATAC Library Prep,16237R1,
```

Demultiplex by changing a few subtle options in the type argument and the trimming arguments as instructed by [the 10x Genomics online docs for
ATAC libraries](https://support.10xgenomics.com/single-cell-atac/sequencing/doc/specifications-sequencing-requirements-for-single-cell-atac).

```bash
python3 alademux/alademux.py -r 190710_A00421_0084_AHC7K3DRXX  \
  -t 10x-atac \
  -b "--use-bases-mask='Y50,I8,Y16,Y49'"
```

### Long ranger

If you have long ranger libraries,

1. Change the demultiplexing type:

  `-t='10x-long'`

2. Per the [10x Genomics website](https://support.10xgenomics.com/genome-exome/sequencing/doc/specifications-sequencing-requirements-for-genome-and-exome) trim the reads:

 `-b "--use-bases-mask='Y150,I8,n*,Y150'"`

As of this writing, QC metrics are not generated for the long ranger results
and will not be reported to the lab.

## Patch PCR

A MiSeq run contains KT-Varley's Patch PCR libraries. Demultiplexing follows similar steps as before. First preview the run so we can trim  reads appropriately.

```bash
python3 alademux/preview_demux.py 190917_M00736_0336_MS8428226-300V2
```

```
2019-09-25 17:52 Flow cell configuration:
2019-09-25 17:52 Number | NumCycles | IndexRead
2019-09-25 17:52 1      | 151       | N
2019-09-25 17:52 2      | 9       | Y
2019-09-25 17:52 3      | 9       | Y
2019-09-25 17:52 4      | 151       | N


Lane,Application,Sample_Project
1,Varley Patch PCR Target Enrichment,17400R,
```

In this case, the script must have a 'y' for the **R3** read, to capture the UMI
sequence. This run is from a paired-end configuration, but sometimes Patch PCR
runs have a single end configuration.

```bash
python3 alademux/alademux.py -r 190917_M00736_0336_MS8428226-300V2 \
  -t patchpcr \
  -b --use-bases-mask='y150n,i8n,y9,y150n'
```
