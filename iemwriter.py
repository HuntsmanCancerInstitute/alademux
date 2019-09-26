#!/usr/bin/env python
# create_samplesheet.py - given the name of a run folder, create the sample
# sheet for the run.

import os
import sys
import datetime as dt
import gnomex
import pipelineparams as params

class MissingBarcodeError(LookupError):
    '''Barcode missing in lane where required for demultiplexing.'''

class IEMWriter:
    def __init__(self, run_id, file_path, lanes= None):
        self.run_id = run_id
        self.file_path = '/'.join([file_path,'SampleSheet.csv'])
        self.lanes = lanes
        self.query = """select flowcell.barcode,
        flowcellchannel.number,
        sample.number,
        genomebuild.genomebuildname,
        sample.barcodesequence,
        appuser.firstname,
        appuser.lastname,
        request.number,
        sample.barcodesequenceb
        from flowcell
        join flowcellchannel on flowcellchannel.idflowcell=flowcell.idflowcell
        join sequencelane on sequencelane.idflowcellchannel =
        flowcellchannel.idflowcellchannel
        join sample on sequencelane.idsample = sample.idsample
        left outer join genomebuild on sequencelane.idgenomebuildalignto =
        genomebuild.idgenomebuild
        join request on sequencelane.idrequest = request.idrequest
        join appuser on request.idappuser = appuser.idappuser
        where flowcellchannel.filename = '%s'\n""" % self.run_id
        if self.lanes:
            self.query += "and flowcellchannel.number in (%s)\n" % \
            str(','.join(list(map(str,self.lanes))))
            self.query += "order by flowcellchannel.number, sample.number;"
    def execute_query(self):
        g = gnomex.GNomExConnection()
        connection = g.GnConnect(params.db_user,params.db_password,asdict=False)
        c = connection.cursor()
        # if self.lanes:
        #    # Append lane selection to where clause.
        #    self.query += "and flowcellchannel.number in (%s)\n" % str(','.join(map(str,self.lanes)))
        #    self.query += "order by flowcellchannel.number, sample.number;"
        try:
            c.execute(self.query)
        except pymssql.OperationError:
            sys.stderr.write( self.query )
            raise
        records = c.fetchall()
        connection.close()
        return(records)
    def transform_records(self, records):
        # Count how many samples are in each lane.
        sample_count = {}
        for rec in records:
            try:
                lane = rec[1]
            except KeyError:
                print(rec)
                raise
            try:
                sample_count[lane]+=1
            except KeyError:
                sample_count[lane] = 1
        # Write the records.
        csv_rows = ''
        for i,rec in enumerate(records):
            row = self.clean_row(rec, sample_count)
            csv_rows += ','.join(row) + '\n'
        return(csv_rows)
    def erase_commas(self, field):
        """Removes all commas."""
        if field:
            return ''.join(field.split(','))
        else:
            return field
    def clean_row(self, rec, sample_count):
        """
        Handle edge cases in records.
        """
        dr = {'Lane' : self.erase_commas(str(rec[1])),
        'Sample_ID' : self.erase_commas(rec[2]),
        'Sample_Project' : self.erase_commas(rec[7]),
        'index' : self.erase_commas(rec[4]),
        'index2' : self.erase_commas(rec[8])
        }
        # replace nones with emptry string
        # no trailing white space
        for k,v in dr.items():
            if v is None:
                dr[k] = ''
            dr[k] = dr[k].strip()
        # If lane has a single sample, don't write its barcode.
        lane = rec[1]
        if sample_count[lane] == 1:
            dr['index'] = ''
            dr['index2'] = ''
        elif sample_count[lane] > 1 and not len(dr['index']):
            raise MissingBarcodeError('''Multiple samples in lane but missing
            barcode for sample: ''' + dr['Sample_ID'])
        # Drop potential trailing integers on sequence request number
        # e.g. 1480R1 -> 14806R
        # still works on clean sequence request number
        dr['Sample_Project'] = dr['Sample_Project'].split('R')[0] + 'R'
        # build the row
        row = [dr['Lane'], dr['Sample_ID'],
               '_'.join([dr['Sample_ID'], self.run_id]),
               dr['Sample_Project'],
               dr['index'], dr['index2']]
        return(row)
    def generate_header(self, iem_header=True):
        """ Write IEM header (optional) and column names """
        data_cols = ['Lane','Sample_ID','Sample_Name','Sample_Project','index','index2']
        if iem_header:
            header = '[Header]\n'
            raw_date = dt.date.today()
            todays_date = '/'.join([str(raw_date.month),str(raw_date.day),str(raw_date.year)])
            h1 = [('Date', todays_date),
                  ('Workflow', 'GenerateFASTQ'),
                  ('Application', 'FASTQ Only'),
                  ]
            for k,v in h1:
                header+=str(k) + ',' + str(v) + '\n'
            header+='\n[Settings]\n'+'\n[Data]\n'
            header+=','.join(data_cols) + '\n'
        else:
            header = ','.join(data_cols) + '\n'
        return(header)
    def write_sample_sheet(self, iem_header=True):
        """
        1. Query GNomEx DB.
        2. Write out a config file (i.e. sample sheet) in the style of
        Illumina Experiment Manager.
        """
        db_records = self.execute_query()
        if not len(db_records):
            print("No GNomEx records found for run id: %s\n" % self.run_id)
            return(False)

        records = self.transform_records(db_records)
        if iem_header:
            header = self.generate_header()
        else:
            header = self.generate_header(iem_header=False)
        sample_sheet_content = header + records
        with open(self.file_path, "w") as file:
            file.write(sample_sheet_content)
        return(True)
