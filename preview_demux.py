#!/usr/bin/env python
# execute_pipeline.py - top level code to run hci_demux pipelines.
import sys
import os
import xml.dom.minidom

from logger import Logger
import pipelineparams as params
import gnomex

def SummarizeRunInfo(run_full_path):
  """
  Given the full path name of a run folder, this function returns the
  number of data reads for the sequencing run by parsing the RunInfo.xml
  file.
  """
  # Parse the RunInfo.xml file.
  configfile = os.path.join(run_full_path, "RunInfo.xml")

  if not os.path.exists(configfile):
    Log().Log('File %s does not exist!' %s)
    return None
  doc = xml.dom.minidom.parse(configfile)


  Logger().Log("Flow cell configuration: ")
  # Get the Read elements from the document.
  reads=doc.getElementsByTagName("Read")
  keys = ['Number', 'NumCycles', 'IsIndexedRead']
  Logger().Log('Number | NumCycles | IndexRead')
  for read in reads:
    vals = [read.getAttribute(key) for key in keys]
    Logger().Log("%s      | %s       | %s" %  (vals[0],vals[1],vals[2]))

  return None

def IdentifyPipeline(run_id):
  """
  IdentifyPipeline(run_id) - given id of a sequencing run (which is the
  directory name of the run folder, not full path) returns
  the class that is to process the run, or None of the type of run
  isn't recognized.
  """
  # Connect to GNomex.
  connection = gnomex.GNomExConnection().GnConnect(params.db_user,params.db_password)
  c = connection.cursor()
  # Select the sequencing application name(s) in use on the flow cell.
  query="""select distinct fcc.number as Lane, app.application as Application, req.number as Sample_Project
                from application as app
                join seqlibprotocolapplication as slpapp
                    on app.codeapplication = slpapp.codeapplication
                join seqlibprotocol as slp
                    on slp.idseqlibprotocol = slpapp.idseqlibprotocol
                join sample
                    on sample.idseqlibprotocol = slp.idseqlibprotocol
                join request as req
                    on sample.idRequest = req.idRequest
                join sequencelane as seqlane
                    on seqlane.idsample = sample.idsample
                join flowcellchannel as fcc
                    on seqlane.idflowcellchannel = fcc.idflowcellchannel
                    where fcc.filename = '%s'\n""" % run_id
  c.execute(query)
  recs=c.fetchall()
  if len(recs) == 0:
    Logger().Log("No records found for run id: %s" % run_id)
    return None

  # recs is a list of tuples each containing a unicode string. Convert
  # this into a new-line delimited string.
  csv_string = """Lane,Application,Sample_Project"""
  print(csv_string)
  for rec in recs:
    out = ''
    for r in rec:
      out+= str(r) + ","
    print(out)
    csv_string+=out+'\n'

  return None

def main():
  """
  main() function for previewing demultiplexing requirements.
  """
  run_id=sys.argv[1]
  run_full_path='/Repository/SeqStore04/IlluminaRuns/' + run_id
  if os.path.isdir(run_full_path):
    SummarizeRunInfo(run_full_path)
  else:
    Logger().Log("Directory %s does not exist." % run_full_path)
  print("\n")
  IdentifyPipeline(run_id)

if __name__ == "__main__":
  main()
