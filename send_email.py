import os
import sys
import smtplib
import string
import argparse
import pipelineparams as params
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.application import MIMEApplication

class Emailer:
  """Base class for objects that need to email stuff."""

  def send_msg( self, from_addr, to_addr, subject, message, attachments=[], msgtype="plain" ):
    msg = MIMEMultipart()
    msg['Subject'] = subject
    msg['From'] = from_addr
    msg['To'] = ', '.join(to_addr)
    msg.attach(MIMEText(message))

    for f in attachments:
      fp=open(f,"rb")
      msg.attach(MIMEApplication(fp.read(),
        Content_Disposition='attachment; filename="%s"' % os.path.basename(f),
        Name=os.path.basename(f)))

    s = smtplib.SMTP('hci-mail.hci.utah.edu')
    # Send the message.
    retval = s.sendmail( from_addr, to_addr, msg.as_string() )
    s.quit()
    return retval

# Example of use:
def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("runid", help="Illumina run ID",
      type=str)
  parser.add_argument('attachments', nargs='*', type=str)
  args = parser.parse_args()
  flowcell = args.runid
  attachments = args.attachments
  for file in attachments:
    if not os.path.exists(file):
      print("File '%s' not found." % file)
      sys.exit(1)
  from_addr = params.from_addr
  subject = "[Demux complete] %s" % flowcell
  message = "Demultiplexed output was copied to GNomEx."
  to_addr = params.to_addr 
  unsubscribe = 'If you no longer wish to receive this automated email message'
  unsubscribe+= ' please email ' + params.unsubscribe_addr
  message = message +  "\n\n" + unsubscribe 
  Emailer().send_msg( from_addr, to_addr, subject, message, attachments )

if __name__ == "__main__":
  main()
