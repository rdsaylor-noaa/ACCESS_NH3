#!/usr/bin/env python
#
# pltobs.py - plot NH3 flux observations
#
# Rick D. Saylor, April 2015
#
import os
import sys
import matplotlib.pylab as plt
import matplotlib.dates as mdates
import numpy as np
import csv
from datetime import datetime
import seaborn as sns

def plotall(astartdt, aenddt, otype):

   # extra settings from seaborn
   sns.set_context("talk")
   sns.set_style("ticks") 

   # set colors
   colors = ['black', 'peru', 'royalblue', 'red', 'green', 'gray', 'darkviolet']

   # set formatting parameters
   tfsize = 18    # plot title font size
   lfsize = 14    # legend font size
   msize  = 4     # marker size
   yfsize = 18    # y-axis title font size
   xfsize = 18    # x-axis title font size
   tlbsize = 16   # tick label size

   startdt = datetime.strptime(astartdt, '%Y-%m-%d-%H:%M')
   enddt = datetime.strptime(aenddt, '%Y-%m-%d-%H:%M')

   # get flux gradient fluxes
   fname1 = './Illinois2014_30MinDatav02.csv'
   fh1 = open(fname1, 'rU')
   reader = csv.DictReader(fh1)
   odts = []
   fgobs = []
   npts=0
   for row in reader:
      dt = datetime.strptime(row['TIMESTAMP'], '%m/%d/%y %H:%M')
      fg = float(row['NH3 Flux (ug/m2s)'])
      if (dt >= startdt and dt <= enddt and fg > -6999.):
         odts.append(dt)
         fgobs.append(1000.*fg)   # convert ug to ng
         npts+=1

   # get REA fluxes
   fname2 = './Illinois-2014-REA.csv'
   fh2 = open(fname2, 'rU')
   reader = csv.DictReader(fh2)
   rdts = []
   fobs = []
   nrea=0
   for row in reader:
      ssdt = row['Date']+' '+row['Start']
      esdt = row['Date']+' '+row['End']
      sdt = datetime.strptime(ssdt, '%m/%d/%y %H:%M')
      edt = datetime.strptime(esdt, '%m/%d/%y %H:%M')
      flx = row['Flux (ng/m2s)']
      if (sdt >= startdt and edt <= enddt):
         rdts.append(sdt)
         rdts.append(edt)
         fobs.append(flx)
         fobs.append(flx)
         nrea+=1 

   print 'nrea=',nrea

   fig, ax = plt.subplots(figsize=(22,6))

   ax.plot(odts, fgobs, color=colors[6], linestyle='None', marker='o', markersize=msize)

   for n in range(nrea):
      m=2*n
      print m, m+1
      print rdts[m:m+2]
      print fobs[m:m+2]
      ax.plot(rdts[m:m+2], fobs[m:m+2], color=colors[3], linestyle='-', linewidth=2.0, marker='None')

   # details
   days = mdates.DayLocator()
   ax.xaxis.set_major_locator(days)
   ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
   if (npts > 49):
      hours = mdates.HourLocator(byhour=range(24), interval=4)
   else:
      hours = mdates.HourLocator(byhour=range(24), interval=1)
   ax.xaxis.set_minor_locator(hours)

   plt.grid(b=True, which='major', color='gray', linewidth=1.0, alpha=0.5)
   plt.grid(b=True, which='minor', color='gray', linewidth=0.5, alpha=0.25)

   ax.tick_params(which="both", direction="out")
   ax.tick_params(which="major", length=8)
   ax.tick_params(which="minor", length=5)
   ax.tick_params(which="both", labelsize=tlbsize)

   plt.ylabel('ng m$^{-2}$ s$^{-1}$', fontsize=yfsize, labelpad=0)
   plt.title('Measured NH$_3$ Flux', fontsize=tfsize)

   if (otype == 'pdf'):
      plt.savefig('./pltobs.pdf')
   elif (otype == 'png'):
      plt.savefig('./pltobs.png')
   else:
      plt.show()


# main
def main(argv=None):
   if argv is None:
      argv = sys.argv

   # enforce proper number of arguments
   if len(argv) != 4:
      print "usage: %s STARTDT ENDDT OUTTYPE" % os.path.basename(sys.argv[0])
      print "  where,"
      print "          STARTDT as YYYY-mm-dd-HH:MM"
      print "          OUTDT as YYYY-mm-dd-HH:MM"
      print "          OUTTYPE = 'x11', 'png', or 'pdf'"
      return 2

   # get filename
   astartdt = argv[1]
   aenddt   = argv[2]
   otype    = argv[3]

   # call routine
   plotall(astartdt, aenddt, otype)

   # all done!
   return 0

if __name__ == "__main__":
    sys.exit(main())
