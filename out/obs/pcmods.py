#!/usr/bin/env python
#
# pcmods.py ... plot comparison of model flux results against observations
#
# Rick D. Saylor, September 2017
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
   colors = ['black', 'peru', 'royalblue', 'red', 'magenta', 'green', 'orange', 'darkviolet']

   # set formatting parameters
   tfsize = 18    # plot title font size
   lfsize = 10    # legend font size
   msize  = 4     # marker size
   yfsize = 18    # y-axis title font size
   xfsize = 18    # x-axis title font size
   tlbsize = 16   # tick label size

   startdt = datetime.strptime(astartdt, '%Y-%m-%d-%H:%M')
   enddt = datetime.strptime(aenddt, '%Y-%m-%d-%H:%M')

   # get all data
   fn = './AllModelFluxes.csv'
   fh = open(fn, 'rU')
   reader = csv.DictReader(fh)
   dts = []
   freas = []
   ffgs = []
   fsurf1s = []
   fsurf2s = []
   fsurf3s = []
   faccs = []
   famms = []
   npts = 0
   for row in reader:
      date = row['Date']
      time = row['Time']
      sdt = date+' '+time
      dt = datetime.strptime(sdt, '%m/%d/%y %H:%M:%S')
      fsurf1 = row['SURFATM-NH3-1 (ng/m2/s)']
      fsurf2 = row['SURFATM-NH3-2 (ng/m2/s)'] 
      fsurf3 = row['SURFATM-NH3-3 (ng/m2/s)']
      ffg  = row['FG (ng/m2/s)']
      frea = row['REA (ng/m2/s)']
      facc   = row['ACCESS-NH3-d (ng/m2/s)']
      famm   = row['AMM2L-a (ng/m2/s)']
      if (dt >= startdt and dt <= enddt):
         dts.append(dt)
         fsurf1s.append(float(fsurf1))
         fsurf2s.append(float(fsurf2))
         fsurf3s.append(float(fsurf3))
         freas.append(frea)
         ffgs.append(ffg)
         faccs.append(facc)
         famms.append(famm)
         npts+=1 

   afreas = np.zeros(npts)
   affgs  = np.zeros(npts)
   afaccs = np.zeros(npts)
   afamms = np.zeros(npts)

   for i in range(npts):
      if (faccs[i] == ''):
         faccs[i] = -6999.
      if (famms[i] == ''):
         famms[i] = -6999. 
      if (freas[i] == ''):
         freas[i] = -6999.
      if (ffgs[i] == ''):
         ffgs[i] = -6999. 

   for i in range(npts):
      if (faccs[i] <= -6999.):
         afaccs[i] = np.nan
      else:
         afaccs[i] = float(faccs[i])

   for i in range(npts):
      if (famms[i] <= -6999.):
         afamms[i] = np.nan
      else:
         afamms[i] = float(famms[i])

   for i in range(npts):
      if (ffgs[i] <= -6999.):
         affgs[i] = np.nan
      else:
         affgs[i] = float(ffgs[i])

   for i in range(npts):
      if (freas[i] <= -6999.):
         afreas[i] = np.nan
      else:
         afreas[i] = float(freas[i])

   fig, ax = plt.subplots(figsize=(12,6))

   ax.plot(dts, fsurf1s, color=colors[0], linestyle='-', marker='None', linewidth=2.0, label='SURFATM-NH3 constant')
   ax.plot(dts, fsurf2s, color=colors[1], linestyle='-', marker='None', linewidth=2.0, label='SURFATM-NH3 no inhibitor')
   ax.plot(dts, fsurf3s, color=colors[2], linestyle='-', marker='None', linewidth=2.0, label='SURFATM-NH3 w/ inhibitor')

   ax.plot(dts, afaccs, color=colors[5], linestyle='-', marker='None', linewidth=2.0, label='ACCESS-NH3')
   ax.plot(dts, afamms, color=colors[6], linestyle='-', marker='None', linewidth=2.0, label='AMM2L')
   ax.plot(dts, afreas, color=colors[3], linestyle='-', marker='None', linewidth=2.0, label='REA')
   ax.plot(dts, affgs, color=colors[4], linestyle='None', marker='o', markersize=msize, label='Flux Grad')

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
   plt.legend(loc=1, fontsize=lfsize, bbox_to_anchor=(0.275, 0.98))

   plt.ylabel('ng m$^{-2}$ s$^{-1}$', fontsize=yfsize, labelpad=0)
   plt.title('NH$_3$ Flux', fontsize=tfsize)

   if (otype == 'pdf'):
      plt.savefig('./pcmods.pdf')
   elif (otype == 'png'):
      plt.savefig('./pcmods.png')
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
