#!/usr/bin/env python
#
# ppsrc.py - plot NH3 source profile results from Coweeta ACCESS_NH3 simulations
#
# Rick D. Saylor, November 2017
#
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import math
from datetime import datetime
import csv
import seaborn as sns

def plotprofs(simname, otype, profno, hmax, hc, avgyn):

   # extra settings from seaborn
   sns.set_context("talk")
   sns.set_style("ticks")

   dirout='../'

   fig = plt.figure(figsize=(8, 10))
   ax  = fig.add_subplot(111)
   
   # set colors
   colors = ['gray', 'peru', 'brown', 'red', 'royalblue', 'green', 'violet', 'magenta', 'cyan']

   # set profile number sim hours
   profhrs = {
   '69':  [8, 9, 10, 11, 12],
   '70':  [13, 14, 15, 16, 17],
   '71':  [32, 33, 34, 35, 36],
   '72':  [37, 38, 39, 40, 41],
   '73':  [56, 57, 58, 59, 60],
   '74':  [61, 62, 63, 64, 65],
   '76':  [10, 11, 12, 13, 14],
   '77':  [15, 16, 17, 18, 19],
   '78':  [32, 33, 34, 35, 36],
   '79':  [37, 38, 39, 40, 41],
   '80':  [56, 57, 58, 59, 60],
   '81':  [61, 62, 63, 64, 65] 
   }

   # set formatting parameters
   tfsize = 18    # plot title font size
   tyloc  = 1.02  # title y location
   lfsize = 18    # legend font size
   yfsize = 18    # y-axis title font size
   xfsize = 18    # x-axis title font size
   tlmaj  = 6     # major tick length
   tlmin  = 4     # minor tick length
   tlbsize = 17   # tick label font size
   tlbpad = 3     # tick label padding
   lnwdth = 1.5   # linewidth

   # read elapsed hour/datetime key file
   indt = dirout+simname+'/ACCESS_timekey.dat'
   f_dt = open(indt)
   lines = f_dt.readlines()
   lines = lines[1:]    # ignore the header line
   f_dt.close()
   hrs = []
   dts = []
   for line in lines:
      data = line.split()
      hr   = data[0]
      date = data[1]
      time = data[2]
      dt = date+' '+time
      dts.append(datetime.strptime(dt, "%Y-%m-%d %H:%M:%S"))
      hrs.append(time[0:5])

   # sp
   fname = dirout+simname+'/s/s.dat'
   f_fname = open(fname)
   lines = f_fname.readlines()
   hdr   = lines[0]
   cols  = hdr.split()
   nhr   = len(cols) - 1
   lines = lines[1:]     # ignore the header line
   f_fname.close()
   nz = len(lines)

   z   = np.zeros(nz)
   nh3 = np.zeros( (nz, nhr) )
   iz = 0
   for line in lines:
      data = line.split()
      z[iz] = float(data[0])
      ihr = 0
      data = data[1:]
      for value in data:
         nh3[iz, ihr] = float(value)
         ihr+=1
      iz+=1 

   # average model profiles or not?
   nh3avg = np.zeros(nz)
   nh3sd = np.zeros(nz)
   phrs = profhrs[profno] 
   if (avgyn.upper() == 'Y'): 
      for iz in range(nz):
         nh3sum = 0.0
         for i in phrs: 
            nh3sum+= nh3[iz, i] 
         nh3avg[iz] = nh3sum/float(len(phrs)) 
         nh3sdsum = 0.0
         for i in phrs:
            nh3sdsum+= (nh3[iz, i] - nh3avg[iz])**2.0
         nh3sd[iz] = (nh3sdsum/float(len(phrs)-1))**0.5
      xerr=nh3sd
      hrstt=hrs[phrs[0]]
      hrend=hrs[phrs[len(phrs)-1]]
      labstr=datetime.strftime(dts[phrs[0]], "%Y-%m-%d") + ' ' + hrstt + '-' + hrend
      plt.errorbar(nh3avg, z, color='0.15', ecolor='0.7', xerr=xerr, linestyle='-', linewidth=lnwdth, label=labstr)
   else:
      ic=0
      for i in phrs:
         labstr = dts[i]
         plt.plot(nh3[:, i], z, color=colors[ic], linestyle='-', linewidth=lnwdth, label=labstr)
         ic+=1
         if (ic > len(colors)-1):
           ic=0

   plt.minorticks_on()
   plt.ylim(-0.1, hmax)

   # draw line showing canopy height
   xbnds = list(ax.get_xlim())
   xbnds[0] = round(xbnds[0], 2)
   print 'xbnds=', xbnds
   ahc = [hc, hc]
   plt.plot(xbnds, ahc, color='0.25', linestyle='--', linewidth=lnwdth)
   ax.set_xlim(xmin=xbnds[0], xmax=xbnds[1])

   # make it tidy
   plt.grid(b=True, which='major', color='gray', linewidth=1.0, alpha=0.5)
   plt.grid(b=True, which='minor', color='gray', linewidth=0.5, alpha=0.25)

   ax.tick_params(which="both", direction="out")
   ax.tick_params(which="major", length=tlmaj)
   ax.tick_params(which="minor", length=tlmin)
   ax.tick_params(which="both", labelsize=tlbsize, pad=tlbpad)

   plt.legend(loc=1, fontsize=lfsize, bbox_to_anchor=(0.65, 0.93))

   plt.xlabel('ng m$^{-3}$ s$^{-1}$', fontsize=xfsize, labelpad=10)
   plt.ylabel('z (m)', fontsize=yfsize, labelpad=10)
   plt.title('NH$_3$ Canopy SourceProfile - '+simname, fontsize=tfsize, y=tyloc)

   oname = './img/'+simname+'_'+str(profno)+'_ppsrc'
   if (otype == 'pdf'):
      plt.savefig(oname+'.pdf')
   elif (otype == 'png'):
      plt.savefig(oname+'.png')
   elif (otype == 'tiff'):
      plt.savefig(oname+'.tiff')
   else:
      plt.show()

   return 0

# main
def main(argv=None):
   if argv is None:
      argv = sys.argv

   # enforce proper number of arguments passed
   if len(argv) !=5:
      print "usage: %s PROFILE# SIMNAME AVGYN OUTTYPE" % os.path.basename(sys.argv[0])
      return 2

   # extract argument
   profno  = argv[1]
   simname = argv[2]
   avgyn   = argv[3]
   otype   = argv[4]
   hmax    = 43.0
   hc      = 33.0
#  if (int(profno) > 74):
#     simname = 'cow0728a'
#  else:
#     simname = 'cow0722a'

   # call routine
   plotprofs(simname, otype, profno, hmax, hc, avgyn)
 
   # all done!
   return 0

if __name__ == "__main__":
    sys.exit(main())
