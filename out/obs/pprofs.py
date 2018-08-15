#!/usr/bin/env python
#
# pprofs.py - plot NH3 profile results from Coweeta ACCESS_NH3 simulations along with
#             observed concentration profiles
#
# Rick D. Saylor, November 2017
#
import os
import sys
import matplotlib.pylab as plt
import numpy as np
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
   fname = dirout+simname+'/sp/sp.dat'
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
      plt.errorbar(nh3avg, z, color='0.15', ecolor='0.70', xerr=xerr, linestyle='-', linewidth=lnwdth, label=labstr)
   else:
      ic=0
      for i in phrs:
         labstr = dts[i]
         plt.plot(nh3[:, i], z, color=colors[ic], linestyle='-', linewidth=lnwdth, label=labstr)
         ic+=1
         if (ic > len(colors)-1):
           ic=0

   # plot observations
   # obs data
   profobs = {
   '69':  [0.355362592, 0.302686033, 0.239429972, 0.182263993, 0.219309514, 0.295367658, 0.314590512, 0.266674335, 0.305643530, 0.310166473],
   '70':  [0.383545679, 0.211495924, 0.232042935, 0.171461270, 0.336749725, 0.475650441, 0.286076686, 0.302310689, 0.361011620, 0.377921136],
   '71':  [0.388907842, 0.206305755, 0.158095811, 0.292218679, 0.368670186, 0.453537286, 0.340720874, 0.382774370, 0.294886894, 0.419157592],
   '72':  [0.425279607, 0.216278872, 0.213988599, 0.197393076, 0.208246146, 0.225368456, 0.262928398, 0.374167279, 0.268627733, 0.268828453],
   '73':  [0.418028325, 0.188282357, 0.173292096, 0.160811452, 0.192939524, 0.187102373, 0.199750627, 0.284963607, 0.184367368, 0.209509278],
   '74':  [0.570353358, 0.229000126, 0.226928396, 0.170792180, 0.205778623, 0.209353014, 0.307837410, 0.311107130, 0.307648029, 0.246976095],
   '76':  [0.379063715, 0.212519415, 0.180330861, 0.121151765, 0.131353927, 0.157234207, 0.166498651, 0.169351934, 0.169305337, 0.259780897, 0.185796233],
   '77':  [0.353647819, 0.301414286, 0.102101025, 0.071822805, 0.096016834, 0.101985090, 0.096160549, 0.107225510, 0.144769632, 0.114614153, 0.124090133],
   '78':  [0.327849217, 0.184746984, 0.094580087, 0.103543718, 0.144313511, 0.160624812, 0.124451435, 0.174464558, 0.145672245, 0.194103598, 0.192180395],
   '79':  [0.437244662, 0.166550936, 0.123435180, 0.132661561, 0.149211208, 0.160423296, 0.186431359, 0.210612531, 0.257595165, 0.218268974, 0.194438350],
   '80':  [0.162693591, 0.157785512, 0.163158502, 0.091850650, 0.248539238, 0.326772894, 0.287265835, 0.290968669, 0.331810930, 0.359841864, 0.362765229],
   '81':  [0.237861735, 0.113895179, 0.116711305, 0.100937322, 0.141518093, 0.150007216, 0.134403519, 0.138395179, 0.147266864, 0.203266624, 0.185594405] 
   }

   zobs1 = [0.53, 1.78, 4.42, 9.91, 17.23, 24.85, 28.20, 32.01, 34.60, 43.00]
   zobs2 = [0.53, 1.78, 4.42, 9.91, 17.23, 24.85, 28.20, 32.01, 34.60, 39.33, 43.00]

   nh3obs = profobs[profno]
   if (int(profno) > 74):
      zobs = zobs2
   else:
      zobs = zobs1

   labstr = 'Profile #' + str(profno)
   plt.plot(nh3obs, zobs, linestyle='None', marker='o', markersize=8.0, mfc=colors[4], label=labstr)

   plt.minorticks_on()
   plt.ylim(-0.1, hmax)

   # draw line showing canopy height
   xbnds = list(ax.get_xlim())
   xbnds[0] = 0.0
   xbnds[1] = round(xbnds[1], 2)
   print('xbnds=', xbnds)
   hc = [33.0, 33.0]
   plt.plot(xbnds, hc, color='0.25', linestyle='--', linewidth=lnwdth)
   ax.set_xlim(xmin=xbnds[0], xmax=xbnds[1])

   # make it tidy
   plt.grid(b=True, which='major', color='gray', linewidth=1.0, alpha=0.5)
   plt.grid(b=True, which='minor', color='gray', linewidth=0.5, alpha=0.25)

   ax.tick_params(which="both", direction="out")
   ax.tick_params(which="major", length=tlmaj)
   ax.tick_params(which="minor", length=tlmin)
   ax.tick_params(which="both", labelsize=tlbsize, pad=tlbpad)

   plt.legend(loc=1, fontsize=lfsize, bbox_to_anchor=(0.59, 0.94))
#  plt.legend(loc=1, fontsize=lfsize, bbox_to_anchor=(0.99, 0.94))

   plt.xlabel('NH$_3$ ($\mu$g m$^{-3}$)', fontsize=xfsize, labelpad=10)
   plt.ylabel('z (m)', fontsize=yfsize, labelpad=10)
   plt.title('NH$_3$ Canopy Profile - '+simname, fontsize=tfsize, y=tyloc)

   oname = './img/'+simname+'-'+str(profno)+'_pprof'
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
      print("usage: %s PROFILE# SIMNAME AVGYN OUTTYPE" % os.path.basename(sys.argv[0]))
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
