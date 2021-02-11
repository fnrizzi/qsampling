#!/usr/bin/env python

import numpy as np
import sys, re, os
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from numpy import linalg as la
from scipy import linalg

# radius of the earth
earthRadius = 6371. #km
mantleThickness = 2891. # km
cmbRadius = earthRadius - mantleThickness

def plotCMB(ax):
  cmbTh = np.linspace(0, 2*np.pi, 100)
  cmbRa = cmbRadius*np.ones(100)
  ax.plot(cmbTh, cmbRa, c='k', linestyle='-', linewidth=0.5, zorder=2)

def plotEarthSurf(ax):
  surfTh = np.linspace(0, 2*np.pi, 100)
  surfRa = earthRadius*np.ones(100)
  ax.plot(surfTh, surfRa, c='k', linewidth=0.5)

def doPlot(ax, th, r, z, bd):
  cm1 = plt.cm.get_cmap('cividis')

  h1=ax.pcolormesh(th, r, z, cmap=cm1, shading = "flat",
                    vmin=bd[0], vmax=bd[1], zorder=1)

  ax.set_ylim([cmbRadius, earthRadius])
  ax.set_yticks([]) #[3480, 5701, 6371])
  #plt.yticks(fontsize=13)
  ax.set_thetamin(-90)
  ax.set_thetamax(90)
  ax.set_xticks(np.pi/180. * np.linspace(-90, 90., 7, endpoint=True))
  ax.set_xticklabels([r'$\pi$', r'$5\pi/6$', r'$4\pi/6$',
                       r'$\pi/2$', r'$2\pi/6$', r'$\pi/6$', r'$0$'],
                      fontsize=11)

  ax.set_rorigin(-1)
  plotEarthSurf(ax)
  plotCMB(ax)

  mycolor = 'w'
  ax.xaxis.label.set_color(mycolor);
  ax.tick_params(axis='x', colors=mycolor)
  ax.yaxis.label.set_color(mycolor);
  ax.tick_params(axis='y', colors=mycolor)
  #fig1.savefig(outName, format="png",bbox_inches='tight', dpi=300, transparent=True)
  #plt.show()

if __name__== "__main__":
  nr, nth = 200, 1000
  cc  = np.loadtxt("./coords_vp.txt")
  th0, r0 = -cc[:,0]+np.pi/2., cc[:, 1]/1000. #m to km
  th, r = th0.reshape((nr,nth)), r0.reshape((nr,nth))

  d = np.loadtxt("./snaps_vp_0", skiprows=1)
  # print(np.min(d0), np.max(d0))
  # d = 2.*(d0 - np.min(d0))/np.ptp(d0)-1.
  # print(np.min(d), np.max(d))
  # print(d.shape)

  U, lam, VH = np.linalg.svd(d, full_matrices=False)

  s = 100

  fig1 = plt.figure(0)
  ax = fig1.add_subplot(111, projection='polar')
  doPlot(ax, th, r, d[:,-1].reshape((nr, nth)), [-5e-9, 5e-9])

  q4, r4, p4 = linalg.qr(U.T, pivoting=True)
  ax.scatter(th0[p4[0:s]], r0[p4[0:s]], c='b', s=10, zorder=2)

  plt.tight_layout()
  plt.show()



  #doPlot1(th, r, d[:,-1].reshape((nr, nth)), 1,[-5e-9, 5e-9])
