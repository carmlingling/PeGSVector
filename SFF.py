#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 11:08:23 2023

@author: cllee3
"""

#Import required modules


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib as mpl
import matplotlib.collections as mcoll
import matplotlib.path as mpath
import math

import scipy as sc
#from descartes import PolygonPatch
#import alphashape
import os


###########File names################
#topdir = '/eno/cllee3/DATA/esbilili/shear/img1/'
topdir = '/Users/carmenlee/Documents/GitHub/20150731reprocessed/'
adjacency_list_name = 'Adjacency_list.txt'
position_data_name = 'centers_tracked.txt'
reset = 'reset.txt'
pix2m = 0.00019853

#######################################user defined functions########################################
def sorting(contactangles, angle, forcemag,forcenor, forcetan, l0, mv):
    normals = []
    tangents = []
    mag = []
    contactdens = []
    contactlength = []
    for a in range(len(contactangles)-1):
        idx = np.where((angle >= contactangles[a]) & (angle < contactangles[a+1]))
        ind = np.asarray(idx)[0]
        fnorm = np.mean(forcenor[ind])
        ftan = np.mean(forcetan[ind])
       # fabs = np.mean(forcemag[ind])
      
        normals.append(fnorm)
        tangents.append(ftan)
        mag.append(len(ind))
        contactdens.append(np.mean(mv[ind]))
        contactlength.append(np.mean(l0[ind]))
        #print('tangents=', tangents)
    return(mag, normals, tangents, contactdens, contactlength)

def contactdis(thetas, a, phi):
    return((1+a*np.cos(2*(thetas-phi))))

def normaldis(thetas, a, phi):
    return(1+a*np.cos(2*(thetas-phi)))

def tangentdis(thetas, a, phi):
    return(-a*np.sin(2*(thetas-phi)))


################################Import data###############################

# Extract some information from the data set:
try:
    adj_data = np.genfromtxt(topdir+adjacency_list_name, delimiter=',')
    #print(adj_data)
    # nsteps is the maximal frame number
    pos_data = np.genfromtxt(topdir+position_data_name,
                             delimiter=',', skip_header=1)
    #print(pos_data)
    
    reset_data = np.genfromtxt(topdir+reset, delimiter = ',')
    
    reset_data = reset_data.astype(int)
    reset_data=np.reshape(reset_data, (-1, 2))
    print(reset_data)
    nsteps = np.unique(pos_data[:, 0]).astype(int)
    
    #print(nsteps)
except:
    print('No valid data found for experiment ')
    # if no data is found for an experiment 

###############################From reset data note which frame number is part of a full cycle#####################    #here we define the image names that correspond to images that are part of a compression cycle
valid_data = []

for chunk in range(len(reset_data)):
    points = np.arange(reset_data[chunk,0], reset_data[chunk, 1]+1)
    valid_data.append(points)
valid = np.asarray([item for sub_list in valid_data for item in sub_list])
#valid = np.arange(reset_data[0], reset_data[1]+1)
#############################Collect data by frame number###################

anglemaster = []
forceanglemaster = []
weights = []
strainx = []
strainy = []
fnormal = []
stress = []
ftangential = []
stresstensor = []
forces = []


l0 = []

mv = []

for i in range(len(nsteps)):
    print(i)
   
    #importing force data
    ind = np.argwhere(adj_data[:, 0].astype(int) == nsteps[i])
    forcedata = adj_data[ind, 1:]
    
    #importing position data
    indpos = np.argwhere(pos_data[:, 0].astype(int) == nsteps[i])
    positiondata = pos_data[indpos, :]
    #print(positiondata)
    
    newpositiondata = positiondata.reshape(positiondata.shape[0], -1)
    y = 2848 - newpositiondata[:,3] #convert to RH coordinate
    newpositiondata[:,3] = y
    
    Sigma = np.asarray([[0,0],[0,0]]) #stress tensor set up 
    
    
    #find the start of the compression and note the initial chamber size
    if nsteps[i] in reset_data[:,0]:
        deltay0 = max(newpositiondata[:, 3]+newpositiondata[:, 4]/pix2m) - \
            min(newpositiondata[:, 3]-newpositiondata[:, 4]/pix2m)
        deltax0 = max(newpositiondata[:, 2]+newpositiondata[:, 4]/pix2m) - \
            min(newpositiondata[:, 2]-newpositiondata[:, 4]/pix2m)
        print('deltay0', deltay0)
    
    if len(forcedata) == 0: #ignore empty force data
        pass
    elif nsteps[i] in valid:

        neighbourpair1 = forcedata[:, :, 0] #retrieve particle ID
        neighbourpair2 = forcedata[:, :, 1]
        #print(forcedata)
        N = len(neighbourpair1)
        #Z= 
        #print(neighbourpair1)
        conversion1 = np.argwhere(newpositiondata[:, 1] == neighbourpair1)[:, 1]
        conversion2 = np.argwhere(newpositiondata[:, 1] == neighbourpair2)[:, 1]

        xy1 = newpositiondata[conversion1][:, 2:4] #positions
        xy2 = newpositiondata[conversion2][:, 2:4]
        r = newpositiondata[conversion1][:,4]/pix2m #radii
        l0.append(np.average(r))#(*pix2m))
        
        
        angle = np.arctan2(xy2[:, 1]-xy1[:, 1], xy2[:, 0]-xy1[:, 0])
        
        ly = r*np.sin(angle)
        lx = r*np.cos(angle)
        ''' if nsteps[i] ==7:
            print(angle, ly, lx)'''
        anglemaster.append(angle)

        forcedata = forcedata.reshape(forcedata.shape[0], -1)
        
        forceabs = np.sqrt(forcedata[:, 3]**2+forcedata[:, 2]**2)
        
        if len(forcedata[:, 3]) != len(angle):
            pass
        else:
            forcex = forcedata[:, 3] * np.cos(angle)-forcedata[:, 2]*np.sin(angle)
            forcey = forcedata[:, 3] * np.sin(angle)-forcedata[:, 2]*np.cos(angle)
            forceangle = np.arctan2(forcey, forcex)
            forceanglemaster.append(forceangle)
            fnormal.append(forcedata[:, 3])
            ftangential.append(forcedata[:, 2])
            weights.append(forceabs)
            
            sigmaxx = sum(lx*forcex)
            
            sigmayy = sum(ly*forcey)
            sigmaxy = sum(ly*forcex)
            sigmayx = sum(lx*forcey)
            #print(lx, forcex)
            Sigma[0,0] = sigmaxx
            Sigma[1, 0] =sigmayx
            Sigma[0,1] = sigmaxy
            Sigma[1,1]= sigmayy
            
        deltay = max(newpositiondata[:, 3]+newpositiondata[:, 4]/pix2m) - \
            min(newpositiondata[:, 3]-newpositiondata[:, 4]/pix2m)
        deltax = max(newpositiondata[:, 2]+newpositiondata[:, 4]/pix2m) - \
            min(newpositiondata[:, 2]-newpositiondata[:, 4]/pix2m)
        mv.append(N/(deltay*deltax))
        strainx.append((deltax-deltax0)/deltax0)
        strainy.append((deltay-deltay0)/deltay0)
        #stress.append(stress_data[i, 2])
        forces.append(sum(forceabs))
        stresstensor.append(Sigma/(deltax*deltay))
        
strainx = np.asarray(strainx)
strainy = np.asarray(strainy)
Sigmas  =np.asarray(stresstensor)

fig = plt.figure()
axqp = fig.add_subplot(111)
#aax = axqp.twinx()
c = np.arange(len(strainx))

axqp.scatter([np.sqrt(strainy**2+strainx**2)], (Sigmas[:,0,0]-Sigmas[:,1,1])/(Sigmas[:,1,1]+Sigmas[:,0,0]),c = c)
axqp.set(xlabel = r'|$\sqrt{\left(\frac{y-y_0}{y_0}\right)^2 + \left(\frac{x-x_0}{x_0}\right)^2}$', ylabel = r'$(\sigma_{yy}-\sigma_{xx})/(\sigma_{yy}+\sigma_{xx})$')
plt.show()


angle = np.asarray([item for sub_list in anglemaster for item in sub_list])
radii, theta = np.histogram(angle, bins=36, range=[-np.pi, np.pi])

forceangle = [item for sub_list in forceanglemaster for item in
              sub_list]
fangle = np.asarray(forceangle)
forcemag = np.asarray([item for sub_list in weights for item in
                       sub_list])
forcenor = np.asarray([item for sub_list in fnormal for item in
                       sub_list])

forcetan = np.asarray([item for sub_list in ftangential for item in
                       sub_list])
print(forcetan)
strain = []
l = []
m = []
for k in range(len(strainy)):
    strain.append([np.sqrt(strainy[k]**2+strainx[k]**2)]*len(ftangential[k]))
    l.append([l0[k]]*len(ftangential[k]))
    m.append([mv[k]]*len(ftangential[k]))
    
strains = np.asarray([item for sub_list in strain for item in
                            sub_list]) 
l0s = np.asarray([item for sub_list in l for item in
                            sub_list]) 
mvs = np.asarray([item for sub_list in m for item in
                            sub_list]) 

data = np.vstack((angle, forcemag, forcenor, forcetan, strains, l0s, mvs))
dataint = data[:,np.abs(data[2,:]-np.mean(data[2,:]))<=(4* data[2,:].std())]
datanew = dataint[:,np.abs(dataint[3,:]-np.mean(data[3,:]))<=(4 * data[3,:].std())]
np.savetxt(topdir+'data.txt', np.transpose(datanew))


fig, ax = plt.subplots(ncols = 2)
ax[0].plot(np.sqrt(strainy**2+strainx**2), Sigmas[:, 0,0], label = r'$\sigma_{xx}$')
ax[1].plot(np.sqrt(strainy**2+strainx**2), Sigmas[:, 1,1], label = r'$\sigma_{yy}$')

radii, theta = np.histogram(datanew[0,:], bins=36, range=[-np.pi, np.pi])
strainlimits = np.linspace(0, max(strains), 11)

an = []
at = []
a =[]
strainmid = []
forcescale = []
colors = plt.cm.viridis(np.linspace(0,1,len(strainlimits)))

for m in range(len(strainlimits)-1):
    datanew2 = datanew[:,((datanew[4,:]>strainlimits[m]) & (datanew[4,:]<=strainlimits[m+1]))]
    mag, normals, tangents, contactdens, contactlength= sorting(theta, datanew2[0,:], datanew2[1,:], datanew2[2,:], datanew2[3,:], datanew2[5,:],datanew2[6,:])
    
    contactdens = np.mean(datanew2[5,:])
    contactlength = np.mean(datanew2[6,:])
    f0 = np.mean(normals)
    e0 = np.mean(mag)



    fig5 = plt.figure()
    ax5 = fig5.add_subplot(111, projection='polar')
    print(theta, tangents)
    bars = ax5.bar(theta[0:-1], np.abs(tangents), bottom=0.0, color=colors)
    
    ax5.set_rmax(0.004)
    plt.show()
    
    p0 = [max(mag)-e0/e0, theta[np.argmax(mag)]]
    poptc, pcov = sc.optimize.curve_fit(contactdis, theta[0:-1], np.asarray(mag)/e0, p0 =p0, bounds=([0, -np.pi], [np.inf, np.pi]))
    
    p0 = [(max(normals)-f0)/f0, theta[np.argmax(normals)]]
    poptn, pcov = sc.optimize.curve_fit(normaldis, theta[0:-1], normals/f0, p0 =p0,bounds=([0, -np.pi], [np.inf, np.pi]))
    
    
    p0 = [max(tangents)/f0, theta[np.argmax(tangents)]]
    print(p0)
    poptt, pcov = sc.optimize.curve_fit(tangentdis, theta[0:-1], tangents/f0, p0 =p0,bounds=([0, -np.pi], [np.inf, np.pi]))
    a.append(poptc[0])
    an.append(poptn[0])
    at.append(poptt[0])
    forcescale.append(f0*contactdens*contactlength/2)
    strainmid.append((strainlimits[m]+strainlimits[m+1])/2)
    
    
a = np.asarray(a)
at= np.asarray(at)
an = np.asarray(an)
ax[0].plot(strainmid, forcescale*(1+a*an/2+(a+an+at)/2))
ax[1].plot(strainmid, forcescale*(1+a*an/2-(a+an+at)/2))
plt.show()
