import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib as mpl
import matplotlib.collections as mcoll
import matplotlib.path as mpath
import math
import networkx as nx
import scipy as sc
#from descartes import PolygonPatch
#import alphashape
import os
topdir='/Users/carmenlee/Documents/Github/20150731reprocessed/'
#topdir = '/eno/cllee3/DATA/esbilili/uniaxial/img1/'
adjacency_list_name = 'Adjacency_list.txt'
position_data_name = 'centers_tracked.txt'
stress_name = 'stress.txt'
reset = 'reset.txt'

pix2m = 0.00019853
########################################################user defined functions########################################
def sorting(contactangles, angle, forcemag,forcenor, forcetan):
    normals = []
    tangents = []
    mag = []
    for a in range(len(contactangles)-1):
        idx = np.where((angle >= contactangles[a]) & (angle < contactangles[a+1]))
        ind = np.asarray(idx)[0]
        fnorm = np.mean(forcenor[ind])
        ftan = np.mean(forcetan[ind])
       # fabs = np.mean(forcemag[ind])
        normals.append(fnorm)
        tangents.append(ftan)
        mag.append(len(ind))
        #print('tangents=', tangents)
    return(mag, normals, tangents)

def contactdis(thetas, a, phi):
    return((1+a*np.cos(2*(thetas-phi))))

def normaldis(thetas, a, phi):
    return(1+a*np.cos(2*(thetas-phi)))

def tangentdis(thetas, a, phi):
    return(-a*np.sin(2*(thetas-phi)))

def removeDanglers(forcedata):
    forcedatacopy = forcedata.reshape(forcedata.shape[0], -1)
    #print(forcedatacopy.shape)
    val, cnt = np.unique(forcedatacopy[:,1], return_counts=True)
    val2, cnt2 = np.unique(forcedatacopy[:,0], return_counts=True)
    singles = val[cnt==1]
    singles2 = val2[cnt2==1]
    #print(val[cnt == 1])
    danglers = len(val[cnt == 1])+len(val2[cnt2==1])
    
         
    while danglers > 0:
        indices = np.argwhere(np.isin(forcedatacopy[:,1], singles, invert = True))
        
        
        newind = []
        for m in range(len(indices)):
            #print(indices[m][0])
            newind.append(indices[m][0])
        #print(indices)
        forcedatacopy = forcedatacopy[newind, :]
        indices2 = np.argwhere(np.isin(forcedatacopy[:,0], singles2, invert = True))
        newind = []
        for m in range(len(indices2)):
            #print(indices[m][0])
            newind.append(indices2[m][0])
        #print(indices)
        forcedatacopy = forcedatacopy[newind, :]
        #print(forcedatacopy)
        val, cnt = np.unique(forcedatacopy[:,1], return_counts=True)
        val2, cnt2 = np.unique(forcedatacopy[:,0], return_counts=True)
        singles = val[cnt==1]
        singles2 = val2[cnt2==1]
        #print(val[cnt == 1])
        
        danglers = len(val[cnt == 1])+len(val2[cnt2==1])
        #print(danglers)
        
        
    return(forcedatacopy)

def find_centroid(xy):
    """ Computes the centroid (a.k.a. center of gravity) for a non-self-intersecting polygon.

    Parameters
    ----------
    polygon : list of two-dimensional points (points are array-like with two elements)
        Non-self-intersecting polygon (orientation does not matter).

    Returns
    -------
    center_of_gravity : list with 2 elements
        Coordinates (or vector) to the centroid of the polygon.
    """
    polygon = np.vstack((xy, xy[0]))
    offset = polygon[0]
    #print('o', offset)
    center_of_gravity = [0.0, 0.0]
    double_area = 0.0
    for ii in range(len(polygon)):
        p1 = polygon[ii]
        #print(p1)
        p2 = polygon[ii-1]
        f = (p1[0]-offset[0])*(p2[1]-offset[1]) - (p2[0]-offset[0])*(p1[1]-offset[1])
        double_area += f
        center_of_gravity[0] += (p1[0] + p2[0] - 2*offset[0]) * f
        center_of_gravity[1] += (p1[1] + p2[1] - 2*offset[1]) * f
    center_of_gravity[0] = center_of_gravity[0] / (3*double_area) + offset[0]
    center_of_gravity[1] = center_of_gravity[1] / (3*double_area) + offset[1]
    #return center_of_gravity
    # If you want to return both the CoG and the area, comment the return above
    return center_of_gravity, abs(double_area/2)
###########################################################################################################

# Extract some information from the data set:
try:
    adj_data = np.genfromtxt(topdir+adjacency_list_name, delimiter=',')
    #print(adj_data)
    # nsteps is the maximal frame number
    pos_data = np.genfromtxt(topdir+position_data_name,
                             delimiter=',', skip_header=1)
    #print(pos_data)
    #stress_data = np.genfromtxt(topdir+stress_name, delimiter=',')
    #print(stress_data.shape)
    reset_data = np.genfromtxt(topdir+reset, delimiter = ',')
    print(reset_data)
    reset_data = reset_data.astype(int)
    
    nsteps = np.unique(pos_data[:, 0]).astype(int)
    
    #print(nsteps)
except:
    print('No valid data found for experiment ')
    # if no data is found for an experiment number go to the next experiment
    
########################begin combining adjacency and position data##############################
#reset_data=np.reshape(reset_data, (1,2))
#print(reset_data)
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
maxstress = 0
bigboy =0
l = []
polygonlist = np.empty(len(nsteps))
#here we define the image names that correspond to images that are part of a compression cycle
valid_data = []
for chunk in range(len(reset_data)):
    points = np.arange(reset_data[chunk,0], reset_data[chunk, 1]+1)
    valid_data.append(points)
valid = np.asarray([item for sub_list in valid_data for item in sub_list])



#####################now we go over each image in the data set

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
    y = 2848 - newpositiondata[:,3]
    newpositiondata[:,3] = y
    
    #sanity check
    '''if nsteps[i] ==7:
        figure = plt.figure()
        axis = figure.add_subplot(111)
    
        axis.scatter(newpositiondata[:,2], newpositiondata[:,3], s=2*newpositiondata[:, 4]/pix2m, facecolors = 'none', edgecolors = 'k')
        axis.axis('scaled')
        #plt.show()'''
    #figurec = plt.figure()
    #axc = figurec.add_subplot(111)
    
    Sigma = np.asarray([[0,0],[0,0]]) #stress tensor set up 
    area = []
    #print(Sigma)
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

        neighbourpair1 = forcedata[:, :, 0]
        neighbourpair2 = forcedata[:, :, 1]
        #print(forcedata)
        
        
        #print(neighbourpair1)
        conversion1 = np.argwhere(newpositiondata[:, 1] == neighbourpair1)[:, 1]
        conversion2 = np.argwhere(newpositiondata[:, 1] == neighbourpair2)[:, 1]

        xy1 = newpositiondata[conversion1][:, 2:4]
        xy2 = newpositiondata[conversion2][:, 2:4]
        r = newpositiondata[conversion1][:,4]/pix2m
        lx = xy2[:, 0]-xy1[:, 0]
        ly = xy2[:, 1]-xy1[:, 1]
        '''if nsteps[i] == 7:
            axis.plot([xy1[:,0], xy2[:,0]], [xy1[:,1], xy2[:,1]], 'r-')
            plt.show()'''
        finalforcedatacopy = removeDanglers(forcedata)
        if len(finalforcedatacopy)>0:
            tempneighbourpair1 = finalforcedatacopy[:,0].reshape((-1, 1))
            tempneighbourpair2 = finalforcedatacopy[:,1].reshape((-1, 1))
            
            
            # G = nx.Graph()
            # for k in range(len(tempneighbourpair1)):
            #     G.add_edge(tempneighbourpair1[k][0], tempneighbourpair2[k][0])
           
            # #print(sorted(G))
            # cycles = nx.simple_cycles(G)
            # #print(cycles)
            # loc = sorted(cycles)
            # #print(loc)
            # H = nx.Graph()
            # for l in range(len(loc)):
            #     nx.add_cycle(H, loc[l])
            # polygons = sorted(nx.minimum_cycle_basis(H))
            # print('polygons', polygons)
            # ordered = []
            # for poly in range(len(polygons)):
                
            #     for bigpoly in range(len(loc)):
                    
            #         trackers = True
            #         for k in loc[bigpoly]:
            #             if k not in polygons[poly]:
            #                 trackers = False
            #                 break 
            #         if trackers:
            #             ordered.append(loc[bigpoly])
                        
            # polygonlist[i] = [nsteps[i], ordered]      
            # print('ordered',ordered)
            # for item in range(len(ordered)):
            #     temp = np.asarray(ordered[item]).reshape((-1, 1))
            #     locations = []
            #     xy = []
            #     for k in ordered[item]:
            #         ind = np.argwhere(newpositiondata[:,1] == k)[0][0]
            #         #print(ind)
            #         locations.append(ind)
            #         xy.append(newpositiondata[ind,2:4])
                #locations = np.argwhere(newpositiondata[:,1] == temp)[:,1]
                
                # xy = np.asarray(xy)
                # #print(locations, xy)
                # #xy = newpositiondata[locations][:,2:4]
                # #axc.scatter(newpositiondata[locations,2], newpositiondata[locations,3], s=2*newpositiondata[locations, 4]/pix2m, facecolors = 'none', edgecolors = 'k')
                # #axc.plot(xy[:,0], xy[:,1])
                # #for m in range(len(locations)):
                #    # axc.text(newpositiondata[locations[m],2], newpositiondata[locations[m],3], str(locations[m]))
                # #axis.plot(xy[:,0], xy[:,1], 'yo')
                # centroid, area = find_centroid(xy)
                #print(centroid, area)
                #axis.plot(centroid[0], centroid[1], 'go')
            
            
            conversion1f = np.argwhere(newpositiondata[:, 1] == tempneighbourpair1)[:,1]
            conversion2f= np.argwhere(newpositiondata[:, 1] == tempneighbourpair2)[:,1]

            xy1f = newpositiondata[conversion1f][:, 2:4]
            xy2f = newpositiondata[conversion2f][:, 2:4]
            lxf = xy2[:, 0]-xy1[:, 0]
            lyf= xy2[:, 1]-xy1[:, 1]
        
            #for m in range(len(tempneighbourpair1)):
                #axis.text(xy1f[m,0], xy1f[m,1], str(tempneighbourpair1[m]))
        
            #axis.plot([xy1f[:,0], xy2f[:,0]], [xy1f[:,1], xy2f[:,1]], 'b-')
        
        angle = np.arctan2(xy2[:, 1]-xy1[:, 1], xy2[:, 0]-xy1[:, 0])
        
        ly = r*np.sin(angle)
        lx = r*np.cos(angle)
        ''' if nsteps[i] ==7:
            print(angle, ly, lx)'''
        anglemaster.append(angle)

        forcedata = forcedata.reshape(forcedata.shape[0], -1)
        # print(len(forcedata[:,3]))
        # print(len(angle))
        forceabs = np.sqrt(forcedata[:, 3]**2+forcedata[:, 2]**2)
        # print(forceabs)
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
            print(forcex)
            sigmayy = sum(ly*forcey)
            sigmaxy = sum(ly*forcex)
            sigmayx = sum(lx*forcey)
            #print(lx, forcex)
            Sigma[0,0] = sigmaxx
            Sigma[1, 0] =sigmayx
            Sigma[0,1] = sigmaxy
            Sigma[1,1]= sigmayy
            
            # #comment here until the end of the loop if you want to plot a rose plot for each image
           
            # radii, theta = np.histogram(angle, bins=36, range=[-np.pi, np.pi])
            # data = np.vstack([angle, forceabs, forcedata[:,3],forcedata[:,2]])
            # datanew = data[:,np.abs(data[3,:]-( 0.0006497852319827032))<=(3 * 0.16848932776403405)]
            # mag, normals, tangents = sorting(theta, datanew[0,:], datanew[1,:], datanew[2,:], datanew[3,:])
            # normals = np.asarray([0 if math.isnan(x) else x for x in normals])
            # tangents = np.asarray([0 if math.isnan(x) else x for x in tangents])
            # mag = np.asarray([0 if math.isnan(x) else x for x in mag])
            # # #mag = mag.replace(np.nan,0)
            # # print(normals, tangents, mag)
            # if len(normals)!=0 and len(tangents)!=0 and len(radii)!=0:
            #     #print('pink')
            #     fig2 = plt.figure(2)
            #     ax2 = fig2.add_subplot(111, projection='polar')
            #     bars = ax2.bar(theta[0:-1], radii, bottom=0.0)
            #     ax2.set_rmax(120)
            #     fig2.savefig(topdir+'filteredcontactrose'+str(nsteps[i])+'.png')
                
            #     fig4 = plt.figure(3)
            #     ax4 = fig4.add_subplot(111, projection='polar')
                
            #     bars = ax4.bar(theta[0:-1], normals, bottom=0.0)
            #     ax4.set_rmax(0.3)
            #     fig4.savefig(topdir+'filterednormalrose'+str(nsteps[i])+'.png')
            #     colors = plt.cm.viridis(20*tangents/max(tangents))
            #     fig5 = plt.figure(4)
            #     ax5 = fig5.add_subplot(111, projection='polar')
                
            #     bars = ax5.bar(theta[0:-1], np.abs(tangents), bottom=0.0, color=colors)
            #     ax5.set_rmax(0.01)
            #     fig5.savefig(topdir+'filteredtangentialrose'+str(nsteps[i])+'.png')
            #     plt.close()
                 
        #print(np.concatenate((forcedata, forcex.reshape(len(forcex), 1), forcey.reshape(len(forcex), 1)), axis = 1))
        deltay = max(newpositiondata[:, 3]+newpositiondata[:, 4]/pix2m) - \
            min(newpositiondata[:, 3]-newpositiondata[:, 4]/pix2m)
        deltax = max(newpositiondata[:, 2]+newpositiondata[:, 4]/pix2m) - \
            min(newpositiondata[:, 2]-newpositiondata[:, 4]/pix2m)
        area.append(deltay*deltax)
        strainx.append((deltax-deltax0)/deltax0)
        strainy.append((deltay-deltay0)/deltay0)
        #stress.append(stress_data[i, 2])
        forces.append(sum(forceabs))
        stresstensor.append(Sigma/(deltax*deltay))
        # if stress_data[i,2]>maxstress:
        #     bigboy = nsteps[i]
        #     maxstress = stress_data[i,2]
            
        #     maxangle = angle
        #     maxtangent = forcedata[:,2]
        #     maxnormal = forcedata[:,3]
    
        #colors = plt.cm.viridis(tangents/max(tangents))
   # axc.axis('scaled')
  #  plt.show()
#np.savetxt(topdir+'polygons', polygonlist)
strainx = np.asarray(strainx)
strainy = np.asarray(strainy)


Sigmas  =np.asarray(stresstensor)
#print(Sigmas.shape)



a, b = os.path.split(topdir[:-1])
a, c = os.path.split(a)
angle = np.asarray([item for sub_list in anglemaster for item in sub_list])
radii, theta = np.histogram(angle, bins=36, range=[-np.pi, np.pi])
fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection='polar')
bars = ax2.bar(theta[0:-1], radii, bottom=0.0)
fig2.savefig(topdir+'contactrose'+c+b+'.png')
forceangle = [item for sub_list in forceanglemaster for item in
              sub_list]
fangle = np.asarray(forceangle)
forcemag = np.asarray([item for sub_list in weights for item in
                       sub_list])
forcenor = np.asarray([item for sub_list in fnormal for item in
                       sub_list])

forcetan = np.asarray([item for sub_list in ftangential for item in
                       sub_list])
strain = []
sigma = []

l0 = []
for k in range(len(strainy)):
    strain.append([np.sqrt(strainy[k]**2+strainx[k]**2)]*len(ftangential[k]))
    sigma.append([(Sigmas[k,0,0]-Sigmas[k,1,1])/(Sigmas[k,1,1]+Sigmas[k,0,0])]*len(ftangential[k]))
strains = np.asarray([item for sub_list in strain for item in
                        sub_list])   
qp = np.asarray([item for sub_list in sigma for item in
                        sub_list])   
data = np.vstack((angle, forcemag, forcenor, forcetan, strains))
dataint = data[:,np.abs(data[2,:]-np.mean(data[2,:]))<=(4* data[2,:].std())]
datanew = dataint[:,np.abs(dataint[3,:]-np.mean(data[3,:]))<=(4 * data[3,:].std())]
np.savetxt(topdir+'data.txt', np.transpose(datanew))
fig = plt.figure()
axqp = fig.add_subplot(111)
#aax = axqp.twinx()
c = np.arange(len(strainx))

axqp.scatter([np.sqrt(strainy**2+strainx**2)], (Sigmas[:,0,0]-Sigmas[:,1,1])/(Sigmas[:,1,1]+Sigmas[:,0,0]),c = c)
axqp.set(xlabel = r'|$\sqrt{\left(\frac{y-y_0}{y_0}\right)^2 + \left(\frac{x-x_0}{x_0}\right)^2}$', ylabel = r'$(\sigma_{yy}-\sigma_{xx})/(\sigma_{yy}+\sigma_{xx})$')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(range(len(Sigmas[:,0,0])),Sigmas[:,0,0],label = r'$\sigma_{xx}$')
ax.plot(range(len(Sigmas[:,0,0])),Sigmas[:,1,0],label = r'$\sigma_{yx}$')
ax.plot(range(len(Sigmas[:,0,0])),Sigmas[:,0,1],label = r'$\sigma_{xy}$')
ax.plot(range(len(Sigmas[:,0,0])),Sigmas[:,1,1],label = r'$\sigma_{yy}$')
ax.set(xlabel = 'frame number', ylabel = 'stress tensor')
fig.legend()


print('mean=',np.mean(data[3,:]), 'std=',np.std(data[3,:]))
fig, ax = plt.subplots(ncols = 2)
ax[0].plot(np.sqrt(strainy**2+strainx**2), Sigmas[:, 0,0], label = r'$\sigma_{xx}$')
ax[1].plot(np.sqrt(strainy**2+strainx**2), Sigmas[:, 1,1], label = r'$\sigma_{yy}$')
#ax[0].plot(forcenor,'.')
#ax[0].plot(datanew[2,:], '.')
#ax[0].set(xlabel= 'force id', ylabel ='normal force')
#ax[1].plot(forcetan, '.')
#ax[1].plot(datanew[3,:], '.')
#ax[1].set(xlabel= 'force id', ylabel ='tangential force')



radii, theta = np.histogram(datanew[0,:], bins=36, range=[-np.pi, np.pi])
strainlimits = np.linspace(0, max(strains), 11)

#strainlimits = np.linspace(0, 0.010, 6)
print(strainlimits)

figlintan, axlintan = plt.subplots()
figlinnor, axlinnor = plt.subplots()
figlincon, axlincon = plt.subplots()
axlintan.set(ylabel = 'f_t/f0', xlabel = 'angle')
axlinnor.set(ylabel = 'f_n/f0', xlabel = 'angle')
axlincon.set(ylabel = '# contacts/sum contacts', xlabel = 'angle')
#figtan, axtan = plt.subplots(ncols = 5, subplot_kw={'projection': 'polar'})
an = []
at = []
a =[]
strainmid = []
forcescale = []
colors = plt.cm.viridis(np.linspace(0,1,len(strainlimits)))
print(colors)
for m in range(len(strainlimits)-1):
    datanew2 = datanew[:,((datanew[4,:]>strainlimits[m]) & (datanew[4,:]<=strainlimits[m+1]))]
    mag, normals, tangents= sorting(theta, datanew2[0,:], datanew2[1,:], datanew2[2,:], datanew2[3,:])
    # fig3 = plt.figure()
    # ax3 = fig3.add_subplot(111, projection='polar')
    # bars = ax3.bar(theta[0:-1], datanew2[0,:], bottom=0.0)
    # fig3.savefig(topdir+'strainavgcontactrose'+c+b+'strain'+str(m)+'.png')

    #radii, theta = np.histogram(angle, bins = 360, weights = weights, range = [-np.pi, np.pi])
    f0 = np.mean(normals)
    e0 = np.mean(mag)
    
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
    forcescale.append(f0)
    strainmid.append((strainlimits[m]+strainlimits[m+1])/2)
    print(p0,poptn)
    #norm = mpl.colors.Normalize(vmin=-np.pi, vmax=np.pi)
    #colors = plt.cm.viridis(tangents/max(tangents))
    axlincon.plot(theta[0:-1], mag/e0, '.-', color = colors[m])
    axlincon.plot(theta[0:-1], contactdis(theta[0:-1], poptc[0], poptc[1]), color = colors[m])
    axlinnor.plot(theta[0:-1], normals/f0,'.-', color = colors[m] )
    #axlinnor.plot(theta[0:-1], np.ones(len(normals)))
    axlinnor.plot(theta[0:-1], normaldis(theta[0:-1], poptn[0], poptn[1]), color=colors[m])
    #axlinnor.plot(theta[0:-1], normaldis(theta[0:-1], p0))
    #axtan[m].bar(theta[0:-1], normals, bottom=0.0)
    #axtan[m].set_rmax(0.3)
    

    fig5 = plt.figure()
    ax5 = fig5.add_subplot(111, projection='polar')
    bars = ax5.bar(theta[0:-1], np.abs(tangents), bottom=0.0, color=colors)
    axlintan.plot(theta[0:-1], (tangents)/f0, '.-',color = colors[m])
    axlintan.plot(theta[0:-1], tangentdis(theta[0:-1], poptt[0], poptt[1]), color = colors[m])
    ax5.set_rmax(0.004)
    #fig5.savefig(topdir+'strainavgtangentialrose'+c+b+'strain'+str(m)+'.png')
    plt.close(fig5)
a = np.asarray(a)
at= np.asarray(at)
an = np.asarray(an)
ax[0].plot(strainmid, f0*(1+a*an/2-(a+an+at)/2))
ax[1].plot(strainmid, f0*(1+a*an/2+(a+an+at)/2))
    
fig6, ax6 = plt.subplots(ncols = 3)
ax6[0].plot(strainmid, a, '.')
ax6[1].plot(strainmid, an, '.')
ax6[2].plot(strainmid, at, '.')
ax6[0].set(xlabel= 'strain', ylabel ='a')
ax6[1].set(xlabel= 'strain', ylabel ='an')
ax6[2].set(xlabel= 'strain', ylabel ='at')
axqp.plot(strainmid, (np.asarray(a)+np.asarray(an)+np.asarray(at))/2, label = r'$a+a_n+a_t$')
fig.legend()
#figtan.savefig(topdir+'strainavgnormalrose'+c+b+'all'+'.png')
plt.show()
mag, normals, tangents= sorting(theta, datanew[0,:], datanew[1,:], datanew[2,:], datanew[3,:])
fig3 = plt.figure()
ax3 = fig3.add_subplot(111, projection='polar')
bars = ax3.bar(theta[0:-1], mag, bottom=0.0)
#fig3.savefig(topdir+'forcerose'+c+b+'.png')

#radii, theta = np.histogram(angle, bins = 360, weights = weights, range = [-np.pi, np.pi])
norm = mpl.colors.Normalize(vmin=-np.pi, vmax=np.pi)
colors = plt.cm.viridis(tangents/max(tangents))
fig4 = plt.figure()
ax4 = fig4.add_subplot(111, projection='polar')
bars = ax4.bar(theta[0:-1], normals, bottom=0.0)
#fig4.savefig(topdir+'normalrose'+c+b+'.png')

fig5 = plt.figure()
ax5 = fig5.add_subplot(111, projection='polar')
bars = ax5.bar(theta[0:-1], np.abs(tangents), bottom=0.0, color=colors)
#fig5.savefig(topdir+'tangentialrose'+c+b+'.png')
plt.show()
