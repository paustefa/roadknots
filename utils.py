#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 10:27:17 2018

@author: stefan
"""

import os
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np

#import customized function
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from grid import *
        
def plot_path(m, file=None, elev=None, age='red', lat=[], lng=[], lat2=[], lng2=[], title='',  interpolation='none', description=""):
    
    rcParams['image.composite_image'] = True

    plt.figure(figsize=(10,10), facecolor='w', edgecolor='k')#, figsize=(14, 14))
       
    m.drawcoastlines()
    m.drawcountries()
    #m.drawstates()
    #m.fillcontinents(color='#FFFFFF',lake_color='#FFFFFF', zorder=0)  #04BAE3 gives nice blue
    m.drawmapboundary(fill_color='#FFFFFF')
    #length=round(haversine((m.lonmin, m.latmin), (m.lonmax, m.latmin))*0.1, -2)
    m.drawmapscale(m.boundarylonmin+abs(m.boundarylonmax-m.boundarylonmin)*0.9, m.llcrnrlat+0.1*abs(m.urcrnrlat-m.llcrnrlat),m.boundarylonmin+abs(m.boundarylonmax-m.boundarylonmin)*0.9, m.llcrnrlat+0.1*abs(m.urcrnrlat-m.llcrnrlat),100, barstyle='fancy', zorder=100)
    m.drawparallels(np.arange(-90,90,5), labels=[1,0,0,0],fontsize=10)
    m.drawmeridians(np.arange(0,360,5), labels=[0,0,0,1],fontsize=10)
    
    if lat!=[]:
        x, y = m(lng, lat)
        plt.scatter(x,y, marker='o', c=age, s=3, alpha=1, label='Neolithic Site', zorder=5, cmap='Reds')
        lgnd=plt.legend(scatterpoints=1, loc='lower left')
        lgnd.legendHandles[0]._sizes = [30]
        lgnd.legendHandles[0].set_color('red')

    if lat2!=[]:
        x, y = m(lng2, lat2)
        plt.scatter(x,y,  marker='o', color='grey', s=3, alpha=0.3, label='Modern City', zorder=10)
        lgnd=plt.legend(scatterpoints=1, loc='lower left',)                
        lgnd.legendHandles[0]._sizes = [30]
        lgnd.legendHandles[0].set_color('red')
        lgnd.legendHandles[1]._sizes = [30]        
    
    if interpolation=='Gaussian':
        from scipy.ndimage.filters import gaussian_filter
        filtered_arr=gaussian_filter(elev, 2)
        m.imshow(filtered_arr, alpha=0.4, origin='upper', zorder=10)
    else:
        m.imshow(elev, interpolation=interpolation, alpha=0.6, origin='upper', zorder=8)
    cbar=m.colorbar()
    cbar.set_label('Betweenness')
    
    #for contour lines to work we need to change source code x.shape[0]/2; 
    #https://stackoverflow.com/questions/44519804/indexerror-with-basemap-contour-when-using-certain-projections
    #X,Y = np.meshgrid(df.lon_left+np.array(range(df.nrows)), df.lat_low+np.array(range(df.ncols)))
    #m.contour(x=X,y=Y,data=df.elev)
    
    plt.title(title, fontsize=15)
    plt.annotate(description, (0,0),(0,-20),xycoords='axes fraction', textcoords='offset points', va='top')
    
    #plt.axis([3, 7, 3, 7])
    
    if file!=None:

        plt.savefig(fname=file, dpi=280, bbox_inches='tight', pad_inches=0)
        
    else:
        plt.show()
        
        

