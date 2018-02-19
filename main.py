#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 14:08:38 2018

This computes betweenness centrality for Turkey

as baseline coordinates for turkey we take
(25, 35) (46, 43)

@author: stefan
"""
import os
dir1=os.path.join(os.sep, 'home', 'stefan', 'Dropbox', 'ancient_trade_ra', '180219_roadknots_stefan')
#dir1=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

import sys
sys.path.append(os.path.join(dir1, 'Code'))

from grid import *


x_ll=25
y_ll=35

x_ur=46
y_ur=43

x_dist=abs(x_ll-x_ur)
y_dist=abs(y_ll-y_ur)


file=os.path.join(dir1, 'Input', 'GloElev_5min.asc')

fac=1

df=dist_mat(file, y_upper=y_ur+fac*y_dist, y_lower=y_ur-fac*y_dist, x_left=x_ll-fac*x_dist, x_right=x_ur+x_dist)

dist=df.grid_dist()

heat=df.get_centrality(k=10000, parallel=False, free_cpus=1, dist_max=10)

import pickle
pickle.dump(df, open(os.path.join(dir1, 'Temp', 'df.pkl'), 'wb'))


"""
from utils import plot_path
from mpl_toolkits.basemap import Basemap
m = Basemap(llcrnrlat=df.lat_low,urcrnrlat=df.lat_up, llcrnrlon=df.lon_left,urcrnrlon=df.lon_right,resolution='l', projection='mill')

cent=np.copy(heat)
cent[cent>100]=100
plot_path(m, elev=cent)
"""


