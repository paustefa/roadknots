#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 13:26:06 2018

It seems to be important to have the newest version of scipy. 
There was an issue with pickling dok_matrix in older versions.

To Dos

    - make path computations more efficient -> parallelize + shared memory 
    
    see 
    https://docs.python.org/2/library/multiprocessing.html#multiprocessing.Array
    https://stackoverflow.com/questions/14124588/shared-memory-in-multiprocessing
    
@author: stefan
"""

##############################################################################
# 1. Define function calculating distance
##############################################################################

import numpy as np
from haversine import haversine
from collections import Counter
from scipy.sparse import dok_matrix
from scipy.sparse.csgraph import dijkstra


# parameters governing walking times on flat and sloped surfaces
coeff_flat=0.72
coeff_up=6
coeff_down_mod=2
coeff_down_steep=-2 
# those parameters are from Alessio/Aitken/Langmuir

coeff_sea_premium=0.9
coeff_boat_loading=1
"""
thomas:
% I totally made up those numbers: 
% Travel by sea (flat surface) is X% faster walking the same distance on
% flat land; and I add a Y hours penalty to load/unload from a boat. I
% arbitrarily use X=0.1 and Y=1.

% Note that all distances are measured in meters (not in kms). Travel times
% are in hours (not that we really care...)."""

class grid:
    
    def __init__(self, elev, cellsize, xllcorner, yllcorner, NA_value, snapshot, snap_rowcol):
        self.elev=elev
        self.cellsize=cellsize
        self.lon_left=xllcorner
        self.lat_low=yllcorner
        self.NA=NA_value
        self.nrows=elev.shape[0]
        self.ncols=elev.shape[1]
        self.length=self.nrows*self.ncols #calculate number of grid cells
        self.lat_up=self.lat_low+self.cellsize*self.nrows #calculate upper limit of latitude
        self.lon_right=self.lon_left+self.cellsize*self.ncols
        
        #save info on the snapshot
        self.row_l, self.lat_up_adj=self.to_row(snapshot[0][1])
        self.row_u, self.lat_low_adj=self.to_row(snapshot[1][1])
        self.col_l, self.lon_left_adj=self.to_col(snapshot[0][0])
        self.col_r, self.lon_right_adj=self.to_col(snapshot[1][0])
        
        #self.snapshot=snapshot
        self.snapshot=[(self.lon_left_adj, self.lat_low_adj), (self.lon_right_adj, self.lat_up_adj)]
        self.snap_rowcol=snap_rowcol
        
        #calculate lat and long of grid cell centers
        #self.lat=np.repeat(self.lat_up-np.arange(self.cellsize/2, self.nrows*self.cellsize, self.cellsize),self.ncols)
        self.lat=self.lat_up-np.arange(self.cellsize/2, self.nrows*self.cellsize, self.cellsize)
        #self.lng=np.tile(self.lon_left+np.arange(self.cellsize/2, self.ncols*self.cellsize, self.cellsize), self.nrows)
        self.lng=self.lon_left+np.arange(self.cellsize/2, self.ncols*self.cellsize, self.cellsize)

        #self.coord=np.concatenate((np.matrix(self.lng).T, np.matrix(self.lat).T), axis=1)
        
        #create ids / nodes
        self.nodes=np.arange(0, self.length).reshape(self.nrows, self.ncols)
        self.nodes_land=self.nodes[elev>0]
        
        #initiate graph /distance matrix and centrality measures
        self.dist=None
        self.centrality=None
        self.centrality_draws=None
        
        
        
    # this functions are problematic when the snapshot is larger than the actual grid passed. 
    # (in this case it should just give the first and last row/column
    
    def to_row(self,lat):
        lat_vec=self.lat_up-np.arange(self.cellsize/2, self.nrows*self.cellsize, self.cellsize)
        idx = (np.abs(lat_vec-lat)).argmin() #yuk
        return idx, lat_vec[idx]
    
    
    def to_col(self,lng):
        lng_vec=self.lon_left+np.arange(self.cellsize/2, self.ncols*self.cellsize, self.cellsize)
        idx = (np.abs(lng_vec-lng)).argmin()
        return idx, lng_vec[idx]       
    
    
    def to_node(self,lat,lng):
        return self.nodes[self.to_row(lat)[0], self.to_col(lng)[0]]

    #define a function that calculates the distance
    
    def distance(self, orig_lat, orig_lng, dest_lat, dest_lng):  
        
        #extract elevation of origin and destination
        orig_elev=self.elev[orig_lat, orig_lng]
        dest_elev=self.elev[dest_lat, dest_lng]
        
        #to be exact we need to add / substract half a cellsize to be calculating from center to center of grid cells
        orig=(self.lat_up-orig_lat*self.cellsize-self.cellsize/2, self.lon_left+orig_lng*self.cellsize+self.cellsize/2)
        dest=(self.lat_up-dest_lat*self.cellsize-self.cellsize/2, self.lon_left+dest_lng*self.cellsize+self.cellsize/2)
        
        #if orig_elev!=self.NA and dest_elev!=self.NA:
               
        if orig_elev>=0 and dest_elev>=0: #travel over land
            if dest_elev-orig_elev>0 : #uphill
                return (1/3600)*(coeff_flat*haversine(orig,dest)*1000+coeff_up*(dest_elev-orig_elev))
            elif dest_elev-orig_elev<0 and dest_elev-orig_elev>=-0.2125:
                return (1/3600)*(coeff_flat*haversine(orig,dest)*1000+coeff_down_mod*(dest_elev-orig_elev))
            else:
                return (1/3600)*(coeff_flat*haversine(orig,dest)*1000+coeff_down_steep*(dest_elev-orig_elev))

        elif orig_elev<0 and dest_elev<0: #travel over sea
            return (1/3600)*coeff_sea_premium*coeff_flat*haversine(orig,dest)*1000
            
        else: #in and out of water => extra penalty
            return coeff_boat_loading+(1/3600)*coeff_flat*haversine(orig,dest)*1000
        
        #return value such that this path is never taken
        #else:
        #    return 9999

##############################################################################
# 3. Define function that takes grid and outputs matrix with three columns: orig_id, dest_id, distance
##############################################################################
    
    #compute the distance to neighboring nodes, this is surprisingly fast even for large grids
    def grid_dist(self):
        
        dist=dok_matrix((self.length, self.length))
        
        for i in range(1, self.nrows-1):
            for j in range(1, self.ncols-1):
                id_orig=self.nodes[i,j]
                dist[id_orig, self.nodes[i-1,j]]=self.distance(i, j, i-1, j)
                dist[id_orig, self.nodes[i-1,j+1]]=self.distance(i, j, i-1, j+1)
                dist[id_orig, self.nodes[i,j+1]]=self.distance(i, j, i, j+1)
                dist[id_orig, self.nodes[i+1,j+1]]=self.distance(i, j, i+1, j+1)
                dist[id_orig, self.nodes[i+1,j]]=self.distance(i, j, i+1, j)
                dist[id_orig, self.nodes[i+1,j-1]]=self.distance(i, j, i+1, j-1)
                dist[id_orig, self.nodes[i,j-1]]=self.distance(i, j, i, j-1)
                dist[id_orig, self.nodes[i-1,j-1]]=self.distance(i, j, i-1, j-1)
         
        self.dist=dist
        return dist
    
    #auxiliary function, that extracts shortest paths from predecessor matrix
    #used in get_orig_paths
    def get_path(self, d, pred):
        it=pred[0][d]
        ls=[]
        while it!=-9999:
            ls.append(it)
            it=pred[0][it]
        return ls
    
    #based on origin vector and nodes to be kept as destinations this does the actual computation
    def get_orig_paths(self, origins, keep=[], dist_max=None):
        
        #initiate counter 
        res=Counter(dict.fromkeys(np.arange(0,self.length,1),0))
        
        #loop through origins and depending on whether maximal distance is specified or not compute predecessor matrix
        #with dijkstra algo
        for orig in origins:

            if dist_max!=None:
                trav_time, pred = dijkstra(self.dist, indices=[orig], limit=dist_max, directed=True, unweighted=False, return_predecessors=True)
                if len(keep)>0:
                    dest=np.where((pred!=-9999) & (keep==True))[1]
                else:
                    dest=np.where((pred!=-9999))[1]
            else:
                trav_time, pred = dijkstra(self.dist, indices=[orig], directed=True, unweighted=False, return_predecessors=True)
                dest = np.arange(0, self.length)
                if len(keep)>0:
                    dest = dest[keep]
            
            if dest!=[]:
            
                for d in dest:
                    res.update(self.get_path(d, pred))
            
            else:
                print('No nodes left.')
        return res
    
    #function that computes actual road knot scores, supports parallelized computing
    def get_centrality(self, k=100, dist_max=None, parallel=False, free_cpus=0):
        
        self.centrality_draws=k
        origins=np.random.choice(self.nodes_land, k)
        keep=self.elev.flatten()>0
        scores=Counter(dict.fromkeys(np.arange(0,self.length,1),0))
        
        if parallel==False:
            scores.update(self.get_orig_paths(origins, keep, dist_max=dist_max))
                
        else:
            
            import multiprocessing
            cpus=multiprocessing.cpu_count()
            chunk_ls=self.chunks(origins, cpus-free_cpus)
            
            from joblib import Parallel, delayed
            res=Parallel(n_jobs=cpus-free_cpus, backend="multiprocessing")(delayed(self.get_orig_paths)(chunk, keep, dist_max=dist_max) for chunk in chunk_ls)
        
            for r in res:
                scores.update(r)
    
        self.centrality=np.array(list(scores.values())).reshape(self.nrows, self.ncols)
        return self.centrality
            
    #auxiliary function to divide origins into chunks which can then be fed to parallel processes    
    def chunks(self, vec, n):
        import math
        chunk_size=math.ceil(len(vec)/n)
        for i in range(0, len(vec), chunk_size):
            yield vec[i:i+chunk_size]        
            
##############################################################################
# 5. Define function that takes asc file as input and outputs grid class object
##############################################################################

#iterator, avoids loading entire elevation grid
import csv
from itertools import islice

def iterator(file, row_start, row_end, col_start, col_end):   
    with open(file, 'r') as f:
        r = csv.reader(islice(f, 6 + row_start, 6+row_end), delimiter=' ')
        result=np.array([np.array([int(s) for s in row if s!=''])[col_start:col_end] for row in r])
        return result
    

#allows converting shape files to raster format
import rasterio
from rasterio import features
import geopandas as gpd

def to_raster(shapefile, rasterfile, size=0.1, touched=False):
    
    gpd_file = gpd.read_file(shapefile)
    
    transform = rasterio.transform.from_origin(-180, 90, size, size)
    cols = int(360 / size)
    rows = int(180 / size)
        
    with rasterio.open(rasterfile, 'w', driver='GTiff', height=rows, width=cols, count=1, dtype='uint8',
                       crs='+proj=latlong', transform=transform) as out:
        out_arr = out.read(1)
        
        shape = [(r, 1) for r in gpd_file.geometry]
        
        #all_touched option to be discussed
        burned = features.rasterize(shapes=shape, fill=0, out=out_arr, transform=out.transform, all_touched=touched)
                
        out.write_band(1, burned) 
        
        return burned


#initiates grid() instance based on coordinates and elevation file
def dist_mat(file, x_left, x_right, y_upper, y_lower):
        
    info={}
    #extract relevant information from .asc file
    with open(file, 'r') as f:
        r = csv.reader(islice(f, 0, 6), delimiter=' ')
        for row in r:
            info[row[0]]=float(row[-1])
            
    yulcorner=info['yllcorner']+info['cellsize']*info['nrows']
    share=(yulcorner-y_upper)/(info['cellsize']*info['nrows'])
    row_start=int(round(share*info['nrows']))
    
    share=(yulcorner-y_lower)/(info['cellsize']*info['nrows'])
    row_end=int(round(share*info['nrows']))
    
    share=(x_left-info['xllcorner'])/(info['cellsize']*info['ncols'])
    col_start=int(round(share*info['ncols']))
    
    share=(x_right-info['xllcorner'])/(info['cellsize']*info['ncols'])
    col_end=int(round(share*info['ncols']))
                
    elev=iterator(file, row_start, row_end, col_start, col_end)
        
    #pass the snapshot of interest
    snap=[(x_left, y_lower), (x_right, y_upper)]
    
    snap_rowcol=[(row_start, row_end), (col_start, col_end)]
      
    #notice here that for the sake of being exact, we don't pass arguments x_left und y_lower, but recalculate limits
    return grid(elev, info['cellsize'], info['xllcorner']+col_start*info['cellsize'], 
                info['yllcorner']+(info['nrows']-row_end)*info['cellsize'], info['NODATA_value'], snap, snap_rowcol)
            
    


        