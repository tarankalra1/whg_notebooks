# -*- coding: utf-8 -*-
"""
Created on Fri May 14 14:34:24 2021
make raster from SWAN grid outputs from a STRUCTURED grid with a single horizontal resolution  -- translated from matlab
Makes swan HS raster in meters
Or setup raster in FEET (output setup)
Needs Gdal, Scipy, and numpy installed to enviroment
@author: gmedley
"""

import numpy as np
from osgeo import gdal
from osgeo import osr
import os
import math
import scipy.io as sci

#name rasters
rastername = "Restored_HS_50perc_Ft.tif"

os.chdir(r'X:\2020-0048_NFWF\modeling\NCRF_2019\65760_CoralReefs_Guam\Swan_Modeling\SwanRuns\SouthGrid\Restored_Reef\50perc')
matfi = sci.loadmat('50perc_RestoredReef_outs')
x = matfi['Xp']
y = matfi['Yp']
Hsig = matfi['Hsig']
#Hsig = matfi['Setup']
Hsig = Hsig * 3.28084
gridsize = 15

gridx = 579 #j direction 
gridy = 463 #i direction

xMin = 247840 #this is your x origin in state plane m (or whatever you used)
yMin = 1469140 #this is your y origin in state plane m (or whatever you used)
xMax = gridx*gridsize + xMin # gets your dimension in state plane m by your resolution and grid size
yMax = gridy*gridsize + yMin

Xm = np.linspace(xMin, xMax, math.floor((xMax - xMin)/gridsize))
Ym = np.linspace(yMin, yMax, math.ceil((yMax - yMin)/gridsize))

Ymsiz = np.size(Ym)
Xmsiz = np.size(Xm)

#Turn into a raster
#make geotransform (in state plane m)
os.chdir(r'X:\2020-0048_NFWF\modeling\NCRF_2019\65760_CoralReefs_Guam\Swan_Modeling\SwanRuns\SouthGrid\Restored_Reef\50perc')
geotransform = (xMin,gridsize,0,yMin,0, gridsize)  #stateplanexmin, gridsize, scalefactor, stateplaneymin, scalefactor, gridsize) 

output_raster = gdal.GetDriverByName('GTiff').Create(rastername, Xmsiz, Ymsiz, 1, gdal.GDT_Float32)  # Open the file #+1 adds extra cell to each end of raster object to make the raster not shift
output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
srs = osr.SpatialReference()                 # Establish its coordinate encoding
srs.ImportFromEPSG(32655)  
output_raster.SetProjection( srs.ExportToWkt() )   # Exports the coordinate system 
                                                   # to the file
output_raster.GetRasterBand(1).WriteArray(Hsig)   # Writes the array to the raster

output_raster.FlushCache()