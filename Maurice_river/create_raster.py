# -*- coding: utf-8 -*-
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat 
from osgeo import gdal
from osgeo import osr
import os 
import math 
"""
Created on Mon Apr  3 09:53:50 2023
Toolbox of post-processing and plotting functions 
@author: gmedley - coordinate rotation added by ptaeb
"""
"""
PCS_NAD83_UTM_zone_3N =	26903
   PCS_NAD83_UTM_zone_4N =	26904
   PCS_NAD83_UTM_zone_5N =	26905
   PCS_NAD83_UTM_zone_6N =	26906
   PCS_NAD83_UTM_zone_7N =	26907
   PCS_NAD83_UTM_zone_8N =	26908
   PCS_NAD83_UTM_zone_9N =	26909
   PCS_NAD83_UTM_zone_10N =	26910
   PCS_NAD83_UTM_zone_11N =	26911
   PCS_NAD83_UTM_zone_12N =	26912
   PCS_NAD83_UTM_zone_13N =	26913
   PCS_NAD83_UTM_zone_14N =	26914
   PCS_NAD83_UTM_zone_15N =	26915
   PCS_NAD83_UTM_zone_16N =	26916
   PCS_NAD83_UTM_zone_17N =	26917
   PCS_NAD83_UTM_zone_18N =	26918
   PCS_NAD83_UTM_zone_19N =	26919
   PCS_NAD83_UTM_zone_20N =	26920
   PCS_NAD83_UTM_zone_21N =	26921
   PCS_NAD83_UTM_zone_22N =	26922
   PCS_NAD83_UTM_zone_23N =	26923
   PCS_NAD83_Alabama_East =	26929
   PCS_NAD83_Alabama_West =	26930
   PCS_NAD83_Alaska_zone_1 =	26931
   PCS_NAD83_Alaska_zone_2 =	26932
   PCS_NAD83_Alaska_zone_3 =	26933
   PCS_NAD83_Alaska_zone_4 =	26934
   PCS_NAD83_Alaska_zone_5 =	26935
   PCS_NAD83_Alaska_zone_6 =	26936
   PCS_NAD83_Alaska_zone_7 =	26937
   PCS_NAD83_Alaska_zone_8 =	26938
   PCS_NAD83_Alaska_zone_9 =	26939
   PCS_NAD83_Alaska_zone_10 =	26940
   PCS_NAD83_California_1 =	26941
   PCS_NAD83_California_2 =	26942
   PCS_NAD83_California_3 =	26943
   PCS_NAD83_California_4 =	26944
   PCS_NAD83_California_5 =	26945
   PCS_NAD83_California_6 =	26946
   PCS_NAD83_Arizona_East =	26948
   PCS_NAD83_Arizona_Central =	26949
   PCS_NAD83_Arizona_West =	26950
   PCS_NAD83_Arkansas_North =	26951
   PCS_NAD83_Arkansas_South =	26952
   PCS_NAD83_Colorado_North =	26953
   PCS_NAD83_Colorado_Central =	26954
   PCS_NAD83_Colorado_South =	26955
   PCS_NAD83_Connecticut =	26956
   PCS_NAD83_Delaware =	26957
   PCS_NAD83_Florida_East =	26958
   PCS_NAD83_Florida_West =	26959
   PCS_NAD83_Florida_North =	26960
   PCS_NAD83_Hawaii_zone_1 =	26961
   PCS_NAD83_Hawaii_zone_2 =	26962
   PCS_NAD83_Hawaii_zone_3 =	26963
   PCS_NAD83_Hawaii_zone_4 =	26964
   PCS_NAD83_Hawaii_zone_5 =	26965
   PCS_NAD83_Georgia_East =	26966
   PCS_NAD83_Georgia_West =	26967
   PCS_NAD83_Idaho_East =	26968
   PCS_NAD83_Idaho_Central =	26969
   PCS_NAD83_Idaho_West =	26970
   PCS_NAD83_Illinois_East =	26971
   PCS_NAD83_Illinois_West =	26972
   PCS_NAD83_Indiana_East =	26973
   PCS_NAD83_Indiana_West =	26974
   PCS_NAD83_Iowa_North =	26975
   PCS_NAD83_Iowa_South =	26976
   PCS_NAD83_Kansas_North =	26977
   PCS_NAD83_Kansas_South =	26978
   PCS_NAD83_Kentucky_North =	26979
   PCS_NAD83_Kentucky_South =	26980
   PCS_NAD83_Louisiana_North =	26981
   PCS_NAD83_Louisiana_South =	26982
   PCS_NAD83_Maine_East =	26983
   PCS_NAD83_Maine_West =	26984
   PCS_NAD83_Maryland =	26985
   PCS_NAD83_Massachusetts =	26986
   PCS_NAD83_Massachusetts_Is =	26987
   PCS_NAD83_Michigan_North =	26988
   PCS_NAD83_Michigan_Central =	26989
   PCS_NAD83_Michigan_South =	26990
   PCS_NAD83_Minnesota_North =	26991
   PCS_NAD83_Minnesota_Cent =	26992
   PCS_NAD83_Minnesota_South =	26993
   PCS_NAD83_Mississippi_East =	26994
   PCS_NAD83_Mississippi_West =	26995
   PCS_NAD83_Missouri_East =	26996
   PCS_NAD83_Missouri_Central =	26997
   PCS_NAD83_Missouri_West =	26998
"""


def make_swan_raster_StructuredGrid(rastername,matfile,varname,EPSG,xMin,yMin,gridx,gridy,gridsize,angle_degrees):
    matfile = loadmat(matfile)
    var = matfile[varname]

    yMax = yMin * gridy
    #Turn into a raster

    #make geotransform (in state plane m)
    #geotransform = (xMin,gridsize,0,yMin,0, -gridsize)  #stateplanexmin, gridsize, scalefactor, stateplaneymin, scalefactor, gridsize)

    #calculate the angle in radians
    angle_radians = math.radians(angle_degrees)

    #calculate geotransform by accounting for rotated coordinates
    geotransform = (xMin, gridsize * math.cos(angle_radians), -gridsize * math.sin(angle_radians), yMin,
                   gridsize * math.sin(angle_radians), gridsize * math.cos(angle_radians))

    driver = gdal.GetDriverByName('GTiff')
    output_raster = driver.Create(rastername,
                                  gridx,
                                  gridy,
                                  1,
                                  gdal.GDT_Float32)  # Open the file #+1 adds extra cell to each end of raster object to make the raster not shift

    output_raster.SetGeoTransform(geotransform)
    srs = osr.SpatialReference()                       # Establish its coordinate encoding
    srs.ImportFromEPSG(EPSG)
    output_raster.SetProjection( srs.ExportToWkt() )   # Exports the coordinate system
    output_raster.GetRasterBand(1).WriteArray(var)    # Writes the array to the raster
    output_raster.FlushCache()
    print("raster finished")
