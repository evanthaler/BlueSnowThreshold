# BlueSnowThreshold

 
  This code is designed to run on PlanetScope Analytic Surface Reflectance Imagery (SR) which is a reflectance product scaled by 10,000.
  This means the range of the pixel values is 0-10,000 rather than 0-1, which is typical of other SR products. 
  The 10000 scaling factor is worked into this code, which means adjustments must be made to use this code on other SR products. 
  Future version of the code will include parameters which can be easily modified to designate a scaling factor.
  

'''
Parameters
----------
wd : string
    full path to directory where input rasters are located and where output rasters will be saved
inras : string
    name of input raster
outras : string
    name of output snow raster
plotting : Boolean, optional
    Display kernel density plot of blue reflectance values. The default is False.
sigma : integer, optional
    The window size for gaussian filtering.This is called 'sigma' in the gaussian_filter1d function.
    The default here is 5.
blueband : integer, optional
    The raster band which corresponds to the blue band. The default is 1.
meanbluethresh : float, optional
    The mean blue reflectance threshold, used to assess whether an image has a high percentage of snow cover. 
    The default is 0.70
nosnowthresh : float, optional
    The lower bound on the mean blue reflectance threshold, which is used to determine if the image is snowless
    The default is 1000
Returns
-------
None.
'''
