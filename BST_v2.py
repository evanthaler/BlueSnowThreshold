import numpy as np
import os,glob
from osgeo import gdal
import matplotlib.pyplot as plt
from scipy.signal import argrelmin
import matplotlib.pyplot as plt
from scipy.stats.kde import gaussian_kde
from scipy.ndimage.filters import gaussian_filter1d

##Setting environmental variables for GDAL--might not be necessary for some users
##These lines can probably be commented out on most machines
os.environ['PROJ_LIB'] = r'C:\Users\361045\Anaconda3\envs\pygeo\Library\share\proj'
os.environ['GDAL_DATA'] = r'C:\Users\361045\Anaconda3\envs\pygeo\Library\share'



def BST (wd,inras,outras,plotting=False,sigma=3,blueband=1,meanbluethresh=0.70,nosnowthresh=1000):
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
    #Change to directory where images are located
    os.chdir(wd)
    
   
    #open image
    tif = gdal.Open(inras)
    #get image metadata
    meta = tif.GetMetadata()
    if (tif.RasterCount >1): #here we are making sure that tif is a true color image.
        #read in band 1 (change number in GetRasterBand() to read in different bands)
        blue = tif.GetRasterBand(blueband)
        #get nodata value for the band
        nodat = blue.GetNoDataValue()
        if nodat is None:
            nodat=0.
        #read band values into array
        blue = np.array(blue.ReadAsArray()).astype(float)
        blue[blue==nodat]=np.nan
        ##remove nans from blue array and flatten to 1d
        bluef =blue[~np.isnan(blue)].ravel()
        
        #need to smooth array to remove noise--this will help us find the true peaks
        #here we set the filter window size to the default of 50
        bluef_sm = gaussian_filter1d(bluef,sigma=sigma)
        #bluef_sm = bluef_sm/np.max(bluef_sm) ##uncomment to normalize data to range 0-1
        bluef_sm = bluef_sm
        meanblue = np.mean(bluef_sm)
        print('blue mean is: ', meanblue)
        if ((meanblue <(meanbluethresh*10000)) and (meanblue >(1000))):
            print('likely an image with mixed snow and non-snow pixels')
            export=True
                
    
    
            #First, need to check the order of reflectance values. Does the order go from large to small or small to large?
            #if the order goes large to small, we will index from right to left, else we index from left to right.
            if (bluef_sm[bluef_sm>np.mean(bluef_sm[bluef_sm>0])][0]) > (bluef_sm[bluef_sm>np.mean(bluef_sm[bluef_sm>0])][-1]):
                print('inverted array: So, we need to grab the last value in the array ')
                #Here we get all of the array values greater than the mean and then find the first local minimum greater than the mean
                lmin=argrelmin(bluef_sm[bluef_sm>np.mean(bluef_sm[bluef_sm>0])])[0][-1]
                #Index blue at local min--this will be the threshold for the image.
                #The index value isn't necessarily the trough bewteen the peaks in a bimodal distribution
                bluethresh = bluef_sm[bluef_sm>np.mean(bluef_sm)][lmin]
                print('threshold value is: ',bluethresh)
            else:
                #Here we get all of the array values greater than the mean and then find the first local minimum greater than the mean
                lmin=argrelmin(bluef_sm[bluef_sm>np.mean(bluef_sm[bluef_sm>0])])[0][0] 
                #Index blue at local min--this will be the threshold for the image.
                #The index value isn't necessarily the trough bewteen the peaks in a bimodal distribution
                bluethresh = bluef_sm[bluef_sm>np.mean(bluef_sm)][lmin]    
                print('threshold value is: ',bluethresh)
        elif meanblue >=(meanbluethresh*10000):
            print('likely an image with all snow  pixels')
            bluethresh = (meanbluethresh*10000)
            print('threshold value is: ',bluethresh)
            export=True
            
        elif meanblue <nosnowthresh:
            print('Are you sure this image has snowy pixels? If you think it does, try decreasing the nosnowthresh value')
            export =False
        
        
        if export ==True:
            
            if plotting==True:
                
                ##Generating plots of the blue distribution significantly increases computation time
                bluef_kde = gaussian_kde(bluef_sm)
                bluebins = np.linspace(np.min(bluef_sm),np.max(10000),50)
                bkdepdf = bluef_kde(bluebins)*100
                fs=16
                
                plt.figure(figsize=(6,6))
                plt.plot(bluebins,bkdepdf,'-k',linewidth=3)
                plt.xlabel('Blue reflectance',fontsize=fs);plt.xticks(size=fs,rotation=20)
                plt.ylabel('Probability density',fontsize=fs);plt.yticks(size=fs)
                plt.axvline(bluethresh,linestyle='--',color='b')
                plt.tight_layout()
                plt.savefig(f[:-4]+'_thresholdKDE.pdf')
                plt.show()
                    
                    
                # plt.figure(figsize=(6,6))
                # plt.hist(bluef_sm,bins=50,color='k')
                # plt.xlabel('Blue reflectance',fontsize=fs);plt.xticks(size=fs,rotation=20)
                # plt.ylabel('Pixel count',fontsize=fs);plt.yticks(size=fs)
                # plt.axvline(bluethresh,linestyle='--',color='b')
                # plt.tight_layout()
                # plt.savefig(f[:-4]+'_threshold.pdf')
                # plt.show()
            
            #generate an array of zeros with the input raster dimensions
            snowBLUE=np.zeros(blue.shape)
            #set all pixels predcited to be snow to a value of 1
            snowBLUE=np.where((blue>bluethresh),snowBLUE+1,snowBLUE)
            
        
        
            #######
            #Below are a bunch of gdal functions to create the output file
            #There shouldn't be a need to change any of the lines
            #######
            
            geotransform = tif.GetGeoTransform()
            originX = geotransform[0]
            originY = geotransform[3]
            pixelWidth = geotransform[1]
            pixelHeight = geotransform[5]
            cols = tif.RasterXSize
            rows = tif.RasterYSize
            
        
            driver = gdal.GetDriverByName('GTiff')
            outRaster = driver.Create(outras, cols, rows, 1, gdal.GDT_Float32,["COMPRESS=LZW"])
            outband=outRaster.GetRasterBand(1)
            
            outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
            outRaster.SetGeoTransform(geotransform)
            outRaster.SetProjection(tif.GetProjection())
            #outband = outRaster.GetRasterBand(1)
            outband.SetNoDataValue(0)
            outband.WriteArray(snowBLUE)
            
            
            outband.FlushCache()
            outRaster=None
            tif=None
    else:
        print('this is a single band raster and probably not a true color image--moving on to next raster')
        tif=None
        outRaster=None
        
        
        
        
#Set working directory
wd=r'C:\Users\361045\Documents\projects\ngee\imagery\t47\vegetation'
#Get list of tif files in working directory
flist=glob.glob(wd+'\\'+'*.tif')#for macOS/linux,might need to change '\\' to '/'
for f in flist: ##we're going to loop through the tifs in the directory and calculate snow on the true color images
    
    inras=f
    #generate path for output file...here we are just throwing the new snow rasters in the working directory
    outras=f[:-4]+'_snow2.tif'
    BST(wd,inras,outras,False)



