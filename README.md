# BlueSnowThreshold

BST.py is a Python implementation of the Blue Snow Threshold algorithm (BST) (Thaler, Crumley, Bennett (2023), Remote Sensing of Environment). The BST uses a smoothed distribution of blue reflectance values to predict a threshold value which defines snow vs non-snow pixels. The algorithm makes assumptions about the mean blue reflectance for images that are a) completely snowless,b) a mix of snow and non-snow pixels, and c) nearly completely snow-covered. Each of these assumptions is included as a function argument with default values, but they can be changed by providing values in the function call. 

Be aware that this code is designed to run on PlanetScope Analytic Surface Reflectance Imagery (SR) which is a reflectance product scaled by 10,000.
This means the range of the pixel values is 0-10,000 rather than 0-1, which is typical of other SR products. 
The 10000 scaling factor is worked into this code as a default for a function argument.  Significant modifications to the code will need to made if the range of values pixel values is not altered by a simple multiplicative scaling. 

Python GDAL is required for this code
The diptest Python package is also required

  

