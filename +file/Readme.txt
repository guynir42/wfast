Package for controlling input/output of data products. 
Author: guy.nir@weizmann.ac.il

*AstroData defines all the data products that need to be moved around:
images, timestamps, psfs, cutouts, lightcurves, and so on. 

*Reader gets data from files (typically HDF5 files). 

*BufferWheel keeps a bunch of structs with the data products, 
and can dump them to disk asynchronuously using a mex file. 
It is also used to get data from the camera. 

*Deflator simply connects a Reader and a BufferWheel and allows loading
and resaving data in a compressed format. 

