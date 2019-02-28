Package for dealing with astronomical image metadata. 
Author: guy.nir@weizmann.ac.il

*Parameters object contains all the metadata, 
some of it in sub-objects (described below). 
You should give a handle to this object to every
object that does calculations on the data and needs access
to the metadata (e.g., pixel scale, exposure time). 

*Ephemeris contains the time and coordinates of the observations
and can calculate lots of astronomical angles and times 
like LST, alt/az of the observations, etc. 

*Filter contains the filter response and can be used 
to convert magnitude to counts in the detector. 

*Star keeps track of the properties of stars in the field, 
e.g., magnitude, color, xy position. A vector of Star objects
is essentially like a catalog. 

