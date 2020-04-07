Package for dealing with astronomical image metadata. 
Author: guy.nir@weizmann.ac.il

Header object contains all the metadata, 
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
This is mostly unused now. 
We used it when observing single or double stars
when working on high contrast imaging. 

*WorldCoordinates: keeps the info on the transformation
between sky and image coordinates. 
Can be used to give the field rotation and central coordinates.

*Catalog: takes star positions and matches them, 
using Eran's astrometry, to the GAIA catalog. 
Saves the information in a matlab table that can 
be used or saved to disk. 

*Parameters: inherits all from Header, 
and is used only for backward compatibility. 
Use h = cast(p) to turn Parameters object into Header object. 

