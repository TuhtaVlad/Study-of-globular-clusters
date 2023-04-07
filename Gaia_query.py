from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import Angle
import numpy as np

##Vizier query of Gaia DR2 data'

def Get_Gaia(Ra, Dec, Radius, Max): 
    ##make astropy SkyCoord object
    Coo = SkyCoord(ra=Ra*u.degree, dec=Dec*u.degree, frame='icrs')
    ##set search cone
    Cone = Angle(Radius * u.deg) #0.35
    ##set columns
##    V = Vizier()
    V = Vizier(columns=['RAJ2000', 'e_RAJ2000', 'DEJ2000', 'e_DEJ2000', \
                           'pmRA',    'e_pmRA',    'pmDE',    'e_pmDE', \
                        'Plx', 'e_Plx', 'RPlx', \
                        'Gmag', 'e_Gmag', 'BPmag', 'e_BPmag', 'RPmag', 'e_RPmag', \
                        'RV', 'e_RV', 'Teff'])
    
    ##set limit of rows
    V.ROW_LIMIT=int(Max)
    Vizier_result = V.query_region(Coo, radius=Cone, catalog=['I/345'])
    Vizier_stars = Vizier_result[0]
    
    return Vizier_stars										 

##Cat = Get_Gaia(0, 0, 0.1, 100)
