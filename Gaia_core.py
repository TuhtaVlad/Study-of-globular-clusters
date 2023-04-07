import numpy as np
from Gaia_query import Get_Gaia
from Ra_Dec_Density import Ra_Dec_Density

from Plx_Density import Plx_Density
from PM_Density import PM_Density
from PM_Density import PM_Density_Gauss

import warnings
warnings.simplefilter("ignore")


##get Gaia data
##set query here (Ra(deg), Dec(deg), Cone(deg), Maximum rows (int))
RA = 88.07662 #11.777 #
DEC = 32.55316 #85.246 #
Cone = 0.4
Max = 100000
Cat = Get_Gaia(RA, DEC, Cone*1.42, Max)
##print (Cat.info) ##here you can print table information

##data cleaning
##get indexes of valid data
Index = np.where((np.isfinite(Cat['Plx'])) & (Cat['Plx']>0.1) & (Cat['Plx']<7))[0] ##where parallax is exists and positive
##clear data, delete rows where parallax isn't a real
Cat = Cat[Index]
print (str(len(Index)) + ' objects after Plx cleaning')
Index = np.where((Cat['pmRA']<30) & (Cat['pmRA']>-30) & (Cat['pmDE']<30) & (Cat['pmDE']>-30))[0]
##clear data, delete rows where parallax isn't a real
Cat = Cat[Index]
print (str(len(Index)) + ' objects after PM cleaning')



####draw probability density for Ra-Dec data
Cat, Center = Ra_Dec_Density(Cat, RA, DEC, Cone)    ##add Cat['X'], Cat['Y'], Cat['R_deg'], Cat['RaDecPb']
print('Center Ra = %9.5f' % Center[0])
print('Center Dec = %9.5f' % Center[1])

Pb_Factor = 0.1
Cat, C_Dist, C_Rad = Plx_Density(Cat, Center, Pb_Factor)    ##add Cat['PlxPb']

PMRange = 5
PM_Density(Cat, PMRange)






##PM_Density_Gauss(Cat)

#plot probability for pm
##search max-min
##Size = 200
##R=2000
##X = Cat['pmRA']
##Y = Cat['pmDE']
##e_X = Cat['e_pmRA']
##e_Y = Cat['e_pmDE']
##
##Xmax = np.max(X)
##Xmin = np.min(X)
##Ymax = np.max(Y)
##Ymin = np.min(Y)
##x_Scale = 0.5 #200/(Xmax - Xmin)
##y_Scale = 0.5 #200/(Ymax - Ymin)
##print(x_Scale)
##print(y_Scale)
####Plot2D(X, Y, ('pmRa', 'pmDec'))
##surf = np.zeros(shape=(Size, Size))
##Gauss(X*x_Scale, Y*y_Scale, e_X *x_Scale, e_Y *y_Scale, surf)

###plot probability for pm
####search max-min
##Size = 200
##Xmax = np.max(Cat['RAJ2000'])
##Xmin = np.min(Cat['RAJ2000'])
##Ymax = np.max(Cat['DEJ2000'])
##Ymin = np.min(Cat['DEJ2000'])
##x_Scale = 200/(12)
##y_Scale = 200/(12)
##print(x_Scale)
##print(y_Scale)
####X_index = np.linspace(Xmax, Xmin,  num=10) ##index array
####Y_index = np.linspace(Ymax, Ymin,  num=10) ##index array
####xv, yv = np.meshgrid(X_index, Y_index)
##surf = np.zeros(shape=(Size, Size))
##Gauss(Cat['RAJ2000']*x_Scale, Cat['DEJ2000']*y_Scale, Cat['e_RAJ2000']*x_Scale, Cat['e_DEJ2000']*y_Scale, surf)



##data filtering
##get indexes of valid data
##Low_Limit = 1 ##parallax in mas = 0.001arcsec!
##Up_Limit = 100
##Index = np.where((Cat['Plx']>=Low_Limit) & (Cat['Plx']<=Up_Limit))[0] ##where Low_Limit<=parallax<=Up_Limit mas
##Cat = Cat[Index]
##print (str(len(Index)) + ' objects has valid values')

####plot positions of objects in Ra-Dec axis
##Plot_RaDec(Cat['RAJ2000'], Cat['DEJ2000'], Cat['Gmag'])


##print (Cat['RAJ2000'], Cat['DEJ2000'], Cat['Plx'])
##Plot2D(Cat['pmRA'], Cat['pmDE'], ('pmRa', 'pmDec'))
##Plot2D(Cat['pmRA'], Cat['Plx'], ('pmRa', 'Plx'))
##Plot2D(Cat['pmDE'], Cat['Plx'], ('pmDec', 'Plx'))

##Plot3D(Cat['pmRA'], Cat['pmDE'], Cat['Plx'], ('pmRa', 'pmDec', 'Plx'))
