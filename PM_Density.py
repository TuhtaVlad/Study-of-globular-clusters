import matplotlib.pyplot as plt
from astropy import wcs
import numpy as np
from matplotlib import gridspec
import matplotlib.cm as cm
import matplotlib
import matplotlib.pyplot
from scipy import stats

def Plot_Map(X, Y, Z):
    Z = Z - Z.min()
    Z = 50*Z / Z.max()
    plt.scatter(X, Y, s=Z, c=Z, alpha=0.5)
    plt.show()
    
def Plot_PDF(Z, xmin, xmax, ymin, ymax, pix_x, pix_y, C_C):
    fig = plt.figure(figsize=(9,7))
    ax = fig.add_axes([0.05, 0.08, 0.8, 0.8])
    pdf = ax.imshow(Z, cmap=plt.cm.jet, extent=[xmin, xmax, ymin, ymax]) #, extent=[xmin, xmax, ymin, ymax
    ax.set_xlabel('pmRa (mas/year)')
    ax.set_ylabel('pmDec (mas/year)')
    plt.title('Density')
    plt.grid()
    cbaxes = fig.add_axes([0.8, 0.08, 0.05, 0.8]) 
    fig.colorbar(pdf, cax = cbaxes, label="stars per square mas/year")
    ax.plot(pix_x, pix_y, 'w.', markersize=2, alpha=0.3)
    
    ax.scatter(C_C[0], C_C[1], s=200, edgecolor='black', facecolor='none')
    plt.show()

    
def Draw_Center(ROI, x, y):
    Max = np.unravel_index(ROI.argmax(), ROI.shape) #get maximum
    plt.imshow(ROI, cmap=plt.cm.jet) #cm.gist_earth_r
    plt.plot(x, y, 'ro', fillstyle='none', markersize = 5, label='centroid')
    plt.plot(Max[1], Max[0], 'gx', fillstyle='none', markersize = 5, label='maximum')
    plt.xlabel('pix')
    plt.ylabel('pix')
    plt.legend()
    plt.title('Proper motion cluster center ')
    plt.savefig('PM_cluster_center.png')
    plt.show()
    
def gaussian(B, C, x0, y0):
    """Returns a gaussian function with the given parameters"""
    return lambda y,x: np.exp(-((((x0-x)**2)/(2*(B**2)))+(((y0-y)**2)/(2*(C**2))))) /(2*np.pi*B*C)
##    return lambda y,x: np.exp(-(((x0-x)/B)**2+((y0-y)/C)**2)/2)

def Get_Bkg(Data):
    Bkg_median = np.median(Data)
    print('Median background=', Bkg_median)             #good background level estimate

    ##more robust estimate of background using histogram
    arr_h = np.ravel(Data)                                  #flatten array
    hist_l = Data.min()                                     #low boundaries for histigram
    hist_h = Data.max()                                     #high boundaries for histigram
    n_bin = 50                                              #set number of bin for histogram
    hist = np.histogram(arr_h, n_bin, (hist_l, hist_h))     #calc histogram
    bins = hist[1][0:-1]+(hist[1][1]-hist[1][0])/2          #array with centers of bins
    counts = hist[0]                                        #array with counts in bins
    Bckgrnd_hist = (bins[np.argmax(counts)])                #max element of histogram
    print('Histogram background=', Bckgrnd_hist)
    plt.plot(bins, counts)
    plt.vlines(Bckgrnd_hist, 0, max(hist[0]), color='blue')
    plt.show()
    return Bckgrnd_hist

##calculate center of mass
def centroid(R1, R2, R3, arr):
    total=0
    Ry=arr.shape[0]/2
    Rx=arr.shape[1]/2
    
    #mask
    X_index = np.arange(0, arr.shape[1], 1) ##index array
    Y_index = np.arange(0, arr.shape[0], 1) ##index array
    distance = np.sqrt(np.power(np.ones(arr.shape)*(X_index[None, :]-Rx), 2) + np.power(np.ones(arr.shape)*(Y_index[:, None]-Ry), 2)) ##distance array

    ##mean sky
    annulus_mask = np.copy(distance)
    annulus_mask[annulus_mask < R2]=0
    annulus_mask[annulus_mask > R3]=0
    annulus_mask[annulus_mask > 0 ]=1
    masked = arr*annulus_mask
    MSky=np.median(masked[np.nonzero(masked)])

    ##centroid
    aperture_mask = np.copy(distance)
    distance[distance <= R1]=1
    distance[distance > R1] =0
    masked = arr*distance
    total=np.sum(masked)
    
    X = np.sum(masked*X_index[None, :])/total
    Y = np.sum(masked*Y_index[:, None])/total
    return X, Y

def PM_Density(Cat, Range):
##    Index = np.where((Cat['pmRA']<Range) & (Cat['pmRA']>-Range) & (Cat['pmDE']<Range) & (Cat['pmDE']>-Range))[0]
##    Cat = Cat[Index]
##    print (str(len(Index)) + ' objects after cleaning')
    
    pix_x = np.array(Cat['pmRA'])
    pix_y = np.array(Cat['pmDE'])
    RaDecPb =  Cat['RaDecPb']
    PlxPb = Cat['PlxPb']
    Z = PlxPb*RaDecPb
    Z = Z - Z.min()
    Z = 1 + 50*Z / Z.max()
    plt.scatter(pix_x, pix_y, s=Z, c=Z, )
    plt.xlabel('pmRA')
    plt.ylabel('pmDE')
    plt.title('Full probability')
    plt.show()
    
    Scale = 150/Range ##set scale as pix/mas
    
    Mesh_Size = 301j
    X, Y = np.mgrid[-Range:Range:Mesh_Size, -Range:Range:Mesh_Size]
    
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([pix_x, pix_y])
    kernel = stats.gaussian_kde(values, 'silverman')
    Z = np.reshape(kernel(positions).T, X.shape)
##    Pix_Area = (1/Scale)**2                     ##pixel area
##    Z = Z*Pix_Area*len(pix_x)                   ##probabilty per pixel
##    Z = Z / Pix_Area                            ##probabilty per square mas/y
    Z = Z*len(pix_x)
    Z = np.rot90(Z)
##    Z = np.log(Z)

    ##calc PDF for stars
    Stars_PDF = kernel(values)
    Stars_PDF = Stars_PDF*len(pix_x)
    
    ##get Background
    Bkg = Get_Bkg(Stars_PDF)
    
    ##delete background
    Stars_PDF = (Stars_PDF-Bkg)/Stars_PDF
    Stars_PDF = np.clip(Stars_PDF, 0, None)
    Stars_PDF = Stars_PDF #* RaDecPb

    Plot_Map(pix_x, pix_y, Stars_PDF)

    ##get Background
    Bkg = Get_Bkg(Z)
    
    ##delete background
    Z = (Z-Bkg)/Z
    Z = np.clip(Z, 0, None)
    
    ##get max
    Max = np.unravel_index(Z.argmax(), Z.shape) #get maximum
##    print(Max)
    
    ##fit center R~0.5mas
    R = 0.5*Scale
    y_start = int(Max[0]-R)
    if y_start < 0:
        y_start = 0
    y_stop = int(Max[0]+R)
    if y_stop > Z.shape[0]:
        y_stop = Z.shape[0]
    x_start = int(Max[1]-R)
    if x_start < 0:
        x_start = 0
    x_stop = int(Max[1]+R)
    if x_stop > Z.shape[0]:
        x_stop = Z.shape[0]
        
    ROI = Z[y_start:y_stop, x_start:x_stop]
    x, y = centroid(int(min(ROI.shape)/4), int(min(ROI.shape)/2), min(ROI.shape), ROI)
    
    ##draw center
    Draw_Center(ROI, x, y)
    
    C_C = np.array([-150+x + x_start, 150-(y + y_start)]) / Scale
    print('Center pmRa = %7.3f' % C_C[0])
    print('Center pmDec = %7.3f' % C_C[1])

    Plot_PDF(Z, -Range-0.5/Scale, Range+0.5/Scale, -Range-0.5/Scale, Range+0.5/Scale, pix_x, pix_y, C_C )


def PM_Density_Gauss(Cat, Range):
    Index = np.where((Cat['pmRA']<Range) & (Cat['pmRA']>-Range) & (Cat['pmDE']<Range) & (Cat['pmDE']>-Range))[0]
    Cat = Cat[Index]
    print (str(len(Index)) + ' objects after cleaning')

    Size = 200 ##size of picture
    Scale = Size / (2* Range) ##scale factor = pix/mas
    PS = np.zeros(shape=(Size, Size)) ##probability surface
    
    pmRa = Scale * np.array(Cat['pmRA'])
    pmDec = Scale * np.array(Cat['pmDE'])
    e_pmRa = Scale * np.array(Cat['e_pmRA'])
    e_pmDec = Scale * np.array(Cat['e_pmDE'])

    for ii in range (0, pmRa.shape[0]):
        p = [e_pmRa[ii], e_pmDec[ii], pmRa[ii]+PS.shape[1]/2, pmDec[ii]+PS.shape[0]/2]
        Gauss_model = gaussian(*p)(*np.indices(PS.shape))
        PS = PS + Gauss_model

    draw_pic(PS)
