import matplotlib.pyplot as plt
from matplotlib import gridspec
from astropy import wcs
import numpy as np
from scipy  import stats
##from scipy.optimize import curve_fit
##import scipy.ndimage.filters as filters


def Plot_Probabilty(X, Y, Z, w):
    Z = Z - Z.min()
    Z = 50*Z / Z.max()
##    plt.scatter(X, Y, s=Z, c=Z, alpha=0.5)
##    plt.show()

    fig = plt.figure(figsize=(8,6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[10, 1])
    ax = fig.add_subplot(gs[0], projection=w) 
    
    ax.scatter(X, Y, s=Z, c=Z, alpha=0.5)
    ax.invert_yaxis()
    ax.set_xlabel('Ra')
    ax.set_ylabel('Dec')
    plt.title('Probability map')
    plt.grid()
    plt.savefig('Probability_RaDec_map.png')
    plt.show()

    
def Plot_Map(X, Y, Z, w):
    Z_ = Z
    Z = 1 - (min(Z)-Z)/(min(Z)-max(Z))
    Z = 3+30*Z
    
    scale_y = np.arange(round(min(Z_), 0), round(max(Z_), 0), 1)
    scale_z = 1 - (min(scale_y)-scale_y)/(min(scale_y)-max(scale_y))
    scale_z = 3+30*scale_z
    scale_x = np.zeros(scale_y.shape)
    
    fig = plt.figure(figsize=(8,6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[10, 1])
    ax = fig.add_subplot(gs[0], projection=w) #, projection=w
    
    ax.scatter(X, Y, s=Z, c='red', alpha=0.5)
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 300)
    ax.invert_yaxis()
    ax.set_xlabel('Ra')
    ax.set_ylabel('Dec')
    plt.title('Cluster map')
    plt.grid()
    
    ax2 = fig.add_subplot(gs[1])
    ax2.scatter(scale_x, scale_y, s=scale_z, c='red', alpha=0.5)
    ax2.set_ylim(ax2.get_ylim()[::-1])
    ax2.yaxis.tick_right()
    ax2.set_xticks([])
    ax2.set_ylabel ('G mag')
    ax2.yaxis.set_label_position("right")
    plt.savefig('cluster_map.png')
    plt.show()

def Draw_Center(ROI, x, y):
    Max = np.unravel_index(ROI.argmax(), ROI.shape) #get maximum
    plt.imshow(ROI, cmap=plt.cm.jet) #cm.gist_earth_r
    plt.plot(x, y, 'ro', fillstyle='none', markersize = 5, label='centroid')
    plt.plot(Max[1], Max[0], 'gx', fillstyle='none', markersize = 5, label='maximum')
    plt.xlabel('pix')
    plt.ylabel('pix')
    plt.legend()
    plt.title('Cluster center')
    plt.savefig('cluster_center.png')
    plt.show()

def Plot_PDF(Z, pix_x, pix_y, Pbl, w, C_Ra, C_De): #, center, s):
    Pbl = Pbl - Pbl.min()
    Pbl = 50*Pbl / Pbl.max()
    
    fig = plt.figure(figsize=(9,7))
    ax = fig.add_axes([0.05, 0.08, 0.8, 0.8], projection=w) #, projection=w

    ##draw PDF and stars
    pdf = ax.imshow(Z, cmap=plt.cm.jet) # cm.gist_earth_r
    ax.set_xlim(0, 300)
    ax.set_ylim(300, 0)
    ax.set_xlabel('Ra')
    ax.set_ylabel('Dec')
    plt.title('Probability')
    plt.grid()
    cbaxes = fig.add_axes([0.8, 0.08, 0.05, 0.8]) 
    fig.colorbar(pdf, cax = cbaxes, label="stars per square arcmin")
##    ax.plot(pix_x, pix_y, 'k.', markersize=2, alpha=0.3)
    ax.scatter(pix_x, pix_y, s=10, c='black', alpha=0.3)
    ax.scatter(C_Ra, C_De, transform=ax.get_transform('icrs'), s=300, edgecolor='white', facecolor='none')
##    ax.scatter(C_x, C_y, s=R_pix, edgecolor='white', facecolor='none')

    plt.savefig('cluster_PDF.png')
    plt.show()

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
    plt.plot(bins, counts, label = 'Background')
    plt.title('Histogram for PDF')
    plt.vlines(Bckgrnd_hist, 0, max(hist[0]), color='blue', label = 'Mode')
    plt.xlabel('Density')
    plt.ylabel('Bins')
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
    
def Ra_Dec_Density(Cat, RA, DEC, Cone): ##Dataset, center, radii
    ##create WCS solution for RaDec->XY projection
    Scale = Cone/150.                           ##set scale for fig size 800x800
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [150, 150]                    ##refernce pixel of WCS at center of fig
    w.wcs.crval = [RA, DEC]                     ##RaDec of reference pixel (deg)
    w.wcs.cdelt = np.array([-Scale, -Scale])     ##scale deg/pix
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]      ##type of projection
    w.wcs.cunit = ['deg', 'deg']                ##WCS units

    pix_x, pix_y = w.wcs_world2pix(Cat['RAJ2000'], Cat['DEJ2000'], 0)   ##transform RaDec to XY


    
    ##plot cluster map
    Plot_Map(pix_x, pix_y, Cat['Gmag'], w)

    ##calculation of PDF
    Mesh_Size = 300j
    X, Y = np.mgrid[0:300:Mesh_Size, 0:300:Mesh_Size] ##300
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([pix_x, pix_y])
    kernel = stats.gaussian_kde(values, 'silverman')
    Z = np.reshape(kernel(positions).T, X.shape)
    Pix_Area = (300/(Mesh_Size.imag - 1))*(300/(Mesh_Size.imag - 1)) ##pixel area 
    Z = Z*Pix_Area*len(pix_x)                   ##probabilty per pixel
##    print(Z.sum())                              ##full probability ~ number of stars
    Z = Z / (3600 * Scale * Scale * Pix_Area)   ####probabilty per square arcmin

    Z = np.rot90(Z)
    Z = np.flipud(Z)

    ##calc PDF for stars
    Stars_PDF = kernel(values)
    Stars_PDF = Stars_PDF*Pix_Area*len(pix_x)
    Stars_PDF = Stars_PDF / (3600 * Scale * Scale * Pix_Area)
    
    ##get Background
    Bkg = Get_Bkg(Stars_PDF)
    
    ##delete background
    Stars_PDF = (Stars_PDF-Bkg)/Stars_PDF
    Stars_PDF = np.clip(Stars_PDF, 0, None)
    
    
    ##get max
    Max = np.unravel_index(Z.argmax(), Z.shape) #get maximum

    ##fit center
    R = 50
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

    C_C = np.array([x + x_start, (y + y_start)])
##    print(C_C)
    R_pix = np.sqrt((pix_x - C_C[0])**2 + (pix_y - C_C[1])**2)
    plt.plot(R_pix, Stars_PDF, 'ro')
    plt.title('Probability vs distance from center(pix)')
    plt.show()
    
    R_deg = R_pix*Scale
    Cat['X'] = pix_x
    Cat['Y'] = pix_y
    Cat['R_deg'] = R_deg
    Cat['RaDecPb'] = Stars_PDF
    
    Center = w.wcs_pix2world(C_C[0], C_C[1], 0) 

    Plot_PDF(Z, pix_x, pix_y, Stars_PDF, w, Center[0], Center[1])
    return Cat, Center
    


