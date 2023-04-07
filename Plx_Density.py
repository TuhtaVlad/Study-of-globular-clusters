import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.special import erf

##calculate center of mass
def centroid(arr):
    total=0
    Ry=arr.shape[0]/2
    Y_index = np.arange(0, arr.shape[0], 1) ##index array
    total=np.sum(arr)
    Y = np.sum(arr*Y_index)/total
    return Y

def func(x, a, x0, sigma):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

##calculate center of mass
def Gauss_Center(arr):
    X = np.arange(0, arr.shape[0], 1)
    x0=arr.shape[0]/2
    popt, pcov = curve_fit(func, X, arr, p0=(arr.max(), x0, 1.))
    y_fit = func(X, *popt)
    plt.plot(X, y_fit, color='r')
    plt.plot(X, arr, 'go')
    plt.xlabel('Local distance (pix)')
    plt.ylabel('Density')
    plt.legend()
##    plt.plot(Dist, np.asarray(Cat['RPlx']), 'g.')
    plt.title('Density maximum')
    plt.show()
    return popt[1]

def Plx_Density(Cat, Center, Pb_Factor):
    R_deg = Cat['R_deg']
    RaDecPb = Cat['RaDecPb']
    Ap_Rad = np.max(R_deg[np.where(RaDecPb>Pb_Factor)])
    print('Core raduis (deg) = %5.2f' % Ap_Rad)
    Max_Rad = np.max(R_deg)
    print('Max raduis (deg) = %5.2f' % Max_Rad)
    Ap_Area = np.pi * Ap_Rad**2
    Max_Area = np.pi * Max_Rad**2 - Ap_Area

    Dist = 1000/np.asarray(Cat['Plx'])
    Dist_Err = 1000*np.asarray(Cat['e_Plx'])/(np.asarray(Cat['Plx'])**2)
##    plt.plot(Dist, Dist_Err, 'r.')
##    plt.title('Distance Errors vs Distance')
##    plt.show()
    
    Dist_Core = Dist[np.where(R_deg<Ap_Rad)]
    Err_Core = Dist_Err[np.where(R_deg<Ap_Rad)]
    Dist_Bkg = Dist[np.where(R_deg>=Ap_Rad)]
    Err_Bkg = Dist_Err[np.where(R_deg>=Ap_Rad)]

    plt.plot(Dist_Core,  Err_Core, 'r.', label='Core')
    plt.plot(Dist_Bkg, Err_Bkg, 'b.', label='Background')
    plt.xlabel('Distance (pc)')
    plt.ylabel('Errors (pc)')
    plt.legend()
    plt.title('Distance Errors')
    plt.show()
    
##    plt.plot(Dist_Ap,  Err_Ap*Dist_Ap*Dist_Ap/1000, 'r.')
##    plt.plot(Dist_Bkg, Err_Bkg*Dist_Bkg*Dist_Bkg/1000, 'b.')
##    plt.title('Errors')
##    plt.show()
    
    X = np.linspace(0, 10000, 1000)
    kernel = stats.gaussian_kde(Dist_Core.ravel())
    Core = kernel(X)
    Core = Core * 10 * len(Dist_Core)/Ap_Area

    kernel = stats.gaussian_kde(Dist_Bkg.ravel())
    Bkg = kernel(X)
    Bkg = Bkg * 10 * len(Dist_Bkg)/Max_Area

##    Dist_Best = Dist_Ap[np.where(Err_Ap > 20)]
##    kernel = stats.gaussian_kde(Dist_Best.ravel())
##    Best = kernel(X)
##    Best = Best * 10 * len(Dist_Best)/Max_Area

##    fig, ax = plt.subplots(figsize=(10, 6))
    plt.hist(Dist_Core, bins=1000, color = 'red', alpha=0.5) #, normed=True
    plt.hist(Dist_Bkg, bins=1000, color = 'blue', alpha=0.5) #, normed=True
##    plt.hist(Dist_Best, bins=1000, color = 'green', alpha=0.5) #, normed=True
    plt.plot(X, Core, 'r.', label='Core')
    plt.plot(X, Bkg, 'b.', label='Background')
    plt.plot(X, Core - Bkg, 'g.', label='Difference')
##    plt.plot(X, Best, 'm.')
    plt.xlabel('Distance (pc)')
    plt.ylabel('Stars / 10 pc')
    plt.title('Density vs distance')
    plt.legend()
##    plt.savefig('Distance_PDF.png')
    plt.show()

    ##get max
    Y = Core - Bkg
    Max = Y.argmax() #get maximum

    ##fit center
    R = 15
    start = int(Max-R)
    if start < 0:
        y_start = 0
    stop = int(Max+R)
    if stop > Y.shape[0]:
        stop = Y.shape[0]
    ROI = Y[start:stop]
    Center = Gauss_Center(ROI)
    C_Dist = (Center + start)*10.
##    print(Center)
    C_Rad = C_Dist * (Ap_Rad*3600/206265)
    print('Center distance (pc) = %9.1f' % C_Dist)
    print('Cluster radius (pc) = %9.1f' % C_Rad)

    ##calc probability for center of cluster for every star using Plx and error
    Pbl_Center = np.exp(-((Dist - C_Dist)**2)/(2*Dist_Err**2)) /(Dist_Err * np.sqrt(2*np.pi))
##    Cat['PlxPb'] = Pbl

    ##calc probability as integral for diameter of cluster
    S1 = C_Dist - C_Rad
    S2 = C_Dist + C_Rad
    P1 = 0.5* (1 + erf((S1 - Dist)/np.sqrt(2*Dist_Err)))
    P2 = 0.5* (1 + erf((S2 - Dist)/np.sqrt(2*Dist_Err)))
##    plt.plot(P1, 'r.')
##    plt.plot(P2, 'b.')
##    plt.show()
    
    Pbl_erf = P2-P1
    plt.plot(Dist, Pbl_erf, 'r.')
    plt.xlabel('Dist')
    plt.ylabel('Pbl_erf')
##    plt.title('Distance probability')
    plt.show()
    Cat['PlxPb'] = Pbl_erf
    Pbl = Pbl_erf

    
    G = Cat['Gmag']
    B = Cat['BPmag']
    R = Cat['RPmag']
    
##    plt.plot(B-R, G, 'r.', alpha = 0.5, label='all')
    Index = np.where(Pbl_erf > 0.01)[0]
    plt.plot((B-R)[Index], G[Index], 'g.', alpha = 0.5, label='Pbl_erf > 0.01' )
##    Index = np.where(RaDecPb > 0.1)[0]
##    plt.plot((B-R)[Index], G[Index], 'b.', alpha = 0.5, label='RaDecPb > 0.1' )
    Z = Pbl*RaDecPb
    Index = np.where(Z > 0.005)[0]
    plt.plot((B-R)[Index], G[Index], 'm.', alpha = 0.5, label='Z > 0.005' )
    plt.gca().invert_yaxis()
    plt.legend()
    plt.show()
    
    Pbl_Core = Pbl[np.where(R_deg<Ap_Rad)]
    Pbl_Bkg  = Pbl[np.where(R_deg>=Ap_Rad)]
    
    plt.plot(Dist_Core, Pbl_Core, 'r.')
    plt.plot(Dist_Bkg, Pbl_Bkg, 'b.')
    plt.xlabel('Distance')
    plt.ylabel('Pbl')
    plt.title('Distance probability')
    plt.show()
    
    Pbl_Core = Pbl_Core * RaDecPb[np.where(R_deg<Ap_Rad)]
    Pbl_Bkg  = Pbl_Bkg  * RaDecPb[np.where(R_deg>=Ap_Rad)]
    plt.plot(Dist_Core, Pbl_Core, 'r.')
    plt.plot(Dist_Bkg, Pbl_Bkg, 'b.')
    plt.xlabel('Distance')
    plt.ylabel('Pbl')
    plt.title('Full probability')
    plt.show()
    
    plt.plot(Pbl, RaDecPb, 'r.')
    plt.xlabel('Distance pbl')
    plt.ylabel('RaDec Pbl')
    plt.show()

    ##calc 3D distance from center of cluster
    R2D = R_deg * Dist * 3600/206265
    R3D = np.sqrt(R2D**2 + (Dist - C_Dist)**2)
    plt.plot(R3D[np.where(R_deg<Ap_Rad)], Pbl_Core, 'r.')
    plt.plot(R3D[np.where(R_deg>=Ap_Rad)], Pbl_Bkg, 'b.')
    plt.xlabel('Distance from center')
    plt.ylabel('Pbl')
    plt.title('Full probability vs distance from center')
    plt.show()



    
    Z = Z - Z.min()
    Z = 1 + 100*Z / Z.max()
    plt.scatter(R2D, Dist - C_Dist, s=Z, c=Z, )
    plt.xlabel('R2D')
    plt.ylabel('Z')
    plt.title('Z Map')
    plt.show()
    
    return Cat, C_Dist, C_Rad

    
