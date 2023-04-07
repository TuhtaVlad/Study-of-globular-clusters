import numpy as np
from scipy import stats
from scipy import optimize
from math import sqrt
import matplotlib
import matplotlib.pyplot
import matplotlib.cm as cm

def draw_pic(Data):
##    Data = np.log(Data)
    _min = np.mean(Data)
    med=np.median(Data)
    stdv=np.std(Data)
    fig_pic = matplotlib.pyplot.figure(1)
    matplotlib.pyplot.cla()
    matplotlib.pyplot.imshow(Data, cmap=cm.Greys_r, aspect='equal',
                             norm= matplotlib.colors.Normalize(vmin=_min-2*stdv, vmax=_min+stdv*5), interpolation='nearest', origin='lower')
    matplotlib.pyplot.show()


def gaussian(B, C, x0, y0):
    """Returns a gaussian function with the given parameters"""
    return lambda y,x: np.exp(-((((x0-x)**2)/(2*B**2))+(((y0-y)**2)/(2*C**2)))) /(2*np.pi*B*C) 

def Gauss(X, Y, eX, eY, surface):
    out = np.zeros_like(surface)
    for ii in range (0, X.shape[0]):
        p = [eX[ii], eY[ii], X[ii]+surface.shape[1]/2, Y[ii]+surface.shape[0]/2]
        Gauss_model = gaussian(*p)(*np.indices(surface.shape))
        out = out + Gauss_model

    draw_pic(out)
