import matplotlib.colors as colors
import numpy as np
import matplotlib.pylab as plt

def plot_data(x, y, z, gsize=30, levels=15, grid_method='linear', scatter=False, cmap='jet'):
    """Plotting function
    This plots a rasterscan as an intensity image
    the scatter parameter adds markers to indecate the data points
    x,y,z must be the same length and z is the amplitude"""
    from scipy.interpolate import griddata
    import matplotlib.pyplot as plt
    import numpy as np
    print(type(z))
    if type(z) == np.ma.core.MaskedArray:
        if np.ndim(z) == 1:
            m = ~z.mask
            x = x[m]
            y = y[m]
            z = z[m]
        if np.ndim(z) == 2:
            zz = np.empty(z.shape[0])
            mz = np.zeros(z.shape[0], dtype='bool')[:]
            for i in np.arange(z.shape[0]):
                zz[i] = np.ma.mean(z.data[i, :])
                if np.ma.mean(z.mask[i, :]) == 1.0: mz[i] == True
            print(np.where(mz))
            x = x[~mz]
            y = y[~mz]
            z = zz[~mz]

    # define grid.
    npts = z.shape[0]
    xi = np.linspace(x.min(), x.max(), gsize)
    yi = np.linspace(y.min(), y.max(), gsize)
    # grid the data.
    zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method = grid_method)
    #zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='cubic')
    # print(zi.shape)
    # contour the gridded data, plotting dots at the randomly spaced data points.
    #CS = plt.contour(xi, yi, zi, 15, linewidths=0.5, colors='k')
    #CS = plt.contourf(xi, yi, zi, 255, cmap=plt.cm.jet, linewidth=0.0)
    CS = plt.imshow(zi[::-1, :], extent=(x.min(),x.max(),y.min(),y.max()), cmap=plt.get_cmap(cmap,255), aspect='auto')
    plt.colorbar()
    CS = plt.contour(xi, yi, zi, levels, linewidths=0.5, colors='k')
    #plt.contourf(xi, yi, zi, 255, cmap=plt.cm.jet, linewidth=0.0),
    # CS = plt.contourf(xi, yi, zi, 100, cmap=plt.get_cmap('jet',255))
    #plt.colorbar()  # draw colorbar
    # plot data points.
    plt.xlim(x.min(), x.max())
    plt.ylim(y.min(), y.max())
    #plt.title('griddata test (%d points)' % npts)

def plot_mdata(x, y, z, gsize=90, levels=6, x_mask=1,y_mask=1, grid_method='linear', cmap='jet', scatter=False, vmin=None, vmax=None):
    from scipy.interpolate import griddata
    import matplotlib.pyplot as plt
    import numpy as np
    print(type(z))
    if type(z) == np.ma.core.MaskedArray:
        if np.ndim(z) == 1:
            m = ~z.mask
            #masked area
            mx= x[z.mask]
            my= y[z.mask]
            #plot area
            xx = x[m]
            yy = y[m]
            zz = z[m]
        else:
            print ('shape error')
    # define grid.
    #npts = z.shape[0]
    xi = np.linspace(x.min(), x.max(), gsize)
    yi = np.linspace(y.min(), y.max(), gsize)
    
    x_nan_list=[]
    y_nan_list=[]
    for i in range(len(mx)):      
        ii=np.where(abs(xi-mx[i])==np.min(abs(xi-mx[i])))[0][0]
        jj=np.where(abs(yi-my[i])==np.min(abs(yi-my[i])))[0][0]
        x_nan_list.append(ii)
        y_nan_list.append(jj)
    
    # grid the data.
    zi = griddata((xx, yy), zz, (xi[None, :], yi[:, None]), method = grid_method)
    
    for i in range(len(x_nan_list)):
        zi[y_nan_list[i]-y_mask:y_nan_list[i]+y_mask+1,x_nan_list[i]-x_mask:x_nan_list[i]+x_mask+1]=np.NaN #imshow, (y,x) is confirmed by plot#        
    # contour the gridded data, plotting dots at the randomly spaced data points.
    CS = plt.imshow(zi[::-1, :], vmin=vmin, vmax=vmax, extent=(x.min(),x.max(),y.min(),y.max()), cmap=cmap, aspect='auto' )
    plt.colorbar()
    CS = plt.contour(xi, yi, zi, levels, linewidths=0.5, colors='k')
    plt.xlim(x.min(), x.max())
    plt.ylim(y.min(), y.max())


####cut part of colorbar###
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


####rescaled colorbar###
import matplotlib
#matplotlib.use("TkAgg")
import numpy as np
import matplotlib.cm as cm
import matplotlib.pylab as plt
from matplotlib.colors import Normalize
class DivergingNorm(Normalize):
    def __init__(self, vcenter, vmin=None, vmax=None):
        """
        Normalize data with a set center.

        Useful when mapping data with an unequal rates of change around a
        conceptual center, e.g., data that range from -2 to 4, with 0 as
        the midpoint.

        Parameters
        ----------
        vcenter : float
            The data value that defines ``0.5`` in the normalization.
        vmin : float, optional
            The data value that defines ``0.0`` in the normalization.
            Defaults to the min value of the dataset.
        vmax : float, optional
            The data value that defines ``1.0`` in the normalization.
            Defaults to the the max value of the dataset.

        Examples
        --------
        This maps data value -4000 to 0., 0 to 0.5, and +10000 to 1.0; data
        between is linearly interpolated::

            >>> import matplotlib.colors as mcolors
            >>> offset = mcolors.DivergingNorm(vmin=-4000.,
                                               vcenter=0., vmax=10000)
            >>> data = [-4000., -2000., 0., 2500., 5000., 7500., 10000.]
            >>> offset(data)
            array([0., 0.25, 0.5, 0.625, 0.75, 0.875, 1.0])
        """

        self.vcenter = vcenter
        self.vmin = vmin
        self.vmax = vmax
        if vcenter is not None and vmax is not None and vcenter >= vmax:
            raise ValueError('vmin, vcenter, and vmax must be in '
                             'ascending order')
        if vcenter is not None and vmin is not None and vcenter <= vmin:
            raise ValueError('vmin, vcenter, and vmax must be in '
                             'ascending order')


    def autoscale_None(self, A):
        """
        Get vmin and vmax, and then clip at vcenter
        """
        super(type(self), self).autoscale_None( A)
        if self.vmin > self.vcenter:
            self.vmin = self.vcenter
        if self.vmax < self.vcenter:
            self.vmax = self.vcenter


    def __call__(self, value, clip=None):
        """
        Map value to the interval [0, 1]. The clip argument is unused.
        """
        result, is_scalar = self.process_value(value)
        self.autoscale_None(result)  # sets self.vmin, self.vmax if None

        if not self.vmin <= self.vcenter <= self.vmax:
            raise ValueError("vmin, vcenter, vmax must increase monotonically")
        result = np.ma.masked_array(
            np.interp(result, [self.vmin, self.vcenter, self.vmax],
                      [0, 0.5, 1.]), mask=np.ma.getmask(result))
        if is_scalar:
            result = np.atleast_1d(result)[0]
        return result

def cmap1():
    cmap = plt.cm.viridis
    cmap.set_bad('w', 1.0)
    return cmap

def cmap2():
    cmap= plt.cm.RdBu_r
    cmap.set_bad('w', 1.0)
    return cmap

def cmap3():
    cmap=plt.cm.YlOrRd
    cmap.set_bad('w', 1.0)
    return cmap
