def plot_data(x, y, z, gsize=30, levels=15, grid_method='linear', scatter=False):
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
    CS = plt.imshow(zi[::-1, :], extent=(x.min(),x.max(),y.min(),y.max()), cmap=plt.get_cmap('jet',255), aspect='auto' )
    plt.colorbar()
    CS = plt.contour(xi, yi, zi, levels, linewidths=0.5, colors='k')
    #plt.contourf(xi, yi, zi, 255, cmap=plt.cm.jet, linewidth=0.0),
    # CS = plt.contourf(xi, yi, zi, 100, cmap=plt.get_cmap('jet',255))
    #plt.colorbar()  # draw colorbar
    # plot data points.
    plt.xlim(x.min(), x.max())
    plt.ylim(y.min(), y.max())
    #plt.title('griddata test (%d points)' % npts)

def plot_mdata(x, y, z, gsize=90, levels=6, x_mask=1,y_mask=1, grid_method='linear', scatter=False, vmin=None, vmax=None):
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
            print 'shape error'
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
    CS = plt.imshow(zi[::-1, :], vmin=vmin, vmax=vmax, extent=(x.min(),x.max(),y.min(),y.max()), cmap=plt.get_cmap('jet',255), aspect='auto' )
    plt.colorbar()
    CS = plt.contour(xi, yi, zi, levels, linewidths=0.5, colors='k')
    plt.xlim(x.min(), x.max())
    plt.ylim(y.min(), y.max())

