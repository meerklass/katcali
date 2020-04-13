def plot_data(x, y, z, gsize=30, scatter=False):
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
    zi = griddata((x, y), z, (xi[None, :], yi[:, None]))
    #zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='cubic')
    # print(zi.shape)
    # contour the gridded data, plotting dots at the randomly spaced data points.
    #CS = plt.contour(xi, yi, zi, 15, linewidths=0.5, colors='k')
    #CS = plt.contourf(xi, yi, zi, 255, cmap=plt.cm.jet, linewidth=0.0)
    CS = plt.imshow(zi[::-1, :], extent=(x.min(),x.max(),y.min(),y.max()), cmap=plt.get_cmap('jet',255), aspect='auto' )
    plt.colorbar()
    CS = plt.contour(xi, yi, zi, 15, linewidths=0.5, colors='k')
    #plt.contourf(xi, yi, zi, 255, cmap=plt.cm.jet, linewidth=0.0),
    # CS = plt.contourf(xi, yi, zi, 100, cmap=plt.get_cmap('jet',255))
    #plt.colorbar()  # draw colorbar
    # plot data points.
    plt.xlim(x.min(), x.max())
    plt.ylim(y.min(), y.max())
    #plt.title('griddata test (%d points)' % npts)
