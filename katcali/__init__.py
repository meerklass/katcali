__version__='2.0.0'

from astropy.utils import iers
iers.Conf.iers_auto_url.set('ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all')
iers.Conf.iers_auto_url.set('http://toshi.nofs.navy.mil/ser7/finals2000A.all')

