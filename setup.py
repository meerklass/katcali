#!/usr/bin/env python2
from setuptools import setup

setup(
    name="katcali",
    author="Jingying Wang",
    author_email="astro.jywang@gmail.com",
    version="0.1.3",
    packages=["katcali"],
    zip_safe=True,
)

from astropy.utils import iers
iers.Conf.iers_auto_url.set('ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all')
