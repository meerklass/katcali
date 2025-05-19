#!/usr/bin/env python3
from setuptools import setup, find_packages

setup(
    name="katcali",
    author="Jingying Wang",
    author_email="astro.jywang@gmail.com",
    version="3.0.9",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'katcali': ['config.yaml'],
    },
    zip_safe=True,
)

