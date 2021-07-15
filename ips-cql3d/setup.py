#!/usr/bin/env python3
from setuptools import setup

setup(
    name="ips_cql3d",
    version="1.0.0",
    install_requires=["ipsframework",
                      "ips_component_utilities"],
    packages=['ips_cql3d'],
    package_dir={
        'ips_cql3d': '.',
    },
)
