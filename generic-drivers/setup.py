#!/usr/bin/env python3
from setuptools import setup

setup(
    name="ips_generic_driver",
    version="1.0.0",
    install_requires=["ipsframework",
                      "ips_component_utilities"],
    packages=['ips_generic_driver'],
    package_dir={
        'ips_generic_driver': '.',
    },
)
