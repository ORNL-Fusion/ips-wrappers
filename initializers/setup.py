#!/usr/bin/env python3
from setuptools import setup

setup(
    name="ips_initializers",
    version="1.0.0",
    install_requires=["ipsframework",
                      "ips_component_utilities"],
    packages=['ips_initializers', 'ips_initializers.generic_ps_init'],
    package_dir={
        'ips_initializers': '.',
    },
)
