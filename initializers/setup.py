#!/usr/bin/env python3
from setuptools import setup, find_packages

setup(
    name="ips-initializers",
    version="1.0.0",
    install_requires=["ipsframework",
                      "ips_component_utilities"],
    packages=find_packages()
)
