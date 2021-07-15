#!/usr/bin/env python3
from setuptools import setup

setup(
    name="ips_model_epa",
    version="1.0.0",
    install_requires=["ipsframework",
                      "ips_component_utilities"],
    packages=['ips_model_epa'],
    package_dir={
        'ips_model_epa': '.',
    },
)
