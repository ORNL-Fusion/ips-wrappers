from setuptools import setup, find_packages
from distutils.command.build import build
from distutils.command.install import install
import os

class build_deps(build):
    def run(self):
        os.system("python3 ../utilities/setup.py build")
        build.run(self)

class install_deps(install):
    def run(self):
        os.system("python3 ../utilities/setup.py install")
        install.run(self)

setup(
    name="ips_vmec",
    version="1.0.0",
    install_requires=["ips_component_utilities"],
    packages=find_packages(),
    cmdclass={
        'build'   : build_deps,
        'install' : install_deps
    },
)
