from setuptools import setup, find_packages
from distutils.command.build import build
from distutils.command.install import install
import os

class build_deps(build):
    def run(self):
        pwd = os.getcwd()
        os.chdir("../utilities")
        os.system("python3 setup.py build")
        os.chdir(pwd)
        build.run(self)

class install_deps(install):
    def run(self):
        pwd = os.getcwd()
        os.chdir("../utilities")
        os.system("python3 setup.py install")
        os.chdir(pwd)
        install.run(self)

setup(
    name="massive_serial_runner",
    version="1.0.0",
    install_requires=["ips_component_utilities"],
    packages=find_packages(),
    scripts=["launch"],
    cmdclass={
        'build'   : build_deps,
        'install' : install_deps
    },
)
