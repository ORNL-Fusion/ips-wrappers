from setuptools import setup, find_packages
from distutils.command.build import build
from distutils.command.install import install
import os

class build_deps(build):
    def run(self):
        os.system("python3 ../ips-massive-serial-runner/setup.py build")
        build.run(self)

class install_deps(install):
    def run(self):
        os.system("python3 ../ips-massive-serial-runner/setup.py install")
        install.run(self)

setup(
    name="ml_train",
    version="1.0.0",
    install_requires=["massive_serial_runner"],
    packages=find_packages(),
    cmdclass={
        'build'   : build_deps,
        'install' : install_deps
    },)
