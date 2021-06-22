from setuptools import setup, find_packages
from distutils.command.build import build
from distutils.command.install import install
import os

class build_deps(build):
    def run(self):
        pwd = os.getcwd()
        os.chdir("../ips-v3fit")
        os.system("python3 setup.py build")
        build.run(self)

class install_deps(install):
    def run(self):
        pwd = os.getcwd()
        os.chdir("../ips-v3fit")
        os.system("python3 setup.py install")
        os.chdir(pwd)
        install.run(self)

setup(
    name="ips_cariddi",
    version="1.0.0",
    install_requires=["ips_v3fit"],
    packages=find_packages(),
    scripts=["cariddi_pre.py","cariddi_bin.py"],
    cmdclass={
        'build'   : build_deps,
        'install' : install_deps
    },
)
