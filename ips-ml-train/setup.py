from setuptools import setup, find_packages
from distutils.command.build import build
from distutils.command.install import install
import os

class build_deps(build):
    def run(self):
        pwd = os.getcwd()
        os.chdir("../ips-massive-serial-runner")
        os.system("python3 setup.py build")
        os.chdir(pwd)
        os.chdir("../ips-massive-parallel-runner")
        os.system("python3 setup.py build")
        os.chdir(pwd)
        build.run(self)

class install_deps(install):
    def run(self):
        pwd = os.getcwd()
        os.chdir("../ips-massive-serial-runner")
        os.system("python3 setup.py install")
        os.chdir(pwd)
        os.chdir("../ips-massive-parallel-runner")
        os.system("python3 setup.py install")
        os.chdir(pwd)
        install.run(self)

setup(
    name="ml_train",
    version="1.0.0",
    install_requires=["massive_serial_runner","massive_parallel_runner"],
    packages=find_packages(),
    cmdclass={
        'build'   : build_deps,
        'install' : install_deps
    },)
