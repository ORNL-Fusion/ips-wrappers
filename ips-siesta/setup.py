from setuptools import setup, find_packages

setup(
    name="ips_siesta",
    version="1.0.0",
    install_requires=["ips_vmec"],
    packages=find_packages(),
)
