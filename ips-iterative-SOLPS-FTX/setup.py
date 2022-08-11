from setuptools import setup, find_packages

setup(
    name="ips_solps_ftx",
    version="1.0.0",
    install_requires=["ipsframework","param_handler"],
    scripts=["ips_solps_ftx/python_scripts_for_coupling/binTRIDYN.py"],
    find_packages=find_packages()
)
