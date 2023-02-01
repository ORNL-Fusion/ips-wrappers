from setuptools import setup, find_packages

setup(
    name="ips_solps_ftx",
    packages=["ips_solps_ftx","ips_solps_ftx.python_scripts_for_coupling"],
    version="1.0.0",
    install_requires=["ipsframework","ips_xolotlFT","solps-iter"],
    scripts=["ips_solps_ftx/python_scripts_for_coupling/binTRIDYN.py","ips_solps_ftx/python_scripts_for_coupling/plasmaOut2ftxIn.py"],
    find_packages=find_packages()
)