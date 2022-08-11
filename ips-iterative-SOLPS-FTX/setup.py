from setuptools import setup, find_packages

setup(
    name="ips_solps_ftx",
    version="1.0.0",
    install_requires=["ips_siesta"],
    scripts=["ips_solps_ftx/python_scripts_for_coupling/plasmaOut2ftxIn.py","ips_solps_ftx/python_scripts_for_coupling/binTRIDYN.py"],
    packages=find_pakages()
)
