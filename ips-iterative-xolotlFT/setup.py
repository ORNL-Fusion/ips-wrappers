from setuptools import setup, find_packages

setup(
    name="ips_xolotlFT",
    packages=["ips_xolotlFT","ips_xolotlFT.python_scripts_for_coupling"],
    version="1.0.0",
    install_requires=["ipsframework","numpy","h5py","matplotlib","mpi4py"],
    scripts=["ips_xolotlFT/python_scripts_for_coupling/translate_ftridyn_to_xolotl.py",
             "ips_xolotlFT/python_scripts_for_coupling/get_yields.py"],
    package_data={"ips_xolotlFT.python_scripts_for_coupling": ["table1.txt"]},
    find_packages=find_packages()
)
