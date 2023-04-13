from setuptools import setup, find_packages

setup(
    name="ips_solps_ftx",
    packages=["ips_solps_ftx","ips_solps_ftx.python_scripts_for_coupling"],
    version="1.0.0",
    install_requires=["ipsframework","ips_xolotlFT","solps-iter"],
    scripts=["ips_solps_ftx/python_scripts_for_coupling/binTRIDYN.py","ips_solps_ftx/python_scripts_for_coupling/plasmaOut2ftxIn.py","ips_solps_ftx/python_scripts_for_coupling/write_ftxOut.py", "ips_solps_ftx/python_scripts_for_coupling/writePlasmaOut.py", "ips_solps_ftx/python_scripts_for_coupling/fort44ext.py","ips_solps_ftx/python_scripts_for_coupling/b2fextract.py", "ips_solps_ftx/python_scripts_for_coupling/rizp_extract.py", "ips_solps_ftx/python_scripts_for_coupling/average_SOLPS_input.py"],
    package_data={"ips_solps_ftx.python_scripts_for_coupling": ["ionization_potentials.txt"]},
    find_packages=find_packages()
)
