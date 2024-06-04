from setuptools import setup, find_packages

setup(
    name="ips_solps_ftx",
    packages=["ips_solps_ftx","ips_solps_ftx.python_scripts_for_coupling"],
    version="1.0.0",
    install_requires=["ipsframework","ips_xolotlFT","solps-iter"],
    scripts=["ips_solps_ftx/python_scripts_for_coupling/binTRIDYN.py","ips_solps_ftx/python_scripts_for_coupling/plasmaOut2ftxIn_call.py","ips_solps_ftx/python_scripts_for_coupling/write_ftxOut_call.py", "ips_solps_ftx/python_scripts_for_coupling/SOLPS_outputs_for_FTX_call.py", "ips_solps_ftx/python_scripts_for_coupling/fort44ext.py","ips_solps_ftx/python_scripts_for_coupling/b2fextract.py", "ips_solps_ftx/python_scripts_for_coupling/rizp_extract.py", "ips_solps_ftx/python_scripts_for_coupling/average_SOLPS_input.py", "ips_solps_ftx/python_scripts_for_coupling/updateSOLPSinput.py", "ips_solps_ftx/python_scripts_for_coupling/read_impact_angle.py", "ips_solps_ftx/python_scripts_for_coupling/read_impact_energy.py"], #added "_call" # launch-vs-call scripts April 2024
    package_data={"ips_solps_ftx.python_scripts_for_coupling": ["ionization_potentials.txt"]},
    find_packages=find_packages()
)
