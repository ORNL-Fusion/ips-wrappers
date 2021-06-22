from setuptools import setup, find_packages

setup(
    name="ips_cariddi",
    version="1.0.0",
    install_requires=["ips_v3fit"],
    packages=find_packages(),
    scripts=["cariddi_pre.py","cariddi_bin.py"]
)
