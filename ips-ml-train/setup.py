from setuptools import setup, find_packages

setup(
    name="ml_train",
    version="1.0.0",
    install_requires=["massive_serial_runner"],
    packages=find_packages(),
)
