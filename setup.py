import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="biopal",
    version="1.0.0",
    #author="the BioPAL team",
    #author_email="biopal@esa.int",
    description="BIOMASS Product Algorithm Laboratory",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/BioPAL/BioPAL",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)