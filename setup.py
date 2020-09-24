import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="biopal",
    version="1.0.0",
    #author="Example Author",
    #author_email="author@example.com",
    description="BiomassL2 AGB Prototype Processor",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/BioPAL/BioPAL",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)