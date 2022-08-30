import nox
from pathlib import Path

"""
Run from command line functionalities in isolated python virtual environment

Usage:
    nox -r -fb conda -s FUNCTION
see also nox --help

Details:
    available functions are:
    build_sdist: sdist package for PyPI 
    build_wheel: wheel package for PyPI
    build_doc: generate sphinx documentation from docstrings and rst files
"""


@nox.session()
def build_sdist(session: nox.Session):
    session.install("build")
    session.run(
        "python", "-m", "build", "--sdist", silent=True,
    )


@nox.session()
def build_wheel(session: nox.Session):
    session.install("build")
    session.run(
        "python", "-m", "build", "--wheel", silent=True,
    )


@nox.session(python="3.7")
def build_doc(session: nox.Session):
    session.conda_install("GDAL=3.5", channel="conda-forge")
    session.install("-e", ".")
    session.install("-r", "doc/requirements.txt")
    rst_files = [
        "chains/agb.rst",
        "chains/fd.rst",
        "chains/fh.rst",
        "chains/tomo_fh.rst",
        "utilities/ground_cancellation.rst",
        "utilities/raster.rst",
        "utilities/plot.rst",
    ]
    for rst_file in rst_files:
        session.run(
            "python", "-m", "sphinx.ext.autosummary.generate", str(Path("doc/api").joinpath(rst_file)),
        )
    session.run(
        "python", "-m", "sphinx", "-b=html", "-a", "doc", "doc/_build",
    )
