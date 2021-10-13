"""
biopal quick start

Usage:
  biopal-quickstart FOLDER
  
Arguments:
  FOLDER       Path to the folder to initialize
  
Options:
  -h --help          Show this screen
  --version          Show version
  
Details:
    Initialize the folder with default input and configuration xml files:
    Input_File.xml to be edited before BioPAL run
    Configuration_File.xml optionally to be edited before BioPAL run
"""

import sys
import pkgutil
from biopal import __version__
from pathlib import Path


def main():
    from docopt import docopt

    args = docopt(__doc__, version=__version__)

    folder = Path(args["FOLDER"])

    if folder.exists():
        print("Error: provided folder already exists", file=sys.stderr)
        sys.exit(1)

    folder.mkdir(parents=True)
    folder.joinpath("Input_File.xml").write_bytes(pkgutil.get_data("biopal", "_package_data/inputs/Input_File.xml"))
    folder.joinpath("Configuration_File.xml").write_bytes(
        pkgutil.get_data("biopal", "_package_data/conf/Configuration_File.xml")
    )
    print("Folder intialized '{}'".format(folder.absolute()))


if __name__ == "__main__":
    main()
