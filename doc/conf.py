# -*- coding: utf-8 -*-
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.absolute()))

# -- Project information -----------------------------------------------------
project = "BioPAL"
author = "BioPAL team"
copyright = "2021, BioPAL team"

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "numpydoc",
    "nbsphinx",
]
autosummary_generate = True
exclude_patterns = ["_build", "**.ipynb_checkpoints", "legacy"]

# -- HTML options ------------------------------------------------------------
html_logo = "_static/logo.png"
html_static_path = ["_static"]
html_theme = "furo"
