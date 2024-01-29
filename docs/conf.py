"""Sphinx configuration."""
import os
import sys


sys.path.insert(0, os.path.abspath("../src"))

project = "Proksee Batch"
author = "Lael D. Barlow"
copyright = "2024, Lael D. Barlow"
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_click",
    "myst_parser",
]
autodoc_typehints = "description"
html_theme = "furo"
