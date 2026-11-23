# Sphinx configuration for the FLUME documentation.
#
#     sphinx-build -W --keep-going -b html Docs/sphinx_doc Docs/sphinx_doc/_build/html

project = "FLUME"
author = "FLUME and ERF developers"
copyright = " "
version = ""
release = ""
language = "en"

extensions = [
    "sphinx.ext.mathjax",
    "sphinx_copybutton",
]

master_doc = "index"
source_suffix = ".rst"
exclude_patterns = ["_build"]

# Number figures and tables, so they can be cited with :numref:`label`.
numfig = True
numfig_format = {"figure": "%s", "table": "%s", "code-block": "%s"}

html_theme = "sphinx_rtd_theme"
html_title = "FLUME documentation"
