# Sphinx configuration for the FLUME documentation.
#
#     sphinx-build -W --keep-going -b html Docs/sphinx_doc Docs/sphinx_doc/_build/html

project = "FLUME"
author = "FLUME contributors"
copyright = "2026, FLUME contributors"

extensions = ["sphinx.ext.mathjax"]
templates_path = []
exclude_patterns = ["_build"]

html_theme = "sphinx_rtd_theme"
html_static_path = []

master_doc = "index"

# Number figures and tables, so they can be cited with :numref:`label`.
numfig = True
