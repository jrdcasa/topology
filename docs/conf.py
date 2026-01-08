# Configuration file for the Sphinx documentation builder.

import os
import sys
from datetime import datetime

# -----------------------------------------------------
# Path setup
# -----------------------------------------------------
# Add project root directory to sys.path
sys.path.insert(0, os.path.abspath(".."))


# -----------------------------------------------------
# Project information
# -----------------------------------------------------
project = "topology"
author = "Your Name"
current_year = datetime.now().year
copyright = f"{current_year}, {author}"

# If your project has a __version__
try:
    import topology
    release = topology.__version__
except Exception:
    release = "0.1"


# -----------------------------------------------------
# General configuration
# -----------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

autosummary_generate = True
autodoc_member_order = "bysource"

napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_special_with_doc = True

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -----------------------------------------------------
# HTML output
# -----------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 4,
    "style_external_links": True,
}

# Optional: logo if you have one
# html_logo = "_static/logo.png"

# Optional: add a sidebar
# html_sidebars = {
#     '**': [
#         'globaltoc.html',
#         'sourcelink.html',
#         'searchbox.html'
#     ]
# }

