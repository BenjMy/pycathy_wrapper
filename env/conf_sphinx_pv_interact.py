# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.append(os.path.relpath('..'))
sys.path.insert(0, os.path.abspath('..'))
sys.path.append(os.path.relpath('../pyCATHY/'))
sys.path.insert(0, os.path.abspath('../pyCATHY/'))
#import pyCATHY
import numpy as np

import datetime

import pyvista
from pyvista.plotting.utilities.sphinx_gallery import DynamicScraper
from sphinx_gallery.sorting import FileNameSortKey


import pyvista

# Manage errors
pyvista.set_error_output_file("errors.txt")
# Ensure that offscreen rendering is used for docs generation
pyvista.OFF_SCREEN = True  # Not necessary - simply an insurance policy
# Preferred plotting style for documentation
pyvista.set_plot_theme("document")
pyvista.global_theme.window_size = [1024, 768]
pyvista.global_theme.font.size = 22
pyvista.global_theme.font.label_size = 22
pyvista.global_theme.font.title_size = 22
pyvista.global_theme.return_cpos = False
pyvista.set_jupyter_backend(None)

# Save figures in specified directory
pyvista.FIGURE_PATH = os.path.join(os.path.abspath("./images/"), "auto-generated/")
if not os.path.exists(pyvista.FIGURE_PATH):
    os.makedirs(pyvista.FIGURE_PATH)

# necessary when building the sphinx gallery
pyvista.BUILDING_GALLERY = True
os.environ["PYVISTA_BUILDING_GALLERY"] = "true"



#sys.path.append('.')
#from remove_kernel_metadata import removeK
#removeK()


# Project information
# -----------------------------------------------------------------------------
project = "pyCATHY"
copyright_info = f"2021-{datetime.date.today().year}, The {project} Developers"
#if len(harmonica.__version__.split("+")) > 1 or harmonica.__version__ == "unknown":
#    version = "dev"
#else:
#    version = harmonica.__version__


# -- Project information -----------------------------------------------------

project = 'pyCATHY'
copyright = '2022, B. Mary'
author = 'B. Mary'

# The full version, including alpha/beta/rc tags
release = 'v0.1.1'




# General configuration
# -----------------------------------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.doctest",
    "sphinx.ext.viewcode",
    "sphinx.ext.extlinks",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx_gallery.gen_gallery",
    "sphinx_design",
    "sphinx_copybutton",
    "jupyter_sphinx",
    "pyvista.ext.plot_directive",
    #"myst_nb",
    #"myst_parser",
    #"nbsphinx",
    #"pyvista.ext.viewer_directive",
    "sphinxcontrib.bibtex",
    "sphinx_thebe",
    ]


myst_enable_extensions = [
#     "dollarmath",
#     "amsmath",
#     "deflist",
#     "fieldlist",
#     "html_admonition",
#     "html_image",
#     "colon_fence",
#     "smartquotes",
#     "replacements",
     "linkify",
#     "strikethrough",
#     "substitution",
#     "tasklist",
#     "attrs_inline",
#     "attrs_block",
]

myst_heading_anchors = 3



bibtex_bibfiles = ['content/refs.bib']
bibtex_default_style = 'unsrt'


# Configuration to include links to other project docs when referencing
# functions/classes
intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://docs.scipy.org/doc/numpy/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/reference", None),
    "pandas": ("http://pandas.pydata.org/pandas-docs/stable/", None),
    "matplotlib": ("https://matplotlib.org/", None),
    "pyvista": ("https://docs.pyvista.org", None),
}


# Autosummary pages will be generated by sphinx-autogen instead of sphinx-build
#autosummary_generate = []
#autosummary_generate = True

autosummary_generate = [
                        'pygimli', 
                        'resipy'
]

autosummary_imported_members = True

# Otherwise, the Return parameter list looks different from the Parameters list
napoleon_use_rtype = False
# Otherwise, the Attributes parameter list looks different from the Parameters
# list
napoleon_use_ivar = True


# Always show the source code that generates a plot
plot_include_source = True
plot_formats = ["png"]

# Sphinx project configuration
templates_path = ["_templates"]
exclude_patterns = ["_build", "**.ipynb_checkpoints"]
source_suffix = ['.rst', '.md'] #,'.ipynb']
# The encoding of source files
source_encoding = "utf-8"
master_doc = "index"
pygments_style = "default"
add_function_parentheses = False
jupyter_execute_notebooks_kernel = "pycathy_doc"


# Sphinx-Gallery configuration
# -----------------------------------------------------------------------------
sphinx_gallery_conf = {
    # path to your examples scripts
    #"examples_dirs": "../examples",
    "examples_dirs": ["../examples/SSHydro","../examples/DA"],
    # path where to save gallery generated examples
    #"gallery_dirs": "gallery",
    "gallery_dirs": ["content/SSHydro", "content/DA"],
    #"filename_pattern": "example_.+.ipynb",
    #"filename_pattern": {'py'|'ipynb'}, 
    "example_extensions": {'.py'},
    #"dont_preprocess": [],
    # Remove the "Download all examples" button from the top level gallery
    "download_all_examples": False,
    # Sort gallery example by file name instead of number of lines (default)
    "within_subsection_order": FileNameSortKey,
    # directory where function granular galleries are stored
    "backreferences_dir": "api/generated/backreferences",
    # Modules for which function level galleries are created.  In
    # this case sphinx_gallery and numpy in a tuple of strings.
    "doc_module": "pyCATHY",
    # Insert links to documentation of objects in the examples
    "reference_url": {"pyCATHY": None},
    # Add pyvista to the image scrapers
    "image_scrapers": (DynamicScraper(), "matplotlib"),
    #"image_scrapers": ("pyvista", "matplotlib"),
    'pypandoc': True,
    "first_notebook_cell": "%matplotlib inline\n"
    "from pyvista import set_plot_theme\n"
    'set_plot_theme("document")\n',
}


# Pyvista configurations
# -----------------------------------------------------------------------------
# necessary when building the sphinx gallery
pyvista.BUILDING_GALLERY = True
pyvista.OFF_SCREEN = True

# Optional - set parameters like theme or window size
pyvista.set_plot_theme("document")
pyvista.global_theme.window_size = (1024 * 2, 768 * 2)



# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ["_static"]


# HTML output configuration
# -----------------------------------------------------------------------------
html_title = f'{project} <span class="project-release">{release}</span>'
#html_logo = "_static/sphx_glr_pyCATHY_weilletal_001.png"
#html_favicon = "_static/sphx_glr_pyCATHY_weilletal_thumb.png"
html_last_updated_fmt = "%b %d, %Y"
html_copy_source = True
#html_static_path = ["_static"]
html_extra_path = []
html_show_sourcelink = False
html_show_sphinx = True
html_show_copyright = True
# CSS files are relative to the static path
html_css_files = ["custom.css"]


html_theme_options = {
    "repository_url": "https://github.com/BenjMy/pycathy_wrapper",
    "repository_branch": "master",
    "path_to_docs": "doc",
    "launch_buttons": {
        "binderhub_url": "https://mybinder.org",
        "notebook_interface": "jupyterlab",
        "colab_url": "https://colab.research.google.com",

    },
    "use_edit_page_button": True,
    "use_issues_button": True,
    "use_repository_button": True,
    "use_download_button": True,
    "home_page_in_toc": False,
}
