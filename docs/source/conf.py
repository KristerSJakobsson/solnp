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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

import sphinx_rtd_theme

# -- Project information -----------------------------------------------------

project = 'solnp'
copyright = '2020, Krister S. Jakobsson'
master_doc = 'index'
title = 'Title'
subtitle = 'Thesis'
author = 'Krister S. Jakobsson'
institute = ''
department = ""
source_suffix = '.rst'

# -- General configuration ---------------------------------------------------

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = '%Y-%M-%d, '

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
html_use_smartypants = True

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.ifconfig',
    'sphinx.ext.githubpages',
    "sphinx_rtd_theme",
]

todo_include_todos = True
napoleon_google_docstring = False
napoleon_include_special_with_doc = False


# We do it like this to support multiple sphinx version without having warning.
# Our buildbot consider warning as error.
try:
    from sphinx.ext import imgmath
    extensions.append('sphinx.ext.imgmath')
except ImportError:
    try:
        from sphinx.ext import pngmath
        extensions.append('sphinx.ext.pngmath')
    except ImportError:
        print("Not able to find sphinx module for outputting math formulas in html")

mathjax_config = {
    'extensions': ['tex2jax.js'],
    'jax': ['input/TeX', 'output/HTML-CSS'],
}

latex_packaes = [
    "amsmath",
    "mathtools",
    "amsfonts",
    "amssymb",
    "dsfont"
    # "optidef"
]

default_role = 'math'
pngmath_divpng_args = ['-gamma 1.5','-D 110']
#pngmath_divpng_args = ['-gamma', '1.5', '-D', '110', '-bg', 'Transparent']
imgmath_latex_preamble =  "\n".join([("\\usepackage{%s}" % package) for package in latex_packaes]) +\
                          '\\def\\Z{\\mathbb{Z}}\n'+\
                          '\\def\\R{\\mathbb{R}}\n'+\
                          '\\def\\bX{\\mathbf{X}}\n'+\
                          '\\def\\X{\\mathbf{X}}\n'+\
                          '\\def\\By{\\mathbf{y}}\n'+\
                          '\\def\\Bbeta{\\boldsymbol{\\beta}}\n'+\
                          '\\def\\U{\\mathbf{U}}\n'+\
                          '\\def\\V{\\mathbf{V}}\n'+\
                          '\\def\\V1{\\mathds{1}}\n'+\
                          '\\def\\hU{\\mathbf{\hat{U}}}\n'+\
                          '\\def\\hS{\\mathbf{\hat{\Sigma}}}\n'+\
                          '\\def\\hV{\\mathbf{\hat{V}}}\n'+\
                          '\\def\\E{\\mathbf{E}}\n'+\
                          '\\def\\F{\\mathbf{F}}\n'+\
                          '\\def\\x{\\mathbf{x}}\n'+\
                          '\\def\\h{\\mathbf{h}}\n'+\
                          '\\def\\v{\\mathbf{v}}\n'+\
                          '\\def\\nv{\\mathbf{v^{{\\bf -}}}}\n'+\
                          '\\def\\nh{\\mathbf{h^{{\\bf -}}}}\n'+\
                          '\\def\\s{\\mathbf{s}}\n'+\
                          '\\def\\b{\\mathbf{b}}\n'+\
                          '\\def\\c{\\mathbf{c}}\n'+\
                          '\\def\\W{\\mathbf{W}}\n'+\
                          '\\def\\C{\\mathbf{C}}\n'+\
                          '\\def\\P{\\mathbf{P}}\n'+\
                          '\\def\\T{{\\bf \\mathcal T}}\n'+\
                          '\\def\\B{{\\bf \\mathcal B}}\n'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []
exclude_dirs = ['images', 'scripts', 'sandbox']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# -- Options for LaTeX output ---------------------------------------------

latex_engine = 'xelatex'
latex_show_urls = 'footnote'
latex_additional_files = []
latex_elements = {
    'papersize': 'a4paper',
    'pointsize': '11pt',
    'figure_align': 'htbp',
    'preamble': '',
}

latex_elements['preamble'] += r"""
\usepackage{amsfonts}
\usepackage{parskip}
\usepackage{microtype}
\usepackage{amsmath}
"""

# Latex Titlepage
latex_logo = None
mystyle = 'title'
latex_elements['preamble'] = "\n".join([("\\usepackage{%s}" % package) for package in latex_packaes]) +\
                            '\\def\\Z{\\mathbb{Z}}\n'+\
                            '\\def\\R{\\mathbb{R}}\n'+\
                            '\\def\\bX{\\mathbf{X}}\n'+\
                            '\\def\\X{\\mathbf{X}}\n'+\
                            '\\def\\By{\\mathbf{y}}\n'+\
                            '\\def\\Bbeta{\\boldsymbol{\\beta}}\n'+\
                            '\\def\\bU{\\mathbf{U}}\n'+\
                            '\\def\\bV{\\mathbf{V}}\n'+\
                            '\\def\\V1{\\mathds{1}}\n'+\
                            '\\def\\hU{\\mathbf{\hat{U}}}\n'+\
                            '\\def\\hS{\\mathbf{\hat{\Sigma}}}\n'+\
                            '\\def\\hV{\\mathbf{\hat{V}}}\n'+\
                            '\\def\\E{\\mathbf{E}}\n'+\
                            '\\def\\F{\\mathbf{F}}\n'+\
                            '\\def\\x{\\mathbf{x}}\n'+\
                            '\\def\\h{\\mathbf{h}}\n'+\
                            '\\def\\v{\\mathbf{v}}\n'+\
                            '\\def\\nv{\\mathbf{v^{{\\bf -}}}}\n'+\
                            '\\def\\nh{\\mathbf{h^{{\\bf -}}}}\n'+\
                            '\\def\\s{\\mathbf{s}}\n'+\
                            '\\def\\b{\\mathbf{b}}\n'+\
                            '\\def\\c{\\mathbf{c}}\n'+\
                            '\\def\\W{\\mathbf{W}}\n'+\
                            '\\def\\C{\\mathbf{C}}\n'+\
                            '\\def\\P{\\mathbf{P}}\n'+\
                            '\\def\\T{{\\bf \\mathcal T}}\n'+\
                            '\\def\\B{{\\bf \\mathcal B}}\n'

latex_additional_files += [mystyle + '.sty']


# ('source', 'target', 'title', 'author', 'documentclass')

latex_documents = [
    (master_doc, 'thesis.tex', title, author, 'report'),
]