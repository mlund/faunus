# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

import datetime

def getCurrentYear():
    ''' return current year, 20XX '''
    currentDateTime = datetime.datetime.now()
    date = currentDateTime.date()
    return date.strftime("%Y")

project = 'Faunus'
copyright = '{}, Mikael Lund'.format(getCurrentYear())
author = 'Mikael Lund'
source_suffix = ['.rst', '.md']
master_doc = 'index'
extensions = [
        'recommonmark',
        'sphinx_markdown_tables',
        'sphinx.ext.autosectionlabel',
        'sphinx_rtd_theme']
html_theme = "sphinx_rtd_theme"
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
autosectionlabel_prefix_document = True
