# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'isweep'
copyright = '2025, Seth D. Temple'
author = 'Seth D. Temple'

release = '1.0'

# -- General configuration

import os
import sys
sys.path.insert(0, os.path.abspath('../../src'))

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
]

napoleon_google_docstring = False     # Turn on numpydoc strings
napoleon_numpy_docstring = True     # Turn on numpydoc strings

# modules = ['isweep']
autodoc_mock_imports = ['isweep']

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'

numfig = True
autosectionlabel_prefix_document = True

highlight_language = 'python'
