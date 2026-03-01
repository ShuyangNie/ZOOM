import os
import sys
import sphinx_rtd_theme

sys.path.insert(0, os.path.abspath("."))
sys.path.insert(0, os.path.abspath("../"))
sys.path.insert(0, os.path.abspath("../../"))
sys.path.insert(0, os.path.abspath("../zoom/"))

project = 'ZOOM'
copyright = '2025, Shuyang Nie'
author = 'Shuyang Nie'
release = '0.0.1'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx.ext.viewcode',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinxarg.ext'
]

# Napoleon settings
#napoleon_google_docstring = True
#napoleon_numpy_docstring = True
#napoleon_include_init_with_doc = False
#napoleon_include_private_with_doc = False
#napoleon_include_special_with_doc = True
#napoleon_use_admonition_for_examples = False
#napoleon_use_admonition_for_notes = False
#napoleon_use_admonition_for_references = False
#napoleon_use_ivar = False
#napoleon_use_param = True
#napoleon_use_rtype = True
#napoleon_preprocess_types = True

# Generate the API documentation when building
autosummary_generate = True
autodoc_default_options = {'members': True, 'inherited-members': True}
numpydoc_show_class_members = False
autoclass_content = "class"
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme'
html_logo = 'logo.png'

html_theme_options = {
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

html_static_path = ['_static']

pygments_style = 'sphinx'

html_sidebars = {
    '**': [
        'relations.html',  # needs 'show_related': True theme option to display
        'searchbox.html',
    ]
}

numfig = True
numfig_secnum_depth = 1

autodoc_mock_imports = [
    "numpy", "pandas", "sklearn", "anndata",
    "scanpy", "statsmodels", "nibabel", "scipy",
    "abagen", "neuromaps", "joblib", "tqdm",
    "gseapy"
]
