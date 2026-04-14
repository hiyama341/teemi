from pathlib import Path
import sys


DOCS_DIR = Path(__file__).resolve().parent
REPO_ROOT = DOCS_DIR.parent

sys.path.insert(0, str(REPO_ROOT))

project = "teemi"
copyright = "2024, Technical University of Denmark (DTU)"
author = "Lucas Levassor"

try:
    from teemi import __version__ as release
except Exception:
    release = "0.0.0"

version = release

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

autosummary_generate = True

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}

autodoc_mock_imports = [
    "benchlingapi",
    "Bio",
    "crispr_cas",
    "dnachisel",
    "dotenv",
    "h2o",
    "intermine",
    "matplotlib",
    "matplotlib.patches",
    "matplotlib.pyplot",
    "numpy",
    "openpyxl",
    "pandas",
    "pydna",
    "pylab",
    "requests",
    "scipy",
    "scipy.stats",
    "seaborn",
]

templates_path = ["_templates"]
html_static_path = ["_static"]
html_css_files = ["custom.css"]

exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "LICENSE",
    "MANIFEST.in",
    "Makefile",
    "api.rst",
    "generated",
    "index.rst",
    "installation.rst",
    "modules.rst",
    "readme.rst",
    "teemi.*.rst",
    "teemi.rst",
    "usage.rst",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

master_doc = "index"
language = "en"
pygments_style = "sphinx"

myst_enable_extensions = [
    "colon_fence",
]

html_theme = "pydata_sphinx_theme"
html_logo = str(REPO_ROOT / "pictures" / "teemi_logo.svg")
html_title = "teemi documentation"
html_theme_options = {
    "github_url": "https://github.com/hiyama341/teemi",
    "navigation_depth": 3,
    "show_toc_level": 2,
    "logo": {
        "text": "",
    },
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/hiyama341/teemi",
            "icon": "fa-brands fa-github",
        },
    ],
}
