import subprocess
import sys
from importlib.metadata import version

import geneviewer


def test_version_matches_distribution():
    assert geneviewer.__version__ == version("geneviewer-py")


def test_shiny_is_optional():
    # geneviewer.shiny loads Shiny; the package itself does not
    code = "import sys, geneviewer; assert 'shiny' not in sys.modules"
    subprocess.run([sys.executable, "-c", code], check=True)
