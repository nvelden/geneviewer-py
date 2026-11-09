from importlib.metadata import version

import geneviewer


def test_version_matches_distribution():
    assert geneviewer.__version__ == version("geneviewer-py")
