import os
from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def js_checkout() -> Path:
    """The geneviewer-js checkout named by ``GENEVIEWER_JS``.

    Tests that use it are skipped when the variable is not set.
    """
    path = os.environ.get("GENEVIEWER_JS")
    if not path:
        pytest.skip("GENEVIEWER_JS is not set")
    return Path(path)
