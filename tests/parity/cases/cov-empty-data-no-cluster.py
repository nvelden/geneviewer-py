"""Parity case: the chart of fixture cov-empty-data-no-cluster."""

import warnings

from geneviewer import GC


def build(data):
    # The one row has no start: the chart has no rows left
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return GC(data, height="150px").labels("name")
