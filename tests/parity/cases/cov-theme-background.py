"""Parity case: the chart of fixture cov-theme-background."""

from geneviewer import GC


def build(data):
    return (
        GC(data, group="class", cluster="cluster", height="250px")
        .labels("name")
        .theme("dark", background="#1d2433")
    )
