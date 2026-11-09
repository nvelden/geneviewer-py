"""Chart data processing: links, introns and normalize.

Ported from R in later steps: ``get_links`` (``R/chart_links.R``) and the
gene normalization of ``R/gene_normalization.R`` in M7.6, ``get_introns``
(``R/transcript.R``) in M7.5. Until then, the chart methods that need them
raise ``NotImplementedError``.
"""

from __future__ import annotations

from typing import Any


def _not_yet(what: str) -> NotImplementedError:
    return NotImplementedError(f"{what} is not available yet in this version")


def get_links(state: Any, group: Any, value1: Any, value2: Any, cluster: Any) -> Any:
    """Return links between genes of neighbouring clusters with the same group.

    Parameters
    ----------
    state : Any
        From ``data.chart_state()``.
    group, value1, value2, cluster : Any
        As for ``GC.links()``.

    Returns
    -------
    pandas.DataFrame
        The links.
    """
    raise _not_yet("links(group=...) without data=")


def get_introns(data: Any) -> Any:
    """Return transcript data with an intron between each pair of exons.

    Parameters
    ----------
    data : pandas.DataFrame
        Exons and UTRs of transcripts.

    Returns
    -------
    pandas.DataFrame
        The data with the introns.
    """
    raise _not_yet("transcript()")


def normalize(spec: Any, **kwargs: Any) -> Any:
    """Return the spec with gene positions normalised.

    Parameters
    ----------
    spec : Any
        The spec of the chart.
    **kwargs
        As for ``GC.normalize()``.

    Returns
    -------
    dict
        The new spec.
    """
    raise _not_yet("normalize()")
