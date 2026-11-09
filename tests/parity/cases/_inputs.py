"""Inputs of parity cases whose chart changes the data it is given.

The spec's data of these fixtures is the result, not the input. Until the
datasets are part of the package (M7.8), the input is rebuilt from the result
exactly: the rows in their first order by ``rowID``, the old positions from
``original_start`` and ``original_end``, without the columns R added.
"""

import pandas as pd

TRANSCRIPTS = pd.DataFrame(
    {
        "transcript": ["transcript1"] * 4 + ["transcript2"] * 3,
        "type": ["5_utr", "exon", "exon", "3_utr", "5_utr", "exon", "3_utr"],
        "start": [1, 101, 201, 301, 1, 101, 301],
        "end": [50, 150, 250, 350, 50, 150, 350],
        "strand": ["forward"] * 7,
    }
)


def before_normalize(data: pd.DataFrame) -> pd.DataFrame:
    """The data before normalize(), when it normalised every cluster."""
    data = data.sort_values("rowID", kind="stable").reset_index(drop=True)
    data["start"] = data["original_start"]
    data["end"] = data["original_end"]
    return data.drop(columns=["rowID", "strand", "original_start", "original_end"])


def before_transcript(data: pd.DataFrame, transcript: str) -> pd.DataFrame:
    """The data before transcript(): without the introns R added."""
    data = data[data["type"] != "intron"]
    data = data.sort_values("rowID", kind="stable").reset_index(drop=True)
    columns = ["rowID"] if transcript == "transcript" else ["rowID", "transcript"]
    return data.drop(columns=columns)
