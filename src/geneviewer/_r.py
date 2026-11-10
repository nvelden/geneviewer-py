"""R's rules for values, where the port must follow them exactly.

R compares and names things through ``as.character()``: cluster names,
strand values, ids. These helpers give the same strings for Python values.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd


def is_na(value: Any) -> bool:
    """Return whether a value is R's ``NA``: ``None``, ``NaN``, ``pd.NA`` or ``NaT``.

    Parameters
    ----------
    value : Any
        A single value.

    Returns
    -------
    bool
        True for a missing value.
    """
    if value is None or value is pd.NA or value is pd.NaT:
        return True
    if isinstance(value, (float, np.floating)):
        return math.isnan(value)
    if isinstance(value, np.datetime64):
        return bool(np.isnat(value))
    return False


def is_bool(value: Any) -> bool:
    """Return whether a value is a single logical (R's ``is.logical``, length 1).

    Parameters
    ----------
    value : Any
        The value.

    Returns
    -------
    bool
        True for ``bool`` and numpy's ``bool_``.
    """
    return isinstance(value, (bool, np.bool_))


def is_number(value: Any) -> bool:
    """Return whether a value is a single number (R's ``is.numeric``, length 1).

    Parameters
    ----------
    value : Any
        The value.

    Returns
    -------
    bool
        True for ints and floats, not for logicals.
    """
    return isinstance(value, (int, float, np.integer, np.floating)) and not is_bool(
        value
    )


def as_vector(value: Any) -> list[Any]:
    """Return a value as a list, as R sees a vector.

    Parameters
    ----------
    value : Any
        A scalar, a list or tuple, a numpy array or a pandas Series.

    Returns
    -------
    list
        ``[]`` for ``None``, ``[value]`` for a scalar, else the items.
    """
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return list(value)
    if isinstance(value, (pd.Series, pd.Index, np.ndarray)):
        return value.tolist()
    return [value]


def _number(value: float) -> str:
    # R's as.character() of a double: 15 significant digits, in fixed
    # notation unless the scientific one is shorter
    if math.isinf(value):
        return "Inf" if value > 0 else "-Inf"
    if value == 0:
        return "0"
    sign = "-" if value < 0 else ""
    mantissa, exponent = f"{abs(value):.14e}".split("e")
    digits = mantissa.replace(".", "").rstrip("0") or "0"
    exp = int(exponent)
    n = len(digits)
    if exp >= 0:
        decimals = max(0, n - exp - 1)
        fixed_width = exp + 1 + (decimals + 1 if decimals else 0)
    else:
        decimals = n - exp - 1
        fixed_width = 2 + (-exp - 1) + n
    sci_width = n + (1 if n > 1 else 0) + (4 if abs(exp) < 100 else 5)
    if fixed_width <= sci_width:
        return sign + f"{abs(value):.{decimals}f}"
    sci_mantissa = digits[0] + ("." + digits[1:] if n > 1 else "")
    return f"{sign}{sci_mantissa}e{'-' if exp < 0 else '+'}{abs(exp):02d}"


def as_character(value: Any) -> str | None:
    """Return a value as R's ``as.character()`` writes it.

    Parameters
    ----------
    value : Any
        A single value.

    Returns
    -------
    str or None
        ``None`` for ``NA``; ``"TRUE"``/``"FALSE"`` for logicals; ints as
        they are; floats as R writes doubles (``1`` for ``1.0``, ``1e+05``
        for ``100000.0``); other values with ``str()``.
    """
    if is_na(value):
        return None
    if is_bool(value):
        return "TRUE" if value else "FALSE"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        return _number(float(value))
    return str(value)
