"""Compare a spec built in Python with the one R exported.

The rules of geneviewer-js ``test/compile.mjs``, with one addition:

- one-element arrays are compared as their element (R writes a vector of
  one as a scalar), except under ``clusters`` and ``selection``, which are
  always arrays;
- numbers are equal within 1e-9 relative difference (R writes 16
  significant digits); booleans are not numbers;
- a missing key equals ``null``; ``null`` in an array (R's NA) must match;
- ``[]`` equals ``{}``; key order does not matter, array order does.
"""

from __future__ import annotations

from typing import Any

ALWAYS_ARRAYS = ("clusters", "selection")


def _is_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def _numbers_equal(a: float, b: float) -> bool:
    if a == b:
        return True
    scale = max(abs(a), abs(b))
    if scale < 1e-12:
        return abs(a - b) < 1e-12
    return abs(a - b) / scale <= 1e-9


def _is_empty(x: Any) -> bool:
    return isinstance(x, (list, dict)) and len(x) == 0


def differences(
    expected: Any, actual: Any, path: tuple = (), out: list | None = None
) -> list[tuple[str, Any, Any]]:
    """Return where two specs differ, as (path, expected, actual).

    Parameters
    ----------
    expected : Any
        The spec R exported (``spec.json``).
    actual : Any
        The spec Python built (``to_spec()``).
    path : tuple
        The key path so far.
    out : list, optional
        The differences so far.

    Returns
    -------
    list of tuple
        The differences.
    """
    out = [] if out is None else out
    key = path[-1] if path else None
    if key not in ALWAYS_ARRAYS:
        if isinstance(expected, list) and len(expected) == 1:
            expected = expected[0]
        if isinstance(actual, list) and len(actual) == 1:
            actual = actual[0]
    if expected is None and actual is None:
        return out
    if _is_empty(expected) and _is_empty(actual):
        return out
    where = ".".join(str(p) for p in path)
    if _is_number(expected) and _is_number(actual):
        if not _numbers_equal(expected, actual):
            out.append((where, expected, actual))
        return out
    if isinstance(expected, list) and isinstance(actual, list):
        if len(expected) != len(actual):
            out.append((where, expected, actual))
            return out
        for i, (e, a) in enumerate(zip(expected, actual, strict=True)):
            # In arrays, null is a value (NA)
            if (e is None) != (a is None):
                out.append((f"{where}.{i}", e, a))
            else:
                differences(e, a, (*path, i), out)
        return out
    if isinstance(expected, dict) and isinstance(actual, dict):
        for name in dict.fromkeys([*expected, *actual]):
            differences(expected.get(name), actual.get(name), (*path, name), out)
        return out
    if type(expected) is not type(actual) or expected != actual:
        out.append((where, expected, actual))
    return out
