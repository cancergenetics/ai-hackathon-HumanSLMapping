"""Genetic-interaction label definitions (Costanzo et al. 2016 confidence thresholds).

THE PROJECT LABEL is `stringent_negative`: P < 0.05 AND epsilon < -0.12 (both strict inequalities).
Changing anything in this file requires a PI decision recorded in RESEARCH_TREE.md (AGENTS.md rule 1).

The other thresholds are provided for characterisation only (Costanzo 2016 supplement / Costanzo 2010):
  lenient       : P < 0.05
  intermediate  : P < 0.05 and |eps| > 0.08
  stringent     : P < 0.05 and (eps < -0.12 or eps > 0.16)
"""
from __future__ import annotations

import numpy as np

P_THRESHOLD = 0.05
STRINGENT_NEG_EPS = -0.12
STRINGENT_POS_EPS = 0.16
INTERMEDIATE_ABS_EPS = 0.08


def _arr(x) -> np.ndarray:
    return np.asarray(x, dtype=float)


def stringent_negative(eps, p) -> np.ndarray:
    """Project label. NaN eps or p -> False (but callers should not pass NaNs; the NxN file has none)."""
    eps, p = _arr(eps), _arr(p)
    return (p < P_THRESHOLD) & (eps < STRINGENT_NEG_EPS)


def stringent_positive(eps, p) -> np.ndarray:
    eps, p = _arr(eps), _arr(p)
    return (p < P_THRESHOLD) & (eps > STRINGENT_POS_EPS)


def intermediate_negative(eps, p) -> np.ndarray:
    eps, p = _arr(eps), _arr(p)
    return (p < P_THRESHOLD) & (eps < -INTERMEDIATE_ABS_EPS)


def intermediate_positive(eps, p) -> np.ndarray:
    eps, p = _arr(eps), _arr(p)
    return (p < P_THRESHOLD) & (eps > INTERMEDIATE_ABS_EPS)


def lenient_negative(eps, p) -> np.ndarray:
    eps, p = _arr(eps), _arr(p)
    return (p < P_THRESHOLD) & (eps < 0)


def lenient_positive(eps, p) -> np.ndarray:
    eps, p = _arr(eps), _arr(p)
    return (p < P_THRESHOLD) & (eps > 0)
