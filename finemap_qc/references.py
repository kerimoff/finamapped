"""Strict helpers shared by reference preparation scripts."""
import numpy as np
import pandas as pd

from .io import normalize_variant_ids


def high_pip_variants(df, threshold, variant_col="variant", pip_col="pip"):
    if not np.isfinite(threshold) or not 0 <= threshold <= 1:
        raise ValueError("PIP threshold must be finite and in [0, 1]")
    pip = pd.to_numeric(df[pip_col], errors="raise")
    if not (np.isfinite(pip) & pip.between(0, 1)).all():
        raise ValueError("reference PIPs must be finite and in [0, 1]")
    parsed = normalize_variant_ids(df[variant_col])["variant"]
    if parsed.isna().any():
        raise ValueError("reference has invalid coordinate variant IDs (rsIDs are not supported)")
    return set(parsed[pip > threshold])
