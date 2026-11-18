"""Pairwise colour separations of every map layer under simulated colour-vision deficiency.

Separations are OKLab Euclidean differences multiplied by 100 (the unit used by
tests/test_rendered_colors.py). Simulations follow Machado, Oliveira and Fernandes
(2009) at severity 1.0, applied in linear RGB. Writes a CSV and a Markdown table.

Usage: python paper/cvd_separations.py OUTPUT_DIR
"""
from itertools import combinations
from pathlib import Path
import sys

import numpy as np
import pandas as pd

MACHADO = {
    "normal": np.eye(3),
    "protanopia": np.array([[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216],
                            [-0.003882, -0.048116, 1.051998]]),
    "deuteranopia": np.array([[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413],
                              [-0.011820, 0.042940, 0.968881]]),
    "tritanopia": np.array([[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602],
                            [0.004733, 0.691367, 0.303900]]),
}

# Colours as used by chemgridmap.plotting and web/app.js (0.6.0).
LAYERS = {
    "activity (supported fill)": {"active": "#2a9d8f", "medium": "#e6a23c", "inactive": "#c44e52",
                                  "not measured": "#ffffff"},
    "support (uncertain tint)": {"active": "#b8e0db", "medium": "#f6dfb8", "inactive": "#ebc0c1"},
    "support (supported vs uncertain, same class)": None,  # computed below
    "reliability bins": {"0-0.2": "#cde1f5", "0.2-0.4": "#a0bdda", "0.4-0.6": "#7599bf",
                         "0.6-0.8": "#4a77a4", "0.8-1": "#1d5689"},
    "neighbour evidence (outline)": {"active": "#2a9d8f", "medium": "#e6a23c", "inactive": "#c44e52",
                                     "mixed": "#7a5195", "insufficient": "#8a99a6", "none": "#b8c2cc"},
}
FADED_REFERENCE = {"active": "#b8e0db", "medium": "#f6dfb8", "inactive": "#ebc0c1"}  # before greyscale
DARK_RING = "#1f2933"


def _hex(value):
    return np.array([int(value[i:i + 2], 16) / 255 for i in (1, 3, 5)])


def _linear(rgb):
    return np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)


def _gamma(rgb):
    rgb = np.clip(rgb, 0, 1)
    return np.where(rgb <= 0.0031308, 12.92 * rgb, 1.055 * rgb ** (1 / 2.4) - 0.055)


def simulate(color, condition):
    return _gamma(MACHADO[condition] @ _linear(_hex(color)))


def oklab(rgb):
    r, g, b = _linear(np.asarray(rgb, dtype=float))
    lms = np.cbrt([0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b,
                   0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b,
                   0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b])
    l, m, s = lms
    return np.array([0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
                     1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
                     0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s])


def delta(first, second, condition):
    return float(np.linalg.norm(oklab(simulate(first, condition)) - oklab(simulate(second, condition))) * 100)


def greyscale(color):
    """CSS filter grayscale(1): luminance-weighted grey (Rec. 709 weights on sRGB values)."""
    r, g, b = _hex(color)
    y = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "#%02x%02x%02x" % ((round(y * 255),) * 3)


def pairs():
    rows = []
    for layer, palette in LAYERS.items():
        if palette is None:
            items = [(f"{k} supported", f"{k} uncertain", LAYERS["activity (supported fill)"][k], FADED_REFERENCE[k])
                     for k in ("active", "medium", "inactive")]
        else:
            items = [(a, b, palette[a], palette[b]) for a, b in combinations(palette, 2)]
        for a, b, ca, cb in items:
            rows.append({"layer": layer, "first": a, "second": b, "first_hex": ca, "second_hex": cb,
                         **{c: delta(ca, cb, c) for c in MACHADO}})
    # Neighbour-evidence layer: coloured inner frame against the new compound's own fill, the
    # dark outer ring against the frame, and greyed reference cells against coloured frames.
    fills = LAYERS["activity (supported fill)"]
    for frame in ("active", "medium", "inactive"):
        for fill in ("active", "medium", "inactive", "not measured"):
            ca, cb = LAYERS["neighbour evidence (outline)"][frame], fills[fill]
            rows.append({"layer": "neighbour evidence (frame vs own fill)", "first": f"frame {frame}",
                         "second": f"fill {fill}", "first_hex": ca, "second_hex": cb,
                         **{c: delta(ca, cb, c) for c in MACHADO}})
        rows.append({"layer": "neighbour evidence (dark ring vs frame)", "first": "dark ring", "second": f"frame {frame}",
                     "first_hex": DARK_RING, "second_hex": LAYERS["neighbour evidence (outline)"][frame],
                     **{c: delta(DARK_RING, LAYERS["neighbour evidence (outline)"][frame], c) for c in MACHADO}})
    for frame, colour in LAYERS["neighbour evidence (outline)"].items():
        for ref, tint in FADED_REFERENCE.items():
            grey = greyscale(tint)
            rows.append({"layer": "neighbour evidence (frame vs greyed reference)", "first": f"frame {frame}",
                         "second": f"reference {ref} (grey)", "first_hex": colour, "second_hex": grey,
                         **{c: delta(colour, grey, c) for c in MACHADO}})
    return pd.DataFrame(rows)


def main(out):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    table = pairs()
    table.to_csv(out / "cvd_separations.csv", index=False, float_format="%.1f")
    summary = table.groupby("layer")[list(MACHADO)].min().round(1)
    lines = ["| Layer | Pair | Normal | Protanopia | Deuteranopia | Tritanopia |", "|---|---|---|---|---|---|"]
    for _, r in table.iterrows():
        lines.append(f"| {r.layer} | {r['first']} – {r['second']} | {r.normal:.1f} | {r.protanopia:.1f} | "
                     f"{r.deuteranopia:.1f} | {r.tritanopia:.1f} |")
    (out / "cvd_separations.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(summary.to_string())


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "cvd_out")
