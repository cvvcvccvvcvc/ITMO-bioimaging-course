from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "biological_interpretation"
ASSETS = OUT / "assets"
ASSETS.mkdir(parents=True, exist_ok=True)

CHANNEL_DIR = ROOT / "data" / "channels"
DOWNSAMPLE_NPZ = ROOT / "marker_expectation" / "downsampled_mean_f16.npz"
CELL_QC = (
    ROOT
    / "outputs"
    / "manual-20260611-bioimage"
    / "presentations"
    / "codex-multiplex-analysis-v2"
    / "assets"
    / "single_cell_quantification_54_full_qc.csv"
)
CLUSTERS = (
    ROOT
    / "outputs"
    / "manual-20260611-bioimage"
    / "presentations"
    / "codex-multiplex-analysis-v2"
    / "assets"
    / "single_cell_clusters_full.csv"
)

PIXEL_UM = 0.377
DS = 16


MARKER_GROUPS = {
    "nuclear": ["Hoechst1", "DRAQ5", "Ki67"],
    "epithelial": [
        "Cytokeratin",
        "CK7",
        "CDX2",
        "SOX9",
        "CD49f",
        "CD49a",
        "CD44",
        "CD66",
        "MUC1",
        "MUC2",
        "ITLN1",
        "OLFM4",
        "CHGA",
        "Synaptophysin",
        "aDefensin5",
    ],
    "immune": [
        "CD45",
        "CD45RO",
        "CD3",
        "CD4",
        "CD8",
        "CD7",
        "CD25",
        "CD127",
        "CD69",
        "CD161",
        "CD57",
        "NKG2G",
        "CD56",
        "CD19",
        "CD21",
        "CD38",
        "CD138",
        "BCL2",
        "CD11c",
        "HLADR",
        "CD123",
        "CD15",
        "CD16",
        "CD68",
        "CD163",
        "CD206",
        "CD117",
    ],
    "stromal_vascular": [
        "Vimentin",
        "aSMA",
        "FAP",
        "CollagenIV",
        "CD31",
        "CD34",
        "CD36",
        "CD90",
        "Podoplanin",
    ],
}

MARKER_INTERPRETATION = {
    "Hoechst1": ("DNA dye", "all nuclei; use for tissue geometry and segmentation"),
    "MUC1": ("apical epithelial mucin", "luminal/apical epithelial surface"),
    "CD25": ("IL2RA", "activated/T-reg immune cells; sparse if specific"),
    "CDX2": ("intestinal nuclear transcription factor", "intestinal epithelial nuclei"),
    "Synaptophysin": ("synaptic vesicle protein", "neuroendocrine cells and neural elements"),
    "CD57": ("HNK-1 epitope", "NK/senescent T/neural or rare epithelial contexts"),
    "NKG2G": ("NKG2D/CD314 panel marker", "NK/cytotoxic immune context"),
    "Vimentin": ("intermediate filament", "mesenchymal stroma, endothelium, some immune cells"),
    "CD4": ("T-helper marker", "CD4 T cells; also some APC contexts"),
    "CD19": ("B-cell receptor complex", "B cells/lymphoid aggregates"),
    "CD7": ("T/NK marker", "T/NK cells with CD3 or CD56/NKG2D support"),
    "CD11c": ("integrin ITGAX", "dendritic cells/myeloid APC"),
    "CD161": ("KLRB1", "MAIT/NK/T-cell subsets"),
    "CD15": ("Lewis X/SSEA-1", "granulocytes/neutrophil-rich foci"),
    "CD34": ("endothelial/progenitor marker", "vascular endothelium and some stromal cells"),
    "CD16": ("Fc gamma receptor", "myeloid/NK cells"),
    "ITLN1": ("intelectin-1", "goblet/secretory epithelial signal"),
    "HLADR": ("MHC class II", "antigen-presenting cells; possible activated epithelium"),
    "CD123": ("IL3RA", "plasmacytoid dendritic cells; sparse"),
    "CD66": ("CEACAM family", "mature colonocytes/epithelial and granulocyte contexts"),
    "CD3": ("T-cell receptor complex", "T lymphocytes"),
    "CD45RO": ("PTPRC memory isoform", "memory T cells / activated lymphocytes"),
    "CD38": ("activation/plasma-cell marker", "plasma cells and activated immune cells"),
    "CD90": ("THY1", "fibroblast/stromal and some T-cell contexts"),
    "CK7": ("keratin 7", "patchy epithelial keratin; not dominant in normal colon"),
    "aSMA": ("ACTA2", "smooth muscle, pericytes, myofibroblasts"),
    "CD117": ("KIT", "mast cells and interstitial cells of Cajal"),
    "CD127": ("IL7R", "T-cell/ILC survival receptor"),
    "MUC2": ("gel-forming mucin", "goblet cells/mucus in colon epithelium"),
    "CHGA": ("chromogranin A", "enteroendocrine cells"),
    "FAP": ("fibroblast activation protein", "activated fibroblasts/stromal bands"),
    "CollagenIV": ("basement-membrane ECM", "epithelial basement membrane and vascular ECM"),
    "CD21": ("CR2", "B cells and follicular dendritic-cell networks"),
    "BCL2": ("anti-apoptotic protein", "lymphoid survival; can mark crypt/immune contexts"),
    "CD31": ("PECAM1", "blood vascular endothelium"),
    "Ki67": ("cell-cycle marker", "proliferating nuclei, crypt transit-amplifying cells"),
    "SOX9": ("epithelial progenitor TF", "crypt/progenitor and glandular epithelial nuclei"),
    "CD8": ("cytotoxic T-cell marker", "CD8 T cells / intraepithelial lymphocytes"),
    "CD36": ("scavenger/fatty-acid receptor", "endothelium, macrophage, stromal contexts"),
    "CD138": ("syndecan-1", "plasma cells; also epithelial surface contexts"),
    "CD69": ("early activation marker", "activated/tissue-resident lymphocytes"),
    "CD49f": ("integrin alpha-6", "basal epithelial/basement-membrane interface"),
    "CD49a": ("integrin alpha-1", "collagen-binding stromal/epithelial/TRM contexts"),
    "CD68": ("lysosomal macrophage marker", "macrophages/mononuclear phagocytes"),
    "OLFM4": ("intestinal stem/progenitor marker", "crypt-base/progenitor epithelial cells"),
    "Podoplanin": ("PDPN/D2-40", "lymphatic endothelium and fibroblast-like stroma"),
    "CD45": ("PTPRC pan-leukocyte", "all leukocytes except mature erythroid cells"),
    "CD163": ("scavenger receptor", "M2-like/resident macrophages"),
    "CD44": ("adhesion receptor", "epithelial stem/immune/stromal adhesion contexts"),
    "CD56": ("NCAM1", "NK cells, neural, neuroendocrine/ICC-like contexts"),
    "Cytokeratin": ("pan-keratin", "epithelial cytoskeleton"),
    "CD206": ("mannose receptor", "M2-like macrophages"),
    "aDefensin5": ("DEFA5", "Paneth cells; expected rare/low in colon"),
    "DRAQ5": ("DNA dye", "all nuclei; should match Hoechst"),
}


def marker_to_path(marker: str) -> Path:
    matches = sorted(CHANNEL_DIR.glob(f"??_{marker}.tif"))
    if not matches:
        raise FileNotFoundError(marker)
    return matches[0]


def read_crop(marker: str, box: tuple[int, int, int, int]) -> np.ndarray:
    x0, y0, x1, y1 = box
    arr = tifffile.imread(marker_to_path(marker))
    return arr[y0:y1, x0:x1]


def norm(a: np.ndarray, lo_p: float = 1, hi_p: float = 99.7, gamma: float = 0.75) -> np.ndarray:
    a = np.asarray(a, dtype=np.float32)
    if a.size == 0:
        return np.zeros(a.shape, dtype=np.float32)
    # Log display preserves both sparse puncta and broad structures.
    a = np.log1p(a)
    lo, hi = np.percentile(a[np.isfinite(a)], [lo_p, hi_p])
    if hi <= lo:
        return np.zeros(a.shape, dtype=np.float32)
    out = np.clip((a - lo) / (hi - lo), 0, 1)
    return np.power(out, gamma)


def norm8(a: np.ndarray, lo_p: float = 1, hi_p: float = 99.7, gamma: float = 0.75) -> np.ndarray:
    return (norm(a, lo_p, hi_p, gamma) * 255).astype(np.uint8)


def downsample_rgb(channels: dict[str, int], gamma: float = 0.75) -> np.ndarray:
    z = np.load(DOWNSAMPLE_NPZ)
    stack = z["stack"]
    names = [str(x) for x in z["names"]]
    idx = {n: i for i, n in enumerate(names)}
    rgb = np.zeros((*stack.shape[1:], 3), dtype=np.float32)
    for marker, channel in channels.items():
        rgb[..., channel] = np.maximum(rgb[..., channel], norm(stack[idx[marker]], gamma=gamma))
    return np.clip(rgb * 255, 0, 255).astype(np.uint8)


def add_boxes(img: np.ndarray) -> Image.Image:
    pil = Image.fromarray(img).convert("RGB")
    draw = ImageDraw.Draw(pil)
    font = ImageFont.load_default()

    boxes = {
        "A crypt mucosa": (500, 4900, 5000, 8900),
        "B stromal/smooth muscle fold": (5000, 0, 9800, 4600),
        "C vessel-rich field": (7600, 4700, 9500, 7100),
        "D pan-marker artifact": (3946, 5420, 4714, 6188),
    }
    colors = {
        "A crypt mucosa": (255, 235, 59),
        "B stromal/smooth muscle fold": (0, 188, 212),
        "C vessel-rich field": (255, 87, 34),
        "D pan-marker artifact": (255, 0, 255),
    }
    for label, (x0, y0, x1, y1) in boxes.items():
        ds_box = tuple(int(v / DS) for v in (x0, y0, x1, y1))
        color = colors[label]
        draw.rectangle(ds_box, outline=color, width=3)
        draw.text((ds_box[0] + 5, ds_box[1] + 5), label, fill=color, font=font)
    return pil


def make_global_assets() -> None:
    context = downsample_rgb({"Cytokeratin": 0, "MUC2": 1, "Hoechst1": 2})
    add_boxes(context).resize((1250, 1190), Image.Resampling.LANCZOS).save(
        ASSETS / "01_global_orientation_boxes.png", quality=95
    )

    composites = [
        ("CK / MUC2 / nuclei", downsample_rgb({"Cytokeratin": 0, "MUC2": 1, "Hoechst1": 2})),
        ("MUC2 / SOX9 / CDX2", downsample_rgb({"MUC2": 0, "SOX9": 1, "CDX2": 2})),
        ("CD3 / CD19 / CD68", downsample_rgb({"CD3": 0, "CD19": 1, "CD68": 2})),
        ("aSMA / CollagenIV / CD31", downsample_rgb({"aSMA": 0, "CollagenIV": 1, "CD31": 2})),
        ("Podoplanin / FAP / CD34", downsample_rgb({"Podoplanin": 0, "FAP": 1, "CD34": 2})),
        ("CD45 / CK / CD31", downsample_rgb({"CD45": 0, "Cytokeratin": 1, "CD31": 2})),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(13, 8.4))
    for ax, (title, im) in zip(axes.flat, composites):
        ax.imshow(im)
        ax.set_title(title, fontsize=10)
        ax.set_axis_off()
    fig.tight_layout(pad=0.3)
    fig.savefig(ASSETS / "02_global_composites.png", dpi=180)
    plt.close(fig)


def make_crop_panel(
    name: str,
    box: tuple[int, int, int, int],
    markers: list[str],
    cols: int = 4,
) -> None:
    rows = int(np.ceil(len(markers) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 3.1, rows * 2.9))
    axes = np.atleast_1d(axes).flat
    for ax, marker in zip(axes, markers):
        crop = read_crop(marker, box)
        ax.imshow(norm8(crop), cmap="gray", interpolation="nearest")
        ax.set_title(marker, fontsize=9)
        ax.set_axis_off()
    for ax in list(axes)[len(markers) :]:
        ax.set_axis_off()
    fig.tight_layout(pad=0.45)
    fig.savefig(ASSETS / f"{name}.png", dpi=180)
    plt.close(fig)


def make_crop_rgb(name: str, box: tuple[int, int, int, int], mapping: dict[str, int]) -> None:
    crops = {m: read_crop(m, box) for m in mapping}
    first = next(iter(crops.values()))
    rgb = np.zeros((*first.shape, 3), dtype=np.float32)
    for marker, ch in mapping.items():
        rgb[..., ch] = np.maximum(rgb[..., ch], norm(crops[marker]))
    im = Image.fromarray((np.clip(rgb, 0, 1) * 255).astype(np.uint8))
    im.thumbnail((1400, 1100), Image.Resampling.LANCZOS)
    im.save(ASSETS / f"{name}.png", quality=95)


def make_roi_assets() -> None:
    lower_left = (500, 4900, 5000, 8900)
    upper_right = (5000, 0, 9800, 4600)
    vessel = (7600, 4700, 9500, 7100)
    artifact = (3946, 5420, 4714, 6188)

    make_crop_rgb(
        "03_lower_left_rgb_ck_muc2_nuclei",
        lower_left,
        {"Cytokeratin": 0, "MUC2": 1, "Hoechst1": 2},
    )
    make_crop_panel(
        "04_lower_left_crypt_marker_panels",
        lower_left,
        [
            "Hoechst1",
            "DRAQ5",
            "Cytokeratin",
            "CDX2",
            "SOX9",
            "OLFM4",
            "Ki67",
            "MUC2",
            "ITLN1",
            "MUC1",
            "CD49f",
            "CollagenIV",
            "aSMA",
            "CD31",
            "CD45",
            "CD68",
        ],
    )
    make_crop_rgb(
        "05_upper_right_rgb_asma_coliv_cd31",
        upper_right,
        {"aSMA": 0, "CollagenIV": 1, "CD31": 2},
    )
    make_crop_panel(
        "06_upper_right_stromal_marker_panels",
        upper_right,
        [
            "Hoechst1",
            "Cytokeratin",
            "MUC2",
            "CDX2",
            "SOX9",
            "aSMA",
            "CollagenIV",
            "FAP",
            "Vimentin",
            "CD90",
            "CD49a",
            "Podoplanin",
            "CD31",
            "CD34",
            "CD45",
            "CD68",
        ],
    )
    make_crop_rgb(
        "07_vessel_field_rgb_cd31_coliv_asma",
        vessel,
        {"CD31": 0, "CollagenIV": 1, "aSMA": 2},
    )
    make_crop_panel(
        "08_vessel_field_marker_panels",
        vessel,
        [
            "Hoechst1",
            "CD31",
            "CD34",
            "CollagenIV",
            "aSMA",
            "Podoplanin",
            "FAP",
            "Vimentin",
            "CD45",
            "CD68",
            "HLADR",
            "CD3",
        ],
        cols=4,
    )
    make_crop_panel(
        "09_pan_marker_artifact_panels",
        artifact,
        [
            "Hoechst1",
            "Cytokeratin",
            "CDX2",
            "MUC2",
            "CD45",
            "CD3",
            "CD19",
            "CD68",
            "aSMA",
            "CollagenIV",
            "CD31",
            "FAP",
        ],
        cols=4,
    )


def robust_z(df: pd.DataFrame, markers: list[str]) -> pd.DataFrame:
    x = np.log1p(df[markers].astype(float))
    med = x.median()
    mad = (x - med).abs().median().replace(0, np.nan)
    z = (x - med) / (1.4826 * mad)
    return z.replace([np.inf, -np.inf], np.nan).fillna(0)


def summarize_regions() -> None:
    cells = pd.read_csv(CELL_QC)
    clusters = pd.read_csv(CLUSTERS)[["label", "leiden", "cluster_label"]]
    df = cells.merge(clusters, on="label", how="left")
    markers = [c for c in cells.columns if c not in {"label", "area", "x", "y"}]
    z = robust_z(df, markers)

    regions = {
        "A_lower_left_crypt_mucosa": (df.x.between(500, 5000) & df.y.between(4900, 8900)),
        "B_upper_right_stromal_fold": (df.x.between(5000, 9800) & df.y.between(0, 4600)),
        "C_right_vessel_field": (df.x.between(7600, 9500) & df.y.between(4700, 7100)),
        "whole_qc": pd.Series(True, index=df.index),
    }

    region_rows = []
    top_marker_rows = []
    for region, mask in regions.items():
        sub = df[mask].copy()
        subz = z.loc[mask]
        counts = sub["cluster_label"].value_counts()
        region_rows.append(
            {
                "region": region,
                "n_cells": int(len(sub)),
                "x_min": float(sub.x.min()) if len(sub) else np.nan,
                "x_max": float(sub.x.max()) if len(sub) else np.nan,
                "y_min": float(sub.y.min()) if len(sub) else np.nan,
                "y_max": float(sub.y.max()) if len(sub) else np.nan,
                **{f"label_count__{k}": int(v) for k, v in counts.items()},
            }
        )
        medz = subz.median().sort_values(ascending=False)
        for marker, value in medz.head(15).items():
            top_marker_rows.append(
                {
                    "region": region,
                    "marker": marker,
                    "median_robust_z_log1p": float(value),
                    "raw_median": float(sub[marker].median()) if len(sub) else np.nan,
                    "interpretation": MARKER_INTERPRETATION.get(marker, ("", ""))[1],
                }
            )

    pd.DataFrame(region_rows).to_csv(ASSETS / "region_celltype_counts.csv", index=False)
    pd.DataFrame(top_marker_rows).to_csv(ASSETS / "region_top_markers.csv", index=False)

    broad_rows = []
    for region, mask in regions.items():
        subz = z.loc[mask]
        for group, members in MARKER_GROUPS.items():
            present = [m for m in members if m in subz]
            broad_rows.append(
                {
                    "region": region,
                    "marker_group": group,
                    "median_of_group_median_z": float(subz[present].median().median()),
                    "top_group_markers": "; ".join(subz[present].median().sort_values(ascending=False).head(5).index),
                }
            )
    pd.DataFrame(broad_rows).to_csv(ASSETS / "region_group_scores.csv", index=False)

    marker_rows = []
    inv = pd.read_csv(ROOT / "marker_expectation" / "channel_inventory.csv")
    stats = pd.read_csv(ROOT / "marker_expectation" / "channel_statistics.csv")
    merged = inv.merge(stats[["marker", "mean", "zero_fraction", "p95", "p99", "p999"]], on="marker")
    for _, row in merged.iterrows():
        marker = row["marker"]
        meaning, tells = MARKER_INTERPRETATION.get(marker, ("", ""))
        marker_rows.append(
            {
                "index": int(row["channel_index"]),
                "marker": marker,
                "group": row["category"],
                "meaning": meaning,
                "what_it_tells_us": tells,
                "observed_global_pattern": row["expected_pattern"],
                "mean_intensity": float(row["mean"]),
                "zero_fraction": float(row["zero_fraction"]),
                "p99": float(row["p99"]),
                "p999": float(row["p999"]),
            }
        )
    pd.DataFrame(marker_rows).to_csv(ASSETS / "marker_interpretation_54.csv", index=False)


def main() -> None:
    make_global_assets()
    make_roi_assets()
    summarize_regions()


if __name__ == "__main__":
    main()
