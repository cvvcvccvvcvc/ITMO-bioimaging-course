from __future__ import annotations

import json
import math
import os
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path("marker_expectation")
OUT.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".mplconfig"))

import matplotlib.pyplot as plt
import tifffile
from scipy import ndimage as ndi
from skimage import color, filters, measure, morphology


OME_PATH = Path("data/reg001_expr.ome.tif")
CHANNELS_DIR = Path("data/channels")
DOWNSAMPLE = 16
MAX_WORKERS = 4


MARKER_META = {
    "Hoechst1": ("nuclear", "DNA/nuclei; should outline essentially all nuclei."),
    "DRAQ5": ("nuclear", "DNA/nuclei; should broadly match Hoechst."),
    "Ki67": ("nuclear/proliferation", "Proliferating nuclear cells, especially crypt proliferative zones."),
    "MUC1": ("epithelial/mucin", "Apical epithelial mucin/glycoprotein; patchy epithelial signal is expected."),
    "MUC2": ("epithelial/mucin", "Goblet-cell mucin; strong granular/mucus epithelial signal is expected."),
    "CDX2": ("epithelial/nuclear TF", "Intestinal epithelial nuclear transcription factor."),
    "SOX9": ("epithelial/nuclear TF", "Crypt progenitor/stem epithelial nuclear marker."),
    "Cytokeratin": ("epithelial", "Pan-epithelial cytoskeleton."),
    "CK7": ("epithelial", "Epithelial keratin; limited/patchy in colon is plausible."),
    "CD49f": ("epithelial/basal", "Integrin alpha-6; epithelial basement/basal compartments."),
    "CD49a": ("epithelial/stromal", "Integrin alpha-1; tissue-resident immune/stromal/epithelial contexts."),
    "CD44": ("epithelial/immune/stromal", "Adhesion/stem/immune marker; broad but not pan-everything signal."),
    "CD66": ("epithelial/myeloid", "CEACAM family; epithelial/myeloid localization can occur."),
    "OLFM4": ("epithelial/stem", "Crypt stem/progenitor-associated epithelial marker."),
    "ITLN1": ("epithelial/goblet", "Goblet-cell/epithelial lectin."),
    "CHGA": ("neuroendocrine", "Enteroendocrine cells; sparse punctate epithelial cells expected."),
    "Synaptophysin": ("neuroendocrine/neural", "Sparse neuroendocrine or neural structures expected."),
    "aDefensin5": ("epithelial/Paneth", "Paneth-cell marker; low/rare in colon, strong diffuse signal would be suspect."),
    "CD45": ("immune/pan-leukocyte", "Pan-leukocyte; lamina propria and lymphoid aggregates."),
    "CD45RO": ("immune/T-memory", "Memory T cells; should overlap immune-rich areas."),
    "CD3": ("immune/T", "T cells; should overlap CD45/CD7 and subsets CD4/CD8."),
    "CD4": ("immune/T-helper", "Helper T cells and some myeloid contexts."),
    "CD8": ("immune/T-cytotoxic", "Cytotoxic T cells; should overlap CD3/CD45."),
    "CD7": ("immune/T/NK", "T/NK cells; should overlap CD3/CD45."),
    "CD25": ("immune/T-reg/activation", "Activated/regulatory T cells; sparse immune marker."),
    "CD127": ("immune/T", "IL7R; T-cell subset marker."),
    "CD69": ("immune/activation", "Activated tissue-resident immune cells; sparse/patchy."),
    "CD161": ("immune/T/NK", "MAIT/NK/T subsets; immune-rich areas."),
    "CD57": ("immune/NK/T", "NK/senescent T/neural contexts; sparse immune-like signal."),
    "NKG2G": ("immune/NK", "NK-like marker name in panel; should be immune-associated if specific."),
    "CD56": ("immune/NK/neural", "NK cells/neural/endocrine contexts; sparse."),
    "CD19": ("immune/B", "B cells; should form lymphoid aggregate or scattered lamina propria signal."),
    "CD21": ("immune/B/FDC", "B-cell/follicular dendritic network marker; lymphoid follicles."),
    "CD38": ("immune/plasma", "Activated B/plasma cells; immune-rich lamina propria."),
    "CD138": ("immune/plasma/epithelial", "Plasma cells and epithelial contexts; should not be universal."),
    "BCL2": ("immune/survival", "Anti-apoptotic marker; lymphocytes/crypt contexts possible."),
    "CD11c": ("immune/myeloid/APC", "Dendritic/myeloid antigen-presenting cells."),
    "HLADR": ("immune/APC", "Antigen-presenting cells and some epithelium; immune-rich areas."),
    "CD123": ("immune/pDC", "Plasmacytoid dendritic cells; sparse."),
    "CD15": ("immune/granulocyte", "Granulocytes; sparse or luminal inflammatory signal."),
    "CD16": ("immune/myeloid/NK", "Myeloid/NK marker; sparse/patchy."),
    "CD68": ("immune/macrophage", "Macrophages; lamina propria puncta."),
    "CD163": ("immune/macrophage", "Macrophage subset; should overlap CD68/CD206 partly."),
    "CD206": ("immune/macrophage", "Macrophage/mannose receptor; should overlap macrophage-rich stroma."),
    "CD117": ("immune/mast", "Mast cells/interstitial cells; sparse punctate."),
    "Vimentin": ("stromal", "Mesenchymal/stromal cells and some immune cells."),
    "aSMA": ("stromal/smooth muscle", "Smooth muscle/pericytes/myofibroblasts; fiber-like structures."),
    "FAP": ("stromal/fibroblast", "Activated fibroblasts; stromal bands."),
    "CollagenIV": ("ECM/basement membrane", "Basement membrane/vascular ECM; linear structures."),
    "CD31": ("vascular/endothelial", "Blood endothelium; vessel-like loops/tubes."),
    "CD34": ("vascular/stromal", "Endothelium/stromal progenitor; vessel/stromal signal."),
    "CD36": ("vascular/stromal/immune", "Endothelium/adipocyte/macrophage contexts; patchy."),
    "CD90": ("stromal/T", "Fibroblast/stromal and T-cell contexts."),
    "Podoplanin": ("lymphatic/stromal", "Lymphatic endothelium and fibroblast-like stroma."),
}


PROJECT_SELECTED = {
    p.name.split("_", 1)[1].rsplit(".", 1)[0]
    for p in CHANNELS_DIR.glob("*.tif")
    if "_" in p.name
}


def parse_ome_metadata(path: Path):
    with tifffile.TiffFile(path) as tif:
        root = ET.fromstring(tif.ome_metadata)
        ns = {"ome": root.tag.split("}")[0].strip("{")}
        pixels = root.find(".//ome:Pixels", ns)
        channels = root.findall(".//ome:Channel", ns)
        names = [ch.attrib.get("Name", f"Channel_{i}") for i, ch in enumerate(channels)]
        physical_x = float(pixels.attrib["PhysicalSizeX"])
        physical_y = float(pixels.attrib["PhysicalSizeY"])
        size_x = int(pixels.attrib["SizeX"])
        size_y = int(pixels.attrib["SizeY"])
        dtype = pixels.attrib["Type"]
    return names, physical_x, physical_y, size_x, size_y, dtype


def robust_display(img, low=1.0, high=99.8):
    finite = np.isfinite(img)
    if not finite.any():
        return np.zeros_like(img, dtype=np.float32)
    lo, hi = np.percentile(img[finite], [low, high])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(np.nanmin(img)), float(np.nanmax(img))
    if hi <= lo:
        return np.zeros_like(img, dtype=np.float32)
    out = (img.astype(np.float32) - lo) / (hi - lo)
    return np.clip(out, 0, 1)


def crop_to_factor(arr: np.ndarray, factor: int):
    h = (arr.shape[0] // factor) * factor
    w = (arr.shape[1] // factor) * factor
    return arr[:h, :w], h, w


def block_mean(arr: np.ndarray, factor: int):
    cropped, h, w = crop_to_factor(arr, factor)
    return cropped.reshape(h // factor, factor, w // factor, factor).mean(
        axis=(1, 3), dtype=np.float32
    )


def highpass_profile_metric(profile: np.ndarray, smooth_window: int = 61):
    profile = profile.astype(np.float32)
    smooth = ndi.uniform_filter1d(profile, size=smooth_window, mode="nearest")
    hp = profile - smooth
    denom = np.median(np.abs(profile - np.median(profile))) + 1e-6
    return {
        "hp_mad_ratio": float(np.median(np.abs(hp - np.median(hp))) / denom),
        "hp_p99": float(np.percentile(np.abs(hp), 99)),
        "hp_rms": float(np.sqrt(np.mean(hp * hp))),
    }


def seam_ratio(profile: np.ndarray, period: int):
    diffs = np.abs(np.diff(profile.astype(np.float32)))
    if len(diffs) < period * 2:
        return np.nan
    seam_idx = np.arange(period, len(profile), period)
    seam_idx = seam_idx[(seam_idx > 0) & (seam_idx < len(diffs))]
    if len(seam_idx) == 0:
        return np.nan
    return float((np.mean(diffs[seam_idx]) + 1e-6) / (np.median(diffs) + 1e-6))


def save_contact_sheet(images, titles, out_path, cols=9, cmap="gray"):
    n = len(images)
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 2.1, rows * 2.05))
    axes = np.atleast_1d(axes).ravel()
    for ax, img, title in zip(axes, images, titles):
        ax.imshow(img, cmap=cmap, interpolation="nearest")
        ax.set_title(title, fontsize=7)
        ax.set_axis_off()
    for ax in axes[n:]:
        ax.set_axis_off()
    fig.tight_layout(pad=0.3)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def make_rgb(ds_stack, names, marker_to_rgb, gamma=0.85):
    h, w = ds_stack.shape[1:]
    rgb = np.zeros((h, w, 3), dtype=np.float32)
    name_to_idx = {n: i for i, n in enumerate(names)}
    for marker, channel in marker_to_rgb.items():
        if marker not in name_to_idx:
            continue
        img = np.log1p(ds_stack[name_to_idx[marker]])
        disp = robust_display(img, 1, 99.6) ** gamma
        rgb[..., channel] = np.maximum(rgb[..., channel], disp)
    return np.clip(rgb, 0, 1)


def save_overview_composites(ds_stack, names, out_path):
    composites = [
        (
            "Nuclei / epithelium / immune",
            make_rgb(
                ds_stack,
                names,
                {"Hoechst1": 2, "Cytokeratin": 1, "CD45": 0},
            ),
        ),
        (
            "Epithelial differentiation",
            make_rgb(
                ds_stack,
                names,
                {"MUC2": 0, "SOX9": 1, "CDX2": 2},
            ),
        ),
        (
            "T / B / macrophage",
            make_rgb(
                ds_stack,
                names,
                {"CD3": 0, "CD19": 1, "CD68": 2},
            ),
        ),
        (
            "Stroma / ECM / endothelium",
            make_rgb(
                ds_stack,
                names,
                {"aSMA": 0, "CollagenIV": 1, "CD31": 2},
            ),
        ),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11, 10))
    for ax, (title, img) in zip(axes.ravel(), composites):
        ax.imshow(img, interpolation="nearest")
        ax.set_title(title, fontsize=10)
        ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def save_heatmap(arr, title, out_path, cmap="magma", vmax=None):
    fig, ax = plt.subplots(figsize=(9, 8))
    im = ax.imshow(arr, cmap=cmap, interpolation="nearest", vmax=vmax)
    ax.set_title(title)
    ax.set_axis_off()
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def analyze_multimarker(ds_stack, names):
    log_stack = np.log1p(ds_stack.astype(np.float32))
    non_nuclear = [i for i, n in enumerate(names) if n not in {"Hoechst1", "DRAQ5"}]
    thresholds = np.percentile(log_stack[non_nuclear].reshape(len(non_nuclear), -1), 99.5, axis=1)
    pos = log_stack[non_nuclear] > thresholds[:, None, None]
    pos_count = pos.sum(axis=0).astype(np.uint8)

    threshold = 15
    mask = pos_count >= threshold
    if not mask.any():
        threshold = 12
        mask = pos_count >= threshold
    if not mask.any():
        threshold = 10
        mask = pos_count >= threshold

    labels, n_labels = ndi.label(mask)
    regions = measure.regionprops(labels, intensity_image=pos_count)

    nuclear_high = np.zeros_like(pos_count, dtype=bool)
    for nuc in ["Hoechst1", "DRAQ5"]:
        if nuc in names:
            idx = names.index(nuc)
            nuclear_high |= log_stack[idx] > np.percentile(log_stack[idx], 95)

    rows = []
    name_to_group = {name: MARKER_META.get(name, ("unknown", ""))[0] for name in names}
    categories = sorted(set(name_to_group.values()))
    for reg in regions:
        if reg.area < 3:
            continue
        coords = reg.coords
        y, x = np.average(coords, axis=0, weights=pos_count[coords[:, 0], coords[:, 1]])
        member = labels == reg.label
        marker_hits = []
        group_hits = {cat: 0 for cat in categories}
        for local_j, global_i in enumerate(non_nuclear):
            if pos[local_j][member].mean() > 0.25:
                marker_hits.append(names[global_i])
                group_hits[name_to_group[names[global_i]]] += 1
        rows.append(
            {
                "component": int(reg.label),
                "area_blocks": int(reg.area),
                "area_um2": float(reg.area * (DOWNSAMPLE * 0.3774) ** 2),
                "centroid_y_px": float(y * DOWNSAMPLE + DOWNSAMPLE / 2),
                "centroid_x_px": float(x * DOWNSAMPLE + DOWNSAMPLE / 2),
                "max_non_nuclear_positive_markers": int(reg.max_intensity),
                "mean_non_nuclear_positive_markers": float(pos_count[member].mean()),
                "nuclear_high_fraction": float(nuclear_high[member].mean()),
                "hit_markers": ";".join(marker_hits),
                "hit_marker_count": len(marker_hits),
                **{f"group_hits__{k}": v for k, v in group_hits.items() if v},
            }
        )
    components = pd.DataFrame(rows).sort_values(
        ["max_non_nuclear_positive_markers", "area_blocks"], ascending=False
    )
    return pos_count, threshold, components


def marker_limits_from_stats(stats: pd.DataFrame):
    out = {}
    for _, row in stats.iterrows():
        lo = float(row["p1"])
        hi = float(row["p998"])
        if hi <= lo:
            hi = float(row["max"])
        out[row["marker"]] = (lo, hi)
    return out


def normalize_crop(crop: np.ndarray, limits):
    lo, hi = limits
    if hi <= lo:
        hi = float(crop.max())
    if hi <= lo:
        return np.zeros_like(crop, dtype=np.float32)
    return np.clip((crop.astype(np.float32) - lo) / (hi - lo), 0, 1)


def save_crop_panel(crop_by_marker, markers, title, out_path, limits, cols=5):
    rows = math.ceil(len(markers) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 2.2, rows * 2.25))
    axes = np.atleast_1d(axes).ravel()
    for ax, marker in zip(axes, markers):
        img = crop_by_marker.get(marker)
        if img is None:
            ax.set_axis_off()
            continue
        ax.imshow(normalize_crop(img, limits[marker]), cmap="gray", interpolation="nearest")
        ax.set_title(marker, fontsize=8)
        ax.set_axis_off()
    for ax in axes[len(markers):]:
        ax.set_axis_off()
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def save_spot_crops(ome_path, names, stats, components, out_dir):
    if components.empty:
        return
    out_dir.mkdir(exist_ok=True)
    limits = marker_limits_from_stats(stats)
    name_to_idx = {n: i for i, n in enumerate(names)}
    top = components.head(6).copy()
    crop_size = 768
    half = crop_size // 2
    selected = [
        "Hoechst1",
        "DRAQ5",
        "Cytokeratin",
        "CDX2",
        "SOX9",
        "MUC2",
        "OLFM4",
        "aDefensin5",
        "CD45",
        "CD3",
        "CD4",
        "CD8",
        "CD19",
        "CD21",
        "CD68",
        "CD163",
        "CD206",
        "CD138",
        "Vimentin",
        "CollagenIV",
        "aSMA",
        "FAP",
        "CD31",
        "Podoplanin",
    ]
    selected = [m for m in selected if m in name_to_idx]
    all_first_spot = list(names)

    coords = []
    for _, row in top.iterrows():
        cy, cx = int(row["centroid_y_px"]), int(row["centroid_x_px"])
        y0 = max(0, cy - half)
        x0 = max(0, cx - half)
        y1 = min(9514, y0 + crop_size)
        x1 = min(9990, x0 + crop_size)
        y0 = max(0, y1 - crop_size)
        x0 = max(0, x1 - crop_size)
        coords.append((int(row["component"]), y0, y1, x0, x1))

    selected_crops = {component: {} for component, *_ in coords}
    first_component = coords[0][0]
    all_marker_crops = {}
    needed = sorted(set(selected + all_first_spot), key=lambda n: name_to_idx[n])
    with tifffile.TiffFile(ome_path) as tif:
        for marker in needed:
            arr = tif.pages[name_to_idx[marker]].asarray(maxworkers=MAX_WORKERS)
            for component, y0, y1, x0, x1 in coords:
                if marker in selected:
                    selected_crops[component][marker] = arr[y0:y1, x0:x1].copy()
                if component == first_component:
                    all_marker_crops[marker] = arr[y0:y1, x0:x1].copy()

    for row, (component, y0, y1, x0, x1) in zip(top.itertuples(index=False), coords):
        title = (
            f"component {component}, crop x={x0}:{x1}, y={y0}:{y1}, "
            f"max positive markers={row.max_non_nuclear_positive_markers}, "
            f"nuclear high fraction={row.nuclear_high_fraction:.2f}"
        )
        save_crop_panel(
            selected_crops[component],
            selected,
            title,
            out_dir / f"spot_{component:03d}_selected_markers.png",
            limits,
            cols=6,
        )
    save_crop_panel(
        all_marker_crops,
        all_first_spot,
        f"All 54 markers in top multi-marker component {first_component}",
        out_dir / f"spot_{first_component:03d}_all_54_markers.png",
        limits,
        cols=9,
    )


def analyze_nuclear_halo(ome_path, names, stats, out_path):
    if "Hoechst1" not in names:
        return pd.DataFrame()
    idx = names.index("Hoechst1")
    with tifffile.TiffFile(ome_path) as tif:
        arr = tif.pages[idx].asarray(maxworkers=MAX_WORKERS)
    # Use the same neighborhood as the downstream Cellpose ROI, but enlarge it.
    roi_info_path = Path("data/cellpose_roi.json")
    if roi_info_path.exists():
        roi = json.loads(roi_info_path.read_text())
        cy = int(roi["y0"] + roi["height"] / 2)
        cx = int(roi["x0"] + roi["width"] / 2)
    else:
        cy, cx = arr.shape[0] // 2, arr.shape[1] // 2
    crop_size = 1024
    y0 = max(0, min(arr.shape[0] - crop_size, cy - crop_size // 2))
    x0 = max(0, min(arr.shape[1] - crop_size, cx - crop_size // 2))
    crop = arr[y0 : y0 + crop_size, x0 : x0 + crop_size]
    sm = filters.gaussian(crop.astype(np.float32), sigma=1)
    thresh = np.percentile(sm, 97)
    nuclei = sm > thresh
    nuclei = morphology.remove_small_objects(nuclei, 20)
    core = morphology.binary_dilation(nuclei, morphology.disk(2))
    ring = morphology.binary_dilation(nuclei, morphology.disk(10)) & ~core
    far = morphology.binary_dilation(nuclei, morphology.disk(26)) & ~morphology.binary_dilation(
        nuclei, morphology.disk(14)
    )
    far &= crop > np.percentile(crop, 5)

    ring_med = float(np.median(crop[ring])) if ring.any() else np.nan
    far_med = float(np.median(crop[far])) if far.any() else np.nan
    bg_med = float(np.median(crop[crop < np.percentile(crop, 40)]))
    core_med = float(np.median(crop[core])) if core.any() else np.nan
    halo_ratio_far = float(ring_med / far_med) if np.isfinite(far_med) and far_med > 0 else np.nan
    halo_ratio_bg = float(ring_med / bg_med) if bg_med > 0 else np.nan

    overlay = color.gray2rgb(robust_display(crop, 1, 99.8))
    overlay[ring, 0] = 1
    overlay[ring, 1:] *= 0.2
    overlay[core, 1] = 1
    overlay[core, [0, 2]] *= 0.2

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    axes[0].imshow(robust_display(crop, 1, 99.8), cmap="gray")
    axes[0].set_title("Hoechst ROI")
    axes[1].imshow(overlay)
    axes[1].set_title("nuclei cores green, rings red")
    axes[2].hist(crop[ring].ravel(), bins=80, alpha=0.6, label="ring")
    if far.any():
        axes[2].hist(crop[far].ravel(), bins=80, alpha=0.6, label="far local")
    axes[2].set_title("local background distributions")
    axes[2].legend(fontsize=8)
    for ax in axes[:2]:
        ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)

    return pd.DataFrame(
        [
            {
                "marker": "Hoechst1",
                "roi_y0": y0,
                "roi_x0": x0,
                "core_median": core_med,
                "ring_median": ring_med,
                "far_local_median": far_med,
                "low_background_median": bg_med,
                "ring_to_far_ratio": halo_ratio_far,
                "ring_to_low_background_ratio": halo_ratio_bg,
                "nucleus_core_pixels": int(core.sum()),
                "ring_pixels": int(ring.sum()),
                "far_pixels": int(far.sum()),
            }
        ]
    )


def main():
    names, physical_x_nm, physical_y_nm, size_x, size_y, dtype = parse_ome_metadata(OME_PATH)
    inventory = []
    thumbnails = []
    downsampled = []
    stats_rows = []
    stripe_rows = []

    print(f"OME: {len(names)} channels, {size_x}x{size_y}, {dtype}, {physical_x_nm} nm/px")
    tile_period_ds = 512 // 4

    with tifffile.TiffFile(OME_PATH) as tif:
        for i, name in enumerate(names):
            page = tif.pages[i]
            print(f"[{i + 1:02d}/{len(names)}] {name}", flush=True)
            arr = page.asarray(maxworkers=MAX_WORKERS)
            sample = arr[::8, ::8].astype(np.float32)
            pct_names = ["p0_1", "p1", "p5", "p10", "p50", "p95", "p99", "p995", "p998", "p999", "p9999"]
            pct_vals = np.percentile(sample, [0.1, 1, 5, 10, 50, 95, 99, 99.5, 99.8, 99.9, 99.99])
            row = {
                "channel_index": i,
                "marker": name,
                "category": MARKER_META.get(name, ("unknown", ""))[0],
                "expected_pattern": MARKER_META.get(name, ("unknown", "Not annotated."))[1],
                "used_in_current_work": name in PROJECT_SELECTED,
                "min": int(arr.min()),
                "max": int(arr.max()),
                "mean": float(arr.mean()),
                "std": float(arr.std()),
                "zero_fraction": float((arr == 0).mean()),
                "near_saturation_fraction": float((arr >= 65500).mean()),
            }
            row.update({k: float(v) for k, v in zip(pct_names, pct_vals)})
            stats_rows.append(row)

            ds = block_mean(arr, DOWNSAMPLE)
            downsampled.append(ds)
            thumbnails.append(robust_display(np.log1p(ds), 1, 99.8))

            log_sample = np.log1p(arr[::4, ::4].astype(np.float32))
            row_profile = np.median(log_sample, axis=1)
            col_profile = np.median(log_sample, axis=0)
            row_metrics = highpass_profile_metric(row_profile)
            col_metrics = highpass_profile_metric(col_profile)
            stripe_rows.append(
                {
                    "channel_index": i,
                    "marker": name,
                    "horizontal_profile_hp_mad_ratio": row_metrics["hp_mad_ratio"],
                    "horizontal_profile_hp_p99": row_metrics["hp_p99"],
                    "horizontal_profile_hp_rms": row_metrics["hp_rms"],
                    "vertical_profile_hp_mad_ratio": col_metrics["hp_mad_ratio"],
                    "vertical_profile_hp_p99": col_metrics["hp_p99"],
                    "vertical_profile_hp_rms": col_metrics["hp_rms"],
                    "horizontal_tile_seam_ratio_512px": seam_ratio(row_profile, tile_period_ds),
                    "vertical_tile_seam_ratio_512px": seam_ratio(col_profile, tile_period_ds),
                }
            )

            inventory.append(
                {
                    "channel_index": i,
                    "marker": name,
                    "category": MARKER_META.get(name, ("unknown", ""))[0],
                    "expected_pattern": MARKER_META.get(name, ("unknown", "Not annotated."))[1],
                    "used_in_current_work": name in PROJECT_SELECTED,
                }
            )

    stats = pd.DataFrame(stats_rows)
    stripes = pd.DataFrame(stripe_rows)
    inventory = pd.DataFrame(inventory)
    ds_stack = np.stack(downsampled, axis=0).astype(np.float32)

    inventory.to_csv(OUT / "channel_inventory.csv", index=False)
    stats.to_csv(OUT / "channel_statistics.csv", index=False)
    stripes.to_csv(OUT / "stripe_metrics.csv", index=False)
    np.savez_compressed(OUT / f"downsampled_mean_f{DOWNSAMPLE}.npz", stack=ds_stack, names=np.array(names))

    titles = [f"{i:02d} {n}" for i, n in enumerate(names)]
    save_contact_sheet(thumbnails, titles, OUT / "all_54_marker_overview_log_scaled.png", cols=9)

    for category, idxs in inventory.groupby("category").groups.items():
        idxs = list(idxs)
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", category)
        save_contact_sheet(
            [thumbnails[i] for i in idxs],
            [titles[i] for i in idxs],
            OUT / f"group_{safe}.png",
            cols=min(6, max(1, len(idxs))),
        )

    save_overview_composites(ds_stack, names, OUT / "biological_group_rgb_overviews.png")

    pos_count, positive_threshold, components = analyze_multimarker(ds_stack, names)
    components.to_csv(OUT / "multimarker_positive_components.csv", index=False)
    save_heatmap(
        pos_count,
        f"Non-nuclear markers above per-marker p99.5; components threshold >= {positive_threshold}",
        OUT / "multimarker_positive_count_heatmap.png",
        cmap="magma",
        vmax=max(5, int(pos_count.max())),
    )

    corr_x = np.log1p(ds_stack.reshape(ds_stack.shape[0], -1).T)
    corr = np.corrcoef(corr_x, rowvar=False)
    corr_df = pd.DataFrame(corr, index=names, columns=names)
    corr_df.to_csv(OUT / "downsampled_log_marker_correlation.csv")
    fig, ax = plt.subplots(figsize=(11, 10))
    im = ax.imshow(corr, vmin=-0.3, vmax=1, cmap="coolwarm")
    ax.set_xticks(range(len(names)))
    ax.set_yticks(range(len(names)))
    ax.set_xticklabels(names, rotation=90, fontsize=5)
    ax.set_yticklabels(names, fontsize=5)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    fig.tight_layout()
    fig.savefig(OUT / "downsampled_log_marker_correlation_heatmap.png", dpi=220)
    plt.close(fig)

    # Save profile views for the channels where profile striping is strongest.
    stripe_score = (
        stripes["horizontal_profile_hp_mad_ratio"].fillna(0)
        + stripes["vertical_profile_hp_mad_ratio"].fillna(0)
    )
    top_stripe = stripes.assign(stripe_score=stripe_score).sort_values("stripe_score", ascending=False).head(8)
    save_contact_sheet(
        [thumbnails[int(i)] for i in top_stripe["channel_index"]],
        [
            f"{int(row.channel_index):02d} {row.marker}\nH {row.horizontal_profile_hp_mad_ratio:.2f} V {row.vertical_profile_hp_mad_ratio:.2f}"
            for row in top_stripe.itertuples(index=False)
        ],
        OUT / "top_stripe_metric_channels.png",
        cols=4,
    )

    halo = analyze_nuclear_halo(OME_PATH, names, stats, OUT / "hoechst_nuclear_halo_roi.png")
    halo.to_csv(OUT / "nuclear_halo_metrics.csv", index=False)

    save_spot_crops(OME_PATH, names, stats, components, OUT / "multimarker_spot_crops")

    summary = {
        "ome_path": str(OME_PATH),
        "channel_count": len(names),
        "shape_yx": [size_y, size_x],
        "physical_size_nm": [physical_y_nm, physical_x_nm],
        "downsample_factor": DOWNSAMPLE,
        "multimarker_component_threshold": positive_threshold,
        "max_non_nuclear_positive_markers_in_block": int(pos_count.max()),
        "multimarker_component_count": int(len(components)),
        "top_multimarker_components": components.head(10).to_dict(orient="records"),
        "top_stripe_channels": top_stripe.to_dict(orient="records"),
        "halo_metrics": halo.to_dict(orient="records"),
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2))
    print("Wrote outputs to", OUT.resolve())


if __name__ == "__main__":
    main()
