from __future__ import annotations

import json
import math
import re
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import tifffile
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage as ndi
from skimage import filters, measure, morphology


OUT = Path("marker_expectation")
OME_PATH = Path("data/reg001_expr.ome.tif")
DOWNSAMPLE = 16
MAX_WORKERS = 4
PHYSICAL_SIZE_UM = 0.3774


def robust_display(img, low=1.0, high=99.8):
    finite = np.isfinite(img)
    if not finite.any():
        return np.zeros_like(img, dtype=np.uint8)
    lo, hi = np.percentile(img[finite], [low, high])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(np.nanmin(img)), float(np.nanmax(img))
    if hi <= lo:
        return np.zeros_like(img, dtype=np.uint8)
    out = np.clip((img.astype(np.float32) - lo) / (hi - lo), 0, 1)
    return (out * 255).astype(np.uint8)


def get_font(size=12):
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
        "/Library/Fonts/Arial.ttf",
    ]
    for path in candidates:
        if Path(path).exists():
            try:
                return ImageFont.truetype(path, size=size)
            except Exception:
                pass
    return ImageFont.load_default()


def paste_centered(canvas, tile, box):
    x0, y0, x1, y1 = box
    w, h = x1 - x0, y1 - y0
    tile = Image.fromarray(tile) if isinstance(tile, np.ndarray) else tile
    tile.thumbnail((w, h), Image.Resampling.LANCZOS)
    canvas.paste(tile, (x0 + (w - tile.width) // 2, y0 + (h - tile.height) // 2))


def save_sheet(tiles, titles, out_path, cols=6, tile_w=210, tile_h=190, title_h=34):
    rows = math.ceil(len(tiles) / cols)
    canvas = Image.new("RGB", (cols * tile_w, rows * (tile_h + title_h)), "white")
    draw = ImageDraw.Draw(canvas)
    font = get_font(11)
    for k, (tile, title) in enumerate(zip(tiles, titles)):
        r, c = divmod(k, cols)
        x = c * tile_w
        y = r * (tile_h + title_h)
        if isinstance(tile, np.ndarray):
            if tile.ndim == 2:
                img = Image.fromarray(tile).convert("RGB")
            else:
                img = Image.fromarray(tile)
        else:
            img = tile
        paste_centered(canvas, img, (x + 4, y + title_h, x + tile_w - 4, y + title_h + tile_h - 4))
        draw.multiline_text((x + 6, y + 4), str(title), fill=(0, 0, 0), font=font, spacing=1)
    canvas.save(out_path)


def make_rgb(ds_stack, names, marker_to_rgb):
    h, w = ds_stack.shape[1:]
    rgb = np.zeros((h, w, 3), dtype=np.float32)
    name_to_idx = {n: i for i, n in enumerate(names)}
    for marker, channel in marker_to_rgb.items():
        if marker not in name_to_idx:
            continue
        disp = robust_display(np.log1p(ds_stack[name_to_idx[marker]]), 1, 99.6).astype(np.float32) / 255
        rgb[..., channel] = np.maximum(rgb[..., channel], disp ** 0.85)
    return (np.clip(rgb, 0, 1) * 255).astype(np.uint8)


def save_colormap(arr, out_path, vmin=None, vmax=None, cmap=cv2.COLORMAP_MAGMA):
    arr = arr.astype(np.float32)
    if vmin is None:
        vmin = float(np.nanmin(arr))
    if vmax is None:
        vmax = float(np.nanmax(arr))
    if vmax <= vmin:
        vmax = vmin + 1
    gray = np.clip((arr - vmin) / (vmax - vmin), 0, 1)
    colored = cv2.applyColorMap((gray * 255).astype(np.uint8), cmap)
    colored = cv2.cvtColor(colored, cv2.COLOR_BGR2RGB)
    Image.fromarray(colored).save(out_path)


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
    labels, _ = ndi.label(mask)
    regions = measure.regionprops(labels, intensity_image=pos_count)

    nuclear_high = np.zeros_like(pos_count, dtype=bool)
    for marker in ["Hoechst1", "DRAQ5"]:
        if marker in names:
            idx = names.index(marker)
            nuclear_high |= log_stack[idx] > np.percentile(log_stack[idx], 95)

    category = pd.read_csv(OUT / "channel_inventory.csv").set_index("marker")["category"].to_dict()
    rows = []
    for reg in regions:
        if reg.area < 3:
            continue
        member = labels == reg.label
        coords = reg.coords
        y, x = np.average(coords, axis=0, weights=pos_count[coords[:, 0], coords[:, 1]])
        marker_hits = []
        group_hits = {}
        for local_j, global_i in enumerate(non_nuclear):
            if pos[local_j][member].mean() > 0.25:
                marker = names[global_i]
                marker_hits.append(marker)
                group = category.get(marker, "unknown")
                group_hits[group] = group_hits.get(group, 0) + 1
        row = {
            "component": int(reg.label),
            "area_blocks": int(reg.area),
            "area_um2": float(reg.area * (DOWNSAMPLE * PHYSICAL_SIZE_UM) ** 2),
            "centroid_y_px": float(y * DOWNSAMPLE + DOWNSAMPLE / 2),
            "centroid_x_px": float(x * DOWNSAMPLE + DOWNSAMPLE / 2),
            "max_non_nuclear_positive_markers": int(reg.max_intensity),
            "mean_non_nuclear_positive_markers": float(pos_count[member].mean()),
            "nuclear_high_fraction": float(nuclear_high[member].mean()),
            "hit_marker_count": len(marker_hits),
            "hit_markers": ";".join(marker_hits),
        }
        row.update({f"group_hits__{k}": v for k, v in group_hits.items()})
        rows.append(row)
    components = pd.DataFrame(rows)
    if not components.empty:
        components = components.sort_values(
            ["max_non_nuclear_positive_markers", "area_blocks"], ascending=False
        )
    return pos_count, threshold, components


def normalize_crop(crop, lo, hi):
    if hi <= lo:
        hi = float(crop.max())
    if hi <= lo:
        return np.zeros_like(crop, dtype=np.uint8)
    return (np.clip((crop.astype(np.float32) - lo) / (hi - lo), 0, 1) * 255).astype(np.uint8)


def save_spot_crops(names, stats, components):
    if components.empty:
        return
    out_dir = OUT / "multimarker_spot_crops"
    out_dir.mkdir(exist_ok=True)
    limits = {r.marker: (float(r.p1), float(r.p998)) for r in stats.itertuples(index=False)}
    name_to_idx = {n: i for i, n in enumerate(names)}
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
    top = components.head(5)
    crop_size = 768
    half = crop_size // 2
    coords = []
    for row in top.itertuples(index=False):
        cy, cx = int(row.centroid_y_px), int(row.centroid_x_px)
        y0 = max(0, min(9514 - crop_size, cy - half))
        x0 = max(0, min(9990 - crop_size, cx - half))
        coords.append((int(row.component), y0, y0 + crop_size, x0, x0 + crop_size))

    selected_crops = {component: {} for component, *_ in coords}
    first_component = coords[0][0]
    all_marker_crops = {}
    needed = sorted(set(selected + list(names)), key=lambda n: name_to_idx[n])
    with tifffile.TiffFile(OME_PATH) as tif:
        for marker in needed:
            print(f"  crop read {marker}", flush=True)
            arr = tif.pages[name_to_idx[marker]].asarray(maxworkers=MAX_WORKERS)
            lo, hi = limits[marker]
            for component, y0, y1, x0, x1 in coords:
                crop = normalize_crop(arr[y0:y1, x0:x1], lo, hi)
                if marker in selected:
                    selected_crops[component][marker] = crop
                if component == first_component:
                    all_marker_crops[marker] = crop

    for row, (component, y0, y1, x0, x1) in zip(top.itertuples(index=False), coords):
        tiles = [selected_crops[component][m] for m in selected]
        titles = [m for m in selected]
        save_sheet(
            tiles,
            titles,
            out_dir / f"spot_{component:03d}_selected_markers.png",
            cols=6,
            tile_w=170,
            tile_h=160,
            title_h=24,
        )
        meta = {
            "component": component,
            "crop_yx": [y0, y1, x0, x1],
            "max_non_nuclear_positive_markers": int(row.max_non_nuclear_positive_markers),
            "mean_non_nuclear_positive_markers": float(row.mean_non_nuclear_positive_markers),
            "nuclear_high_fraction": float(row.nuclear_high_fraction),
            "hit_markers": row.hit_markers,
        }
        (out_dir / f"spot_{component:03d}_metadata.json").write_text(json.dumps(meta, indent=2))
    save_sheet(
        [all_marker_crops[m] for m in names],
        [f"{i:02d} {m}" for i, m in enumerate(names)],
        out_dir / f"spot_{first_component:03d}_all_54_markers.png",
        cols=9,
        tile_w=150,
        tile_h=135,
        title_h=28,
    )


def analyze_halo(names):
    if "Hoechst1" not in names:
        return pd.DataFrame()
    idx = names.index("Hoechst1")
    with tifffile.TiffFile(OME_PATH) as tif:
        arr = tif.pages[idx].asarray(maxworkers=MAX_WORKERS)
    roi_path = Path("data/cellpose_roi.json")
    if roi_path.exists():
        roi = json.loads(roi_path.read_text())
        cy = int(roi["y0"] + roi["height"] / 2)
        cx = int(roi["x0"] + roi["width"] / 2)
    else:
        cy, cx = arr.shape[0] // 2, arr.shape[1] // 2
    crop_size = 1024
    y0 = max(0, min(arr.shape[0] - crop_size, cy - crop_size // 2))
    x0 = max(0, min(arr.shape[1] - crop_size, cx - crop_size // 2))
    crop = arr[y0 : y0 + crop_size, x0 : x0 + crop_size]
    sm = filters.gaussian(crop.astype(np.float32), sigma=1)
    nuclei = sm > np.percentile(sm, 97)
    nuclei = morphology.remove_small_objects(nuclei, 20)
    core = morphology.binary_dilation(nuclei, morphology.disk(2))
    ring = morphology.binary_dilation(nuclei, morphology.disk(10)) & ~core
    far = morphology.binary_dilation(nuclei, morphology.disk(26)) & ~morphology.binary_dilation(
        nuclei, morphology.disk(14)
    )
    far &= crop > np.percentile(crop, 5)
    core_med = float(np.median(crop[core])) if core.any() else np.nan
    ring_med = float(np.median(crop[ring])) if ring.any() else np.nan
    far_med = float(np.median(crop[far])) if far.any() else np.nan
    bg_med = float(np.median(crop[crop < np.percentile(crop, 40)]))

    base = Image.fromarray(robust_display(crop, 1, 99.8)).convert("RGB")
    overlay = np.array(base).copy()
    overlay[ring] = [255, 35, 35]
    overlay[core] = [45, 255, 45]
    overlay = Image.fromarray(overlay)
    save_sheet(
        [base, overlay],
        ["Hoechst ROI", "core green / ring red"],
        OUT / "hoechst_nuclear_halo_roi.png",
        cols=2,
        tile_w=430,
        tile_h=420,
        title_h=26,
    )
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
                "ring_to_far_ratio": ring_med / far_med if far_med and far_med > 0 else np.nan,
                "ring_to_low_background_ratio": ring_med / bg_med if bg_med > 0 else np.nan,
                "nucleus_core_pixels": int(core.sum()),
                "ring_pixels": int(ring.sum()),
                "far_pixels": int(far.sum()),
            }
        ]
    )


def main():
    print("load compact data", flush=True)
    inventory = pd.read_csv(OUT / "channel_inventory.csv")
    stats = pd.read_csv(OUT / "channel_statistics.csv")
    stripes = pd.read_csv(OUT / "stripe_metrics.csv")
    npz = np.load(OUT / f"downsampled_mean_f{DOWNSAMPLE}.npz")
    ds_stack = npz["stack"]
    names = [str(x) for x in npz["names"]]

    print("write marker sheets", flush=True)
    thumbnails = [robust_display(np.log1p(ds_stack[i]), 1, 99.8) for i in range(len(names))]
    save_sheet(thumbnails, [f"{i:02d} {n}" for i, n in enumerate(names)], OUT / "all_54_marker_overview_log_scaled.png", cols=9)
    for category, idxs in inventory.groupby("category").groups.items():
        idxs = list(idxs)
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(category))
        save_sheet(
            [thumbnails[i] for i in idxs],
            [f"{i:02d} {names[i]}" for i in idxs],
            OUT / f"group_{safe}.png",
            cols=min(6, max(1, len(idxs))),
        )

    print("write RGB composites", flush=True)
    composites = [
        make_rgb(ds_stack, names, {"Hoechst1": 2, "Cytokeratin": 1, "CD45": 0}),
        make_rgb(ds_stack, names, {"MUC2": 0, "SOX9": 1, "CDX2": 2}),
        make_rgb(ds_stack, names, {"CD3": 0, "CD19": 1, "CD68": 2}),
        make_rgb(ds_stack, names, {"aSMA": 0, "CollagenIV": 1, "CD31": 2}),
    ]
    save_sheet(
        composites,
        [
            "blue Hoechst / green CK / red CD45",
            "red MUC2 / green SOX9 / blue CDX2",
            "red CD3 / green CD19 / blue CD68",
            "red aSMA / green CollagenIV / blue CD31",
        ],
        OUT / "biological_group_rgb_overviews.png",
        cols=2,
        tile_w=430,
        tile_h=410,
        title_h=34,
    )

    print("detect multimarker components", flush=True)
    pos_count, threshold, components = analyze_multimarker(ds_stack, names)
    components.to_csv(OUT / "multimarker_positive_components.csv", index=False)
    save_colormap(pos_count, OUT / "multimarker_positive_count_heatmap.png", vmin=0, vmax=max(5, int(pos_count.max())))

    print("compute correlations", flush=True)
    flat = np.log1p(ds_stack.reshape(ds_stack.shape[0], -1))
    sample_count = min(90000, flat.shape[1])
    sample_idx = np.linspace(0, flat.shape[1] - 1, sample_count).astype(int)
    corr = np.corrcoef(flat[:, sample_idx])
    pd.DataFrame(corr, index=names, columns=names).to_csv(OUT / "downsampled_log_marker_correlation.csv")
    save_colormap(corr, OUT / "downsampled_log_marker_correlation_heatmap.png", vmin=-0.3, vmax=1.0, cmap=cv2.COLORMAP_COOL)

    print("summarize stripes", flush=True)
    stripe_score = (
        stripes["horizontal_profile_hp_mad_ratio"].fillna(0)
        + stripes["vertical_profile_hp_mad_ratio"].fillna(0)
    )
    top_stripe = stripes.assign(stripe_score=stripe_score).sort_values("stripe_score", ascending=False).head(8)
    save_sheet(
        [thumbnails[int(i)] for i in top_stripe["channel_index"]],
        [
            f"{int(row.channel_index):02d} {row.marker}\nH {row.horizontal_profile_hp_mad_ratio:.2f} V {row.vertical_profile_hp_mad_ratio:.2f}"
            for row in top_stripe.itertuples(index=False)
        ],
        OUT / "top_stripe_metric_channels.png",
        cols=4,
    )

    print("analyze halo", flush=True)
    halo = analyze_halo(names)
    halo.to_csv(OUT / "nuclear_halo_metrics.csv", index=False)

    print("write spot crops", flush=True)
    save_spot_crops(names, stats, components)

    summary = {
        "channel_count": len(names),
        "downsample_factor": DOWNSAMPLE,
        "multimarker_component_threshold": int(threshold),
        "max_non_nuclear_positive_markers_in_block": int(pos_count.max()),
        "multimarker_component_count": int(len(components)),
        "top_multimarker_components": components.head(10).to_dict(orient="records"),
        "top_stripe_channels": top_stripe.to_dict(orient="records"),
        "halo_metrics": halo.to_dict(orient="records"),
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2))
    print("done", flush=True)


if __name__ == "__main__":
    main()
