from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from analyze_codex_markers import (
    DOWNSAMPLE,
    MARKER_META,
    OME_PATH,
    OUT,
    analyze_multimarker,
    analyze_nuclear_halo,
    robust_display,
    save_contact_sheet,
    save_heatmap,
    save_overview_composites,
    save_spot_crops,
)


def main():
    print("load tables", flush=True)
    stats = pd.read_csv(OUT / "channel_statistics.csv")
    stripes = pd.read_csv(OUT / "stripe_metrics.csv")
    inventory = pd.read_csv(OUT / "channel_inventory.csv")
    npz = np.load(OUT / f"downsampled_mean_f{DOWNSAMPLE}.npz")
    ds_stack = npz["stack"]
    names = [str(x) for x in npz["names"]]

    print("make thumbnails/contact sheets", flush=True)
    thumbnails = [
        (robust_display(np.log1p(ds_stack[i]), 1, 99.8) * 255).astype(np.uint8)
        for i in range(len(names))
    ]
    titles = [f"{i:02d} {n}" for i, n in enumerate(names)]
    save_contact_sheet(thumbnails, titles, OUT / "all_54_marker_overview_log_scaled.png", cols=9)

    for category, idxs in inventory.groupby("category").groups.items():
        idxs = list(idxs)
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(category))
        save_contact_sheet(
            [thumbnails[i] for i in idxs],
            [titles[i] for i in idxs],
            OUT / f"group_{safe}.png",
            cols=min(6, max(1, len(idxs))),
        )

    print("make biological composites", flush=True)
    save_overview_composites(ds_stack, names, OUT / "biological_group_rgb_overviews.png")

    print("detect multimarker components", flush=True)
    pos_count, positive_threshold, components = analyze_multimarker(ds_stack, names)
    components.to_csv(OUT / "multimarker_positive_components.csv", index=False)
    save_heatmap(
        pos_count,
        f"Non-nuclear markers above per-marker p99.5; components threshold >= {positive_threshold}",
        OUT / "multimarker_positive_count_heatmap.png",
        cmap="magma",
        vmax=max(5, int(pos_count.max())),
    )

    print("compute marker correlations", flush=True)
    flat = np.log1p(ds_stack.reshape(ds_stack.shape[0], -1))
    # A deterministic spatial sample is enough for marker-level correlation and keeps memory low.
    sample_count = min(90000, flat.shape[1])
    sample_idx = np.linspace(0, flat.shape[1] - 1, sample_count).astype(int)
    corr = np.corrcoef(flat[:, sample_idx])
    corr_df = pd.DataFrame(corr, index=names, columns=names)
    corr_df.to_csv(OUT / "downsampled_log_marker_correlation.csv")

    import matplotlib.pyplot as plt

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

    print("summarize stripe metrics", flush=True)
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

    print("analyze nuclear halo", flush=True)
    halo = analyze_nuclear_halo(OME_PATH, names, stats, OUT / "hoechst_nuclear_halo_roi.png")
    halo.to_csv(OUT / "nuclear_halo_metrics.csv", index=False)

    print("save multimarker spot crops", flush=True)
    save_spot_crops(OME_PATH, names, stats, components, OUT / "multimarker_spot_crops")

    print("write summary", flush=True)
    summary = {
        "ome_path": str(OME_PATH),
        "channel_count": len(names),
        "downsample_factor": DOWNSAMPLE,
        "multimarker_component_threshold": int(positive_threshold),
        "max_non_nuclear_positive_markers_in_block": int(pos_count.max()),
        "multimarker_component_count": int(len(components)),
        "top_multimarker_components": components.head(10).to_dict(orient="records"),
        "top_stripe_channels": top_stripe.to_dict(orient="records"),
        "halo_metrics": halo.to_dict(orient="records"),
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2))
    print("Postprocessing wrote outputs to", OUT.resolve())


if __name__ == "__main__":
    main()
