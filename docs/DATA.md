# Data Setup

Raw microscopy data are intentionally not committed to this repository. The analysis expects local files under `final_work/data/`, and all derived channel TIFFs, masks, Cellpose models, and `outputs/` can be regenerated.

## Source Dataset

- HuBMAP dataset page: <https://portal.hubmapconsortium.org/browse/dataset/a0946b9a99b0940c5e9eb7587deafee5>
- Dataset DOI: <https://doi.org/10.35079/HBM792.FFJT.499>
- HuBMAP ID: `HBM792.FFJT.499`
- Dataset UUID: `a0946b9a99b0940c5e9eb7587deafee5`
- Globus file browser: <https://app.globus.org/file-manager?origin_id=af603d86-eab9-4eec-bb1d-9d26556741bb&origin_path=%2Fa0946b9a99b0940c5e9eb7587deafee5%2F>

Download the processed CODEX OME-TIFF from the dataset files and place it here:

```text
final_work/data/reg001_expr.ome.tif
```

If the downloaded file has a different name, rename it to `reg001_expr.ome.tif`, or update `ome_path` in `final_work/00_extract_selected_channels.ipynb` and the scripts under `final_work/marker_expectation/`.

## Expected Local Layout

After downloading the source file, the local tree should start like this:

```text
final_work/
  data/
    reg001_expr.ome.tif
  00_extract_selected_channels.ipynb
  02_image_preprocessing_hoechst.ipynb
  03_morphological_analysis_nuclei_tissue.ipynb
  04_cell_segmentation_cellpose.ipynb
  05_06_single_cell_quantification_analysis.ipynb
```

The following files and directories are generated and do not need to be downloaded:

```text
final_work/data/channels/*.tif
final_work/data/preprocessed_hoechst_roi.tif
final_work/data/cellpose_models/
final_work/data/cellpose_hoechst_mask.tif
final_work/data/cellpose_contours_overlay.tif
final_work/outputs/
```

## Reproducing the Analysis

Run from the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r final_work/requirements.txt
cd final_work
mkdir -p data outputs
jupyter lab
```

Then run the notebooks in order:

1. `00_extract_selected_channels.ipynb` reads `data/reg001_expr.ome.tif` and writes `data/channels/*.tif`.
2. `02_image_preprocessing_hoechst.ipynb` writes `data/preprocessed_hoechst_roi.tif`.
3. `03_morphological_analysis_nuclei_tissue.ipynb` performs tissue/nuclei morphology checks.
4. `04_cell_segmentation_cellpose.ipynb` downloads the Cellpose-SAM model into `data/cellpose_models/`, writes `data/cellpose_hoechst_mask.tif`, and saves quick-look images in `outputs/`.
5. `05_06_single_cell_quantification_analysis.ipynb` reads `data/channels/*.tif` plus `data/cellpose_hoechst_mask.tif`, then writes QuPath-compatible GeoJSON files into `outputs/`.

The reports in `final_work/marker_expectation/` and `final_work/biological_interpretation/` document QC and interpretation results. Their scripts also assume the working directory is `final_work`.

## Disk Usage

The source and regenerated local artifacts are large. On the original workstation, `final_work/data/` used about 14 GB and `final_work/outputs/` used about 159 MB. They are ignored by git and can be deleted after analysis; rerun the notebooks to recreate them.
