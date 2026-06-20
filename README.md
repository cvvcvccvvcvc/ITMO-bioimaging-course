# Bioimaging Project

Проект по анализу мультиплексных CODEX-изображений HuBMAP для ткани толстого кишечника. Основной пайплайн находится в `final_work` и включает извлечение каналов, preprocessing Hoechst, морфологический анализ, сегментацию Cellpose, single-cell quantification и биологическую интерпретацию.

## Структура

- `final_work/00_extract_selected_channels.ipynb` - извлечение выбранных каналов из OME-TIFF.
- `final_work/02_image_preprocessing_hoechst.ipynb` - preprocessing Hoechst.
- `final_work/03_morphological_analysis_nuclei_tissue.ipynb` - морфология ядер и ткани.
- `final_work/04_cell_segmentation_cellpose.ipynb` - сегментация клеток Cellpose.
- `final_work/05_06_single_cell_quantification_analysis.ipynb` - single-cell quantification, clustering и downstream analysis.
- `final_work/marker_expectation/` - QC маркеров и проверка технических артефактов.
- `final_work/biological_interpretation/` - итоговая биологическая интерпретация и компактные фигуры.
- `final_work/requirements.txt` - Python-зависимости.

## Данные

Сырые microscopy data, Zarr/TIFF, модели Cellpose и сгенерированные `outputs` не хранятся в git из-за размера. Для воспроизведения положите локальные данные в:

```text
final_work/data/reg001_expr.ome.tif
final_work/data/channels/*.tif
```

Использованный HuBMAP dataset: <https://portal.hubmapconsortium.org/browse/dataset/a0946b9a99b0940c5e9eb7587deafee5>

## Запуск

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r final_work/requirements.txt
jupyter lab
```

После этого ноутбуки можно запускать по порядку из `final_work`.
