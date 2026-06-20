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

Сырые microscopy data, Zarr/TIFF, модели Cellpose и сгенерированные `outputs` не хранятся в git из-за размера. Для воспроизведения нужен исходный OME-TIFF из HuBMAP:

```text
final_work/data/reg001_expr.ome.tif
```

Каналы `final_work/data/channels/*.tif`, маска `final_work/data/cellpose_hoechst_mask.tif` и папка `final_work/outputs/` создаются ноутбуками заново. Подробная инструкция по получению данных: [docs/DATA.md](docs/DATA.md).

## Запуск

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r final_work/requirements.txt
jupyter lab
```

После этого откройте Jupyter из папки `final_work` и запускайте ноутбуки по порядку.
