# Reproduce Figure 4 assembly

Run `python montar_figura4.py`. Requires Pillow, numpy, PyMuPDF and Arial or Liberation Sans. Outputs go to `output/` (300 dpi TIFF/PNG and raster PDF).

This reproduces the submitted arrangement from supplied model-output PNGs. It separates the crowded 1-year survival labels in the nomogram and relocates the calibration footer. Assertions protect the nomogram outside its label rectangle and preserve every plotted calibration pixel above its footer. Curves, markers, confidence limits and reported statistics are not refitted or changed. For new model fitting, use the numbered TCGA scripts in the repository; those scripts generate the component plots.

Verification in the preparation environment: the regenerated TIFF was identical pixel-for-pixel to the submitted Fig4_RSOS.tiff.
