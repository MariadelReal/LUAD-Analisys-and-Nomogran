# Reproduce Figure 3

Run `python montar_figura3.py` from any working directory. Requires Python 3, Pillow and reportlab (`python -m pip install pillow reportlab`). Arial, Liberation Sans or DejaVu Sans is required for panel labels. Results are written to `output/`.

Panel mapping: a = <60 MCODE Cluster 1; b = ≥60 MCODE Cluster 1; c = <60 MCC top 35; d = ≥60 MCC top 35. The supplied PNG exports are composited on white and scaled proportionally; visible network content is not redrawn. KIF4A is already present in the supplied corrected labels. The original Cytoscape sessions and identity audit are provided in Online Resource 3.

TIFF/PNG are 5669 × 6142 pixels at 600 ppi (nominal 240 × 260 mm). The PDF embeds the same raster composition and is not a vector network export. This resolution describes the output canvas; resampling cannot restore detail absent from the original PNGs. If the journal requests a smaller printed width, inspect gene-label readability at that width; a vector re-export from the supplied Cytoscape sessions is preferable to further upsampling.
