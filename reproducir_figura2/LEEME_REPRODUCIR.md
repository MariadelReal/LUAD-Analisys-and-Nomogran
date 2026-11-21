# Reproducir la Figura 2 y el suplemento de enriquecimiento

Estos dos scripts generan las imágenes incluidas en esta entrega a partir de las tablas finales. Se ejecutan uno por uno en R o RStudio. No sustituyen ni renumeran los pasos 09 y 10 del análisis GEO.

## Ejecutar

Extrae el ZIP completo. Conserva juntos los scripts y las carpetas `input` y `output`. En RStudio, abre la carpeta `reproducir_figuras` y usa Session > Set Working Directory > Choose Directory. También puedes usar `setwd()` con su ruta completa, usando barras `/` en Windows.

```r
# Ejemplo: adapta solamente la ruta a donde extrajiste el paquete.
setwd("C:/Users/delre/OneDrive/LUAD Royal Society/Scripts/Check run/reproducir_figuras")
stopifnot(dir.exists("input/ORA"), dir.exists("input/GSEA"))

# Una sola vez, si faltan estos paquetes:
faltan <- c("pdftools", "png")
faltan <- faltan[!vapply(faltan, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltan)) install.packages(faltan)

# Primero: regenerar ORA desde los 20 TSV suministrados.
source("01_regenerar_ORA.R", encoding = "UTF-8")

# Después: montar Figura 2 y generar el suplemento.
source("02_montar_figura2_y_suplemento.R", encoding = "UTF-8")
```

El primer script requiere únicamente R base con soporte Cairo. El segundo añade `pdftools` y `png`. Ambos se ejecutaron en R 4.5.0 sobre Linux; los registros enviados por la autora corresponden a R 4.5.1 en Windows y se conservan sin modificar en `input/logs`. Los registros de generación de estas figuras están en `output/sessionInfo_plotting.txt` y `output/sessionInfo_assembly.txt`. No se ha ejecutado aquí la réplica completa del análisis biológico en Windows.

## Salidas y alcance

| Archivo en output | Contenido |
|---|---|
| ORA_UP.pdf / png / tiff | Paneles a–b de ORA para la Figura 2 |
| ORA_DOWN.pdf / png / tiff | Paneles a–b de ORA para Online Resource 1 |
| Fig2_RSOS.pdf / png / tiff | Figura 2 completa, paneles a–d |
| ESM_RSOS_actualizado.pdf | Dos páginas: ORA DOWN y GSEA positivo/negativo |
| ORA_plotted_terms.tsv | Los 80 términos dibujados, FDR, Count y posición numérica |
| ORA_run_summary.tsv | Total de términos y número con FDR ≤0.05 para las 20 tablas |
| ORA_input_md5.tsv | Identificación de los TSV utilizados |
| Online_Resource2_plotted_pathways.tsv | Los 40 Hallmark representados en el suplemento |

Se conservan las 20 tablas ORA sin cambiar sus valores; únicamente se retiraron los sufijos de descarga `(1)` y `(2)` de sus nombres. Los gráficos muestran cinco términos significativos por categoría. Los empates se resuelven conservando el orden del TSV. KEGG está incluido en las tablas, pero no en los paneles, tal como indica la leyenda. Una sustitución de términos empatados no significa que haya cambiado la significancia del análisis.

El montaje conserva los paneles GSEA c–d del PDF `input/GSEA/Fig2cd_GSEA.pdf`. Los TSV fgsea permiten reconstruir las barras y comprobar las anotaciones, pero no contienen el ranking completo necesario para volver a calcular las curvas de enriquecimiento. Para regenerar esas curvas desde el análisis, usa el paso 10 original con sus tablas DEA y el GMT de Hallmark. Aquí no se volvió a ejecutar fgsea ni se modificaron las curvas.

Los PDF ORA son vectoriales. El montaje completo integra los paneles a 300 dpi; las exportaciones PNG/TIFF también son de 300 dpi. Diferencias de fuentes o renderizador entre Windows y Linux pueden modificar unos píxeles, sin cambiar los valores representados.

## Versiones y procedencia

`input/code_reference` conserva dos scripts del commit `21f5598a473f01828d96768dacf6ba02a4d164c0` de `MariadelReal/LUAD-Analisys-and-Nomogran`, como evidencia del método; no son pasos adicionales que debas ejecutar para montar las figuras.

- GO: el paso 09 usa `org.Hs.eg.db` y `clusterProfiler::enrichGO`.
- ORA Reactome/KEGG Legacy: el paso 09 declara MSigDB `2026.1.Hs`.
- GSEA Hallmark: el paso 10 utiliza `h.all.v2025.1.Hs.symbols.gmt`.
- `msigdbr 26.1.0` es una versión de paquete registrada; no convierte todos los análisis a MSigDB 2026.1.

Las versiones anteriores están sustentadas por los scripts y por los registros recibidos, respectivamente. Este paquete reproduce las figuras desde los resultados exportados; no contiene las instantáneas originales de todas las colecciones ni demuestra por sí solo la identidad de los archivos usados en cada ejecución histórica.
