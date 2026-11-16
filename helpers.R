if (!file.exists("TCGA.Rproj")) stop("Abre TCGA.Rproj.")
library(survival)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
read_tab <- function(name, folder = "results") {
  read.delim(file.path(folder, paste0(name, ".tsv")), check.names = FALSE, na.strings = c("", "NA"))
}
save_tab <- function(x, name) {
  write.table(x, file.path("results", paste0(name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}
panels <- function(name, draw, width, height) {
  for (ext in c("pdf","tiff","png")) {
    path <- file.path("figures", paste0(name, ".", ext))
    if (ext == "pdf") pdf(path, width, height, useDingbats = FALSE) else if (ext == "png")
      png(path, width=width,height=height,units="in",res=200) else
      tiff(path, width=width, height=height, units="in", res=300, compression="lzw")
    tryCatch(draw(), finally = dev.off())
  }
}
