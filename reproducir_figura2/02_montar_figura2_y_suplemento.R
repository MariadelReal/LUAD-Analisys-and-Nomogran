# Ejecutar DESPUES de 01_regenerar_ORA.R desde reproducir_figuras.
# Instalar una sola vez, si faltan: install.packages(c("pdftools", "png"))
# Las curvas GSEA se conservan desde el PDF fuente adjunto: NO se reconstruyen
# desde el TSV fgsea (ese TSV no contiene el ranking ni la curva completa).
local({
  for (p in c("pdftools","png")) if(!requireNamespace(p,quietly=TRUE))
    stop("Instala el paquete: install.packages('",p,"')")
  root<-normalizePath(getwd(),winslash="/",mustWork=TRUE)
  out<-file.path(root,"output"); inp<-file.path(root,"input","GSEA")
  required<-c(file.path(out,"ORA_UP.png"),file.path(out,"ORA_DOWN.png"),file.path(inp,"Fig2cd_GSEA.pdf"))
  if(!all(file.exists(required)))stop("Faltan entradas; corre primero el script 01 desde reproducir_figuras.")
  tmp<-tempfile();dir.create(tmp);on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
  rendered<-pdftools::pdf_convert(required[3],format="png",pages=1,dpi=300,
    filenames=file.path(tmp,"page_%d.%s"),verbose=FALSE)
  top<-png::readPNG(required[1]);bottom<-png::readPNG(rendered[1])
  width<-11; htop<-width*dim(top)[1]/dim(top)[2]
  hbottom<-width*dim(bottom)[1]/dim(bottom)[2];gap<-.12;height<-htop+gap+hbottom
  draw_main<-function(){
    grid::grid.newpage()
    grid::grid.raster(top,x=.5,y=1-htop/(2*height),width=1,height=htop/height,interpolate=FALSE)
    grid::grid.raster(bottom,x=.5,y=hbottom/(2*height),width=1,height=hbottom/height,interpolate=FALSE)
  }
  # PDF del montaje a 300 dpi; ORA_UP.pdf conserva por separado los vectores.
  cairo_pdf(file.path(out,"Fig2_RSOS.pdf"),width=width,height=height);draw_main();dev.off()
  png(file.path(out,"Fig2_RSOS.png"),width=round(width*300),height=round(height*300),res=300,type="cairo")
  draw_main();dev.off()
  tiff(file.path(out,"Fig2_RSOS.tiff"),width=round(width*300),height=round(height*300),res=300,type="cairo",compression="lzw")
  draw_main();dev.off()

  # Dos paginas suplementarias: ORA DOWN (a,b), Hallmark GSEA (a-d).
  cap1<-paste("Online Resource 1. ORA of genes commonly downregulated in tumours across GSE75037 and GSE115002.",
    "(a) <60 years; (b) \u226560 years. The five terms with the smallest FDR \u2264 0.05 per category are shown; ties retain the source-table order.",
    "The x-axis shows -log10(FDR); the dashed line marks FDR = 0.05, and n is the number of overlapping genes.",
    "Categories: GO Biological Process, Cellular Component, Molecular Function and Reactome. KEGG results are supplied as tables.")
  cap2<-paste("Online Resource 2. Hallmark preranked GSEA. The ten largest positive or most negative NES among sets with FDR \u2264 0.05 are shown.",
    "(a,b) <60 years; (c,d) \u226560 years. Bar length shows |NES|; labels retain the sign.",
    "Positive/negative enrichment refers to higher/lower positions in the combined tumour-normal ranking.",
    "Selections use the supplied fgsea tables; no enrichment analysis was rerun.")
  cairo_pdf(file.path(out,"ESM_RSOS_actualizado.pdf"),width=11,height=8.7,pointsize=10)
  grid::grid.newpage(); down<-png::readPNG(required[2])
  grid::grid.raster(down,x=.5,y=1-7.6/(2*8.7),width=1,height=7.6/8.7,interpolate=FALSE)
  grid::grid.text(paste(strwrap(cap1,width=155),collapse="\n"),x=.035,y=.105,
    just=c("left","top"),gp=grid::gpar(fontsize=8.5,lineheight=1.15))
  par(mfrow=c(2,2),oma=c(5.3,0,0,0),mar=c(3.4,12.3,3,2.2),mgp=c(2,.5,0),family="sans",cex=1)
  selected<-list()
  for(age in c("young","old")) {
    d<-read.delim(file.path(inp,paste0("fgsea_H_common_",age,".tsv")),stringsAsFactors=FALSE)
    stopifnot(all(c("pathway","NES","padj") %in% names(d)))
    for(direction in c("positive","negative")) {
      s<-d[!is.na(d$padj)&d$padj<=.05 & if(direction=="positive") d$NES>0 else d$NES<0,,drop=FALSE]
      s<-head(s[order(if(direction=="positive") -s$NES else s$NES,s$pathway),],10)
      s$age_stratum<-age;s$direction<-direction;selected[[length(selected)+1]]<-s
      lab<-gsub("_"," ",sub("^HALLMARK_","",s$pathway))
      lab<-vapply(lab,function(x)paste(strwrap(x,width=27),collapse="\n"),character(1))
      yy<-barplot(rev(abs(s$NES)),horiz=TRUE,col=if(direction=="positive") "#C6373C" else "#22629A",
        border=NA,xlim=c(0,3.7),axes=FALSE,space=.38)
      axis(1,cex.axis=.85,tck=-.02);mtext("|Normalised enrichment score|",side=1,line=2,cex=.85)
      text(-.1,yy,rev(lab),adj=1,xpd=NA,cex=.75)
      text(rev(abs(s$NES))+.05,yy,sprintf("%+.2f",rev(s$NES)),adj=0,cex=.8)
      title(main=paste0(letters[length(selected)],"   ",if(age=="young")"<60" else "\u226560"," years | ",direction),adj=0,cex.main=1.05)
    }
  }
  mtext(paste(strwrap(cap2,width=150),collapse="\n"),side=1,outer=TRUE,line=1,adj=.03,cex=.82)
  dev.off()
  write.table(do.call(rbind,selected),file.path(out,"Online_Resource2_plotted_pathways.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  capture.output(sessionInfo(),file=file.path(out,"sessionInfo_assembly.txt"))
  message("[OK] Figura 2 completa (PDF/PNG/TIFF), suplemento de dos paginas y seleccion GSEA en output/.")
})
