# Ejecutar desde la carpeta reproducir_figuras: source("01_regenerar_ORA.R")
# Solo representa los TSV existentes: NO recalcula ORA, FDR ni el universo.
# R base; no requiere instalar paquetes. Conserva los 20 TSV de entrada.
local({
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  inp <- file.path(root, "input", "ORA")
  out <- file.path(root, "output")
  if (!dir.exists(inp)) stop("Abre reproducir_figuras como directorio de trabajo.")
  dir.create(out, showWarnings = FALSE)
  cats <- c("GO_BP", "GO_CC", "GO_MF", "PATHWAY")
  cols <- c(GO_BP="#D76542", GO_CC="#669C35", GO_MF="#3C8EAE", PATHWAY="#8B70B3")
  short <- c(GO_BP="GO BP",GO_CC="GO CC",GO_MF="GO MF",PATHWAY="Reactome")
  all_rows <- list(); summary <- list(); selections <- list()
  clean_label <- function(x) {
    x <- sub("^REACTOME_", "", x)
    x <- gsub("_", " ", x)
    x <- tolower(x)
    paste0(toupper(substr(x,1,1)),substring(x,2))
  }
  for (age in c("young","old")) for (direction in c("UP","DOWN")) {
    key <- paste(age,direction,sep="_"); chosen <- list()
    for (category in c(cats,"KEGG")) {
      path <- file.path(inp,paste0(key,"_ORA_common_",category,".tsv"))
      if (!file.exists(path)) stop("Falta: ",path)
      d <- read.delim(path,check.names=FALSE,stringsAsFactors=FALSE)
      stopifnot(all(c("ID","Description","p.adjust","Count") %in% names(d)),
                !anyDuplicated(d$ID), !anyNA(d$p.adjust),
                all(d$p.adjust > 0 & d$p.adjust <= 1),
                !anyNA(d$Count), all(d$Count >= 0))
      d$age_stratum <- age; d$direction <- direction; d$category <- category
      d$source_file <- basename(path); d$source_row <- seq_len(nrow(d))
      summary[[length(summary)+1L]] <- data.frame(age_stratum=age,direction=direction,
          category=category,n_tested=nrow(d),n_FDR_le_005=sum(d$p.adjust<=.05))
      all_rows[[length(all_rows)+1L]] <- d
      if (category %in% cats) {
        s <- d[d$p.adjust<=.05,,drop=FALSE]
        # FDR ascendente; empates conservan el orden original del TSV.
        s <- head(s[order(s$p.adjust,s$source_row),,drop=FALSE],5L)
        if (category=="PATHWAY" && any(!grepl("^REACTOME_",s$ID)))
          stop("PATHWAY contiene terminos no Reactome: revisar la leyenda.")
        s$minus_log10_FDR <- -log10(s$p.adjust)
        s$label <- clean_label(s$Description)
        chosen[[category]] <- s
      }
    }
    s <- do.call(rbind,chosen)
    s <- s[order(-s$minus_log10_FDR,match(s$category,cats),s$source_row),]
    rownames(s) <- NULL; s$plot_row <- seq_len(nrow(s))
    selections[[key]] <- s
  }
  write.table(do.call(rbind,summary),file.path(out,"ORA_run_summary.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  audit <- do.call(rbind,selections);rownames(audit)<-NULL
  write.table(audit,file.path(out,"ORA_plotted_terms.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  files <- sort(list.files(inp,pattern="\\.tsv$",full.names=TRUE))
  write.table(data.frame(file=basename(files),md5=unname(tools::md5sum(files))),
    file.path(out,"ORA_input_md5.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  plot_panel <- function(d,letter,age,direction,xmax) {
    par(mar=c(4.8,12.7,3.3,2.6),mgp=c(2.1,.6,0),xaxs="i",yaxs="i",family="sans")
    labels <- vapply(d$label,function(x)paste(strwrap(x,width=34),collapse="\n"),character(1))
    nlines <- lengths(strsplit(labels,"\n",fixed=TRUE))
    rowheight <- .8*nlines+.6
    totalheight <- sum(rowheight)
    yy <- totalheight-cumsum(rowheight)+rowheight/2
    plot(NA,xlim=c(0,xmax),ylim=c(0,totalheight),axes=FALSE,xlab="",ylab="")
    ticks <- pretty(c(0,xmax),n=4);ticks <- ticks[ticks>=0 & ticks<=xmax]
    abline(v=ticks,col="#E6E6E6",lwd=.65)
    abline(v=-log10(.05),col="#777777",lty=2,lwd=.8)
    segments(0,yy,d$minus_log10_FDR,yy,col=cols[d$category],lwd=1.3)
    points(d$minus_log10_FDR,yy,pch=16,col=cols[d$category],cex=.72)
    text(d$minus_log10_FDR+xmax*.025,yy,paste0("n=",d$Count),pos=4,offset=0,cex=.70)
    text(-xmax*.04,yy,labels,adj=1,xpd=NA,cex=.78)
    axis(1,at=ticks,cex.axis=.82,tck=-.015,lwd=.65)
    mtext(expression(-log[10](FDR)),side=1,line=2,cex=.9)
    title(main=paste0(letter,"   ",if(age=="young") "<60" else "\u226560"," years | ",direction),
          adj=0,line=1.4,cex.main=1.15)
    legend("bottom",inset=c(0,-.135),legend=unname(short),col=cols,pch=16,
           ncol=2,cex=.72,bty="n",xpd=NA)
  }
  draw_pair <- function(direction) {
    layout(matrix(1:2,nrow=1));par(cex=1)
    mx <- max(vapply(selections[paste(c("young","old"),direction,sep="_")],
                     function(d)max(d$minus_log10_FDR),numeric(1)))
    xmax <- ceiling(mx*1.14)
    for (j in 1:2) plot_panel(selections[[paste(c("young","old")[j],direction,sep="_")]],
       letters[j],c("young","old")[j],direction,xmax)
  }
  for(direction in c("UP","DOWN")) {
    cairo_pdf(file.path(out,paste0("ORA_",direction,".pdf")),width=11,height=7.6,pointsize=11)
    draw_pair(direction);dev.off()
    png(file.path(out,paste0("ORA_",direction,".png")),width=3300,height=2280,res=300,pointsize=11,type="cairo")
    draw_pair(direction);dev.off()
    tiff(file.path(out,paste0("ORA_",direction,".tiff")),width=3300,height=2280,res=300,
         pointsize=11,type="cairo",compression="lzw")
    draw_pair(direction);dev.off()
  }
  capture.output(sessionInfo(),file=file.path(out,"sessionInfo_plotting.txt"))
  message("[OK] ORA: 4 paneles, ",nrow(audit)," terminos; PDF, PNG, TIFF y tablas en output/.")
})
