if (!file.exists("TCGA.Rproj")) stop("Abre TCGA.Rproj antes de ejecutar el paso.")
base_dir <- getwd()
result_dir <- file.path(base_dir,"results")
dir.create(result_dir,showWarnings=FALSE)
## Modelo clínico y modelo integrado en los mismos pacientes.
library(survival)
d <- read.delim(file.path(result_dir,"Nomogram_input.tsv"))
stopifnot(!anyDuplicated(d$barcode), all(d$OS.event %in% 0:1), all(d$OS.time > 0))
d <- d[complete.cases(d[, c("OS.time","OS.event","LP_raw","age","gender","stage","smoking")]), ]
stopifnot(all(d$gender %in% c("Female","Male")), all(d$stage %in% c("I","II","III","IV")),
          all(d$smoking %in% c("Never","Former","Current","Unknown")))
d$gender <- factor(d$gender,c("Female","Male"))
d$stage <- factor(d$stage,c("I","II","III","IV"))
d$smoking <- factor(d$smoking,c("Never","Former","Current","Unknown"))
score_sd <- sd(d$LP_raw)
d$signature <- as.numeric(scale(d$LP_raw))
d$decade <- d$age/10
fit <- coxph(Surv(OS.time,OS.event) ~ signature + decade + gender + stage + smoking,
             data=d, ties="efron", x=TRUE)
s <- summary(fit)
tab <- data.frame(term=rownames(s$coefficients), HR=s$conf.int[,1],
                  lower=s$conf.int[,3], upper=s$conf.int[,4], p=s$coefficients[,5])
write.table(tab,file.path(result_dir,"adjusted_HR.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
ph <- cox.zph(fit)
write.table(ph$table,file.path(result_dir,"PH_tests.tsv"),sep="\t",col.names=NA)
capture.output(print(ph),print(s),cat("Signature SD:",score_sd,"\n"),sessionInfo(),file=file.path(result_dir,"model_summary.txt"))
pdf(file.path(result_dir,"PH_diagnostics.pdf"),width=7,height=5)
plot(ph);dev.off()
print(tab); print(ph); print(score_sd)
clin <- coxph(Surv(OS.time,OS.event) ~ decade + gender + stage + smoking,data=d)
write.table(data.frame(model=c("clinical","integrated"),C=c(concordance(clin)$concordance,concordance(fit)$concordance)), file.path(result_dir,"concordance.tsv"),sep="\t",row.names=FALSE)
write.table(data.frame(barcode=d$barcode,time=d$OS.time,event=d$OS.event,risk=predict(fit),signature=d$signature), file.path(result_dir,"predictions.tsv"),sep="\t",row.names=FALSE)
saveRDS(fit,file.path(result_dir,"cox_fit.rds"))

pdf(file.path(result_dir,"Forest_signature_adjusted.pdf"),width=8,height=5)
y <- rev(seq_len(nrow(tab)))
par(mar=c(5,10,3,2))
plot(tab$HR,y,log="x",xlim=range(tab$lower,tab$upper),yaxt="n",ylab="",xlab="Adjusted hazard ratio (95% CI)",pch=19)
axis(2,at=y,labels=tab$term,las=1)
abline(v=1,lty=2,col="grey60")
segments(tab$lower,y,tab$upper,y)
dev.off()

# Auditar si los pacientes y covariables de la corrida coinciden con la réplica archivada.
current <- read.delim(file.path(result_dir,"Nomogram_input.tsv"),na.strings=c("","NA"))
ref <- read.delim("referencia/clinical.tsv",na.strings=c("","NA"))
comp <- merge(current,ref,by="barcode",all=TRUE,suffixes=c("_run","_ref"))
cols <- c("OS.time","OS.event","age","gender","stage","smoking")
different <- function(a,b) {
  both <- !is.na(a) & !is.na(b)
  out <- xor(is.na(a),is.na(b))
  out[both] <- if(is.numeric(a) && is.numeric(b)) abs(a[both]-b[both])>1e-6 else as.character(a[both])!=as.character(b[both])
  out
}
delta <- do.call(rbind,lapply(cols,function(nm) {
  a <- comp[[paste0(nm,"_run")]]
  b <- comp[[paste0(nm,"_ref")]]
  i <- which(different(a,b))
  data.frame(barcode=comp$barcode[i],field=rep(nm,length(i)),run=as.character(a[i]),reference=as.character(b[i]))
}))
write.table(delta,file.path(result_dir,"clinical_differences_vs_reference.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
if(nrow(delta)) warning("Hay diferencias clínicas con la referencia; revisa clinical_differences_vs_reference.tsv antes de usar las cifras del manuscrito.")
