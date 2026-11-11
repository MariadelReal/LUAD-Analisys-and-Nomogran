source("helpers.R")
library(glmnet)
# Internal validation conditional on the 30 candidates obtained in GEO.
library(survival)
d <- readRDS("data/clinical_for_modeling.rds")
d$decade <- d$age/10
X <- readRDS("data/X_genes_signature_scaled.rds")
e <- readRDS("data/expr_tumor_SYMBOL_filtered.rds")
stopifnot(!anyDuplicated(e$gene), all(colnames(X) %in% e$gene),
          identical(d$barcode,rownames(X)),all(d$barcode %in% names(e)))
raw <- t(as.matrix(e[match(colnames(X),e$gene),d$barcode,drop=FALSE]))
colnames(raw) <- colnames(X)
stopifnot(all(is.finite(raw)),max(abs(scale(raw)-X))<1e-6)
endpoint <- with(d,Surv(OS.time,OS.event))
clinical_formula <- Surv(OS.time,OS.event) ~ decade + gender + stage + smoking
B <- 200L
ids <- seq_len(nrow(d))
score_C <- function(fit,data) {
  risk <- predict(fit,newdata=data,type="lp")
  stopifnot(all(is.finite(risk)))
  concordance(Surv(OS.time,OS.event) ~ risk,data=data,reverse=TRUE)$concordance
}
one_boot <- function(i) {
  set.seed(73100+i)
  idx <- sample(ids,replace=TRUE)
  people <- unique(idx)
  folds <- sample(rep(1:10,length.out=length(people)))[match(idx,people)]
  train_x <- scale(raw[idx,,drop=FALSE])
  test_x <- scale(raw,center=attr(train_x,"scaled:center"),scale=attr(train_x,"scaled:scale"))
  fit <- cv.glmnet(train_x,endpoint[idx,],family="cox",alpha=1,foldid=folds,standardize=TRUE)
  active <- as.numeric(coef(fit,s="lambda.min"))!=0
  test <- d
  test$score <- drop(predict(fit,test_x,s="lambda.min",type="link"))
  train <- test[idx,,drop=FALSE]
  vars <- c("decade","gender","stage","smoking","score")
  train <- train[complete.cases(train[,vars]),]
  test <- test[complete.cases(test[,vars]),]
  clinical <- coxph(clinical_formula,data=train,singular.ok=FALSE)
  integrated <- if(any(active)) coxph(update(clinical_formula,. ~ . + score),data=train,singular.ok=FALSE) else clinical
  stats <- c(clin_train=score_C(clinical,train),clin_test=score_C(clinical,test),
             int_train=score_C(integrated,train),int_test=score_C(integrated,test))
  stopifnot(all(is.finite(stats)))
  c(stats,setNames(as.integer(active),colnames(raw)))
}
res <- vector("list",B)
for(i in seq_len(B)) {
  # Record numerical failures explicitly. A zero-gene selection is a valid replicate.
  res[[i]] <- tryCatch(withCallingHandlers(one_boot(i),warning=function(w) {
    if (grepl("infinite|converg|numerical",conditionMessage(w),ignore.case=TRUE)) stop(conditionMessage(w))
  }),error=function(e) conditionMessage(e))
  if(i %% 10==0) message(i,"/",B)
}
ok <- vapply(res,is.numeric,logical(1))
save_tab(data.frame(iteration=which(!ok),reason=as.character(unlist(res[!ok]))),"bootstrap_failures")
if(!any(ok)) stop("Ninguna réplica válida. Revisa results/bootstrap_failures.tsv.")
z <- as.data.frame(do.call(rbind,res[ok]))
save_tab(data.frame(iteration=which(ok),z),"bootstrap_replicates")
ap <- read_tab("concordance")
op <- c(mean(z$clin_train-z$clin_test),mean(z$int_train-z$int_test))
stopifnot(identical(ap$model,c("clinical","integrated")))
ap$optimism <- op
ap$corrected <- ap$C-op
delta <- data.frame(model="integrated_minus_clinical",C=diff(ap$C),optimism=diff(op),corrected=diff(ap$corrected))
summary_tab <- rbind(ap,delta)
summary_tab$valid <- sum(ok)
summary_tab$attempted <- B
save_tab(summary_tab,"bootstrap_concordance")
freq <- sort(colMeans(z[colnames(raw)]),decreasing=TRUE)
save_tab(data.frame(gene=names(freq),frequency=unname(freq)),"bootstrap_selection")
panels("Supplement_selection_stability",function(){
  par(mar=c(4,6,3,1))
  barplot(rev(freq),horiz=TRUE,las=1,cex.names=.75,xlim=c(0,1),col="#205776",
          xlab="Fraction of successful bootstrap fits",main="LASSO selection stability")
},7,8)
writeLines(capture.output(sessionInfo()),"results/bootstrap_session.txt")
print(summary_tab)
message("Modelos sin genes en lambda.min: ",sum(rowSums(z[colnames(raw)])==0))
if(!all(ok)) warning("Hubo réplicas fallidas. Revisa sus causas antes de reportar la estimación.")
if(all(ok)) {
  methods <- paste("Internal validation used 200 patient-level bootstrap samples.",
    "In each sample, scaling, ten-fold LASSO tuning, gene selection, and clinical and integrated model fitting were repeated.",
    "Repeated copies of a patient shared a cross-validation fold.",
    "Optimism was estimated as bootstrap-sample concordance minus concordance in the original cohort.")
  result <- sprintf(paste("Across %d bootstrap samples, optimism-corrected concordance was %.3f for",
    "the clinical model and %.3f for the integrated model (paired difference %.3f)."),
    sum(ok),ap$corrected[1],ap$corrected[2],diff(ap$corrected))
  writeLines(c(methods,"",result),"results/Texto_bootstrap_EN.txt")
} else if(file.exists("results/Texto_bootstrap_EN.txt")) file.remove("results/Texto_bootstrap_EN.txt")
