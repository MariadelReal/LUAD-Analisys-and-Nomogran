source("helpers.R")
library(rms)
x <- read_tab("Nomogram_input")
x$gender <- factor(x$gender, levels=c("Female","Male"))
x$stage <- factor(x$stage, levels=c("I","II","III","IV"))
x$smoking <- factor(x$smoking, levels=c("Never","Former","Current","Unknown"))
vars <- c("OS.time","OS.event","age","gender","stage","smoking","LP_pos")
d <- x[complete.cases(x[,vars]), vars]
label(d$LP_pos) <- sprintf("%d-gene signature score", nrow(read_tab("lambda.min")))
label(d$age) <- "Age (years)"
label(d$gender) <- "Sex"
label(d$stage) <- "Stage"
label(d$smoking) <- "Smoking status"
dd <- datadist(d); options(datadist="dd")
full <- cph(Surv(OS.time,OS.event) ~ LP_pos + gender + stage + smoking + age,
            data=d, method="efron", x=TRUE, y=TRUE, surv=TRUE, time.inc=1825)
clinical <- coxph(Surv(OS.time,OS.event) ~ gender + stage + smoking + age, data=d)
stopifnot(!full$fail, nrow(d)==length(predict(full,type="lp")))
h <- c(365,1095,1825)
risk <- as.numeric(predict(full,type="lp"))
S <- Survival(full)
nom <- nomogram(full,fun=list(function(lp) S(365,lp),function(lp) S(1095,lp),
  function(lp) S(1825,lp)),funlabel=c("1-year OS","3-year OS","5-year OS"),lp=FALSE)
panels("Fig4c_Nomogram",function() plot(nom,xfrac=.28),11,7)
set.seed(2025)
cal <- calibrate(full,cmethod="KM",method="boot",u=1825,m=80,B=200)
panels("Fig4d_Calibration_5y",function() plot(cal,xlab="Predicted 5-year OS",ylab="Observed 5-year OS"),6.5,5.5)
saveRDS(list(model=full,clinical=clinical,calibration=cal),"results/nomogram_fit.rds")
writeLines(capture.output(sessionInfo()), "results/R_session_nomogram.txt")
