source("helpers.R")
p <- read_tab("predictions")
# Reverse Kaplan-Meier. At tied times deaths precede censoring.
G <- 1
w <- numeric(nrow(p))
for (day in sort(unique(p$time))) {
  at <- p$time == day
  w[at] <- 1/G
  risk_for_censoring <- sum(p$time>=day)-sum(at & p$event==1)
  if (risk_for_censoring>0) G <- G*(1-sum(at & p$event==0)/risk_for_censoring)
}
ord <- order(p$risk,decreasing=TRUE)
last <- !duplicated(p$risk[ord],fromLast=TRUE)
roc <- function(day) {
  a <- ifelse(p$event==1 & p$time<day,w,0)
  b <- as.numeric(p$time>day)
  stopifnot(all(is.finite(a)),sum(a)>0,sum(b)>0)
  tp <- c(0,cumsum(a[ord])[last]/sum(a))
  fp <- c(0,cumsum(b[ord])[last]/sum(b))
  list(tp=tp,fp=fp,AUC=sum(diff(fp)*(head(tp,-1)+tail(tp,-1))/2))
}
days <- c(365,1095,1825)
curves <- lapply(days,roc)
AUC <- vapply(curves,function(x) x$AUC,numeric(1))
save_tab(data.frame(days,AUC),"AUC")
panels("Fig4e_Time_ROC",function() {
  cols <- c("#145B75","#C05B31","#397848")
  plot(0:1,0:1,type="n",xlab="1 - specificity",ylab="Sensitivity")
  abline(0,1,lty=3,col="grey50")
  for(i in seq_along(curves)) lines(curves[[i]]$fp,curves[[i]]$tp,col=cols[i],lwd=2)
  legend("bottomright",sprintf("%d y   AUC %.3f",days/365,AUC),col=cols,lty=1,bty="n")
},6.5,5.5)
print(read_tab("AUC"))
