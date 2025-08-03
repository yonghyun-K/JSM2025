## ----setup, include=FALSE, warning = FALSE--------------------------------------
library(xtable)
options(xtable.comment = FALSE)
knitr::opts_chunk$set(echo = TRUE, warning = FALSE)

## ----sim0, message=FALSE--------------------------------------------------------
# install.packages("GECal", type = "source", # Uncomment this if using macOS
#     repos = "https://cran.r-project.org")
library(GECal)
# Sampled study variable
y=c(5, 4, 7, 9, 11, 10, 13, 12, 15, 15)
# Sampled auxiliary variables
Xs=cbind(
  c(1,1,1,1,1,1,1,1,1,1),
  c(1,1,1,1,1,0,0,0,0,0),
  c(1,3,5,7,9,6,7,8,9,10)
) 
# vector of population totals
total=c(160,124,700)
# base weights before calibration
d = rep(1, 10)


## ----tab, echo=FALSE, include=FALSE, results="asis"-----------------------------
data <- cbind(1:length(d), d, Xs, y)
colnames(data) <- c("$i$", "$\\delta_i$", "$x_{i0}$", "$x_{i1}$", "$x_{i2}$", "$y_i$")
print(xtable(data,
             digits = 0,
             caption = "Simulation results for point estimation (Bias and RMSE)"),
      include.rownames = FALSE,
      sanitize.colnames.function = identity,
      table.placement = "!htb")

## ----tab2, echo=FALSE, results="asis"-------------------------------------------
data <- cbind(1:length(d), d, Xs, y)
colnames(data) <- c("$i$", "$\\delta_i$", "$x_{i0}$", "$x_{i1}$", "$x_{i2}$", "$y_i$")
cat("\\begin{center}")
cat("\\begin{tabular}{@{}rrrrrr@{}}\n")
cat("\\toprule\n")
cat("$i$ & $\\delta_i$ & $x_{i0}$ & $x_{i1}$ & $x_{i2}$ & $y_i$ \\\\ \\midrule\n")

for (i in 1:10) {
  cat(sprintf("%d & %d & %d & %d & %d & %d \\\\\n", 
              i, d[i], Xs[i, 1], Xs[i, 2], Xs[i, 3], y[i]))
}

cat("11 & 0 & NA & NA & NA & NA \\\\\n")
cat("$\\vdots$ & $\\vdots$ & NA & NA & NA & NA \\\\\n")
cat(sprintf("%d & 0 & NA & NA & NA & NA \\\\\n", total[1]))
cat("\\midrule\n")

cat(sprintf("Total & %d & %d & %d & %d & $\\theta = ?$ \\\\\n",
            length(y), total[1], total[2], total[3]))
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\end{center}")


## ----cal1-----------------------------------------------------------------------
# GEC estimator using ET(exponential tilting) divergence
cal_ET <- GEcalib(~ 0 + Xs, dweight = d, const = total, 
            method = "GEC0", entropy = "ET")
head(cal_ET$w)

GECal::estimate(y ~ 1, calibration = cal_ET)$estimate


## ----cal2-----------------------------------------------------------------------
# GEC estimator using EL(empirical likelihood) divergence
cal_EL <- GEcalib(~ 0 + Xs, dweight = d, const = total, 
        method = "GEC0", entropy = "EL")
head(cal_EL$w)

GECal::estimate(y ~ 1, calibration = cal_EL)$estimate


## ----cal3-----------------------------------------------------------------------
# # GEC estimator using CE(cross entropy or shifted KL) divergence
cal_CE <- GEcalib(~ 0 + Xs, dweight = d, const = total,
        method = "GEC0", entropy = "CE", weight.scale = 2)
# design weights should be greater than 1 in CE
head(cal_CE$w)

GECal::estimate(y ~ 1, calibration = cal_CE)$estimate


## ----tab3, echo=FALSE, results="asis"-------------------------------------------
cat("\\begin{table}[]\n")
cat("    \\centering\n")
cat("    \\begin{tabular}{ccccrrrr}\n")  # Now 7 columns: i, y, x1, x2, x3, ET, EL, CE
cat("        \\hline\n")
cat("        & \\multicolumn{4}{c}{} & \\multicolumn{3}{c}{\\textbf{weights}} \\\\\n")
cat("        \\hline\n")
cat("        $i$ & $y_i$ & $x_{i0}$ & $x_{i1}$ & $x_{i2}$ & \\textbf{ET} & \\textbf{EL} & \\textbf{CE} \\\\\n")
cat("        \\hline\n")

for (i in 1:10) {
  cat(sprintf("        %d & %d & %d & %d & %d & %.2f & %.2f & %.2f \\\\\n",
              i,
              y[i],
              Xs[i, 1], Xs[i, 2], Xs[i, 3],
              cal_ET$w[i], cal_EL$w[i], cal_CE$w[i]))
}

cat("        \\hline\n")
cat("    \\end{tabular}\n")
cat("    \\caption{Comparison of ET, EL, and CE(shifted KL) weights.}\n")
cat("\\end{table}\n")


## ----eval=FALSE, include=FALSE--------------------------------------------------
## load("../data/nhis.Rdata")
## rownames(nhis) <- NULL
## N <- nrow(nhis)
## 
## tab1 = table(nhis$AgeGroup, nhis$REGION1, nhis$SEX)
## tab1 = ifelse(tab1 > 15, 5, round(tab1 / 3))
## 
## n = sum(tab1)
## 
## index = c()
## pi_S = c()
## pi = rep(0, N)
## for(AgeGroup in unique(nhis$AgeGroup)){
##   for(REGION1 in unique(nhis$REGION1)){
##     for(SEX in unique(nhis$SEX)){
##       idx_tmp = which(nhis$AgeGroup == AgeGroup &
##                         nhis$REGION1 == REGION1 &
##                         nhis$SEX == SEX)
##       n_h = tab1[AgeGroup,REGION1,SEX]
##       index = c(index, sample(idx_tmp, n_h, replace = FALSE))
##       pi_S = c(pi_S, rep(n_h / length(idx_tmp), n_h))
##       pi[idx_tmp] = n_h / length(idx_tmp)
##     }
##   }
## }
## 
## delta = as.integer(1:N %in% index)
## nhis.samp <- nhis[index, ]
## 
## nhis <- nhis[,c("AgeGroup", "SEX", "REGION1")]
## 
## save(nhis, nhis.samp, file = "nhis.Rdata")


## -------------------------------------------------------------------------------
load("nhis.Rdata")

head(nhis.samp[,c("AgeGroup", "SEX", "REGION1", "Hemo", "Alcohol", "OralExam")])
fortmp <- formula(~ AgeGroup + SEX + REGION1)
const = colMeans(model.matrix(fortmp, nhis))
# const = colSums(model.matrix(fortmp, nhis)) # Used for population total

calibration <- GEcalib(
  fortmp,
  dweight = rep(1, nrow(nhis.samp)),
  data = nhis.samp,
  const = const,
  entropy = "ET",
  method = "GEC0"
)

estimate(Hemo + Alcohol + OralExam ~ 1, data = nhis.samp, 
         calibration = calibration)$estimate

