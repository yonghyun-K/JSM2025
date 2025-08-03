## ----setup, include=FALSE, warning = FALSE--------------------------------------
knitr::opts_chunk$set(warning = FALSE, echo = TRUE)
SIMNUM = 5
set.seed(2025)


## ----sim0, message=FALSE--------------------------------------------------------
# install all the necessary libraries using install.packages(...)
library(CVXR)
library(ggplot2)
library(GGally)

N = 1000 # N is the population size
p = 3; # p is the number of covariates
x = matrix(rnorm(N * p, 2, 1), nc= p) # auxiliary variables
mu = 1 + x[,1] + 2 * x[,2] # E(y | x)
e = rnorm(N, 0, 1) # error
y = mu + e # study variable of interest
pi = 1 / (1 + exp(-(-0.5 - 0.25 * x[,2] + 0.5 * x[,3]))) # 1st order inclusion prob.
# pi does not depend on y conditioning on x -> Missing at Random(MAR) 
delta = rbinom(N, 1, pi) # Sample indicator variable
x_OR = cbind(1, x[,c(1,2)]); x_RP = cbind(1, x[,c(2,3)])

Index_S = (delta == 1)
y_S = y[Index_S]
x_OR_S = x_OR[Index_S,]; x_RP_S = x_RP[Index_S,]


## ----step1----------------------------------------------------------------------
PSmodel = glm(delta ~ 0 + x_RP, family = binomial)
PSmodel$coefficients # Estimated PS model parameters


## ----pihat----------------------------------------------------------------------
pihat = predict.glm(PSmodel, type = "response") # Estimated propensity score
dhat = 1 / pihat; dhat_S = dhat[Index_S]; pihat_S = pihat[Index_S]


## ----step2----------------------------------------------------------------------
w = CVXR::Variable(length(y_S))

# Option 1 ####
# Minimize \sum (\omega_i - \hat d_i)^2 
# s.t. \sum \delta_i \omega_i x_i = \sum x_i
###############

constraints <- list(t(x_OR_S) %*% w == colSums(x_OR))
Phi_R <- CVXR::Minimize(sum((w - dhat_S)^2))

prob <- CVXR::Problem(Phi_R, constraints)
res <- CVXR::solve(prob)
w_S = drop(res$getValue(w))


## ----calibrate------------------------------------------------------------------
sum(w_S * y_S) # Estimated population total of y

betahat = solve(t(x_OR_S) %*% (x_OR_S), t(x_OR_S) %*% (y_S))
sum(x_OR %*% betahat) + sum(dhat_S * (y_S - drop(x_OR_S %*% betahat)))


## ----plot, out.width='50%', fig.align='center', fig.cap = "A scatter plot matrix of \\(\\pi_i^{-1}\\), \\(\\hat\\pi_i^{-1}\\), and \\(\\hat\\omega_i\\)"----
GGally::ggpairs(data.frame(true.inv.prob = 1 / pi[Index_S], 
    fitted.inv.prob = dhat_S, calib.weight = w_S))


## ----opt2, echo=FALSE, eval=FALSE-----------------------------------------------
## # Option 2 ####
## # Minimize \sum \hat d_i (\omega_i / \hat d_i - 1)^2
## # s.t. \sum \delta_i \omega_i x_i = \sum x_i
## ###############
## 
## constraints <- list(t(x_OR_S) %*% w == colSums(x_OR))
## Phi_R <- CVXR::Minimize(sum(dhat_S * (w * pihat_S - 1)^2))


## ----opt3, echo=FALSE, eval=FALSE-----------------------------------------------
## # Option 3 ####
## # Minimize \sum \omega_i^2
## # s.t. \sum \delta_i \omega_i (x_i, \hat d_i) = \sum (x_i, \hat d_i)
## ###############
## 
## x_OR = cbind(x_OR, dhat); x_OR_S = cbind(x_OR_S, dhat_S)
## constraints <- list(t(x_OR_S) %*% w == colSums(x_OR))
## Phi_R <- CVXR::Minimize(sum(w^2))


## ----step3----------------------------------------------------------------------
sum(w_S * (w_S - 1) * (y_S - drop(x_OR_S %*% betahat))^2) # Estimated variance


## ----MCsim, echo=FALSE, warning = FALSE-----------------------------------------
library(statmod)
library(xtable)
options(xtable.comment = FALSE)
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(doRNG))

# Simulation setup
timenow1 = Sys.time()
timenow0 = gsub(' ', '_', gsub('[-:]', '', timenow1))
timenow = paste(timenow0, ".txt", sep = "")

# setup parallel backend to use many processors
cores = min(detectCores() - 3, 101)
print(paste("# of cores in the MC simulation =", cores))
# cl <- makeCluster(cores, outfile = timenow) #not to overload your computer
cl <- makeCluster(cores)
registerDoParallel(cl)
registerDoRNG(seed = 11)

modelcases = expand.grid(c(T,F),c(T,F))
rnames <- apply(apply(modelcases, 2, ifelse, "C", "M"), 1, paste, collapse = "")

gh <- gauss.quad(100, "hermite")
nodes <- sqrt(2) * gh$nodes + 2
weights <- gh$weights / sqrt(base::pi)

grid <- setNames(expand.grid(nodes, nodes), c("x1", "x2"))
w_grid <- expand.grid(weights, weights)

final_res <- vector("list", 4)
theta0_vec = NULL; Epi_vec = NULL

for(cnt in 1:nrow(modelcases)){
  RM = modelcases[cnt, 1]
  OM = modelcases[cnt, 2]
  
  set.seed(11)
  N = 1000 # N is the population size
  p = 3; # p is the number of covariates
  x = matrix(rnorm(N * p, 2, 1), nc= p) # auxiliary variables
  x_OR = x[,c(1,2)]; x_PS = x[,c(2,3)]
  if(OM) f1 <- function(x1, x2) 1 + x1 + 2 * x2
  if(!OM) f1 <- function(x1, x2) 1 + sin(x1) + 0.5 * x2^2
  mu = f1(x_OR[,1], x_OR[,2])
  theta0 <- sum(do.call(f1, grid) * Reduce(`*`, w_grid)) * N
  if(RM) theta0_vec <- c(theta0_vec, theta0)
  e = rnorm(N, 0, 1) # error
  y = mu + e # study variable of interest
  theta = sum(y)
  # if(RM) f2 <- function(x1, x2) 1 / (1 + exp(-(-1 - 0.5 * x1 + x2)))
  # if(!RM) f2 <- function(x1, x2) 1 / (1 + exp(-(-1.25 + 0.3 * x1 * x2 + 0.05 * (x2 - 1)^2)))
  
  if(RM) f2 <- function(x1, x2) 1 / (1 + exp(-(-0.5 - 0.25 * x1 + 0.5 * x2)))
  if(!RM) f2 <- function(x1, x2) 1 / (1 + exp(-(-1 + 0.1 * x1 * x2 + 0.3 * (x2 - 1)^2)))
  # if(!RM) f2 <- function(x1, x2) 1 / (1 + exp(-(-1 + 0.25 * x1 * x2 + 0.05 * (x2 - 1)^2)))
  # if(!RM) f2 <- function(x1, x2) 1 / (1 + exp(-(-0.5 + 0.1 * x1 * x2 + 0.05 * (x2 - 1)^2)))
  # if(!RM) f2 <- function(x1, x2) {
  #   tmp = 1 / (1 + exp(-(-1.25 + 0.3 * x1 * x2 + 0.05 * (x2 - 1)^2)))
  #   ifelse(tmp < 0.9, tmp, 0.9)
  # }
  # if(!RM) f2 <- function(x1, x2) 1 / (1 + exp(-(-1 + ifelse(x1 < 1, -1, 1) + 0.1 * (x2 - 1)^2)))
  
  # if(RM) f2 <- function(x1, x2) 1 / (1 + exp(-(-2.75 - 0.5 * x1 + x2)))
  # if(!RM) f2 <- function(x1, x2) 1 / (1 + exp(-(-3 + 0.3 * x1 * x2 + 0.05 * (x2 - 1)^2)))
  pi = f2(x_PS[,1], x_PS[,2])
  Epi <- sum(do.call(f2, grid) * Reduce(`*`, w_grid)) # E(\pi)
  if(OM) Epi_vec <- c(Epi_vec, Epi)
  # pi does not depend on y conditioning on x -> Missing at Random(MAR)
  x_OR = cbind(1, x_OR)
  x_PS = cbind(1, x_PS)
  
  
  res_foreach <- foreach(
    simnum = 1:SIMNUM,
    .packages = c("nleqslv", "CVXR"),
    .errorhandling="pass") %dopar% {
      delta = rbinom(N, 1, pi) # Sample indicator variable
      Index_S = (delta == 1)
      y_S = y[Index_S]
      x_OR_S = x_OR[Index_S,]; x_PS_S = x_PS[Index_S,]
      
      theta_res = NULL
      var_res = NULL
      CR_res = NULL
      
      PSmodel = glm(delta ~ x_PS, family = binomial)
      PSmodel$coefficients # PS model parameters
      pihat = predict.glm(PSmodel, type = "response") # Propensity score
      dhat = 1 / pihat; dhat_S = dhat[Index_S]; pihat_S = pihat[Index_S]
      
      y_IPW = sum(y_S * dhat_S) / sum(dhat_S)
      theta_res = c(theta_res, IPW = y_IPW * N)
      # var_res = c(var_res, IPW = sum(dhat_S * (dhat_S - 1) * (y_S - y_IPW)^2))
      
      hhat = pihat * x_PS # MLE
      kappa = solve(t(hhat[Index_S,]) %*% (hhat[Index_S,] *
                                             (1 / pihat[Index_S] - 1) / pihat[Index_S]),
                    t(hhat[Index_S,]) %*% ((y_S - y_IPW) *
                                             (1 / pihat[Index_S] - 1) / pihat[Index_S]))
      
      eta = drop(hhat %*% kappa)
      var_res = c(var_res, IPW = sum(dhat_S * (dhat_S - 1) *
                                       (y_S - y_IPW - eta[Index_S])^2))
      # eta[Index_S] = eta[Index_S] + (y_S - y_IPW - eta[Index_S]) / pihat[Index_S]
      # var_res = c(var_res, IPW = var(eta) * N)
      
      for(cnt2 in 1:3){
        w = CVXR::Variable(length(y_S))
        
        z = x_OR; z_S = x_OR_S
        if(cnt2 == 1){
          constraints <- list(t(z_S) %*% w == colSums(z))
          Phi_R <- CVXR::Minimize(sum((w - dhat_S)^2))
          
          # betahat = solve(t(z_S) %*% (z_S), t(z_S) %*% (y_S))
          # etahat = drop(z %*% betahat)
          # etahat[Index_S] = etahat[Index_S] + dhat_S * (y_S - drop(z_S %*% betahat))
          # sum(etahat)
        }else if(cnt2 == 2){
          constraints <- list(t(z_S) %*% w == colSums(z))
          # Phi_R <- CVXR::Minimize(sum(dhat_S * (w * pihat_S - 1)^2))
          Phi_R <- CVXR::Minimize(sum(((w - dhat_S)*dhat_S^{-1})^2))
          
          # betahat = solve(t(z_S) %*% (z_S * dhat_S), t(z_S) %*% (y_S * dhat_S))
          # etahat = drop(z %*% betahat)
          # etahat[Index_S] = etahat[Index_S] + dhat_S * (y_S - drop(z_S %*% betahat))
          # sum(etahat)
        }else if(cnt2 == 3){
          z = cbind(z, dhat); z_S = cbind(z_S, dhat_S)
          constraints <- list(t(z_S) %*% w == colSums(z))
          Phi_R <- CVXR::Minimize(sum(w^2))
          
          # betahat = solve(t(z_S) %*% (z_S), t(z_S) %*% (y_S))
          # etahat = drop(z %*% betahat)
          # sum(etahat)
        }
        
        # theta_res = c(theta_res, setNames(sum(etahat), paste("Cal", cnt2, sep = "")))
        # var_res = c(var_res, setNames(NA, paste("Cal", cnt2, sep = "")))
        
        prob <- CVXR::Problem(Phi_R, constraints)
        res <- CVXR::solve(prob)
        w_S = drop(res$getValue(w))
        
        betahat = solve(t(z_S) %*% (z_S * w_S), t(z_S) %*% (y_S * w_S))
        etahat = drop(z %*% betahat)
        etahat[Index_S] = etahat[Index_S] + w_S * (y_S - drop(z_S %*% betahat))
        
        theta_res = c(theta_res, setNames(sum(y_S * w_S), paste("Cal", cnt2, sep = "")))
        # var_res = c(var_res, Cal = var(etahat) * N)
        var_res = c(var_res, setNames(sum(w_S * (w_S - 1) *
                                            (y_S - drop(z_S %*% betahat))^2), paste("Cal", cnt2, sep = "")))
      }
      CR_res = ifelse(abs(theta_res - theta) < qnorm(0.975) * sqrt(var_res), 1, 0)
      
      list(theta_res - theta, var_res, CR_res)
    }
  final_res[[cnt]] <- res_foreach
}
stopCluster(cl)
timenow2 = Sys.time()
print(timenow2 - timenow1)
paste("# of failure:", sum(!sapply(final_res, function(x) is.numeric(unlist(x)))));

BIAS = sapply(final_res, function(x) rowMeans(do.call(cbind, lapply(x, `[[`, 1))))
colnames(BIAS) = rnames

SE = sapply(final_res, function(x) apply(do.call(cbind, lapply(x, `[[`, 1)), 1,
                                         function(y) sqrt(var(y, na.rm = TRUE) * (length(y)-1) / (length(y)))))
colnames(SE) = rnames

RMSE = sapply(final_res, function(x) sqrt(rowMeans(do.call(cbind, lapply(x, `[[`, 1))^2)))
colnames(RMSE) = rnames

# round(cbind(BIAS, RMSE), 1)

# gsub("\\\\addlinespace", "", kable(cbind(BIAS, RMSE), "latex", booktabs = TRUE, digits = 1) %>%
#        kable_styling())

res_var = sapply(final_res, function(x) rowMeans(do.call(cbind, lapply(x, `[[`, 2))))
colnames(res_var) = rnames

res_CR = sapply(final_res, function(x) rowMeans(do.call(cbind, lapply(x, `[[`, 3))))
colnames(res_CR) = rnames
# round(cbind((res_var - SE^2) / SE^2, res_CR), 3)

# gsub("\\\\addlinespace", "", kable(cbind((res_var - SE^2) / SE^2, res_CR), "latex", booktabs = TRUE, digits = 2) %>%
#        kable_styling())


## ----Point2, echo=FALSE, results="asis"-----------------------------------------
bias_tbl <- xtable(BIAS, digits = 2)
rmse_tbl <- xtable(RMSE, digits = 2)

# Output LaTeX directly
cat("
\\begin{table}[ht]
\\centering
\\begin{minipage}{0.48\\linewidth}
\\centering
")

print(bias_tbl, floating = FALSE, comment = FALSE)

cat("
\\caption{Bias of point estimators}
\\end{minipage}
\\hfill
\\begin{minipage}{0.48\\linewidth}
\\centering
")

print(rmse_tbl, floating = FALSE, comment = FALSE)

cat("
\\caption{RMSE of point estimators}
\\end{minipage}
\\end{table}
")


## ----Var2, echo=FALSE, results="asis"-------------------------------------------
rb_tbl <- xtable((res_var - SE^2) / SE^2, digits = 2)
cr_tbl <- xtable(res_CR, digits = 2)

# Output LaTeX directly
cat("
\\begin{table}[ht]
\\centering
\\begin{minipage}{0.48\\linewidth}
\\centering
")

print(rb_tbl, floating = FALSE, comment = FALSE)

cat("
\\caption{Relative bias of variance estimators}
\\end{minipage}
\\hfill
\\begin{minipage}{0.48\\linewidth}
\\centering
")

print(cr_tbl, floating = FALSE, comment = FALSE)

cat("
\\caption{Coverage rate of 95\\% CI}
\\end{minipage}
\\end{table}
")


## ----Point, eval=FALSE, include=FALSE, results="asis"---------------------------
## xtable::xtable(cbind(BIAS, RMSE), digits = 2,
##                caption = "Simulation results for point estimation(Bias and RMSE)", table.placement = "!htb")


## ----Var, echo=FALSE, include=FALSE, results="asis"-----------------------------
xtable::xtable(cbind((res_var - SE^2) / SE^2, res_CR), digits = 2,
               caption = "Simulation results for variance estimation(Bias and coverage rage)", table.placement = "!htb")


## ----boxplots, echo=FALSE, fig.align='center', out.height='35%', out.width='100%', fig.cap='Performance of the point estimators under four scenarios'----
res_all <- t(do.call(rbind, lapply(final_res, function(x) do.call(cbind, lapply(x, `[[`, 1)))))
colnames(res_all) <- as.vector(outer(rownames(RMSE), rnames, paste, sep = "_"))

# Convert to long format
df_long <- as.data.frame(res_all) %>%
  mutate(row = row_number()) %>%
  pivot_longer(-row, names_to = "Method_Scenario", values_to = "Value") %>%
  separate(Method_Scenario, into = c("Method", "Scenario"), sep = "_")

# For grouped x-axis display
df_long <- df_long %>%
  mutate(Scenario = factor(Scenario, levels = rnames),
         Method = factor(Method, levels = rownames(RMSE)))

# Generate distinct colors
colors <- rainbow(length(rownames(RMSE)))
names(colors) <- rownames(RMSE)

ggplot(df_long, aes(x = Scenario, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.5) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(
    x = "Scenario",
    y = NULL,  # leave y-axis label blank
    fill = "Estimator"
  ) +
  ggtitle(expression(hat(theta) - theta[N])) +  # move to main title
  # geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "solid", color = "red") +
  theme(
    axis.text.x = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 12)
  )

