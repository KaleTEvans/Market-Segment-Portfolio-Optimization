#############################
# Project 1
#############################

####################################################################
# Using code below from provided solution to ensure accuracy of data
####################################################################

# clear environment
rm(list = ls())

# libraries and packages
library(quantmod)
library(tidyverse)
library(moments)
library(CVXR)
library(PerformanceAnalytics)
#library(glpkAPI)

# Cache data once computed correctly
cached_file <- "cache/cached_project1_data.rda"
if ( file.exists(cached_file) ) {
  load(cached_file)
} else {

# Get data
from = as.Date("01/01/2005", format = "%m/%d/%Y")
to = as.Date("04/01/2024", format = "%m/%d/%Y")

yahooSymbols_df = data.frame(symbols = c("IWB", "IWM",
                                         "EFA", "EEM",
                                         "VNQ", "LQD","SHY"),
                             label = c("Large Cap US","Small Cap US",
                                       "Dev. Mkts", "Emer Mkts",
                                       "Global REIT", "Corp. Bonds","Short Tsy"))


invisible(getSymbols(yahooSymbols_df$symbols, src = "yahoo", from = from , to = to))

yahooReturnData_xts <- do.call("cbind",lapply(yahooSymbols_df$symbols, FUN = function(symb){
  temp <- eval(parse(text = symb))
  #  colnames(temp) <-  paste0(symb,c(".Open",".High",".Low",".Close",".Volume",".Adjusted"))
  temp <- temp[,6]
  #  colnames(temp) <- c(paste0(symb,""))
}))
colnames(yahooReturnData_xts) <- gsub(".Adjusted","",colnames(yahooReturnData_xts))

temp <- log(yahooReturnData_xts/stats::lag(yahooReturnData_xts))

temp <- temp[-1,]

financialDataReturns_xts <- temp

dim(financialDataReturns_xts)

# Aggregate returns to 5, 21, 62 and 252 day intervals

return_list = list()
return_list_df = list()
for(numDays in c(1, 5, 21, 62, 252)){
  return_xts = rollapply(financialDataReturns_xts,
                         width = numDays, 
                         sum,
                         by = numDays,
                         align = "right")
  item_name =  paste0("days_",numDays)
  # add non-NA rows of returns_xts to list
  valid_returns_xts = return_xts[rowSums(is.na(return_xts)) == 0,]
  return_list[[item_name]] = valid_returns_xts
  
  # Convert to data frame and set row names as Date
  temp_df <- as.data.frame(valid_returns_xts)
  rownames(temp_df) <- index(valid_returns_xts)  # Setting index as row names
  return_list_df[[item_name]] = temp_df
}

# check number of observations for each time period length
lapply(return_list, dim)

######################################################################
# End copied code
######################################################################

if (!require("ROI")) install.packages("ROI"); library("ROI")
if (!require("ROI.plugin.glpk")) install.packages("ROI.plugin.glpk"); library("ROI.plugin.glpk")
if (!require("ROI.plugin.quadprog")) install.packages("ROI.plugin.quadprog"); library("ROI.plugin.quadprog")
if (!require("ROI.plugin.alabama")) install.packages("ROI.plugin.alabama"); library("ROI.plugin.alabama")

###########################
# Min Variance Objective
###########################
markowitz_objective <- function(r_mat) {
  objective <- Q_objective(Q = 2 * cov(r_mat), 
                           L = rep(0, NCOL(r_mat)))
  list(objective = objective)
}

###########################
# Min MAD Objective
###########################
mad_objective <- function(r_mat){
  x.names <- colnames(r_mat)
  N <- NCOL(r_mat)
  S <- nrow(r_mat)
  mu <- colMeans(r_mat)
  Amat <-  cbind(sweep(as.matrix(r_mat), 2,  mu), - diag(S), diag(S))
  var.names <- c(x.names, 
                 paste0("y_mad_aux", seq_len(S)),
                 paste0("z_mad_aux", seq_len(S)))
  
  constraint <-  L_constraint(L = Amat, dir = rep("==", S), 
                              rhs = rep(0, S), 
                              names = var.names)
  
  objective <- L_objective(L = c(rep(0, N), rep(1/S, 2 * S)))
  
  list(objective = objective, constraint = constraint)
}

#############################
# Min CVAR Objective
#############################

cvar_objective <- function(r_mat, alpha, probs = NULL) {
  x.names <- colnames(r_mat)
  N <- NCOL(r_mat)
  S <- NROW(r_mat)
  mu <- colMeans(r_mat)
  if (is.null(probs)) probs <- rep(1/S, S)
  if (alpha < 0.5) alpha <- 1 - alpha
  
  Amat <- cbind(as.matrix(r_mat),  diag(S), 1)
  var.names <- c(x.names, paste0("z_cvar_aux", seq_len(S)), "gamma")
  
  ## set bounds for gamma (-Inf, Inf) 
  bnds <- ROI::V_bound(li = c(N + S + 1), lb = c( -Inf),
                       ui = c(N + S + 1), ub = c(  Inf))
  
  constraint <- L_constraint(L = Amat, dir = rep(">=", S), 
                             rhs = rep(0, S), 
                             names = var.names)
  
  objective <- L_objective(c(rep(0, N), probs/(1 - alpha), 1))
  
  list(objective = objective, constraint = constraint, bounds = bnds)
}

############################
# Full Investment Constraint
############################
budget_constraint <- function(r_mat, dir = "==", rhs = 1) {
  x.names <- colnames(r_mat)
  L_constraint(L = rep(1, NCOL(r_mat)), 
               dir = dir,  rhs = rhs, names = x.names)
}

#############################
# Target Return Constraint
#############################
reward_constraint <- function(r_mat, dir = ">=", rhs = target) {
  x.names <- colnames(r_mat)
  L_constraint(L = colMeans(r_mat), dir = dir,  
               rhs = rhs, names = x.names)
}

#################################################
# Function for Computing MAD, Min Var, and CVAR
#################################################

compute_risk_measures <- function(returns, target) {
  
  ## Min Variance
  lp  <- OP(objective  =  markowitz_objective(returns)$objective,
            constraints = rbind(budget_constraint(returns),
                                reward_constraint(returns, rhs = target)))
  (sol <- ROI_solve(lp, solver = "quadprog"))
  min_var_solution <- round(solution(sol), 4)
  
  ## MAD
  lp  <- OP(objective  =  mad_objective(returns)$objective,
            constraints = rbind(mad_objective(returns)$constraint,
                                budget_constraint(returns),
                                reward_constraint(returns, rhs = target),
                                use.names = TRUE))
  (sol <- ROI_solve(lp, solver = "glpk"))
  mad_solution <- round(solution(sol)[seq_len(NCOL(returns))], 4)
  
  ## CVAR
  tmp <- cvar_objective(returns, alpha = 0.90)
  lp  <- OP(objective  =  tmp$objective,
            constraints = rbind(tmp$constraint,
                                budget_constraint(returns),
                                reward_constraint(returns, rhs = target),
                                use.names = TRUE),
            bounds = tmp$bounds)
  (sol <- ROI_solve(lp, solver = "glpk"))
  cvar_solution <- round(solution(sol)[seq_len(NCOL(returns))], 4)
  
  return(cbind(min_var_solution, mad_solution, cvar_solution))
}

#############################
# Daily Interval Portfolios

daily_multiplier <- 1/252
weekly_multiplier <- 5/252
monthly_multiplier <- 21/252
quarterly_multiplier <- 62/252
annual_multiplier <- 1

daily_2pct_target <- compute_risk_measures(return_list_df$days_1, (0.02 * daily_multiplier))
daily_4pct_target <- compute_risk_measures(return_list_df$days_1, (0.04 * daily_multiplier))
daily_6pct_target <- compute_risk_measures(return_list_df$days_1, (0.06 * daily_multiplier))

weekly_2pct_target <- compute_risk_measures(return_list_df$days_5, (0.02 * weekly_multiplier))
weekly_4pct_target <- compute_risk_measures(return_list_df$days_5, (0.04 * weekly_multiplier))
weekly_6pct_target <- compute_risk_measures(return_list_df$days_5, (0.06 * weekly_multiplier))

monthly_2pct_target <- compute_risk_measures(return_list_df$days_21, (0.02 * monthly_multiplier))
monthly_4pct_target <- compute_risk_measures(return_list_df$days_21, (0.04 * monthly_multiplier))
monthly_6pct_target <- compute_risk_measures(return_list_df$days_21, (0.06 * monthly_multiplier))

quarterly_2pct_target <- compute_risk_measures(return_list_df$days_62, (0.02 * quarterly_multiplier))
quarterly_4pct_target <- compute_risk_measures(return_list_df$days_62, (0.04 * quarterly_multiplier))
quarterly_6pct_target <- compute_risk_measures(return_list_df$days_62, (0.06 * quarterly_multiplier))

annual_2pct_target <- compute_risk_measures(return_list_df$days_252, (0.02 * annual_multiplier))
annual_4pct_target <- compute_risk_measures(return_list_df$days_252, (0.04 * annual_multiplier))
annual_6pct_target <- compute_risk_measures(return_list_df$days_252, (0.06 * annual_multiplier))

# Save cached variables
save(return_list, return_list_df, 
     daily_2pct_target, daily_4pct_target, daily_6pct_target,
     weekly_2pct_target, weekly_4pct_target, weekly_6pct_target,
     monthly_2pct_target, monthly_4pct_target, monthly_6pct_target,
     quarterly_2pct_target, quarterly_4pct_target, quarterly_6pct_target,
     annual_2pct_target, annual_4pct_target, annual_6pct_target, file = cached_file)

}

convert_to_df <- function(target_matrix) {
  df <- as.data.frame(target_matrix)
  colnames(df) <- c("Min Variance", "Min MAD", "Min CVaR")
  return(df)
}

target_list <- list(
  daily_2pct = daily_2pct_target,
  daily_4pct = daily_4pct_target,
  daily_6pct = daily_6pct_target,
  weekly_2pct = weekly_2pct_target,
  weekly_4pct = weekly_4pct_target,
  weekly_6pct = weekly_6pct_target,
  monthly_2pct = monthly_2pct_target,
  monthly_4pct = monthly_4pct_target,
  monthly_6pct = monthly_6pct_target,
  quarterly_2pct = quarterly_2pct_target,
  quarterly_4pct = quarterly_4pct_target,
  quarterly_6pct = quarterly_6pct_target,
  annual_2pct = annual_2pct_target,
  annual_4pct = annual_4pct_target,
  annual_6pct = annual_6pct_target
)

# Convert each portfolio to a data frame with uniform column names
portfolio_weight_list <- lapply(target_list, convert_to_df)

options(scipen = 999)

# For some reason I was getting errors when attempting to calculate CVaR in get_portfolio_statistics
# The safe_cvar_calc() function will ensure the calculation is performed correctly
safe_cvar_calc <- function(returns, var_value) {
  # Ensure var_value is numeric and not an array/matrix
  var_value <- as.numeric(var_value)
  
  # Subset returns below the negative VaR
  sub_returns <- returns[returns <= -var_value]
  
  # Calculate mean if subset is not empty
  if (length(sub_returns) > 0) {
    return(-mean(sub_returns))
  } else {
    return(NA)
  }
}

get_portfolio_statistics <- function(weights, returns) {
  min_var_weights <- matrix(weights[,1], nrow = 1, ncol = length(daily_2pct_target[,1]), byrow = TRUE)
  mad_weights <- matrix(weights[,2], nrow = 1, ncol = length(daily_2pct_target[,1]), byrow = TRUE)
  cvar_weights <- matrix(weights[,3], nrow = 1, ncol = length(daily_2pct_target[,1]), byrow = TRUE)
  
  min_var_returns <- min_var_weights %*% t(returns)
  mad_returns <- mad_weights %*% t(returns)
  cvar_returns <- cvar_weights %*% t(returns)
  
  min_var_returns <- as.vector(min_var_returns)
  mad_returns <- as.vector(mad_returns)
  cvar_returns <- as.vector(cvar_returns)
  
  min_var_mean <- mean(min_var_returns)
  mad_mean <- mean(mad_returns)
  cvar_mean <- mean(cvar_returns)
  
  min_var_sd <- sd(min_var_returns)
  mad_sd <- sd(mad_returns)
  cvar_sd <- sd(cvar_returns)
  
  min_var_mad <- MeanAbsoluteDeviation(min_var_returns)
  mad_mad <- MeanAbsoluteDeviation(mad_returns)
  cvar_mad <- MeanAbsoluteDeviation(cvar_returns)
  
  min_var_var <- quantile(min_var_returns, 0.10) * -1
  mad_var <- quantile(mad_returns, 0.10) * -1
  cvar_var <- quantile(cvar_returns, 0.10) * -1

  min_var_cvar <- safe_cvar_calc(min_var_returns, min_var_var)
  mad_cvar <- safe_cvar_calc(mad_returns, mad_var)
  cvar_cvar <- safe_cvar_calc(cvar_returns, cvar_var)
  
  stats_df = data.frame(
    "Min Variance" = c(min_var_mean, min_var_sd, min_var_mad, min_var_var, min_var_cvar) * 100,
    "Min MAD" = c(mad_mean, mad_sd, mad_mad, mad_var, mad_cvar) * 100,
    "Min CVaR" = c(cvar_mean, cvar_sd, cvar_mad, cvar_var, cvar_cvar) * 100
  )
  
  row.names(stats_df) <- c("Expected Return", "Standard Deviation", "Mean Absolute Deviation", "Value-at-Risk", "Conditional Value-at-Risk")
  return(stats_df)
}

portfolio_stats_list <- list(
  stats_2pct_daily = get_portfolio_statistics(daily_2pct_target, return_list$days_1),
  stats_4pct_daily <- get_portfolio_statistics(daily_4pct_target, return_list$days_1),
  stats_6pct_daily <- get_portfolio_statistics(daily_6pct_target, return_list$days_1),
  
  stats_2pct_weekly <- get_portfolio_statistics(weekly_2pct_target, return_list$days_5),
  stats_4pct_weekly <- get_portfolio_statistics(weekly_4pct_target, return_list$days_5),
  stats_6pct_weekly <- get_portfolio_statistics(weekly_6pct_target, return_list$days_5),
  
  stats_2pct_monthly <- get_portfolio_statistics(monthly_2pct_target, return_list$days_21),
  stats_4pct_monthly <- get_portfolio_statistics(monthly_4pct_target, return_list$days_21),
  stats_6pct_monthly <- get_portfolio_statistics(monthly_6pct_target, return_list$days_21),
  
  stats_2pct_quarterly <- get_portfolio_statistics(quarterly_2pct_target, return_list$days_62),
  stats_4pct_quarterly <- get_portfolio_statistics(quarterly_4pct_target, return_list$days_62),
  stats_6pct_quarterly <- get_portfolio_statistics(quarterly_6pct_target, return_list$days_62),
  
  stats_2pct_annual <- get_portfolio_statistics(annual_2pct_target, return_list$days_252),
  stats_4pct_annual <- get_portfolio_statistics(annual_4pct_target, return_list$days_252),
  stats_6pct_annual <- get_portfolio_statistics(annual_6pct_target, return_list$days_252)
)

print(portfolio_stats_list$stats_4pct_daily)