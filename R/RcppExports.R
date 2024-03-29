# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

pdf <- function(x, para) {
    .Call(`_EstimWHD_pdf`, x, para)
}

cdf <- function(x, para) {
    .Call(`_EstimWHD_cdf`, x, para)
}

sur <- function(x, para) {
    .Call(`_EstimWHD_sur`, x, para)
}

quan <- function(init, p, para) {
    .Call(`_EstimWHD_quan`, init, p, para)
}

cpc_fun <- function(para, L, U, P0) {
    .Call(`_EstimWHD_cpc_fun`, para, L, U, P0)
}

cpc_grad <- function(para, L, U, P0) {
    .Call(`_EstimWHD_cpc_grad`, para, L, U, P0)
}

GenDataRcpp <- function(init, para, R, k) {
    .Call(`_EstimWHD_GenDataRcpp`, init, para, R, k)
}

loglike <- function(para, X, R, k) {
    .Call(`_EstimWHD_loglike`, para, X, R, k)
}

logMPS <- function(para, X, R, k) {
    .Call(`_EstimWHD_logMPS`, para, X, R, k)
}

prior <- function(para) {
    .Call(`_EstimWHD_prior`, para)
}

Estim <- function(para, X, R, L, U, P0, k, type) {
    .Call(`_EstimWHD_Estim`, para, X, R, L, U, P0, k, type)
}

TK <- function(para, X, R, q, c, L, U, P0, k, type) {
    .Call(`_EstimWHD_TK`, para, X, R, q, c, L, U, P0, k, type)
}

MH_sample <- function(type, para, se, R, X, k, L, U, P0, MC_size, MC_burn, q, c, verbose = 0L, display_progress = TRUE) {
    .Call(`_EstimWHD_MH_sample`, type, para, se, R, X, k, L, U, P0, MC_size, MC_burn, q, c, verbose, display_progress)
}

