ccor = function(x,y,M=3){
    
    if(! "XIM" %in% installed.packages()){
        stop("Must first install the XIM package:
             devtools::install_github('jeskowagner/XIM')")
    }
    if(M > 1000){
        warning("M is very large and computation might take longer.")
    }
    res = list()
    res$estimate = XIM::XIMcalculate(x, y, M=M)
    res$p.value = "NA"
    res$conf.int = "NA"
    res
}

.cor_return = function(cor, what=c("all","cor","conf","pval")){
    # helper to decide what to return from cor_wrap
    if(what=="cor"){
        return(cor$estimate)
    } else if (what=="conf"){
        return(as.numeric(cor$conf.int))
    } else if (what=="pval") {
        return(cor$p.value)
    } else if (what=="all") {
        return(list(cor=cor$estimate,conf=cor$conf.int,pval=cor$p.value))
    } else {
        stop("what must be one of 'all', 'cor', 'conf', or 'pval'")
    }
}
#' @export
#' @title cor_wrap
#'
#' @description Wrapper for cor.test with controlled returns
#'
#' @param x Numeric vector
#' @param y Numeric vector
#' @param method Which correlation to compute.
#' Takes as special value 'chatterjee' to compute revised
#' Chatterjee correlation coefficient (no p values returned in that case)
#' @param return Whether to return only the correlation coefficient, the confidence interval, the p value, or all
#' @param ... Other parameters passed to \code{cor.test}
#' @return what If all, a named list, otherwise a numeric vector
#'
#' @author Jesko Wagner
cor_wrap = function(x, y, method="pearson",
                what=c("all","cor","conf","pval"),
                ...){
    if(method == "chatterjee"){
        res = ccor(x, y, ...)      
    } else {
        res = cor.test(x, y, method=method)
    }
    .cor_return(res, what)
}
