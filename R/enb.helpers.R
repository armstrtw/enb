mldivide <- function(A,B) {
    solve(A) %*% B
}

mrdivide <- function(A, B) {
    stopifnot(is.numeric(A) || is.complex(A), is.numeric(B) || is.complex(B))
    if (is.vector(A))
        A <- t(A)
    if (is.vector(B))
        B <- t(B)
    if (ncol(A) != ncol(B))
        stop("Matrices 'A' and 'B' must have the same number of columns.")
    t(solve(t(B)) %*% t(A))
}


max.cor <- function(factor.returns) {
    x <- round(abs(cor(factor.returns)),2)
    max(x[row(x) < col(x)])
}

risk.remap <- function(bad.pair,wgts,factor.rets,verbose=TRUE) {
    stopifnot(length(bad.pair)==2L)
    stopifnot(length(wgts)==ncol(factor.rets))
    stopifnot(all(bad.pair %in% colnames(factor.rets)))
    stopifnot(all(bad.pair %in% names(wgts)))

    factor.sds <- apply(factor.rets[,bad.pair],2,sd)

    fsd1 <- factor.sds[bad.pair[1]]
    fsd2 <- factor.sds[bad.pair[2]]

    sd1 <- sd(wgts[bad.pair[1]] * factor.rets[,bad.pair[1]])
    sd2 <- sd(wgts[bad.pair[2]] * factor.rets[,bad.pair[2]])

    ## sd1 gets full weight
    if(sd1 >= sd2) {
        wgts[bad.pair[1]] <- wgts[bad.pair[1]] + wgts[bad.pair[2]] * fsd2/fsd1
        wgts <- wgts[-match(bad.pair[2],names(wgts))]
        if(verbose) cat("dropping: ",bad.pair[2],"\n")
    }
    ## sd2 gets full weight
    else {
        wgts[bad.pair[2]] <- wgts[bad.pair[2]] + wgts[bad.pair[1]] * fsd1/fsd2
        wgts <- wgts[-match(bad.pair[1],names(wgts))]
        if(verbose) cat("dropping: ",bad.pair[1],"\n")
    }
    wgts
}

cull.dimension <- function(wgts,factor.rets) {
    x <- abs(cor(factor.rets))
    ## poor man's way to find max position
    ## set lower to 0
    x[row(x) >= col(x)] <- 0
    max.ele <- which.max(x)
    ele.row <- rep(rownames(x),each=ncol(x))[max.ele]
    ele.col <- rep(colnames(x),nrow(x))[max.ele]
    risk.remap(c(ele.row,ele.col),wgts,factor.rets)
}

fix.defects <- function(wgts,factor.rets,thresh=.99) {
    stopifnot(all(names(wgts) == colnames(factor.rets)))
    while(max.cor(factor.rets) > thresh) {
        wgts <- cull.dimension(wgts,factor.rets)
        ## redimension to new wgts
        factor.rets <- factor.rets[,names(wgts)]
    }
    list(wgts=wgts,factor.rets=factor.rets)
}

marginal.risk.contribution <- function(b, Sigma) {
    mrdivide(b * (Sigma %*% b),(t(b) %*% Sigma %*% b) );
}
