## Copyright (c) 2014, Attilio Meucci
## All rights reserved.

## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:

## * Redistributions of source code must retain the above copyright notice, this
##   list of conditions and the following disclaimer.

## * Redistributions in binary form must reproduce the above copyright notice,
##   this list of conditions and the following disclaimer in the documentation
##   and/or other materials provided with the distribution.

## * Neither the name of [project] nor the names of its
##   contributors may be used to endorse or promote products derived from
##   this software without specific prior written permission.

## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
## AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
## FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
## DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
## SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
## CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
## OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

require(expm) ## for sqrtm

## this file is governed by the BSD license, please see:
## http://www.mathworks.com/matlabcentral/fileexchange/43245-portfolio-diversi%C2%85cation-based-on-optimized-uncorrelated-factors

effective.bets <- function(b, Sigma, Tm) {

    ## this funciton computes the Effective Number of Bets and the Diversification distribution
    ## see A. Meucci, A. Santangelo, R. Deguest - "Measuring Portfolio Diversification Based on Optimized Uncorrelated Factors" to appear (2013)
    ##
    ## Last version of code and article available at http://symmys.com/node/599

    ##
    ## INPUTS
    ##  b     : [vector] (n_ x 1) exposures
    ##  Sigma : [matrix] (n_ x n_) covariance matrix
    ##  t     : [matrix] (n_ x n_) torsion matrix
    ##
    ## OUTPUTS
    ##  enb : [scalar] Effetive Number of Bets
    ##  p   : [vector] (n_ x 1) diversification distribution

    p = mrdivide(mldivide(t(Tm),b) * (Tm %*% Sigma %*% b) , ( t(b) %*% Sigma %*% b ) );
    enb = exp(-sum(p*log(1+(p-1)*(p>1e-5))));

    list(enb=enb,p=p)
}


torsion <- function(Sigma, model="pca", method="approximate", max_niter=10000L) {

    ## this funciton computes the Principal Components torsion and the Minimum Torsion for diversification analysis
    ## see A. Meucci, A. Santangelo, R. Deguest - "Measuring Portfolio Diversification Based on Optimized Uncorrelated Factors" to appear (2013)

    ## Last version of code and article available at http://symmys.com/node/599

    ## INPUTS
    ##  Sigma     : [matrix] (n_ x n_) covariance matrix
    ##  model     : [string] choose between 'pca' and 'minimum-torsion' model
    ##  method    : [string] choose between 'approximate' and 'exact' method for 'minimum-torsion' model
    ##  max_niter : [scalar] choose number of iterations of numerical algorithm 
    ## OUTPUTS
    ## t : [matrix] (n_ x n_) torsion matrix

    if(model=="pca") {
        ## R already sorts these
        Tm <- prcomp(Sigma,center=FALSE,scale=FALSE)$rotation
        flip <- Tm[1,] < 0
        Tm[,flip] <- -Tm[,flip]

        ## PCA torsion, transpose
        Tm = t(Tm)
    } else if(model=='minimum-torsion') {
        ## alt
        rho <- cov2cor(Sigma)

        ## Correlation matrix
        sigma = diag(Sigma)^(1/2)
        rho = diag(1/sigma) %*% Sigma %*% diag(1/sigma)
        rho.norm = sqrtm(rho) ## Riccati root of C
            
        if(method=='approximate') {
            Tm = mrdivide(diag(sigma),rho.norm) %*% diag(1/sigma)
        } else if(method=='exact') {
            n_ = ncol(Sigma)
            
            ## initialize
            d = rep(1, n_)
            f = rep(0, max_niter)
            for(i in 1:max_niter) {
                U = diag(d) %*% rho.norm %*% rho.norm %*% diag(d)
                u = sqrtm(U)
                q = solve(u) %*% (diag(d) %*% rho.norm)
                d = diag(q %*% rho.norm)
                pi_ = diag(d) %*% q ## perturbation
                f[i] = norm(rho.norm - pi_, 'F')
        
                if (i > 1 && abs(f[i] - f[i-1])/f[i]/n_ <= 10^-8) {
                    f = f[1:i]
                    break
                } else if (i == max_niter && abs(f[max_niter] - f[max_niter-1])/f[max_niter]/n_ > 10^-8) {
                    warning("number of max iterations reached: n_iter = ", max_niter)
                }
            }
            ##cat(i)
            x = mrdivide(pi_,rho.norm)
            Tm = diag(sigma) %*% x %*% diag(1/sigma)
        }
    }
    Tm
}
