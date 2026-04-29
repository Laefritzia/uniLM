# Installing our package

Use the following R-commands, to install the package directly from Github (with the 'remotes' package):

if ( !("remotes" %in% installed.packages()) ) {
    install.packages("remotes")
}

remotes::install_github("Laefritzia/uniLM")

# uniLM

This package provides a unified framework for robust linear regression centered around the minmax 
Median-of-Means method proposed by Lecué and Lerasle (2020) [1] , to provide both a robust point estimate with bootstrapped variance, as well as an 
outlier-detection diagnostic tool. Further diagnostic tools based on the Ordinary-Least-Squares (OLS)-fit, as well as the Huber-loss M-estimator
are provided.

This is a project developed for a bachelor's course at the University of Vienna, 2026SS. It is currently publicly available for our university-colleagues. The package is functional 
in its first draft, but extensive Software Testing has yet to be performed. Please feel free to message us at leablahout@gmail.com if you have any concerns or questions.

The folder Seminararbeit2026/ contains the material used in the Seminarpaper, which we wrote in named course.

[1] https://projecteuclid.org/journalArticle/Download?urlId=10.1214%2F19-AOS1828
, 
