# Installing our package

- 1. install the 'remotes' package 

if ( !("remotes" %in% installed.packages()) ) {
    install.packages("remotes")
}

- 2. install our package directly from Github:

remotes::install_github("Laefritzia/uniLM")


# uniLM

This is a package providing material and a first approach to implement
the described minmax MOM estimator in 
https://projecteuclid.org/journalArticle/Download?urlId=10.1214%2F19-AOS1828
, developed for a bachelor's course at the University of Vienna, 2026SS.
It is currently publicly available for our university-colleagues. The package is functional 
in its first draft, but extensive Software testing has yet to be performed.
Please feel free to message us at leablahout@gmail.com if you have any concerns or questions.
