r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r) 

install.packages('BiocManager')
BiocManager::install('DESeq2')
