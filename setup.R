install.packages(c('BiocManager','Hmisc', 'RcppEigen','RcppNumerical'), 
                 dependencies='Depends',
                repo = "http://cran.us.r-project.org")
BiocManager::install(c('DESeq2','apeglm', 'rhdf5', 'tximport'))


library(DESeq2)
library(apeglm)
