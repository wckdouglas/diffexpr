install.packages(c('BiocManager','Hmisc'), 
                 dependencies='Depends',
                repo = "http://cran.us.r-project.org")
BiocManager::install(c('DESeq2','apeglm'))
