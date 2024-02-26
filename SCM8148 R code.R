
## APPENDIX I : CODE


setwd("C:/Users/swmit/OneDrive - Queen's University Belfast/SCM8148 Health and Biomedical Informatics and the Exposome/Written Report")

# QUESTION 6

# install required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")

BiocManager::install(c("meshr", "AnnotationHub", "MeSHDbi"))

      # Load required libraries
      library(meshr)
      library(AnnotationHub)
      library(MeSHDbi)
      library(BiocManager)

  # Read the gene set file
  genes <- read.delim("assignment_meshr.csv", header = TRUE)
  
  # Database preparation
    ah <- AnnotationHub()
    d <- display(ah)
  
    dbfile1 <- query(ah, c("MeSHDb", "MeSH.db", "v002"))[[1]]
    dbfileQ <- query(ah, c("MeSHDb", "Homo sapiens", "v002"))[[1]]
    dbfile2 <- query(ah, c("MeSHDb", "MeSH.AOR.db", "v002"))[[1]]
    dbfile3 <- query(ah, c("MeSHDb", "MeSH.PCR.db", "v002"))[[1]]
    MeSH.Hs.db <- MeSHDbi::MeSHDb(dbfileQ)
    MeSH.db <- MeSHDbi::MeSHDb(dbfile1)
    MeSH.AOR.db <- MeSHDbi::MeSHDb(dbfile2)
    MeSH.PCR.db <- MeSHDbi::MeSHDb(dbfile3)
  
    # ORA analysis - defining the relevant parameters
    datameshR <- read.csv("assignment_meshr.csv")
    meshParams <- new("MeSHHyperGParams",
                      geneIds = datameshR$Query,
                      universeGeneIds = datameshR$Background,
                      annotation = "MeSH.Hs.db",
                      meshdb = "MeSH.db", category = "C",
                      database = "gene2pubmed",
                      pvalueCutoff = 1e-6, pAdjust = "BH")
    
    # ORA analysis - Running the test
    meshR <- meshHyperGTest(meshParams)
    summary(meshR)
    
    
# QUESTION 10 
  
  # Install “devtools” package.
    install.packages("devtools")
  # Install phexpo package from github
    devtools::install_github("GHLCLab/phexpo")
  # Load phexpo package
    library(phexpo)
    # run test    
      pent <- perfFishTestChemSingle("3,4,5,3',4'-pentachlorobiphenyl",
                                     enrich_1S= TRUE)
    # Filter results to only those with bonf corrected p-value less than or equal to 0.05
      pent_filtered <- pent[pent$bonf <= 0.05,]
    # print
      print(pent_filtered)
    # write results to CSV file
      write.csv(pent, file = "pentachlorobiphenyl_results.csv",
                row.names = FALSE)
      
