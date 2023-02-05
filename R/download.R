downloadTCGA <- function() {
  library(data.table)
  library(tidyverse)
  library(foreach)
  library(openxlsx)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(HGNChelper)
  library(rapiclient)

  #
  # Read all gene lists and map to correct gene names
  #

  classificationGenes <- fread('data/classificationGenes.csv')

  nanoString <- c('data/LBL-10025-02_PanCancer-Pathways-Gene-List.xlsx',
                  'data/LBL-10043-08_nCounter_PanCancer_Immune_Profiling_Panel_Gene_List.xlsx',
                  'data/LBL-10498-02_IO_360_Gene_List.xlsx')

  genes <- foreach(file = nanoString, .combine = append) %do% {
    excelData <- read.xlsx(file, sheet = 2, startRow = 2) %>%
      setDT
    colnames(excelData)[ 1 ] <- 'Official.Symbol'
    dropFrom <- excelData[ , which(grepl('internal reference', tolower(Official.Symbol))) ]

    excelData[ 1:(dropFrom - 1) ] %>%
      .[ , 1 ] %>%
      deframe
  } %>%
    append(classificationGenes[ , Gene ]) %>%
    HGNChelper::checkGeneSymbols() %>%
    .[ , 'Suggested.Symbol' ] %>%
    unique %>%
    sort %>%
    na.omit %>%
    AnnotationDbi::mapIds(org.Hs.eg.db, ., 'ENTREZID', 'SYMBOL') %>%
    as.list %>%
    .[ !(duplicated(.) | is.na(names(.)) | is.na(.)) ] %>%
    data.table(Entrez = ., Symbol = names(.))


  #
  # Download data from cBioPortal
  #

  # Setup cBioPortal client
  client <- rapiclient::get_api(url = 'https://www.cbioportal.org/api/api-docs')
  operations <- rapiclient::get_operations(client)
  schemas <- rapiclient::get_schemas(client)

  message('Downloading data...')

  entrezIds <- genes[ , Entrez ]

  resp <- operations$fetchAllMolecularDataInMolecularProfileUsingPOST(
    projection = 'SUMMARY',
    molecularProfileId = 'skcm_tcga_rna_seq_v2_mrna_median_all_sample_Zscores',
    entrezGeneIds = entrezIds,
    sampleListId = 'skcm_tcga_all'
  )

  content <- rapiclient::content_or_stop(resp)

  message('Combining data...')
  counts <- as.data.table(content) %>%
    t %>%
    as.data.table %>%
    setnames(names(content[[ 1 ]])) %>%
    .[ , list(Patient = unlist(patientId), Entrez = unlist(entrezGeneId), Score = unlist(value)) ] %>%
    reshape2::dcast(Patient ~ Entrez, fun.aggregate = mean, value.var = 'Score')

  # Rename count columns and drop ids from matrix itself
  rownames(counts) <- counts[ , 1 ]
  counts <- counts[ , -1 ]

  # Set symbols as colnames
  colnames(counts) <- genes[ match(colnames(counts), Entrez), Symbol ]

  dir.create('output', showWarnings = FALSE)
  list(x = counts, map = genes) %>%
    saveRDS('output/raw.RDS')
}
