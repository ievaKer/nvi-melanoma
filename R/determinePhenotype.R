determinePhenotype <- function(x, classifyGenes) {
  inflamed <- function(angiogenesis, immune, stroma) {
    scores <- c(angiogenesis, immune, stroma)
    (which.max(scores) == 2) * 1 + (immune > 0) * 1 + (angiogenesis < 0) * 0.3
  }

  excluded <- function(angiogenesis, immune, stroma) {
    scores <- c(angiogenesis, immune, stroma)
    (max(scores) == immune) * 0.6 + (max(scores) == stroma) * 0.7 + (stroma > 0) * 0.7 + (angiogenesis < 0) * 0.3
  }

  deserted <- function(angiogenesis, immune, stroma) {
    scores <- c(angiogenesis, immune, stroma)
    (immune < 0 & stroma < 0) * 1.7 + (angiogenesis > 0) * 0.3 + (max(scores) == angiogenesis) * 0.3
  }

  whichPhenotype <- function(angiogenesis, immune, stroma) {
    vals <- c('Inf' = inflamed(angiogenesis, immune, stroma),
              'Exc' = excluded(angiogenesis, immune, stroma),
              'Des' = deserted(angiogenesis, immune, stroma))

    names(vals)[ which.max(vals) ]
  }

  x[ , classifyGenes[ , Gene ] ] %>%
    t %>%
    as.data.table(keep.rownames = 'Gene') %>%
    reshape2::melt(id.vars = 'Gene', value.name = 'Value', variable.name = 'SID') %>%
    setDT %>%
    merge(classifyGenes, by = 'Gene') %>%
    .[ , list(Mean = mean(Value)), by = list(SID, Type) ] %>%
    reshape2::dcast(SID ~ Type, value.var = 'Mean') %>%
    setDT %>%
    .[ , Phenotype := whichPhenotype(ANGIOGENESIS, IMMUNE, STROMA), by = SID ] %>%
    .[ , list(SID, Phenotype) ] %>%
    .[ , Phenotype := factor(Phenotype, levels = c('Des', 'Exc', 'Inf')) ]
}
