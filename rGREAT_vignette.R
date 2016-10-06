## ---- echo = FALSE, message = FALSE--------------------------------------
library(knitr)
knitr::opts_chunk$set(
  error = FALSE,
  tidy  = FALSE,
  message = FALSE,
  fig.align = "center")
options(width = 100)
options(markdown.HTML.stylesheet = "custom.css")

## ---- echo = 2-----------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(rGREAT)))
library(rGREAT)

## ------------------------------------------------------------------------
set.seed(123)
bed = circlize::generateRandomBed(nr = 1000, nc = 0)
bed[1:2, ]

## ------------------------------------------------------------------------
job = submitGreatJob(bed)

## ------------------------------------------------------------------------
job

## ---- eval = FALSE-------------------------------------------------------
## job = submitGreatJob(bed, species = "mm9")
## job = submitGreatJob(bed, bg, species = "mm9")
## job = submitGreatJob(bed, adv_upstream = 10, adv_downstream = 2, adv_span = 2000)
## job = submitGreatJob(bed, rule = "twoClosest", adv_twoDistance = 2000)
## job = submitGreatJob(bed, rule = "oneClosest", adv_oneDistance = 2000)

## ---- eval = FALSE-------------------------------------------------------
## job = submitGreatJob(bed, version = "3.0")
## job = submitGreatJob(bed, version = "2.0")

## ------------------------------------------------------------------------
tb = getEnrichmentTables(job)
names(tb)
tb[[1]][1:2, ]

## ------------------------------------------------------------------------
job

## ---- eval = FALSE-------------------------------------------------------
## tb = getEnrichmentTables(job, ontology = c("GO Molecular Function", "BioCyc Pathway"))
## tb = getEnrichmentTables(job, category = c("GO"))

## ------------------------------------------------------------------------
availableOntologies(job)
availableCategories(job)
availableOntologies(job, category = "Pathway Data")

## ---- fig.width = 12, fig.height = 4, fig.align = 'center'---------------
par(mfrow = c(1, 3))
res = plotRegionGeneAssociationGraphs(job)
res[1:2, ]

## ---- eval = FALSE-------------------------------------------------------
## plotRegionGeneAssociationGraphs(job, type = 1)

## ---- fig.width = 12, fig.height = 4-------------------------------------
par(mfrow = c(1, 3))
res = plotRegionGeneAssociationGraphs(job, ontology = "GO Molecular Function",
                                      termID = "GO:0004984")
res[1:2, ]

## ------------------------------------------------------------------------
sessionInfo()

