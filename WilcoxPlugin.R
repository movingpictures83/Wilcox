library(FELLA)
library(org.Mm.eg.db)
library(KEGGREST)

library(igraph)
library(magrittr)

set.seed(1)

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")
source("RIO.R")

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
    pfix = prefix()
  if (length(pfix) != 0) {
     pfix <<- paste(pfix, "/", sep="")
  }
}

run <- function() {}

output <- function(outputfile) {
pdf(outputfile)

  fella.data <<- loadKEGGdata(
    databaseDir = paste(pfix, parameters["database", 2], sep="/"),#tmpdir,
    internalDir = FALSE,
    loadMatrix = "none"
)


analysis.nafld <- readRDS(paste(pfix, parameters["enriched", 2], sep="/"))
cpd.new.fig6 <- readSequential(paste(pfix, parameters["compounds", 2], sep="/"))

pscores.nafld <- getPscores(
    object = analysis.nafld,
    method = "diffusion")

## -----------------------------------------------------------------------------
id.cpd <- getCom(fella.data, level = 5, format = "id") %>% names
wilcox.test(
    x = pscores.nafld[cpd.new.fig6], # metabolites from fig6
    y = pscores.nafld[setdiff(id.cpd, cpd.new.fig6)], # rest of metabolites
    alternative = "less")

}
