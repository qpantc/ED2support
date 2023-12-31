
# library(PEcAn.photosynthesis)
# library(PEcAn.MA)
# library(PEcAn.emulator)




parse.history <- function(historyfile, outfile = "") {
  hist <- XML::xmlToList(XML::xmlParse(historyfile))
  keys <- names(hist$pft)

  cat(paste(keys, collapse = ","), sep = "\n", file = outfile, append = FALSE)
  for (pft in hist) {
    cat(paste(lapply(pft[keys], stringr::str_trim), collapse = ","), sep = "\n", file = outfile, append = TRUE)
  }
} # parse.history


parse.history('./example/histo.xml',outfile = './example/histo.csv')



# meta-analysis ---------------------------------------------------------------------------------------------------

library(PEcAn.MA)

con <- 'ebifarm.pavi'
pft <- "temperate.Early_Hardwood"
pft_id <- PEcAn.DB::db.query("SELECT id FROM pfts WHERE name = $1", con,
                             values = list(pft))[[1]]
traits <- c("SLA", "Vcmax")
trait_string <- paste(shQuote(traits), collapse = ",")

# Load traits and priors from BETY
species <- PEcAn.DB::query.pft_species(pft, con = con)
trait.data <- PEcAn.DB::query.traits(species[["id"]], c("SLA", "Vcmax"), con = con)
prior.distns <- PEcAn.DB::query.priors(pft_id, trait_string, con = con)

# Pre-process data
jagged.data <- lapply(trait.data, PEcAn.MA::jagify)
taupriors <- list(tauA = 0.01,
                  tauB = c(SLA = 1000, Vcmax = 1000))
result <- pecan.ma(jagged.data, prior.distns, taupriors,
                   j.iter = 5000, outdir = tempdir())
