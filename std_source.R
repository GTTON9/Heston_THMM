source("Visualization/std_Viterbi_Visual.R")
source("Visualization/std_BW_Visual.R")

r_files <- list.files(
  "Heston_Standard",
  pattern = "\\.R$",
  recursive = TRUE,
  full.names = TRUE
)


r_files <- sort(r_files)
r_files
invisible(lapply(r_files, function(f) {
  message("Sourcing: ", f)
  source(f, local = .GlobalEnv)
}))
