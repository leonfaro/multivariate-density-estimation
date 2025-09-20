if (!exists("initialize_repo")) {
  source(file.path(getwd(), "R", "loader.R"))
}
root_path <- initialize_repo()
source(file.path(root_path, "experiments", "halfmoon", "halfmoon_pipeline.R"))

if (sys.nframe() == 0L) {
  invisible(NULL)
}
