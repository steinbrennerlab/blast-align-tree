# install_r_deps.R
# ----------------
# Cross-platform installer for BAT R dependencies
# Compatible with macOS, Windows, Linux

required_cran <- c("ggplot2", "optparse", "tidytree", "ape", "broom")
required_bioc <- c("treeio", "ggtree", "Biostrings", "phytools")

# ---- Check R installation ----
if (!requireNamespace("utils", quietly = TRUE)) {
  stop(
    "\nR does not appear to be installed or accessible on your PATH.\n",
    "Please install R from https://cran.r-project.org/ and try again.\n"
  )
}

# ---- Install BiocManager if missing ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# ---- Ensure ggplot2 < 4.0.0 ----
install_compatible_ggplot <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2", repos = "https://cloud.r-project.org", dependencies = TRUE)
  }
  ver <- utils::packageVersion("ggplot2")
  if (ver >= "4.0.0") {
    message("ggplot2 version too new, reinstalling compatible version <4.0.0...")
    install.packages("ggplot2", repos = "https://cloud.r-project.org", dependencies = TRUE)
  }
}
install_compatible_ggplot()

# ---- Ensure ggtree < 4.0.0 ----
install_compatible_ggtree <- function() {
  if (!requireNamespace("ggtree", quietly = TRUE)) {
    suppressMessages(BiocManager::install("ggtree", ask = FALSE))
  }
  ver <- utils::packageVersion("ggtree")
  if (ver >= "4.0.0") {
    message("ggtree version too new, reinstalling compatible version <4.0.0 via BiocManager...")
    suppressMessages(BiocManager::install("ggtree", ask = FALSE))
  }
}
install_compatible_ggtree()

# ---- Install remaining CRAN packages ----
for (pkg in setdiff(required_cran, "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", dependencies = TRUE)
  }
}

# ---- Install remaining Bioconductor packages ----
for (pkg in setdiff(required_bioc, "ggtree")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    suppressMessages(BiocManager::install(pkg, ask = FALSE))
  }
}

# ---- Final message ----
message(sprintf(
  "\nAll required R packages are installed and compatible.\n  ggplot2: %s\n  ggtree: %s\n",
  utils::packageVersion("ggplot2"),
  utils::packageVersion("ggtree")
))
