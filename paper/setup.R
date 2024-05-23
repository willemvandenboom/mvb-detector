## Install the required packages.
for (tmp in c("abind", "brea", "nnet", "pbapply", "remotes")) {
  if(!tmp %in% rownames(installed.packages())) install.packages(
    pkgs = tmp, repos = "https://cloud.r-project.org", dependencies = TRUE
  )
}

# Install `mvb.detector` from GitHub.
remotes::install_github("willemvandenboom/mvb-detector")
