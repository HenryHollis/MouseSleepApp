# Package names
packages <- c("tidyverse", "readxl", "pracma", "shiny", "reader")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

install.packages("broom", type = "binary")

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

path = ""
runApp(path)
