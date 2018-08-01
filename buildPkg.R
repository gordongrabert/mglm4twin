#=======================================================================
# Script to Check, Build and Distribute the `mglm4twin` Package
#
#                                                   mglm4twin Core Team
#=======================================================================
#-----------------------------------------------------------------------
# Packages.

library(devtools)
library(withr)
library(knitr)

# Load the package (to make functions available).
load_all()

# Show all exported objects.
ls("package:mglm4twin")
packageVersion("mglm4twin")

# How many objects in each class.
table(sapply(ls("package:mglm4twin"),
             function(x) class(eval(parse(text=x)))))

#-----------------------------------------------------------------------
# Check.

load_all()

# Create/update NAMESPACE, *.Rd files.
document()

# Check documentation.
check_man()

# Check functions, datasets, run examples, etc. With check_dir = "../",
# it will create a directory named mcglm.Rcheck (right beside this root
# directory) with all the logs, manuals, figures from examples, etc.
check(manual = TRUE, vignettes = TRUE, check_dir = "../",
      cran = TRUE)

#-----------------------------------------------------------------------
# Build the package (it will be one directory up).

build(manual = TRUE, vignettes = TRUE)
# build the binary version for windows (not used)
# build_win()

#-----------------------------------------------------------------------
# Package vignette (use only if necessary)
# Based on: http://r-pkgs.had.co.nz/vignettes.html

# Create the vignette template. Do just once.
# use_vignette("UniModels")

#-----------------------------------------------------------------------
# Generate the README.md to update the GitHub initial page

knit(input = "README.Rmd")

#-----------------------------------------------------------------------
# Examples.

# Run examples from all functions of the package
# run_examples()
# Run examples from a specific function
# dev_example("yscale.components.right")

#-----------------------------------------------------------------------
# Test installation 1: Install from the local .tar.gz.

libTest <- path.expand("~/R-test/")
if (file.exists(libTest)) {
    unlink(libTest, recursive = TRUE)
}
dir.create(path = libTest)

# Install with install.packages() from the .tar.gz. created by build().
pkg <- paste0("../mglm4twin_", packageVersion("mglm4twin"), ".tar.gz")

# Install in a temporary directory.
install.packages(pkg, repos = NULL, lib = libTest)
library(package = "mglm4twin", lib.loc = libTest)
packageVersion("mglm4twin")
ls("package:mglm4twin")
## Before removing it, detach it
detach("package:mglm4twin")

#-----------------------------------------------------------------------

## Create package tarballs
load_all()
pkg <- paste0("../mglm4twin_", packageVersion("mglm4twin"), ".tar.gz")
pkg.win <- paste0("../mglm4twin_", packageVersion("mglm4twin"), ".zip")

## Build the *.zip (Windows version)
cmd.win <- paste("cd ../mglm4twin.Rcheck && zip -r", pkg.win, "mglm4twin")
system(cmd.win)

## PDF manual and network graph
#ntw <- "./data-raw/mcglm_network.html"
man <- "../mglm4twin.Rcheck/mglm4twin-manual.pdf"

##----------------------------------------------------------------------
## Sending package tarballs and manual to remote server to be
## downloadable.
## URL: http://www.leg.ufpr.br/~leg/mglm4twin/

## Send to LEG server
## NOTE: "PATAXO" and "PATAXOP" are exported names in .bashrc (with IP
## and port, respectivelly)
cmd <- paste("scp", pkg, man, pkg.win,
             "leg@leg.ufpr.br:~/public_html/mglm4twin/source")
system(cmd)
#browseURL("http://www.leg.ufpr.br/~leg/mglm4twin/")
### PRECISO VER ISSO

##----------------------------------------------------------------------
## Send to downloads/ folder, so it stays hosted on GitHub
#dest <- "downloads/"
#file.copy(c(pkg, pkg.win, man), dest, overwrite = TRUE)

#-----------------------------------------------------------------------
