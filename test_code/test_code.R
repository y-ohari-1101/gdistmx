install.packages("devtools")
library(devtools)
install.packages("roxygen2")
library(roxygen2)
install.packages("ape")
library(ape)
install.packages("pacman")
library(pacman)

dna <- ape::read.dna("./example/input_fasta/COX1.fasta", format = "fasta")
dna

install.packages("dplyr")
library(dplyr)

gpr_data <- read.csv("./example/sample_group.csv", header = TRUE)
gpr_data

single_d_w_grp <- dist_within_group(single_dna,gpr_data,model = "raw")
single_d_w_grp

class(single_d_w_grp)
class(single_d_w_grp$summary)
class(single_d_w_grp$pairwise)
str(single_d_w_grp)


search()
detach(package:ape, unload=TRUE)
detach(package:dplyr, unload=TRUE)





#' @importFrom dplyr bind_rows

packageVersion("dplyr")
version()

#' @examples
#' dna <- ape::read.dna("./example/input_fasta/COX1.fasta", format = "fasta")
#' group_def <- read.csv("./example/sample_group.csv", header = TRUE)
#' dist_within_group(dna, group_def, model = "raw")
#'
