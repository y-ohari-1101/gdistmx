fasta <- ape::read.dna("./path/to/fasta_file.fasta",format = "fasta")
group <- read.csv("./path/to/csv_file")

group


bw_g_dist <- gdistmx::dist_between_group(fasta, group)
