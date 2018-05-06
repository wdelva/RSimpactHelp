
#' Subset sequences in fasta file based on sequence ID or Name
#'
#' @param pool.seq.file File of sequence data (in fasta format) within which we want to subset
#' @param select.vec Vector of ID or name of sequences we want to select
#' @param name.file Name of fasta file which containes selected sequences
#' @return Files of elcted sequences data
#' @importFrom Rsamtools FaFile
#' @importFrom Rsamtools getSeq
#' @importFrom Rsamtools indexFa
#' @importFrom Rsamtools seqinfo
#' @importFrom ape write.dna
#' @export


choose.sequence.ind <- function(pool.seq.file = "pool.seq.fasta",
                                select.vec = select.vec,
                                name.file = "selected.ind.fasta"){

  sample.name.seq <- select.vec

  Rsamtools::indexFa(paste0(pool.seq.file)) # create an index of file 'file.fasta'

  fa <- Rsamtools::FaFile(paste0(pool.seq.file))  # reference the fasta file and it's index

  gr <- as(Rsamtools::seqinfo(fa), "GRanges") # disover coordinate for each sequence

  seq <- Rsamtools::getSeq(fa, gr[sample.name.seq])

  ape::write.dna(seq, file = paste0(name.file), format = "fasta")

  print(paste("Sequence of ", paste(sample.name.seq), "was considered"))

}
