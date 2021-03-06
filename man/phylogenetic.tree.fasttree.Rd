% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylogenetic.tree.fasttree.R
\name{phylogenetic.tree.fasttree}
\alias{phylogenetic.tree.fasttree}
\title{Construct and calibrate a phylogenetic tree from sequence data}
\usage{
phylogenetic.tree.fasttree(dir.tree = dirfasttree, dir.seq = dirseqgen,
  fasttree.tool = "FastTree", calendar.dates = "dates.csv",
  simseqfile = "C.Epidemic.sequence.fasta", count.start = 1977,
  endsim = 40)
}
\arguments{
\item{dir.tree}{Direcotry where we find the FastTree tool}

\item{dir.seq}{Direcotry where we find sequences and sampling times files, it may be in the folder where simulation of sequences by \code{\link{sequence.simulation.seqgen()}} taken place}

\item{fasttree.tool}{Name of the compiled file of FastTree}

\item{calendar.dates}{File containing named sampling calendar times for each sequence}

\item{simseqfile}{File containing sequences}

\item{count.start}{Calendar year when the simulation started}

\item{endsim}{Number of years when the simulation was done}
}
\value{
A time-stamped phylogenetic  tree with annotation of internal nodes dates
}
\description{
Construct and calibrate a phylogenetic tree from sequence data
}
