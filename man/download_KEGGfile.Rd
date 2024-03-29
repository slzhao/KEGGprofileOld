% Generated by roxygen2 (4.0.2): do not edit by hand
\name{download_KEGGfile}
\alias{download_KEGGfile}
\title{download_KEGGfile}
\usage{
download_KEGGfile(pathway_id = "00010", species = "hsa",
  target_dir = getwd())
}
\arguments{
\item{pathway_id}{the KEGG pathway id, such as '00010'}

\item{species}{the species id in KEGG database, 'hsa' means human, 'mmu' means mouse, 'rno' means rat, etc}

\item{target_dir}{the local directory where the downloaded files are saved}
}
\description{
The function download XML files and png files from KEGG website to local disk
}
\details{
If pathway_id is set as 'all', all KEGG pathway ids in KEGG.db package will be used and downloaded from KEGG website
}
\examples{
download_KEGGfile(pathway_id="00010",species='hsa')
}

