#'
#' Learning dataset
#'
#' Contains factorization model results calculated using \code{MOFA2} and
#' initialized for transfer learning: \code{Fctrzn} and \code{Fctrzn_init}.
#'
#' \describe{
#'  \item{Fctrzn}{output of the factorization analysis of the learning
#'  dataset}
#'  \item{Fctrzn_init}{output of the factorization analysis of the
#'  learning dataset and the initialized values for the transfer learning}
#' }
#'
#' @format List of two \code{MOFA} objects: \code{Fctrzn} and
#' \code{Fctrzn_init}.
#'
#' MOFA2 object is a S4 class and contains the following main slots (the ones
#' relevant and used for the transfer learning):
#' \describe{
#'     \item{data}{input data used for the factorization analysis (mRNA, miRNA, DNAme and SNV)}
#'     \item{samples_metadata}{sample metadata (i.e. sample names and group)}
#'     \item{features_metadata}{features metadata (feature identifiers and views)}
#'     \item{expectations}{expected values of the factors and the loadings}
#'     \item{training_stats}{model training statistics}
#'     \item{data_options}{data processing options}
#'     \item{model_options}{model options}
#'     \item{training_options}{model training options}
#'     \item{dimensions}{dimensions of the model (e.g. \code{M} number of views, \code{N} number of samples, \code{D} number of features)}
#' }
#'
#' For more information about the structure of MOFA object, see the
#' [\code{MOFA2}](https://bioconductor.org/packages/release/bioc/manuals/MOFA2/man/MOFA2.pdf) vignette.
#'
#' @usage data("Lrn")
"Lrn"

#'
#' Target datasets
#'
#' Contains a list of target datasets used in different step in the transfer
#' learning workflow and there associated metadata.
#'
#' @format
#' \describe{
#' \item{YTrg_list}{list of the target dataset - Samples are in columns and
#' features are in rows.
#'     \describe{
#'     \item{\code{mRNA}: data stored into a \code{\link{SummarizedExperiment}} object}
#'     \item{\code{miRNA}: data stored into a \code{\link{SummarizedExperiment}} object}
#'     \item{\code{DNAme}: data stored into a \code{\link{SummarizedExperiment}} object}
#'     \item{\code{SNV}: data stored into a \code{matrix}}
#'     }}
#' \item{Trg_meta}{list of metadata:
#'     \describe{
#'     \item{Five \code{character} - (smpls, ftrs_mRNA, ftrs_miRNA, ftrs_DNAme, ftrs_SNV)}
#'     \item{Four \code{data.frame} with three variables (brcds, submittor, prjct) - (brcds_mRNA, brcds_miRNA, brcds_DNAme, brcds_SNV)}
#'     \item{Six \code{integer} - (Seed, ElbowK_Total, ElbowK_mRNA, ElbowK_miRNA, ElbowK_DNAme, ElbowK_SNV)}
#'     \item{One \code{logical} - (if_vst)}
#'     \item{Four \code{numeric} - (PCVarPrcnt_mRNA, PCVarPrcnt_miRNA, PCVarPrcnt_DNAme, PCVarPrcnt_SNV)}
#'     }}
#' \item{brcds_SS}{a list of 4 list of \code{data.frame} with the sample names.}
#' \item{YTrg_prep}{list of the prepared input target dataset (mRNA,
#' miRNA, DNAme, SNV). Samples are in columns and features are in rows.}
#' }
#'
#' @usage data("Trg")
"Trg"

#'
#' Transfer learning parameters
#'
#' Contains a list of input variables used for the transfer learning.
#'
#' @format
#' \describe{
#'  \item{YTrg}{list of the prepared input target dataset (mRNA,
#' miRNA, DNAme, SNV). Samples are in columns and features are in rows.}
#'  \item{Fctrzn_Lrn_W0}{list of 4 variables (mRNA, miRNA, DNAme and SNV).
#'  Each contains a W0 vector named with the corresponding feature names.}
#'  \item{Fctrzn_Lrn_W}{list of 4 \code{data.frame} (mRNA, miRNA, DNAme and
#' SNV). Each contains factors in columns and features in rows.}
#'  \item{Fctrzn_Lrn_WSq}{list of 4 \code{data.frame} (mRNA, miRNA, DNAme and
#' SNV). Each contains factors in columns and features in rows.}
#'  \item{Tau}{list of 4 \code{data.frame} (mRNA, miRNA, DNAme and
#' SNV). Each contains features in columns.}
#'  \item{TauLn}{list of 4 \code{data.frame} (mRNA, miRNA, DNAme and
#' SNV). Each contains features in columns.}
#'  \item{ZVar}{\code{data.frame} with factors in columns.}
#'  \item{ZMu}{\code{data.frame} with factors in columns and samples in
#' rows.}
#'  \item{ZMu_0}{vector of numerics}
#'  \item{ZMuSq}{\code{data.frame} with factors in columns.}
#' }
#'
#' @usage data("TL_param")
#'
"TL_param"
