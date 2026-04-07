#'
#' Learning dataset
#'
#' Contains factorization model results calculated using MOFA2 and init for TL
#'
#' \code{Fctrzn}: output of the factorization analysis of the learning
#'      dataset
#' \code{Fctrzn_init}: output of the factorization analysis of the learning
#'      dataset and the initialized values for the transfer learning.
#'  MOFA2 object is a S4 class and contains these main slots (relevant for the
#'  transfer learning):
#'
#' @format List of two MOFA object: \code{Fctrzn} and \code{Fctrzn_init}.
#' \describe{
#' \item{\code{data}}{input data used for the factorization analysis (mRNA, miRNA, DNAme and SNV)}
#' \item{\code{samples_metadata}}{sample metadata (sample names and group)}
#' \item{\code{features_metadata}}{features metadata (feature identifier and view)}
#' \item{\code{expectations}}{expected values of the factors and the loadings}
#' \item{\code{training_stats}}{model training statistics}
#' \item{\code{data_options}}{data processing options}
#' \item{\code{model_options}}{model options}
#' \item{\code{training_options}}{model training options}
#' \item{\code{dimensions}}{dimensions of the model (e.g. \code{M} number of views, \code{N} number of samples, \code{D} number of features)}
#' }
#'  For more information about the structure of MOFA object, see the
#'  \link{\code{MOFA}} vignette.
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
#' \code{YTrg_list}: list of the target dataset - \code{mRNA}, \code{miRNA}
#'  and \code{DNAme} data are \code{\link{SummarizedExperiment}} object and
#' \code{SNV} data are stored in a \code{matrix}. Samples are in columns and
#'  features are in rows.
#' \code{Trg_meta}: list of metadata - Five \code{character} (smpls,
#' ftrs_mRNA, ftrs_miRNA, ftrs_DNAme, ftrs_SNV), four \code{data.frame}
#' (brcds_mRNA, brcds_miRNA, brcds_DNAme, brcds_SNV) with three variables
#' (brcds, submittor, prjct),  six \code{integer} (Seed, ElbowK_Total,
#' ElbowK_mRNA, ElbowK_miRNA, ElbowK_DNAme, ElbowK_SNV), one \code{logical}
#' variable (if_vst) and four \code{numeric} (PCVarPrcnt_mRNA, PCVarPrcnt_miRNA,
#' PCVarPrcnt_DNAme, PCVarPrcnt_SNV).
#' \code{brcds_SS}: a list of 4 list of \code{data.frame} with the sample
#' names.
#' \code{YTrg_prep}: list of the prepared input target dataset (mRNA,
#' miRNA, DNAme, SNV). Samples are in columns and features are in rows.
#'
#' @usage data("Trg")
"Trg"

#'
#' Transfer learning parameters
#'
#' Contains a list of input variables used for the transfer learning.
#'
#' @format
#' \code{YTrg}: list of the prepared input target dataset (mRNA,
#' miRNA, DNAme, SNV). Samples are in columns and features are in rows.
#' \code{Fctrzn_Lrn_W0}: list of 4 variables (mRNA, miRNA, DNAme and SNV).
#' Each contains a W0 vector named with the corresponding feature names.
#' \code{Fctrzn_Lrn_W}: list of 4 \code{data.frame} (mRNA, miRNA, DNAme and
#' SNV). Each contains factors in columns and features in rows.
#' \code{Fctrzn_Lrn_WSq}: list of 4 \code{data.frame} (mRNA, miRNA, DNAme and
#' SNV). Each contains factors in columns and features in rows.
#' \code{Tau}: list of 4 \code{data.frame} (mRNA, miRNA, DNAme and
#' SNV). Each contains features in columns.
#' \code{TauLn}: list of 4 \code{data.frame} (mRNA, miRNA, DNAme and
#' SNV). Each contains features in columns.
#' \code{ZVar}: \code{data.frame} with factors in columns.
#' \code{ZMu}:  \code{data.frame} with factors in columns and samples in
#' rows.
#' \code{ZMu_0}: vector of numeric.
#' \code{ZMuSq}: \code{data.frame} with factors in columns.
#'
#' @usage data("TL_param")
#'
"TL_param"
