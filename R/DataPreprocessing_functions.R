## --------------------- TARGET DATASET PREPARATION FUNCTIONS ---

GeoMeans_Lrn_init <- function(view, expdat_meta_Lrn, YTrgFtrs) {
    #'
    #' Retrieve the Geometric means calculated for the learning dataset during
    #' counts normalization
    #'
    #' @param view current view data name
    #' @param expdat_meta_Lrn list of learning dataset factorization metadata
    #' @param YTrgFtrs feature names of the current view
    #'
    #' @returns precalculated Geometric means of the learning dataset
    #'
    #' @examples
    #'
    #' data("Lrn", package = "MOTL")
    #' data("Trg", package = "MOTL")
    #'
    #' expdat_meta_Lrn <- Lrn$Fctrzn@data
    #' YTrgFtrs <- Trg$Trg_meta$ftrs_mRNA
    #'
    #' GeoMeans_Lrn <-
    #'         GeoMeans_Lrn_init(view = "mRNA", expdat_meta_Lrn, YTrgFtrs)
    #'
    #' @export

    if (is.element(view, c("mRNA", "miRNA"))) {
        GeoMeans_Lrn <- expdat_meta_Lrn[[paste0("GeoMeans_", view)]]
        FtrsKeep <- is.element(names(GeoMeans_Lrn), YTrgFtrs)
        GeoMeans_Lrn <- GeoMeans_Lrn[FtrsKeep]
        GeoMeans_Lrn <- GeoMeans_Lrn[match(YTrgFtrs, names(GeoMeans_Lrn))]
    } else {
        GeoMeans_Lrn <- numeric()
    }
    return(GeoMeans_Lrn)
}

GeoMeanFun <- function(x) {
    #'
    #' Calculate the Geometric mean of a vector
    #'
    #' @param x vector of numeric values
    #'
    #' @return Geometric mean of vector x
    #'
    #' @examples
    #' x <- c(125,12,4545,7878,6777,454545,88979)
    #' GeoMeans <- GeoMeanFun(x)
    #'
    #' @export

    GeoMeans <- exp(sum(log(x[x > 0])) / length(x))
    return(GeoMeans)
}

countsNormalization <- function(expdat, GeoMeans) {
    #' Normalize counts data
    #'
    #' Normalize counts data using DESeq2 normalization.
    #' Two ways of normalization:
    #' 1. Use the pre-calculated Geometric means of the learning dataset
    #' 2. Use calculated Geometric means of the \code{expdat} dataset given
    #' in input
    #'
    #' If \code{is.numeric(GeoMeans) == TRUE}, input data are normalized with
    #' pre-calculated Geometric means (from learning dataset).
    #' If non values are provided, Geometric means is calculated based on the
    #' input dataset using \code{\link{GeoMeanFun}} function.
    #' Then, the input dataset is normalized using these Geometric means.
    #'
    #' @param expdat SE object of experimental data (could be miRNA or mRNA)
    #' @param GeoMeans if it's a character, Geometric means will be calculated
    #' for the \code{expdat} variable given in input (learning or target
    #' dataset).
    #' If it's a numerical vector, given Geometric means will be used for
    #' the normalization (the ones pre-calculated from the learning dataset)
    #'
    #' @returns list of data.frame of the counts normalized and Geometric means
    #' calculated
    #'
    #' @import DESeq2
    #'
    #' @examples
    #'
    #' ## Create a matrix with "counts" data
    #' ## Then, create a summarized experiment object
    #' expdat <- matrix(rexp(200, rate = .1), ncol = 20)
    #' expdat <- apply(expdat, MARGIN = 2, round)
    #' expdat <- SummarizedExperiment::SummarizedExperiment(expdat)
    #'
    #' ## With "newGeoMeans", geometric means will be calculated based on the
    #' ## input matrix
    #' GeoMeans <- "newGeoMeans"
    #'
    #' expdat_counts_norm <- countsNormalization(expdat, GeoMeans)
    #'
    #' @export

    ## Create DESeq object
    expdat_dds <- DESeqDataSet(expdat, design = ~1)

    if(is.numeric(GeoMeans)){
        expdat_dds_norm <- estimateSizeFactors(expdat_dds, geoMeans = GeoMeans)
        GeoMeans <- NULL
    }else{
        GeoMeans <- as.vector(apply(counts(expdat_dds), 1, GeoMeanFun))
        expdat_dds_norm <- estimateSizeFactors(expdat_dds, geoMeans = GeoMeans)
    }

    ## Extract normalized counts
    expdat_counts_norm <-
        list("counts" = counts(expdat_dds_norm, normalized = TRUE))

    ## Save Geometric means
    expdat_counts_norm$GeoMeans <- GeoMeans

    return(expdat_counts_norm)
}

countsTransformation <- function(expdat_count, TopD) {
    #' Log2 transform and select top data based on variance
    #'
    #' @param expdat_count data.frame of the counts
    #' @param TopD number of features to keep
    #'
    #' @returns data.frame of the log2 transformed and filtered data
    #'
    #' @examples
    #'
    #' ##
    #' expdat_count <- matrix(rexp(200, rate = .1), ncol = 20)
    #' expdat_count <- apply(expdat_count, MARGIN = 2, round)
    #'
    #' ## input matrix
    #' TopD <- 20
    #'
    #' expdat_counts_fltr <- countsTransformation(expdat_count, TopD)
    #'
    #'
    #' @export

    ## log transform and filter to keep only most variable
    expdat_counts_log <- log2(expdat_count + 1)
    FtrsKeep <-
        base::rank(-matrixStats::rowVars(expdat_counts_log,
                                         na.rm = TRUE,
                                         useNames = FALSE),
                   ties.method = "first"
        ) <= TopD
    expdat_counts_fltr <- expdat_counts_log[FtrsKeep, ]
    return(expdat_counts_fltr)
}

preprocessCountsData <- function(view,
                                 YTrg_list,
                                 normalization = FALSE,
                                 expdat_meta_Lrn,
                                 transformation = FALSE) {
    #'
    #' Preprocess counts data
    #'
    #' Counts data (i.e. mRNA and miRNA) can be normalized and/or transformed.
    #'
    #' Normalization is performed using the \code{\link{countsNormalization}}
    #' function with pre-calculated or new calculated Geometric means.
    #'
    #' Transformation is performed using the \code{\link{countsTransformation}}
    #' function with log2.
    #'
    #' @param view a data view name vector (i.e. mRNA or miRNA)
    #' @param YTrg_list a named list of target data. Names correspond to the
    #' defined views. The list contains matrices.
    #' @param normalization if FALSE, no normalization. If "LrnGeoMeans",
    #' normalization using the pre-calculated Geometric means. If "newGeoMeans",
    #' normalization using Geometric means from dataset. By default, it's set
    #' to FALSE.
    #' @param expdat_meta_Lrn the list of learning set factorization metadata
    #' @param transformation if FALSE, no transformation. If TRUE, log2
    #' normalization.
    #'
    #' @returns Preprocessed counts data for the current view
    #'
    #' @examples
    #'
    #' data("Lrn", package = "MOTL")
    #' data("Trg", package = "MOTL")
    #'
    #' expdat_meta_Lrn <- Lrn$Fctrzn@data
    #' YTrg_list <- Trg$YTrg_list
    #'
    #' mRNA <- preprocessCountsData(view = "mRNA", YTrg_list = YTrg_list,
    #'                            normalization = FALSE,
    #'                            expdat_meta_Lrn = expdat_meta_Lrn,
    #'                            transformation = FALSE)
    #'
    #'
    #' @export

    ## Select current data
    YTrg <- YTrg_list[[view]]

    if (view %in% c("mRNA", "miRNA")) {
        ## Normalization
        ## if "LrnGeoMeans", retrieve the learning set geometric means and
        ## normalize with it
        ## if "newGeoMeans", calculate the Geometric means of the input data
        if (is.character(normalization)) {
            if (normalization == "LrnGeoMeans") {
                message("Normalize using the pre-calculated learning dataset
                        Geometric means")
                GeoMeans <-
                    GeoMeans_Lrn_init(view, expdat_meta_Lrn, rownames(YTrg))
            }
            if (normalization == "newGeoMeans") {
                message("Calculate Geometric means during the normalization")
                GeoMeans <- normalization
            }
            YTrg <- SummarizedExperiment::SummarizedExperiment(assays = list(YTrg))
            YTrg <- countsNormalization(expdat = YTrg, GeoMeans = GeoMeans)
            YTrg <- YTrg$counts
        } else {
            message("No normalization")
        }

        ## Transformation
        if (transformation) {
            message("Log transform data")
            YTrg <- countsTransformation(expdat_count = YTrg, TopD = nrow(YTrg))
        } else {
            message("No transformation")
        }
    }
    return(YTrg)
}

TargetDataPreparation <- function(views,
                                  YTrg_list,
                                  Fctrzn,
                                  smpls,
                                  expdat_meta_Lrn,
                                  normalization = FALSE,
                                  transformation = FALSE) {
    #'
    #' Target data preparation for transfer learning
    #'
    #' The function follows these steps:
    #' 1. Prepare target data for each view
    #' 2. Normalize and/or transform counts data
    #'
    #'
    #' Preparation of data consists on
    #'
    #' - removing features with variance equal to zero
    #' - features harmonization between target and learning data
    #' - column ordering between views.
    #'
    #' Preparation is perform using the
    #' \code{\link{TargetDataPrefiltering}} function. See the documentation
    #' for more details.
    #'
    #' It could be possible to normalize and or transform counts data using the
    #' \code{\link{preprocessCountsData}} function. Normalization can be done
    #' using Geometric means (from learning or target dataset) and transformation
    #' is a log2 transformation of the counts.
    #'
    #' @param views a list of target data views (e.g. \code{c("mRNA", "miRNA")})
    #' @param YTrg_list a named list of target set data. Names correspond to the
    #' defined views. The list contains matrices.
    #' @param Fctrzn the learning factorization model object (from \code{MOFA})
    #' @param smpls a vector of sample names (i.e. column names of the
    #' \code{YTrg_list})
    #' @param normalization if FALSE, no normalization. If "LrnGeoMeans",
    #' normalization using the learning Geometric means. If "newGeoMeans",
    #' Geometric means are calculated using the input data. By default it's
    #' set to FALSE.
    #' @param expdat_meta_Lrn the list of learning set factorization metadata
    #' @param transformation if FALSE, no transformation. If TRUE, log2
    #' transformation of counts data. By default it's set to FALSE.
    #'
    #' @returns a list of prepared target data
    #'
    #' @examples
    #'
    #' data("Lrn", package = "MOTL")
    #' data("Trg", package = "MOTL")
    #'
    #' YTrg_list <- Trg$YTrg_prep
    #' Fctrzn <- Lrn$Fctrzn
    #' smpls <- colnames(YTrg_list$mRNA)
    #' expdat_meta_Lrn <- Lrn$Lrn_meta
    #'
    #'
    #' YTrg_prep <- TargetDataPreparation(views = c("mRNA", "miRNA", "DNAme"),
    #'                                     YTrg_list = YTrg_list,
    #'                                     Fctrzn = Fctrzn, smpls = smpls,
    #'                                     expdat_meta_Lrn = expdat_meta_Lrn,
    #'                                     normalization = FALSE,
    #'                                     transformation = FALSE)
    #'
    #'
    #' @export

    names(views) <- views

    ## Feature variance prefiltering and feature harmonization
    message("Feature prefiltering")
    YTrg <- lapply(views, TargetDataPrefiltering, YTrg_list, Fctrzn, smpls)

    ## Normalization and transformation
    message("Normalization and transformation")
    YTrg <- lapply(
        views,
        preprocessCountsData,
        YTrg,
        normalization,
        expdat_meta_Lrn,
        transformation
    )

    ## Return
    return(YTrg)
}

initTransferLearningParamaters <- function(YTrg,
                                           views,
                                           Fctrzn,
                                           likelihoods) {
    #'
    #' Transfer learning parameters and data objects initialization
    #'
    #' The function performs the following steps:
    #' 1. Extract the factorized learning set weight intercepts \code{W0}
    #' 2. Extract the factorized learning set weights \code{W}
    #' 3. Extract the factorized learning set squared weights \code{Wsq}
    #' 4. Extract the learning set \code{Tau} and log(Tau) \code{TauLn}
    #' For each extracted parameter, common features between learning dataset
    #' and target dataset are kept. Then target data \code{YTrg}, \code{Tau}
    #' and \code{TauLn} are transposed.
    #'
    #' Each parameter are extracted from the \code{Fctrzn} model created using
    #' \code{MOFA2}.
    #'
    #' \code{YTrg} matrices should have the same columns order.
    #'
    #' @param YTrg a named list of target dataset matrices. Names correspond to
    #' the defined views.
    #' @param views a vector of target dataset view names (e.g.
    #' \code{c("mRNA", "miRNA")})
    #' @param Fctrzn the learning dataset factorization model object
    #' (from \code{MOFA})
    #' @param likelihoods a named list of data types. The list can contain
    #' \code{gaussian}, \code{poisson} or \code{bernoulli} depending of the data
    #' type. Names are the view names.
    #'
    #' @returns a list of initialized parameters for transfer learning
    #' 1. \code{YTrg} - the transposed named list of target data
    #' 2. \code{Fctrzn_Lrn_W0} - the Factorized learning set weight intercepts
    #' with same features as YTrg
    #' 3. \code{Fctrzn_Lrn_W} - Factorized learning set weights with same
    #' features as YTrg
    #' 4. \code{Fctrzn_Lrn_WSq} - Factorized learning set squared weights with
    #' same features as YTrg
    #' 5. \code{Tau} - the transposed Tau matrix with same features as YTrg
    #' 6. \code{TauLn} - the transposed log2(Tau) matrix with same features as
    #' YTrg
    #'
    #' @examples
    #'
    #' data("Lrn", package = "MOTL")
    #' data("Trg", package = "MOTL")
    #'
    #' views <- c("mRNA", "miRNA", "DNAme", "SNV")
    #' likelihoods <- c("mRNA" = "gaussian", "miRNA" = "gaussian",
    #'                 "DNAme" = "gaussian", "SNV" = "bernoulli")
    #' Fctrzn <- Lrn$Fctrzn_init
    #' YTrg <- Trg$YTrg_prep
    #'
    #' TLparameter <-
    #'     initTransferLearningParamaters(YTrg = YTrg,
    #'                                     views = views,
    #'                                     Fctrzn = Fctrzn,
    #'                                     likelihoods = likelihoods)
    #'
    #' @export
    #'

    names(views) <- views

    ## Feature names in each data
    YTrgFtrs <- lapply(YTrg, rownames)

    ## FACTORIZED LEARNING WEIGHTS MATRIX ZERO
    message("Factorized learning set weight intercepts")
    Fctrzn_Lrn_W0 <- lapply(views, function(view, Fctrzn, YTrgFtrs) {
        Fctrzn_Lrn_W0 <- Fctrzn@expectations[["W0"]][[view]]
        FtrsKeep <- is.element(names(Fctrzn_Lrn_W0), YTrgFtrs[[view]])
        Fctrzn_Lrn_W0 <- Fctrzn_Lrn_W0[FtrsKeep]
        Fctrzn_Lrn_W0 <-
            Fctrzn_Lrn_W0[match(YTrgFtrs[[view]], names(Fctrzn_Lrn_W0))]
        return(Fctrzn_Lrn_W0)
    }, Fctrzn, YTrgFtrs)

    ## FACTORIZED LEARNING WEIGHTS MATRIX
    message("Factorized learning set weights")
    Fctrzn_Lrn_W <- lapply(views, function(view, Fctrzn, YTrgFtrs) {
        Fctrzn_Lrn_W <- Fctrzn@expectations[["W"]][[view]]
        FtrsKeep <- is.element(rownames(Fctrzn_Lrn_W), YTrgFtrs[[view]])
        Fctrzn_Lrn_W <- Fctrzn_Lrn_W[FtrsKeep, ]
        Fctrzn_Lrn_W <-
            Fctrzn_Lrn_W[match(YTrgFtrs[[view]], rownames(Fctrzn_Lrn_W)), ]
        return(Fctrzn_Lrn_W)
    }, Fctrzn, YTrgFtrs)

    ## FACTORIZED LEARNING WEIGHTS MATRIX SQUARED
    message("Factorized learning set squared weights")
    Fctrzn_Lrn_WSq <- lapply(views, function(view, Fctrzn, YTrgFtrs) {
        Fctrzn_Lrn_WSq <- Fctrzn@expectations[["WSq"]][[view]]
        FtrsKeep <- is.element(rownames(Fctrzn_Lrn_WSq), YTrgFtrs[[view]])
        Fctrzn_Lrn_WSq <- Fctrzn_Lrn_WSq[FtrsKeep, ]
        Fctrzn_Lrn_WSq <-
            Fctrzn_Lrn_WSq[match(YTrgFtrs[[view]], rownames(Fctrzn_Lrn_WSq)), ]
        return(Fctrzn_Lrn_WSq)
    }, Fctrzn, YTrgFtrs)

    ## TAU PARAMETER
    message("Tau")
    Tau <- lapply(views, function(view, Fctrzn, YTrg, YTrgFtrs) {
        Tau <- Fctrzn@expectations[["Tau"]][[view]]$group0
        FtrsKeep <- is.element(rownames(Tau), YTrgFtrs[[view]])
        Tau <- Tau[FtrsKeep, ]
        Tau <- Tau[match(YTrgFtrs[[view]], rownames(Tau)), ]
        Tau <- matrix(
            rowMeans(Tau, na.rm = TRUE),
            nrow = dim(YTrg[[view]])[1],
            ncol = dim(YTrg[[view]])[2],
            byrow = FALSE
        )
        rownames(Tau) <- rownames(YTrg[[view]])
        return(Tau)
    }, Fctrzn, YTrg, YTrgFtrs)

    ## LOG TAU PARAMETER
    message("LOG Tau")
    TauLn <- lapply(views, function(view,
                                    likelihoods,
                                    Fctrzn,
                                    YTrg,
                                    YTrgFtrs) {
        if (likelihoods[view] == "gaussian") {
            TauLn <- Fctrzn@expectations[["TauLn"]][[view]]
            FtrsKeep <- is.element(names(TauLn), YTrgFtrs[[view]])
            TauLn <- TauLn[FtrsKeep]
            TauLn <- TauLn[match(YTrgFtrs[[view]], names(TauLn))]
            TauLn <- matrix(
                TauLn,
                nrow = dim(YTrg[[view]])[1],
                ncol = dim(YTrg[[view]])[2],
                byrow = FALSE
            )
            rownames(TauLn) <- rownames(YTrg[[view]])
        } else {
            TauLn <- numeric()
        }
        return(TauLn)
    }, likelihoods, Fctrzn, YTrg, YTrgFtrs)

    ## Transpose matrices where necessary - to make them samples x features
    YTrg <- lapply(YTrg, t)
    Tau <- lapply(Tau, t)
    TauLn <- lapply(TauLn, t)

    ## Return
    TL_param <- list(
        "YTrg" = YTrg,
        "Fctrzn_Lrn_W0" = Fctrzn_Lrn_W0,
        "Fctrzn_Lrn_W" = Fctrzn_Lrn_W,
        "Fctrzn_Lrn_WSq" = Fctrzn_Lrn_WSq,
        "Tau" = Tau,
        "TauLn" = TauLn
    )
    return(TL_param)
}

## --------------------------------------------------------------


## ----------------------------------- TCGA RELATED FUNCTIONS ---

TCGATargetDataPrefiltering <- function(view, brcds_SS, SS, YTrgFull, Fctrzn) {
    #'
    #' Filter out TCGA target subset data according variance
    #'
    #' This function performs a prefiltering analysis through the steps:
    #' 1. Extract data for the given sample names (\code{brcds_SS})
    #' 2. Remove features with variance equal to zero
    #' 3. Match features between target data set and learning data set
    #'
    #' The function return a pre-filtered target dataframe for the current view.
    #'
    #' @param view a character of current view name data
    #' @param brcds_SS a list of sample names for each view. The list is named
    #' according views (e.g. "brcds_mRNA_SS") and contains list of dataframes
    #' for each view. Each dataframe contains at least one column named "brcds"
    #' (the one used).
    #' @param SS current subset number
    #' @param YTrgFull a named list of target set data. Names correspond to the
    #' defined views. The list contains \code{SummarizedExperiment} for miRNA,
    #' mRNA and DNAme and \code{matrix} for SNV.
    #' @param Fctrzn learning factorization model object (from \code{MOFA})
    #'
    #' @returns the subset data for current view and SS number
    #'
    #' @examples
    #'
    #' # In the paper, several target datasets were created as subset of the
    #' # reference dataset R. This function was used to generate them
    #' # automatically.
    #' # If you are not doing the paper analysis, you can create the brcds_SS
    #' # using the following command line. Replace the "nameX" with the sample
    #' # names of your data.
    #' # You can as much as you want add dataframe on each view.
    #'
    #' brcds_SS_ex <-
    #'         list("brcds_mRNA_SS" =
    #'             list(data.frame("brcds" = c("name01", "name02"))),
    #'         "brcds_miRNA_SS" =
    #'             list(data.frame("brcds" = c("name10", "name11"))))
    #'
    #' # See the doc to create the input parameter
    #'
    #' data("Trg", package = "MOTL")
    #' data("Lrn", package = "MOTL")
    #'
    #' brcds_SS <- Trg$brcds_SS
    #' YTrg_list <- Trg$YTrg_list
    #' Lrn_Fctrzn <- Lrn$Fctrzn
    #'
    #'
    #' expdat_mRNA <- TCGATargetDataPrefiltering(view = "mRNA",
    #' brcds_SS = brcds_SS, SS = 1, YTrgFull = YTrg_list, Fctrzn = Lrn_Fctrzn)
    #' expdat_mRNA[c(1:5), c(1:5)]
    #'
    #'
    #' @export

    ## select samples and subset the YTrg
    brcds <- brcds_SS[[paste0("brcds_", view, "_SS")]][[SS]]
    SmplsKeep <- is.element(colnames(YTrgFull[[view]]), brcds$brcds)
    YTrgSS <- YTrgFull[[view]][, SmplsKeep]

    ## prefiltering, only condition is that variance >0
    if (is.element(view, c("mRNA", "DNAme", "miRNA"))) {
        FtrsKeep <- matrixStats::rowVars(SummarizedExperiment::assay(YTrgSS),
                                         na.rm = TRUE) > 0
    } else {
        FtrsKeep <- matrixStats::rowVars(YTrgSS, na.rm = TRUE) > 0
    }
    FtrsKeep[is.na(FtrsKeep)] <- FALSE
    YTrgSS <- YTrgSS[FtrsKeep, ]

    ## harmonize features between Trg SS and Lrn data
    features_metadata <- Fctrzn@features_metadata
    FtrsLrn <- features_metadata$feature[features_metadata$view == view]
    FtrsCommon <- FtrsLrn[is.element(FtrsLrn, rownames(YTrgSS))]
    FtrsKeep <- is.element(rownames(YTrgSS), FtrsCommon)
    YTrgSS <- YTrgSS[FtrsKeep, ]
    YTrgSS <- YTrgSS[match(FtrsCommon, rownames(YTrgSS)), ]

    return(YTrgSS)
}

TCGATargetDataPreparation <- function(views,
                                      YTrgFull,
                                      brcds_SS,
                                      SS,
                                      Fctrzn,
                                      smpls,
                                      normalization = FALSE,
                                      expdat_meta_Lrn,
                                      transformation = TRUE) {
    #'
    #' Prepare TCGA target dataset for transfer learning
    #'
    #' This function follows these steps:
    #' 1. Filter out features according variance
    #' 2. Reshape data into matrices
    #' 3. Order samples to have the same columns order between different views
    #' 4. Normalize and/or transform counts data
    #'
    #'
    #' First, samples included in the brcds_SS list are selected. Then features
    #' with variance equal to zero are removed. They are also remove if they are
    #' not retrieved in the learning dataset. These steps are performed using
    #' \code{\link{TCGATargetDataPrefiltering}} function.
    #'
    #' The mRNA, miRNA et DNAme data are stored into \code{SummarizedExperiment}
    #' object. For the next step, data have to be stored into a \code{matrix}.
    #' SNV data are already a matrix.
    #' Then, samples are ordered in the same way between views.
    #'
    #' Finally, counts data (e.g. mRNA and miRNA) can be normalized and/or
    #' transformed using \code{\link{preprocessCountsData}} function.
    #' - if \code{normalization = FALSE}: counts data are not normalized
    #' - if \code{normalization = "LrnGeoMeans"}: counts data are normalized
    #' using the learning dataset Geometric means calculated
    #' - if \code{normalization = "newGeoMeans"}: counts data are normalized
    #' using the geometric means calculated on the target dataset.
    #'
    #' Normalization is performed in the \code{\link{countsNormalization}}
    #' function using \code{estimateSizeFactors} from
    #' \code{\link{DESeq2}} package. And transformation is perform using
    #' \code{\link{countsTransformation}} with a log2 transformation.
    #'
    #' Look \code{\link{GeoMeans_Lrn_init}} and \code{\link{GeoMeanFun}} to see
    #' how learning Geometric means are calculated.
    #'
    #' @param views a list of target data views (e.g. \code{c("mRNA", "miRNA")})
    #' @param YTrgFull a named list of target set data. Names correspond to the
    #' defined views. The list contains \code{SummarizedExperiment} for miRNA,
    #' mRNA and DNAme and \code{matrix} for SNV.
    #' @param brcds_SS a list of sample names for each view. The list is named
    #' according views (e.g. "brcds_mRNA_SS") and contains list of dataframes
    #' for each view. Each dataframe contains at least one column named "brcds"
    #' (the one used).
    #' @param SS current subset number
    #' @param Fctrzn the learning factorization model object (from \code{MOFA})
    #' @param smpls a vector of sample names (i.e. column names of the
    #' \code{YTrgFull})
    #' @param normalization if FALSE, no normalization. If "LrnGeoMeans",
    #' normalization using the learning Geometric means. If "newGeoMeans",
    #' normalization with target Geometric means. By default it's set to "FALSE".
    #' @param expdat_meta_Lrn the list of learning set factorization metadata
    #' @param transformation if FALSE, no transformation. If TRUE, log2
    #' transformation of counts data. By default it's set to "FALSE"
    #'
    #' @returns list of prepared subset data for the current subset number
    #'
    #' @examples
    #' # see to create input data
    #'
    #' data("Trg", package = "MOTL")
    #' data("Lrn", package = "MOTL")
    #'
    #' views <- c("mRNA", "miRNA", "DNAme", "SNV")
    #' YTrgFull <- Trg$YTrg_list
    #' brcds_SS <- Trg$brcds_SS
    #' SS <- 1
    #' Fctrzn <- Lrn$Fctrzn
    #' smpls <- colnames(YTrgFull$mRNA)
    #' expdat_meta_Lrn <- Lrn$Lrn_meta
    #'
    #' YTrg_prep <- TCGATargetDataPreparation(views,
    #'                                         YTrgFull,
    #'                                         brcds_SS,
    #'                                         SS,
    #'                                         Fctrzn,
    #'                                         smpls,
    #'                                         normalization = FALSE,
    #'                                         expdat_meta_Lrn,
    #'                                         transformation = FALSE)
    #' YTrg_prep$mRNA[c(1:5), c(1:5)]
    #' YTrg_prep$DNAme[c(1:5), c(1:5)]
    #'
    #' # In the paper, several target datasets were created as subset of the
    #' # reference dataset R. This function was used to generate them
    #' # automatically.
    #' # If you are not doing the paper analysis, you can create the brcds_SS
    #' # using the following command line. Replace the "nameX" with the sample
    #' # names of your data.
    #' # You can as much as you want add dataframe on each view.
    #'
    #' brcds_SS_ex <-
    #'         list("brcds_mRNA_SS" =
    #'             list(data.frame("brcds" = c("name01", "name02"))),
    #'         "brcds_miRNA_SS" =
    #'             list(data.frame("brcds" = c("name10", "name11"))))
    #'
    #' # The SS parameter corresponds to the index of the subset you want to
    #' # prepare for a specific view. It's generated automatically if you used
    #' # the workflow describe in the github paper
    #' # \link{https://github.com/david-hirst/MOTL/blob/main/TCGAStudy/00_TCGAstudy_ReadMe.md}
    #'
    #' @export

    names(views) <- views

    ## Feature variance prefiltering and feature harmonization
    YTrgSSFull <-
        lapply(views, function(view, brcds_SS, SS, YTrgFull, Fctrzn) {
            YTrgSS <- TCGATargetDataPrefiltering(view,
                                                 brcds_SS,
                                                 SS,
                                                 YTrgFull,
                                                 Fctrzn)
            return(YTrgSS)
        }, brcds_SS, SS, YTrgFull, Fctrzn)

    ## Reshape data
    YTrgSSFull <- lapply(views, function(view, YTrgSSFull, smpls) {
        YTrgSS <- YTrgSSFull[[view]]
        if (is.element(view, c("mRNA", "miRNA", "DNAme"))) {
            YTrgSS <- SummarizedExperiment::assay(YTrgSS)
        }

        ## order columns
        colnames(YTrgSS) <- substr(colnames(YTrgSS), 1, 16)
        YTrgSS <- YTrgSS[, match(smpls, colnames(YTrgSS))]

        return(YTrgSS)
    }, YTrgSSFull, smpls)

    ## Normalization and transformation
    YTrgSSFull <- lapply(
        views,
        preprocessCountsData,
        YTrgSSFull,
        normalization,
        expdat_meta_Lrn,
        transformation
    )

    ## Return
    return(YTrgSSFull)
}

## --------------------------------------------------------------


## ------------------------------- NON TCGA RELATED FUNCTIONS ---

mRNA_addVersion <- function(expdat, Lrndat) {
    #'
    #' Format mRNA features to match with learning dataset
    #'
    #'
    #' Get mRNA ensembl ID version from learning dataset
    #' (e.g. ENSG00000122133.17) and attach to the corresponding mRNA ensembl
    #' ID in the target dataset. Feature names need to be similar between
    #' target dataset and learning dataset.
    #'
    #'
    #' @param expdat the mRNA matrix from the target dataset with genes in rows.
    #' Gene names should be in ensembl format and don't contain the version
    #' (e.g. ENSG00000122133). Rownames contain ensembl IDs and colnames sample
    #' names.
    #' @param Lrndat the mRNA W matrix from the learning dataset factorization
    #' with genes in rows. Gene names should be in ensembl format. Rownames
    #' contain ensembl IDs.
    #'
    #' @returns the target mRNA matrix with versions attached
    #'
    #' @import dplyr
    #'
    #' @examples
    #' Lrn_names <-
    #'         c("ENSG00000122133.17", "ENSG00000122194.9", "ENSG00000119411.1")
    #' Lrn_views <- c("mRNA", "mRNA", "mRNA")
    #' expdat_names <-
    #'         c("ENSG00000122133", "ENSG00000122194", "ENSG00000119411")
    #'
    #' Lrndat <- data.frame("view" = Lrn_views, row.names = Lrn_names)
    #' expdat <- data.frame("sample1" = c(1, 52, 4), row.names = expdat_names)
    #' expdat_prep <- mRNA_addVersion(expdat, Lrndat)
    #' expdat
    #' expdat_prep
    #'
    #' @export
    #'
    tmp <- as.data.frame(do.call(rbind, strsplit(rownames(Lrndat), "[.]")))
    # match to stripped ids from target set
    tmp <- data.frame(V1 = rownames(expdat)) %>%
        left_join(tmp, by = c("V1")) %>%
        as.data.frame()
    # rename target dataset features
    rownames(expdat) <- paste0(tmp$V1, ".", tmp$V2)

    # return tidied up matrix
    return(expdat)
}

TargetDataPrefiltering <- function(view, YTrg_list, Fctrzn, smpls) {
    #'
    #' Prepare the target data for a given view
    #'
    #' The function performs the following steps:
    #' 1. Remove the features with variance equal to zero
    #' 2. Harmonize features between the target data and the learning data. Only
    #' the shared features are kept.
    #' 3. Order columns according the order of samples (i.e. \code{smpls})
    #'
    #' @param view current view data name (e.g. "mRNA", or "DNAme")
    #' @param YTrg_list a named list of target data. Names correspond to the
    #' views defined and the corresponding data are saved into \code{matrix}.
    #' @param Fctrzn the learning dataset factorization model object (from
    #' \code{MOFA})
    #' @param smpls an ordered vector of sample names
    #'
    #' @returns a matrix that contains the prepared data for the current view
    #' with the sample ordered.
    #'
    #' @examples
    #'
    #' data("Lrn", package = "MOTL")
    #' data("Trg", package = "MOTL")
    #'
    #' view <- "mRNA"
    #' YTrg_list <- Trg$YTrg_prep
    #' Fctrzn <- Lrn$Fctrzn
    #' smpls <- colnames(YTrg_list$mRNA)
    #'
    #' mRNA_prep <- TargetDataPrefiltering(view, YTrg_list, Fctrzn, smpls)
    #'
    #'
    #' @export

    YTrg <- YTrg_list[[view]]

    ## prefiltering, only condition is that variance >0
    FtrsKeep <- matrixStats::rowVars(YTrg, na.rm = TRUE) > 0
    FtrsKeep[is.na(FtrsKeep)] <- FALSE
    YTrg <- YTrg[FtrsKeep, ]
    ## print(paste0("YTrg dimensions after prefiltering: ", dim(YTrg)))

    ## harmonize features between Trg and Lrn data
    features_metadata <- Fctrzn@features_metadata
    FtrsLrn <- features_metadata$feature[features_metadata$view == view]
    FtrsCommon <- FtrsLrn[is.element(FtrsLrn, rownames(YTrg))]
    FtrsKeep <- is.element(rownames(YTrg), FtrsCommon)
    YTrg <- YTrg[FtrsKeep, ]
    YTrg <- YTrg[match(FtrsCommon, rownames(YTrg)), ]

    ## order columns
    YTrg <- YTrg[, match(smpls, colnames(YTrg))]

    return(YTrg)
}

## --------------------------------------------------------------


## ------------------ LEARNING DATASET PREPARATION FUNCTIONS ----

Tau_init <- function(viewsLrn, Fctrzn, InputModel) {
    #'
    #' Initialization of the Tau values for each view
    #'
    #' Extract the Tau matrix from the MOFA object \code{Fctrzn} for each view.
    #' More explanation about Tau.
    #'
    #' @param viewsLrn the list of learning data views. For TCGA learning data
    #' it will be \code{c("mRNA", "miRNA", "DNAme", "SNV")}).
    #' @param Fctrzn the learning dataset factorization from \code{MOFA2}.
    #' @param InputModel the factorization model object of learning set
    #' \code{MOFA2}
    #'
    #' @returns a named list of Tau matrices. Names correspond to the view
    #' names.
    #'
    #' @examples
    #' \donttest{
    #' viewsLrn <- c("mRNA", "miRNA", "DNAme", "SNV")
    #' InputModel <- "Model.hdf5"
    #' Fctrzn <- load_model(file = InputModel)
    #'
    #' Tau_list <- Tau_init(viewsLrn = viewsLrn,
    #'                     Fctrzn = Fctrzn,
    #'                     InputModel = InputModel)
    #' }
    #' @export

    ## Extract Tau from the factorization of the learning set
    Tau <- rhdf5::h5read(InputModel, "expectations/Tau")
    Tau <- Tau[match(viewsLrn, names(Tau))]

    ## For each view, transfer rownames into the corresponding Tau matrix
    for (i in seq_len(length(viewsLrn))) {
        view <- viewsLrn[i]
        rownames(Tau[[view]]$group0) <-
            rownames(Fctrzn@expectations[["W"]][[view]])
    }

    # Return a named list of Tau matrix
    return(Tau)
}

TauLn_calculation <- function(view,
                              likelihoodsLrn,
                              Fctrzn,
                              LrnFctrnDir,
                              LrnSimple = TRUE) {
    #'
    #' Initialization of the log(Tau) values
    #'
    #' Two ways to initialize the log(Tau) values:
    #' 1. log transformation of the expected Tau (already init in the
    #' \code{Fctrzn}) variable
    #' 2. extract values from a .csv file that saved in \code{LrnFctrnDir}
    #' directory
    #' Tau is initialized only for gaussian data.
    #'
    #' For gaussian data, \eqn{TauLn <- log(Tau)}
    #'
    #' @param view a character of current view name data (e.g. "mRNA")
    #' @param likelihoodsLrn a named list of data types. The list can contain
    #' \code{gaussian}, \code{poisson} or \code{bernoulli} depending of the data
    #' type. Names are the view names.
    #' @param LrnSimple if TRUE, initialization uses the Tau values. If FALSE,
    #' imports values from a .csv file. By default is set to "TRUE".
    #' @param Fctrzn learning dataset factorization model object (from
    #' \code{MOFA})
    #' @param LrnFctrnDir directory where log(Tau) values are saved. Files
    #' should be named like \code{"TauLn_mRNA.csv"}.
    #'
    #' @returns the log(Tau) matrix for the current view
    #'
    #' @examples
    #'
    #' data("Lrn", package = "MOTL")
    #'
    #' Fctrzn <- Lrn$Fctrzn_init
    #' likelihoodsLrn <- get_default_model_options(Fctrzn)$likelihoods
    #'
    #' TauLn_mRNA = TauLn_calculation(view = "mRNA",
    #'                                 likelihoodsLrn = likelihoodsLrn,
    #'                                 Fctrzn = Fctrzn,
    #'                                 LrnSimple = TRUE,
    #'                                 LrnFctrnDir = LrnFctrnDir)
    #'
    #'
    #' @export

    if (likelihoodsLrn[view] == "gaussian") {
        if (LrnSimple) {
            TauLn <- log(Fctrzn@expectations[["Tau"]][[view]]$group0[, 1])
        } else {
            TauLn <- as.vector(utils::read.csv(file.path(
                LrnFctrnDir, paste0("TauLn_", view, ".csv")
            ), header = FALSE)$V1)
            names(TauLn) <- rownames(Fctrzn@expectations[["W"]][[view]])
        }
    } else {
        TauLn <- numeric()
    }
    return(TauLn)
}

WSq_calculation <- function(view, Fctrzn, LrnFctrnDir, LrnSimple = TRUE) {
    #'
    #' Initialization of the squared weight values
    #'
    #' This function load or calculate the squares weight values.
    #'
    #' The squared weight values can be load from a .csv file. This file can be
    #' created during the factorization of the learning data. See the
    #' documentation to learn how to create this file. The file name should
    #' follow this format: \code{WSq_mRNA.csv}.
    #'
    #' The squared weight valued can also be calculated using the weight values
    #' calculated during the factorization of the learning data. These values
    #' are saved in the \code{Fctrzn} variable. See the documentation.
    #'
    #'
    #' @param view a character of current view name data
    #' @param LrnSimple if TRUE, calculates the squared weight values \code{WSq}.
    #' If FALSE, imports values from a .csv file. By default is set to "TRUE".
    #' \eqn{E[W^2]} using the weight values \code{W} \eqn{E[W]}.
    #' If FALSE, load squared weight values from a file. By default, it's set
    #' to TRUE.
    #' @param Fctrzn learning dataset factorization model object (from
    #' \code{MOFA})
    #' @param LrnFctrnDir directory where \code{WSq} values are saved
    #'
    #' @returns the squared weight matrix for the current view
    #'
    #' @examples
    #'
    #' data("Lrn", package = "MOTL")
    #'
    #' Fctrzn <- Lrn$Fctrzn
    #' likelihoodsLrn <- get_default_model_options(Fctrzn)$likelihoods
    #'
    #' WSq_mRNA = WSq_calculation(view = "mRNA",
    #'                             Fctrzn = Fctrzn,
    #'                             LrnFctrnDir = LrnFctrnDir,
    #'                             LrnSimple = TRUE)
    #'
    #'
    #' @export

    if (LrnSimple) {
        WSq <- (Fctrzn@expectations[["W"]][[view]])^2
    } else {
        WSq <-
            utils::read.csv(file.path(LrnFctrnDir, paste0("WSq_", view, ".csv")),
                     header = FALSE)
        WSq <-
            as.matrix(WSq)[, seq_len(dim(Fctrzn@expectations[["W"]][[view]])[2])]
        rownames(WSq) <- rownames(Fctrzn@expectations[["W"]][[view]])
    }
    return(WSq)
}

W0_calculation <- function(view, CenterTrg, Fctrzn, LrnFctrnDir) {
    #'
    #' Initialization of feature weight intercept values
    #'
    #' This function loads or calculates the weight intercept values.
    #'
    #' The weight intercept values can be load from the
    #' \code{EstimatedIntercepts.rds} file. This file can be created using the
    #' \code{\link{intercepts_calculation}} function.
    #'
    #' The weight intercept values can also be initialized using the weight
    #' matrix. The weight matrix is set to zero.
    #'
    #' @param view a character of current view name data
    #' @param CenterTrg if FALSE, use the estimated feature weight intercept
    #' from the \code{EstimatedIntercepts.rds} file. If TRUE, don't use the
    #' estimated feature weight intercept.
    #' @param Fctrzn learning dataset factorization model object (from
    #' \code{MOFA})
    #' @param LrnFctrnDir directory where the extimated intercepts file is.
    #'
    #' @returns a feature weight intercept values matrix for the current data
    #'
    #' @examples
    #'
    #' data("Lrn", package = "MOTL")
    #'
    #' Fctrzn <- Lrn$Fctrzn
    #'
    #' W0_mRNA = W0_calculation(view = "mRNA",
    #'                             CenterTrg = TRUE,
    #'                             Fctrzn = Fctrzn,
    #'                             LrnFctrnDir = LrnFctrnDir)
    #'
    #'
    #' @export

    if (CenterTrg) {
        W0 <- Fctrzn@expectations[["W"]][[view]][, 1] * 0
    } else {
        EstInts <-
            base::readRDS(file.path(LrnFctrnDir, "EstimatedIntercepts.rds"))
        EstInts <- EstInts$Intercepts
        W0 <- EstInts[[view]]
    }
    return(W0)
}

intercepts_calculation <- function(expdat_meta,
                                   Fctrzn,
                                   FctrznDir,
                                   ExpDataDir,
                                   Seed) {
    #'
    #' Intercepts calculation
    #'
    #' Calculate feature weight intercepts for a MOFA factorization based
    #' on MLE. These can be calculated for the learning dataset factorization
    #' and then for the factorization of the target dataset with transfer
    #' learning.
    #'
    #' For Gaussian observed data, weight intercepts are the weight mean for
    #' each feature.
    #' For Poisson and Bernoulli observed data, weight intercepts are calculated
    #' using the maximum likelihood and the \code{\link{mle}} function.
    #'
    #' @param expdat_meta the named list of learning dataset factorization
    #' metadata
    #' @param Fctrzn learning dataset factorization model object (from
    #' \code{MOFA})
    #' @param FctrznDir the learning dataset factorization directory name
    #' @param ExpDataDir the learning dataset directory name
    #' @param Seed random seed number
    #'
    #' @return a file, named EstimatedIntercepts.rds and saved into
    #' \code{FctrznDir} directory.
    #'
    #' @import methods
    #'
    #' @examples
    #' #
    #' \donttest{
    #'
    #' data("Lrn", package = "MOTL")
    #'
    #' expdat_meta <- Lrn$Lrn_meta
    #' Fctrzn <- Lrn$Fctrzn
    #' FctrznDir <- "FctrznDir"
    #' ExpDataDir <- "ExpDataDir"
    #' Seed <- 1234567
    #'
    #' intercepts_calculation(expdat_meta,
    #'                         Fctrzn,
    #'                         FctrznDir,
    #'                         ExpDataDir,
    #'                         Seed)
    #' }
    #' @export
    #'

    message("Estimation of the intercept")

    fit_start_time <- Sys.time()

    ## Extract data from factorization model object
    views <- Fctrzn@data_options$views
    names(views) <- views
    likelihoods <- Fctrzn@model_options$likelihoods
    M <- Fctrzn@dimensions$M
    D <- Fctrzn@dimensions$D

    # loop through the views and estimate the intercept
    # for gaussian data its just the mean
    # for other data will try mle with a naive estimator as backup

    intercepts_list <- lapply(views, function(view,
                                              likelihoods,
                                              D,
                                              expdat_meta,
                                              Fctrzn) {
        ## print(view)

        likelihood <- likelihoods[which(names(likelihoods) == view)]
        DTmp <- D[which(names(D) == view)]

        # YTmp <- read.table(file = file.path(ExpDataDir, paste0(view,'.csv')), sep = ",")
        YTmpFileName <- file.path(ExpDataDir, paste0(view, ".csv"))
        YTmp <- as.data.frame(data.table::fread(file = YTmpFileName, sep = ","))
        YTmp <- t(as.matrix(YTmp))
        rownames(YTmp) <- expdat_meta$smpls
        colnames(YTmp) <-
            expdat_meta[[which(names(expdat_meta) == paste0("ftrs_", view))]]

        expectations <- Fctrzn@expectations
        ZWTmp <- expectations$Z$group0 %*%
            t(expectations$W[[which(names(expectations$W) == view)]])

        invisible(gc())

        if (likelihood == "gaussian") {
            InterceptsNaive <- colMeans(YTmp, na.rm = TRUE)
            Intercepts <- InterceptsNaive
            InterceptsMethod <- rep("Naive", length(Intercepts))
            names(InterceptsMethod) <- names(InterceptsNaive)
        } else if (likelihood == "poisson") {
            ## naive intercept based on approximation to feature means of ZW
            InterceptsNaive <- as.vector(log(-1 + exp(colMeans(YTmp))))

            ## mle estimate of intercept for ZW
            ## if optimiser fails for a feature will return the naive estimate

            Intercepts_df <- do.call(rbind, lapply(seq_len(DTmp), function(d, YTmp, ZWTmp) {
                ## compute for each feature vector
                YTmp_d <- YTmp[, d]
                YTmp_d_keep <- !is.na(YTmp_d)
                YTmp_d <- YTmp_d[YTmp_d_keep]

                ZWTmp_d <- ZWTmp[YTmp_d_keep, d]

                ## NLL function to optimize
                # nLL = function(interceptMLE) -sum(stats::dpois(YLrn[,d], log(1 + exp(ZWLrn[,d]+interceptMLE)), log = TRUE))
                nLL <- function(interceptMLE) {
                    -sum(log(stats::dpois(YTmp_d, log(
                        1 + exp(ZWTmp_d + interceptMLE)
                    ))[stats::dpois(YTmp_d[, d], log(1 + exp(ZWTmp_d + interceptMLE))) != 0]))
                }

                ## try to solve it and use the result otherwise use the naive estimate
                # interceptMLEfit = try(as.vector(stats4::mle(nLL, start=list(interceptMLE=0))@coef[1]))
                interceptMLEfit <- try(as.vector(stats4::mle(nLL, start = list(interceptMLE = InterceptsNaive[d]))@coef[1]))

                if (methods::is(interceptMLEfit, "try-error")) {
                    InterceptsTmp <- InterceptsNaive[d]
                    InterceptsMethodTmp <- "Naive"
                } else {
                    InterceptsTmp <- interceptMLEfit
                    InterceptsMethodTmp <- "MLE"
                }

                intercept <- data.frame(
                    "intercept" = InterceptsTmp,
                    "Method" = InterceptsMethodTmp,
                    row.names = names(InterceptsNaive)[d]
                )

                return(intercept)
            }, YTmp, ZWTmp))

            Intercepts <-
                stats::setNames(Intercepts_df$intercept, row.names(Intercepts_df))
            InterceptsMethod <-
                stats::setNames(Intercepts_df$Method, row.names(Intercepts_df))
        } else if (likelihood == "bernoulli") {
            ## naive intercept based on approximation to feature means of ZW
            InterceptsNaive <-
                log(colMeans(YTmp, na.rm = TRUE) / (1 - colMeans(YTmp, na.rm = TRUE)))

            ## mle estimate of intercept for ZW
            ## if optimiser fails for a feature will return the naive estimate

            # DTmp <- 10

            Intercepts_df <-
                do.call(rbind, lapply(seq_len(DTmp), function(d, YTmp, ZWTmp) {
                    ## compute for each feature vector

                    YTmp_d <- YTmp[, d]
                    YTmp_d_keep <- !is.na(YTmp_d)
                    YTmp_d <- YTmp_d[YTmp_d_keep]

                    ZWTmp_d <- ZWTmp[YTmp_d_keep, d]

                    ## NLL function to optimize
                    nLL <- function(InterceptMLE) {
                        -sum(log(stats::dbinom(
                            YTmp_d,
                            size = 1, stats::plogis(ZWTmp_d + InterceptMLE)
                        )[stats::dbinom(YTmp_d, size = 1, stats::plogis(ZWTmp_d + InterceptMLE)) != 0]))
                    }

                    ## try to solve it and use the result otherwise use the naive estimate

                    interceptMLEfit <-
                        try(stats4::mle(nLL, start = list(InterceptMLE = InterceptsNaive[d]))@coef[1])

                    if (methods::is(interceptMLEfit, "try-error")) {
                        InterceptsTmp <- InterceptsNaive[d]
                        InterceptsMethodTmp <- "Naive"
                    } else {
                        InterceptsTmp <- interceptMLEfit
                        InterceptsMethodTmp <- "MLE"
                    }
                    intercept <- data.frame(
                        "intercept" = InterceptsTmp,
                        "Method" = InterceptsMethodTmp,
                        row.names = names(InterceptsNaive)[d]
                    )

                    return(intercept)
                }, YTmp, ZWTmp))

            Intercepts <- stats::setNames(Intercepts_df$intercept, row.names(Intercepts_df))
            InterceptsMethod <- stats::setNames(Intercepts_df$Method, row.names(Intercepts_df))
        } else {
            InterceptsNaive <- numeric()
            Intercepts <- numeric()
            InterceptsMethod <- character()
        }

        intercepts_list <- list(
            "InterceptsNaive" = InterceptsNaive,
            "Intercepts" = Intercepts,
            "InterceptsMethod" = InterceptsMethod
        )

        return(intercepts_list)
    }, likelihoods, D, expdat_meta, Fctrzn)

    ## save the intercepts in the relevant factorization folder

    fit_end_time <- Sys.time()

    EstimatedIntercepts <- list(
        "Seed" = Seed,
        "InterceptsNaive" = lapply(views, function(v, intercepts_list) {
            return(intercepts_list[[v]][["InterceptsNaive"]])
        }, intercepts_list),
        "Intercepts" = lapply(views, function(v, intercepts_list) {
            return(intercepts_list[[v]][["Intercepts"]])
        }, intercepts_list),
        "InterceptsMethod" = lapply(views, function(v, intercepts_list) {
            return(intercepts_list[[v]][["InterceptsMethod"]])
        }, intercepts_list),
        "fit_start_time" = fit_start_time,
        "fit_end_time" = fit_end_time
    )

    saveRDS(
        EstimatedIntercepts,
        file.path(FctrznDir, "EstimatedIntercepts.rds")
    )

    message("finished")
}

## --------------------------------------------------------------
