## FUNCTIONS

## ---------------------------------- TCGA RELATED FUNCTIONS ----

TCGATargetDataPrefiltering <- function(view, brcds_SS, SS, YTrgFull, Fctrzn) {
    #'
    #' Filter out TCGA target subset data according variance
    #'
    #' This function performs a prefiltering analysis through the steps:
    #' 1. Extract data for the given sample names (\code{brcds_SS})
    #' 2. Remove features with variance equal to zero
    #' 3. Match features between target data set and learning data set
    #'
    #' The function return a prefiltered target dataframe for the current view.
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
    #' @importFrom matrixStats rowVars
    #' @importFrom SummarizedExperiment assay
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
    #' brcds_SS_ex <- list("brcds_mRNA_SS" =
    #'                      list(data.frame("brcds" = c("name01", "name02"))),
    #'                    "brcds_miRNA_SS" =
    #'                      list(data.frame("brcds" = c("name10", "name11"))))
    #'
    #' # See the doc to create the input parameter
    #'
    #' expdat_mRNA <- TCGATargetDataPrefiltering(view = "mRNA",
    #' brcds_SS = brcds_SS, SS = 1, YTrgFull = YTrg_list, Fctrzn = Lrn_Fctrzn)
    #' expdat_mRNA[c(1:5), c(1:5)]
    #'
    #'
    #' @export

    print(view)

    ## select samples and subset the YTrg
    brcds <- brcds_SS[[paste0("brcds_", view, "_SS")]][[SS]]
    SmplsKeep <- is.element(colnames(YTrgFull[[view]]), brcds$brcds)
    YTrgSS <- YTrgFull[[view]][, SmplsKeep]
    print(paste0("YTrgSS dimensions: ", dim(YTrgSS)))

    ## prefiltering, only condition is that variance >0
    if (is.element(view, c("mRNA", "DNAme", "miRNA"))) {
        FtrsKeep <- rowVars(assay(YTrgSS), na.rm = TRUE) > 0
    } else {
        FtrsKeep <- rowVars(YTrgSS, na.rm = TRUE) > 0
    }
    FtrsKeep[is.na(FtrsKeep)] <- FALSE
    YTrgSS <- YTrgSS[FtrsKeep, ]
    print(paste0("YTrgSS dimensions after prefiltering: ", dim(YTrgSS)))

    ## harmonize features between Trg SS and Lrn data
    FtrsLrn <- Fctrzn@features_metadata$feature[Fctrzn@features_metadata$view == view]
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
                                      normalization = "Lrn",
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
    #' object. For the next step, data have to be stored into a \code{matrix}. SNV
    #' data are already a matrix.
    #' Then, samples are ordered in the same way between views.
    #'
    #' Finally, counts data (e.g. mRNA and miRNA) can be normalized and/or
    #' transformed using \code{\link{preprocessCountsData}} function.
    #'  - if \code{normalization = FALSE}: counts data are not normalized
    #'  - if \code{normalization = "Lrn"}: counts data are normalized using
    #'  the learning dataset geomeans calculated
    #'  - if \code{normalization = "Trg"}: counts data are normalized without the
    #'  learning dataset geomeans.
    #'
    #'  Normalization is performed in the \code{\link{countsNormalization}}
    #'  function using \code{\link{estimateSizeFactors}} from \code{\link{DESeq2}}
    #'  package. And transformation is perform using
    #'  \code{\link{countsTransformation}} with a log2 transformation.
    #'
    #'  Look \code{\link{GeoMeans_Lrn_init}} and \code{\link{GeoMeanFun}} to see
    #'  how learning GeoMeans are calculated.
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
    #' @param normalization if FALSE, no normalization. If "Lrn", normalization
    #' using the learning GeoMeans. If "Trg", normalization without learning
    #' GeoMeans. By default it's set to "Lrn".
    #' @param expdat_meta_Lrn the list of learning set factorization metadata
    #' @param transformation if FALSE, no transformation. If TRUE, log2
    #' transformation of counts data. By default it's set to TRUE
    #'
    #' @returns list of prepared subset data for the current subset number
    #'
    #' @importFrom SummarizedExperiment assay
    #'
    #' @examples
    #' # see to create input data
    #'
    #' YTrg_prep <- TCGATargetDataPreparation(views, YTrgFull, brcds_SS, SS, Fctrzn, smpls, normalization = "Lrn", expdat_meta_Lrn, transformation = TRUE)
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
    #' brcds_SS_ex <- list("brcds_mRNA_SS" =
    #'                      list(data.frame("brcds" = c("name01", "name02"))),
    #'                    "brcds_miRNA_SS" =
    #'                      list(data.frame("brcds" = c("name10", "name11"))))
    #'
    #' # The SS parameter corresponds to the index of the subset you want to
    #' # prepare for a specific view. It's generated automatically if you used the
    #' # workflow describe in the github paper
    #' # \link{https://github.com/david-hirst/MOTL/blob/main/TCGAStudy/00_TCGAstudy_ReadMe.md}
    #'
    #' @export

    print("Feature prefiltering")

    ## Feature variance prefiltering and feature harmonization
    YTrgSSFull <- sapply(views, function(view, brcds_SS, SS, YTrgFull, Fctrzn) {
        YTrgSS <- TCGATargetDataPrefiltering(view, brcds_SS, SS, YTrgFull, Fctrzn)
        return(YTrgSS)
    }, brcds_SS, SS, YTrgFull, Fctrzn)

    ## Reshape data
    YTrgSSFull <- sapply(views, function(view, YTrgSSFull, smpls) {
        YTrgSS <- YTrgSSFull[[view]]
        if (is.element(view, c("mRNA", "miRNA", "DNAme"))) {
            YTrgSS <- assay(YTrgSS)
        }

        ## order columns
        colnames(YTrgSS) <- substr(colnames(YTrgSS), 1, 16)
        YTrgSS <- YTrgSS[, match(smpls, colnames(YTrgSS))]

        return(YTrgSS)
    }, YTrgSSFull, smpls)

    ## Normalization and transformation
    print("Normalization and transformation")
    YTrgSSFull <- sapply(
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

## ----------------------------------------- INIT PARAMETERS ----

## Non-TCGA Target DATA

mRNA_addVersion <- function(expdat, Lrndat) {
    #'
    #' Format mRNA features to match with learning dataset
    #'
    #'
    #' Get mRNA ensembl ID version from learning dataset (e.g. ENSG00000122133.17)
    #' and attach to the corresponding mRNA ensembl ID in the target dataset.
    #' Feature names need to be similar between target dataset and learning
    #' dataset.
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
    #' Lrndat <- data.frame("view" = c("mRNA", "mRNA", "mRNA"), row.names = c("ENSG00000122133.17", "ENSG00000122194.19", "ENSG00000119411.11"))
    #' expdat <- data.frame("sample1" = c(1, 52, 4), row.names = c("ENSG00000122133", "ENSG00000122194", "ENSG00000119411"))
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
    #' @param YTrg_list a named list of target data. Names correspond to the views
    #' defined and the corresponding data are saved into \code{matrix}.
    #' @param Fctrzn the learning factorization model object (from \code{MOFA})
    #' @param smpls an ordered vector of sample names
    #'
    #' @returns a matrix that contains the prepared data for the current view with
    #' the sample ordered.
    #'
    #' @importFrom matrixStats rowVars
    #'
    #' @examples
    #' #
    #' TargetDataPrefiltering(view, YTrg_list, Fctrzn, smpls)
    #'
    #'
    #' @export

    YTrg <- YTrg_list[[view]]

    ## prefiltering, only condition is that variance >0
    FtrsKeep <- rowVars(YTrg, na.rm = TRUE) > 0
    FtrsKeep[is.na(FtrsKeep)] <- FALSE
    YTrg <- YTrg[FtrsKeep, ]
    print(paste0("YTrg dimensions after prefiltering: ", dim(YTrg)))

    ## harmonize features between Trg and Lrn data
    FtrsLrn <- Fctrzn@features_metadata$feature[Fctrzn@features_metadata$view == view]
    FtrsCommon <- FtrsLrn[is.element(FtrsLrn, rownames(YTrg))]
    FtrsKeep <- is.element(rownames(YTrg), FtrsCommon)
    YTrg <- YTrg[FtrsKeep, ]
    YTrg <- YTrg[match(FtrsCommon, rownames(YTrg)), ]

    ## order columns
    YTrg <- YTrg[, match(smpls, colnames(YTrg))]

    return(YTrg)
}

GeoMeans_Lrn_init <- function(view, expdat_meta_Lrn, YTrgFtrs) {
    #'
    #' Retrieve the calculated geomeans of the learning set
    #'
    #' @param view current view data name
    #' @param expdat_meta_Lrn list of learning set factorization metadata
    #' @param YTrgFtrs feature names of the current view
    #'
    #' @returns calculated geomeans of the learning set
    #'
    #' @examples
    #'
    #' GeoMeans_Lrn <- GeoMeans_Lrn_init(view = "mRNA",
    #'                                   expdat_meta_Lrn,
    #'                                   YTrgFtrs)
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
    ## ASK_DAVID MT ADD DESCRIPTION
    #' Short description
    #'
    #' Detailed description
    #'
    #' @param x vector of numeric values
    #'
    #' @return mean of non zero values from x vector
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
    #' Calculate geometric means for normalization if data = learning set
    #' Doesn't use geometric means for normalization if data = target set
    #' Use provided geometric means for normalization if transfer learning
    #'
    #' @param expdat SE object of experimental data (could be miRNA or mRNA)
    #' @param GeoMeans "Trg", "Lrn" or vector of numerics
    #'
    #' @returns list of data.frame of the counts normalized and GeoMeans calculated
    #'
    #' @import DESeq2
    #'
    #' @examples
    #'
    #' expdat_counts_norm <- countsNormalization(expdat, GeoMeans)
    #'
    #' @export

    ## create deseq object
    expdat_dds <- DESeqDataSet(expdat, design = ~1)

    if (is.numeric(GeoMeans)) {
        ## normalization
        expdat_dds_norm <- estimateSizeFactors(expdat_dds, geoMeans = GeoMeans)
        GeoMeans <- NULL
    } else if (GeoMeans == "Lrn") {
        ## calculate geometric means to use for normalization of both learning and target sets
        GeoMeans <- apply(counts(expdat_dds), 1, GeoMeanFun)
        ## normalization
        expdat_dds_norm <- estimateSizeFactors(expdat_dds, geoMeans = as.vector(GeoMeans))
    } else if (GeoMeans == "Trg") {
        ## estimate size factors
        expdat_dds_norm <- estimateSizeFactors(expdat_dds)
        GeoMeans <- NULL
    } else {
        print("GeoMeans parameter should be 'Trg' or 'Lrn' or a numeric value")
        stop()
    }

    ## Extract normalized counts
    expdat_counts_norm <- list("counts" = counts(expdat_dds_norm, normalized = TRUE))

    ## Save GeoMeans
    expdat_counts_norm$GeoMeans <- GeoMeans

    return(expdat_counts_norm)
}

countsTransformation <- function(expdat_count, TopD) {
    #' Log2 Transform and select top data based on variance
    #'
    #' @param expdat_count data.frame of the counts
    #' @param TopD number of features to keep
    #'
    #' @returns data.frame of the log2 transformed and filtered data
    #'
    #' @examples
    #' expdat_counts_fltr <- countsTransformation(expdat_count, TopD)
    #'
    #'
    #' @export

    ## log transform and filter to keep only most variable
    expdat_counts_log <- log2(expdat_count + 1)
    FtrsKeep <- base::rank(-rowVars(expdat_counts_log, na.rm = TRUE, useNames = FALSE),
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
    #' function with or without GeoMeans.
    #'
    #' Transformation is performed using the \code{\link{countsTransformation}}
    #' function with log2.
    #'
    #' @param view a data view name vector (i.e. mRNA or miRNA)
    #' @param YTrg_list a named list of target data. Names correspond to the
    #' defined views. The list contains matrices.
    #' @param normalization if FALSE, no normalization. If "Lrn", normalization
    #' using the learning GeoMeans. If "Trg", normalization without learning
    #' GeoMeans. By default, it's set to FALSE.
    #' @param expdat_meta_Lrn the list of learning set factorization metadata
    #' @param transformation if FALSE, no transformation. If TRUE, log2
    #' normalization.
    #'
    #' @returns Preprocessed counts data for the current view
    #'
    #' @importFrom SummarizedExperiment SummarizedExperiment
    #'
    #' @examples
    #' mRNA <- preprocessCountsData(view = "mRNA", YTrg_list = YTrg_list,
    #'                              normalization = "Trg",
    #'                              expdat_meta_Lrn = expdat_meta_Lrn,
    #'                              transformation = TRUE)
    #'
    #'
    #' @export

    ## Select current data
    YTrg <- YTrg_list[[view]]

    if (view %in% c("mRNA", "miRNA")) {
        print(view)

        ## Normalization
        ## if Lrn, retreive the learning set geomeans and normalize with it
        ## if Trg, normalize without geomeans
        if (is.character(normalization)) {
            if (normalization == "Lrn") {
                print("Normalize with the Learning set GeoMeans")
                GeoMeans <- GeoMeans_Lrn_init(view, expdat_meta_Lrn, rownames(YTrg))
            }
            if (normalization == "Trg") {
                print("Normalize without GeoMeans")
                GeoMeans <- normalization
            }
            YTrg <- SummarizedExperiment(assays = list(YTrg))
            YTrg <- countsNormalization(expdat = YTrg, GeoMeans = GeoMeans)
            YTrg <- YTrg$counts
        } else {
            print("No normalization")
        }

        ## Transformation
        if (transformation) {
            print("Log transform data")
            YTrg <- countsTransformation(expdat_count = YTrg, TopD = nrow(YTrg))
        } else {
            print("No transformation")
        }
    }

    ## Return
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
    #' Preparation of data consists on removing features with variance equal to
    #' zero, features harmonization between target and learning data and column
    #' ordering between views. Preparation is perform using the
    #' \code{\link{TargetDataPrefiltering}} function. See the doc for more details.
    #'
    #' It could be possible to normalize and or transform counts data using the
    #' \code{\link{preprocessCountsData}} function. Normalization can be done
    #' using GeoMeans and transformation is a log2 transformation of the counts.
    #'
    #' @param views a list of target data views (e.g. \code{c("mRNA", "miRNA")})
    #' @param YTrg_list a named list of target set data. Names correspond to the
    #' defined views. The list contains matrices.
    #' @param Fctrzn the learning factorization model object (from \code{MOFA})
    #' @param smpls a vector of sample names (i.e. column names of the
    #' \code{YTrgFull})
    #' @param normalization if FALSE, no normalization. If "Lrn", normalization
    #' using the learning GeoMeans. If "Trg", normalization without learning
    #' GeoMeans. By default it's set to FALSE.
    #' @param expdat_meta_Lrn the list of learning set factorization metadata
    #' @param transformation if FALSE, no transformation. If TRUE, log2
    #' transformation of counts data. By default it's set to FALSE.
    #'
    #' @returns a list of prepared target data
    #'
    #' @examples
    #' YTrg_prep <- TargetDataPreparation(views = c("mRNA", "miRNA", "DNAme"),
    #'                                    YTrg_list = YTrg_list,
    #'                                    Fctrzn = Fctrzn, smpls = smpls,
    #'                                    expdat_meta_Lrn = expdat_meta_Lrn,
    #'                                    normalization = FALSE,
    #'                                    transformation = FALSE)
    #'
    #'
    #' @export

    ## Feature variance prefiltering and feature harmonization
    print("Feature prefiltering")
    YTrg <- sapply(views, TargetDataPrefiltering, YTrg_list, Fctrzn, smpls)

    ## Normalization and transformation
    print("Normalization and transformation")
    YTrg <- sapply(
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
                                           expdat_meta_Lrn,
                                           Fctrzn,
                                           likelihoods) {
    #'
    #' Transfer learning parameters initialization
    #'
    #' The function performs the following steps:
    #' 1. Extract the factorized learning set weight intercepts \code{W0}
    #' 2. Extract the factorized learning set weights \code{W}
    #' 3. Extract the factorized learning set squared weights \code{Wsq}
    #' 4. Extract the learning set \code{Tau} and log(Tau) \code{TauLn}
    #' For each extracted parameter, common features between learning and target
    #' data are kept. Then target data \code{YTrg}, \code{Tau} and \code{TauLn}
    #' are transposed.
    #'
    #' Each parameter are extracted from the \code{Fctrzn} model created using
    #' \code{\link{MOFA2}}. More details.
    #'
    #' \code{YTrg} matrices should have the same columns order.
    #'
    #' Define what is each returned parameters?
    #'
    #' @param YTrg a named list of target set data. Names correspond to the
    #' defined views. The list contains \code{matrix}.
    #' @param views a vector of target data views (e.g. \code{c("mRNA", "miRNA")})
    #' @param expdat_meta_Lrn the list of learning set factorization metadata
    #' @param Fctrzn the learning factorization model object (from \code{MOFA})
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
    #' views = c("mRNA", "miRNA", "DNAme", "SNV")
    #' likelihoods = c("mRNA" = "gaussian", "miRNA" = "gaussian",
    #'                 "DNAme" = "gaussian", "SNV" = "bernoulli")
    #'
    #' TLparameter <- initTransferLearningParamaters(YTrg = YTrg,
    #'                                               views = views,
    #'                                               expdat_meta_Lrn = expdat_meta_Lrn,
    #'                                               Fctrzn = Fctrzn,
    #'                                               likelihoods = likelihoods)
    #'
    #' @export
    #'

    ## Feature names in each data
    YTrgFtrs <- lapply(YTrg, rownames)

    ## FACTORIZED LEARNING WEIGHTS MATRIX ZERO
    print("Factorized learning set weight intercepts")
    Fctrzn_Lrn_W0 <- sapply(views, function(view, Fctrzn, YTrgFtrs) {
        Fctrzn_Lrn_W0 <- Fctrzn@expectations[["W0"]][[view]]
        FtrsKeep <- is.element(names(Fctrzn_Lrn_W0), YTrgFtrs[[view]])
        Fctrzn_Lrn_W0 <- Fctrzn_Lrn_W0[FtrsKeep]
        Fctrzn_Lrn_W0 <- Fctrzn_Lrn_W0[match(YTrgFtrs[[view]], names(Fctrzn_Lrn_W0))]
        return(Fctrzn_Lrn_W0)
    }, Fctrzn, YTrgFtrs)

    ## FACTORIZED LEARNING WEIGHTS MATRIX
    print("Factorized learning set weights")
    Fctrzn_Lrn_W <- sapply(views, function(view, Fctrzn, YTrgFtrs) {
        Fctrzn_Lrn_W <- Fctrzn@expectations[["W"]][[view]]
        FtrsKeep <- is.element(rownames(Fctrzn_Lrn_W), YTrgFtrs[[view]])
        Fctrzn_Lrn_W <- Fctrzn_Lrn_W[FtrsKeep, ]
        Fctrzn_Lrn_W <- Fctrzn_Lrn_W[match(YTrgFtrs[[view]], rownames(Fctrzn_Lrn_W)), ]
        return(Fctrzn_Lrn_W)
    }, Fctrzn, YTrgFtrs)

    ## FACTORIZED LEARNING WEIGHTS MATRIX SQUARED
    print("Factorized learning set squared weights")
    Fctrzn_Lrn_WSq <- sapply(views, function(view, Fctrzn, YTrgFtrs) {
        Fctrzn_Lrn_WSq <- Fctrzn@expectations[["WSq"]][[view]]
        FtrsKeep <- is.element(rownames(Fctrzn_Lrn_WSq), YTrgFtrs[[view]])
        Fctrzn_Lrn_WSq <- Fctrzn_Lrn_WSq[FtrsKeep, ]
        Fctrzn_Lrn_WSq <- Fctrzn_Lrn_WSq[match(YTrgFtrs[[view]], rownames(Fctrzn_Lrn_WSq)), ]
        return(Fctrzn_Lrn_WSq)
    }, Fctrzn, YTrgFtrs)

    ## TAU PARAMETER
    print("Tau")
    Tau <- sapply(views, function(view, Fctrzn, YTrg, YTrgFtrs) {
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
    print("LOG Tau")
    TauLn <- sapply(views, function(view,
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

    ## transpose matrices where necessary - to make them samples x features
    YTrg <- lapply(YTrg, t)
    Tau <- lapply(Tau, t)
    TauLn <- lapply(TauLn, t)

    ## return
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

## LEARNING DATA

Tau_init <- function(viewsLrn, Fctrzn, InputModel) {
    #'
    #' Initialization of the Tau values for each view
    #'
    #' Extract the Tau matrix from the MOFA object \code{Fctrzn} for each view.
    #' More explanation about Tau.
    #'
    #' @param viewsLrn the list of learning data views. For TCGA learning data it
    #' will be \code{c("mRNA", "miRNA", "DNAme", "SNV")}).
    #' @param Fctrzn the learning factorization from \code{\link{MOFA2}}.
    #' @param InputModel the factorization model object of learning set
    #' \code{\link{MOFA2}}
    #'
    #' @returns a named list of Tau matrices. Names correspond to the view names.
    #'
    #' @importFrom rhdf5 h5read
    #'
    #' @examples
    #' viewsLrn = c("mRNA", "miRNA", "DNAme", "SNV")
    #' InputModel <- "Model.hdf5"
    #' Fctrzn <- load_model(file = InputModel)
    #'
    #' Tau_list <- Tau_init(viewsLrn = viewsLrn,
    #'                      Fctrzn = Fctrzn,
    #'                      InputModel = InputModel)
    #'
    #' @export

    ## Extract Tau from the factorization of the learning set
    Tau <- h5read(InputModel, "expectations/Tau")
    Tau <- Tau[match(viewsLrn, names(Tau))]

    ## For each view, transfer rownames into the corresponding Tau matrix
    for (i in seq_len(length(viewsLrn))) {
        view <- viewsLrn[i]
        rownames(Tau[[view]]$group0) <- rownames(Fctrzn@expectations[["W"]][[view]])
    }

    # Return a named list of Tau matrix
    return(Tau)
}

## MT - maybe change name into TauLn_init ?
## MT - view = viewsLrn ?
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
    #' imports values from a .csv file.
    #' @param Fctrzn learning factorization model object (from \code{MOFA})
    #' @param LrnFctrnDir directory where log(Tau) values are saved. Files should
    #' be named like \code{"TauLn_mRNA.csv"}.
    #'
    #' @returns the log(Tau) matrix for the current view
    #'
    #' @importFrom utils read.csv
    #'
    #' @examples
    #'
    #' likelihoods = c("mRNA" = "gaussian", "miRNA" = "gaussian",
    #'                 "DNAme" = "gaussian", "SNV" = "bernoulli")
    #'
    #' TauLn_mRNA = TauLn_calculation(view = "mRNA",
    #'                                likelihoodsLrn = likelihoods,
    #'                                Fctrzn = Fctrzn,
    #'                                LrnFctrnDir = LrnFctrnDir)
    #'
    #'
    #' @export

    if (likelihoodsLrn[view] == "gaussian") {
        if (LrnSimple) {
            TauLn <- log(Fctrzn@expectations[["Tau"]][[view]]$group0[, 1])
        } else {
            TauLn <- as.vector(read.csv(file.path(
                LrnFctrnDir, paste0("TauLn_", view, ".csv")
            ), header = FALSE)$V1)
            names(TauLn) <- rownames(Fctrzn@expectations[["W"]][[view]])
        }
    } else {
        TauLn <- numeric()
    }
    return(TauLn)
}

## MT - maybe change name into WSq_init ?
## MT - explain LrnSimple input
WSq_calculation <- function(view, Fctrzn, LrnFctrnDir, LrnSimple = TRUE) {
    #'
    #' Initialization of the squared weight values
    #'
    #' This function load or calculate the squares weight values.
    #'
    #' The squared weight values can be load from a .csv file. This file can be
    #' created during the factorization of the learning data. See the documentation
    #' to learn how to create this file. The file name should follow this format:
    #' \code{WSq_mRNA.csv}.
    #'
    #' The squared weight valued can also be calculated using the weight values
    #' calculated during the factorization of the learning data. These values are
    #' saved in the \code{Fctrzn} variable. See the documentation.
    #'
    #' "factors were ordered in the same way as for other latent variables
    #' if any factors are dropped due to being inactive, they are at the end of the dataset
    #' so can filter based on dimension of W" <-- keep ?
    #'
    #' @param view a character of current view name data
    #' @param LrnSimple if TRUE, calculates the squared weight values \code{WSq}
    #' using the weight values \code{W}. If FALSE, load squared weight values
    #' from a file. By default, it's set to TRUE.
    #' @param Fctrzn learning factorization model object (from \code{MOFA})
    #' @param LrnFctrnDir directory where \code{WSq} values are saved
    #'
    #' @returns the squared weight matrix for the current view
    #'
    #' @importFrom utils read.csv
    #'
    #' @examples
    #'
    #' WSq_mRNA = WSq_calculation(view = "mRNA",
    #'                              Fctrzn = Fctrzn,
    #'                              LrnFctrnDir = LrnFctrnDir,
    #'                              LrnSimple = TRUE)
    #'
    #'
    #' @export

    if (LrnSimple) {
        WSq <- (Fctrzn@expectations[["W"]][[view]])^2
    } else {
        WSq <- read.csv(file.path(LrnFctrnDir, paste0("WSq_", view, ".csv")), header = FALSE)
        WSq <- as.matrix(WSq)[, seq_len(dim(Fctrzn@expectations[["W"]][[view]])[2])]
        rownames(WSq) <- rownames(Fctrzn@expectations[["W"]][[view]])
    }
    return(WSq)
}

## MT - maybe change name into W0_init ?
W0_calculation <- function(view, CenterTrg, Fctrzn, LrnFctrnDir) {
    #'
    #' Initialization of the weight intercept values
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
    #' @param CenterTrg if TRUE, init the weight intercept values using the weight
    #' values. If FALSE, load the estimated weight intercept values from the
    #' \code{EstimatedIntercepts.rds} file.
    #' @param Fctrzn learning factorization model object (from \code{MOFA})
    #' @param LrnFctrnDir directory where the extimated intercepts file is.
    #'
    #' @returns a weight intercepts values matrix of the current data
    #'
    #' @examples
    #'
    #' W0_mRNA = W0_calculation(view = "mRNA",
    #'                          CenterTrg = TRUE,
    #'                          Fctrzn = Fctrzn,
    #'                          LrnFctrnDir = LrnFctrnDir)
    #'
    #'
    #' @export

    if (CenterTrg) {
        W0 <- Fctrzn@expectations[["W"]][[view]][, 1] * 0
    } else {
        EstInts <- base::readRDS(file.path(LrnFctrnDir, "EstimatedIntercepts.rds"))
        EstInts <- EstInts$Intercepts
        W0 <- EstInts[[view]]
    }
    return(W0)
}

intercepts_calculation <- function(expdat_meta,
                                   Fctrzn,
                                   FctrznDir,
                                   ExpDataDir,
                                   Seed,
                                   YTmp) {
    #'
    #' Intercepts calculation
    #'
    #' For Gaussian observed data, intercept is calculated for each feature.
    #' For Poisson and Bernoulli observed data, intercept is calculated using
    #' the maximum likelihood and the \code{\link{mle}} function.
    #'
    #' @param expdat_meta the named list of learning dataset factorization
    #' metadata
    #' @param Fctrzn learning factorization model object (from \code{MOFA})
    #' @param FctrznDir the learning dataset factorization directory name
    #' @param ExpDataDir the learning dataset directory name
    #' @param YTmp ASK_DAVID - REMOVE THIS INPUT
    #' @param Seed random seed number
    #'
    #' @return a file, named EstimatedIntercepts.rds and saved into
    #' \code{FctrznDir} directory.
    #'
    #' @importFrom data.table fread
    #' @importFrom stats dpois, dbinom, plogis, setNames
    #' @importFrom stats4 mle
    #'
    #' @examples
    #' # example code
    #' intercepts_calculation(expdat_meta,
    #'                        Fctrzn,
    #'                        FctrznDir,
    #'                        ExpDataDir,
    #'                        Seed,
    #'                        YTmp)
    #'
    #'
    #'
    #' @export
    #'

    print("Estimation of the intercept")

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
                                              YTmp,
                                              expdat_meta,
                                              Fctrzn) {
        print(view)

        likelihood <- likelihoods[which(names(likelihoods) == view)]
        DTmp <- D[which(names(D) == view)]

        # YTmp <- read.table(file = file.path(ExpDataDir, paste0(view,'.csv')), sep = ",")
        YTmp <- as.data.frame(fread(file = file.path(ExpDataDir, paste0(view, ".csv")), sep = ","))
        YTmp <- t(as.matrix(YTmp))
        rownames(YTmp) <- expdat_meta$smpls
        colnames(YTmp) <- expdat_meta[[which(names(expdat_meta) == paste0("ftrs_", view))]]

        ZWTmp <- Fctrzn@expectations$Z$group0 %*%
            t(Fctrzn@expectations$W[[which(names(Fctrzn@expectations$W) == view)]])

        # mean(colnames(YTmp)==colnames(ZWTmp))
        # mean(rownames(YTmp)==rownames(ZWTmp))

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
                    -sum(log(dpois(YTmp_d, log(
                        1 + exp(ZWTmp_d + interceptMLE)
                    ))[dpois(YTmp_d[, d], log(1 + exp(ZWTmp_d + interceptMLE))) != 0]))
                }

                ## try to solve it and use the result otherwise use the naive estimate
                # interceptMLEfit = try(as.vector(stats4::mle(nLL, start=list(interceptMLE=0))@coef[1]))
                interceptMLEfit <- try(as.vector(mle(nLL, start = list(interceptMLE = InterceptsNaive[d]))@coef[1]))

                if (is(interceptMLEfit, "try-error")) {
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

            Intercepts <- setNames(Intercepts_df$intercept, row.names(Intercepts_df))
            InterceptsMethod <- setNames(Intercepts_df$Method, row.names(Intercepts_df))
        } else if (likelihood == "bernoulli") {
            ## naive intercept based on approximation to feature means of ZW
            InterceptsNaive <- log(colMeans(YTmp, na.rm = TRUE) / (1 - colMeans(YTmp, na.rm = TRUE)))

            ## mle estimate of intercept for ZW
            ## if optimiser fails for a feature will return the naive estimate

            # DTmp <- 10

            Intercepts_df <- do.call(rbind, lapply(seq_len(DTmp), function(d, YTmp, ZWTmp) {
                ## compute for each feature vector

                YTmp_d <- YTmp[, d]
                YTmp_d_keep <- !is.na(YTmp_d)
                YTmp_d <- YTmp_d[YTmp_d_keep]

                ZWTmp_d <- ZWTmp[YTmp_d_keep, d]

                ## NLL function to optimize
                nLL <- function(InterceptMLE) {
                    -sum(log(dbinom(
                        YTmp_d,
                        size = 1, plogis(ZWTmp_d + InterceptMLE)
                    )[dbinom(YTmp_d, size = 1, plogis(ZWTmp_d + InterceptMLE)) != 0]))
                }

                ## try to solve it and use the result otherwise use the naive estimate

                interceptMLEfit <- try(mle(nLL, start = list(InterceptMLE = InterceptsNaive[d]))@coef[1])

                if (is(interceptMLEfit, "try-error")) {
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

            Intercepts <- setNames(Intercepts_df$intercept, row.names(Intercepts_df))
            InterceptsMethod <- setNames(Intercepts_df$Method, row.names(Intercepts_df))
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
    }, likelihoods, D, YTmp, expdat_meta, Fctrzn)

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

    print("finished")
}

## --------------------------------------------------------------


## ----------------------------- TRANSFER LEARNING FUNCTIONS ----

Zeta_calculation <- function(view, likelihoods, E_ZWSq, E_ZE_W) {
    #'
    #' Calculate the Zeta matrix for the current data view
    #'
    #' For the current data view, calculate the Zeta matrix \code{Zeta}.
    #'
    #' For bernoulli data, \eqn{Zeta_nd = sqrt(E[(\sum_{k} z_{n,k} w_{d,k})^2])}.
    #'
    #' For other data type, \eqn{Zeta_nd = E[\sum_{k} z_{n,k} w_{d,k}]}.
    #' So\eqn{Zeta = ZMu %*% t(W)}
    #'
    #' E_ZWSq is calculated using the \code{\link{E_ZWSq_update}} function.
    #'
    #' E_ZE_W is calculated using the \code{\link{E_ZE_W_update}} function.
    #'
    #' Zeta values used for non-gaussian data
    #' for poisson \eqn{Zeta_nd = E[\sum_{k} z_{n,k} w_{d,k}] so Zeta = ZMu %*% t(W)}
    #' for bernoulli \eqn{Zeta_nd = sqrt(E[(\sum_{k} z_{n,k} w_{d,k})^2])}
    #'
    #' @param view a character of current view name data (e.g. \code{mRNA})
    #' @param likelihoods a named list of data types. The list can contain
    #' \code{gaussian}, \code{poisson} or \code{bernoulli} depending of the data
    #' type. Names are the view names.
    #' @param E_ZWSq ASK_DAVID
    #' @param E_ZE_W ASK_DAVID
    #'
    #' @returns Zeta matrix for the current data view
    #'
    #' @examples
    #' view = "mRNA"
    #' likelihoods = c("mRNA" = "gaussian", "miRNA" = "gaussian",
    #'                 "DNAme" = "gaussian", "SNV" = "bernoulli")
    #' E_ZE_W = E_ZE_W_update(view, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
    #' E_Z_SqE_W_Sq = E_Z_SqE_W_Sq_update(view, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
    #' E_ZSqE_WSq = E_ZSqE_WSq_update(view, ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq)
    #' E_ZWSq = E_ZWSq_update(view, E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq)
    #'
    #' Zeta <- Zeta_calculation(view = "mRNA",
    #'                          likelihoods = likelihoods,
    #'                          E_ZWSq = E_ZWSq,
    #'                          E_ZE_W = E_ZE_W)
    #'
    #'
    #' @export

    if (likelihoods[[view]] == "bernoulli") {
        Zeta <- sqrt(E_ZWSq[[view]])
    } else {
        Zeta <- E_ZE_W[[view]]
    }

    return(Zeta)
}

Tau_calculation <- function(view, likelihoods, Zeta, Tau) {
    #'
    #' Update Tau values
    #'
    #' Tau values are updated using Zeta values with the
    #' following equation: \eqn{Tau = (1/2)*(1/Zeta[[view]])*tanh(Zeta[[view]]/2)}
    #' only for bernoulli data.
    #'
    #' @param view a character of current view name data
    #' @param likelihoods a named list of data types. The list can contain
    #' \code{gaussian}, \code{poisson} or \code{bernoulli} depending of the data
    #' type. Names are the view names.
    #' @param Zeta list of Zeta matrix for the current view
    #' @param Tau list of tau matrices
    #'
    #' @returns (updated) Tau matrix for the current view data
    #'
    #' @examples
    #' view = "mRNA"
    #' likelihoods = c("mRNA" = "gaussian", "miRNA" = "gaussian",
    #'                 "DNAme" = "gaussian", "SNV" = "bernoulli")
    #'
    #' Tau <- Tau_calculation(view = view,
    #'                        likelihoods = likelihoods,
    #'                        Zeta = Zeta, Tau = Tau)
    #'
    #'
    #' @export

    if (likelihoods[[view]] == "bernoulli") {
        Tau <- (1 / 2) * (1 / Zeta[[view]]) * tanh(Zeta[[view]] / 2)
    } else {
        Tau <- Tau[[view]]
    }
    return(Tau)
}

YGauss_calculation <- function(view, likelihoods, YTrg, Zeta, Tau, CenterTrg, PoisRateCstnt) {
    #'
    #' Initialize or update pseudo Y values (YGauss)
    #'
    #' For gaussian data this is just the (centered) Y values which are fixed
    #' For non gaussian these are transformed y values that change after each update of z
    #' the y pseudo values are centered at each step if the centering option is selected
    #' For gaussian data this is done for It>=0, for others it is It>0
    #'
    #' @param view a character of current view name data
    #' @param likelihoods a named list of data types. The list can contain
    #' \code{gaussian}, \code{poisson} or \code{bernoulli} depending of the data
    #' type. Names are the view names.
    #' @param YTrg current data matrix
    #' @param Zeta list of Zeta matrices
    #' @param Tau list of Tau matrices
    #' @param CenterTrg if FALSE, use (FALSE) or not (TRUE) use estimated intercepts ??
    #' @param PoisRateCstnt ASK_DAVID
    #'
    #' @returns pseudo Y values for the current view
    #'
    #' @importFrom stats plogis
    #'
    #' @examples
    #' view = "mRNA"
    #' likelihoods = c("mRNA" = "gaussian", "miRNA" = "gaussian",
    #'                 "DNAme" = "gaussian", "SNV" = "bernoulli")
    #'
    #' YGauss <- YGauss_calculation(view = view,
    #'                              likelihoods = likelihoods,
    #'                              YTrg, Zeta, Tau, CenterTrg, PoisRateCstnt)
    #'
    #'
    #' @export

    if (likelihoods[[view]] == "poisson") {
        YGauss <- Zeta[[view]] - plogis(Zeta[[view]]) * (1 - YTrg[[view]] / (log(1 + exp(Zeta[[view]])) + PoisRateCstnt)) / Tau[[view]]
    } else if (likelihoods[[view]] == "bernoulli") {
        YGauss <- (2 * YTrg[[view]] - 1) / (2 * Tau[[view]])
    } else {
        YGauss <- YTrg[[view]]
    }

    if (CenterTrg) {
        YGauss <- sweep(YGauss, 2, as.vector(colMeans(YGauss, na.rm = TRUE)), "-")
    }

    return(YGauss)
}

ZVar_calculation <- function(view, Tau, Fctrzn_Lrn_WSq) {
    #'
    #' Calculation of the Z variances for the current data
    #'
    #' Z variances using initialised / updated tau values and W^2 values
    #' based on the appendix of the mofa paper and Github code
    #'
    #' @param view a character of current view name data
    #' @param Tau list of Tau matrices
    #' @param Fctrzn_Lrn_WSq Factorized learning set squared weights
    #'
    #' @returns calculated Z variances matrix for the current data
    #'
    #' @examples
    #'
    #' ZVar <- ZVar_calculation(view = "mRNA", Tau, Fctrzn_Lrn_WSq)
    #'
    #' @export
    #'
    ZVar_m <- Tau[[view]] %*% Fctrzn_Lrn_WSq[[view]]
    return(ZVar_m)
}

ZMu_calculation <- function(view, k, Fctrzn_Lrn_W, Fctrzn_Lrn_W0, Tau, ZMu_0, ZMu, YGauss) {
    #'
    #' Z mu calculation for the current data
    #'
    #' @param view a character of current view name data
    #' @param k feature index in the current data
    #' @param Fctrzn_Lrn_W list of Factorized learning set weight matrices
    #' @param Fctrzn_Lrn_W0 list of factorized learning set weight intercept matrices
    #' @param Tau list of Tau matrices
    #' @param ZMu_0 list of ZMu intercepts matrices
    #' @param ZMu list of ZMu matrices
    #' @param YGauss list of pseudo Y value matrices
    #'
    #' @returns ZMu values for the current view
    #'
    #' @examples
    #'
    #' ZMu <- ZMu_calculation(view, k, Fctrzn_Lrn_W, Fctrzn_Lrn_W0, Tau, ZMu_0, ZMu, YGauss)
    #'
    #'
    #' @export

    ZMu_tmp1 <- matrix(Fctrzn_Lrn_W[[view]][, k], nrow = dim(Tau[[view]])[1], ncol = dim(Tau[[view]])[2], byrow = TRUE)
    ZMu_tmp1 <- Tau[[view]] * ZMu_tmp1
    ZMu_tmp2 <- cbind(ZMu_0, ZMu[, -k]) %*% t(cbind(Fctrzn_Lrn_W0[[view]], Fctrzn_Lrn_W[[view]][, -k]))
    ZMu_tmp2 <- YGauss[[view]] - ZMu_tmp2
    ZMu_tmp3 <- ZMu_tmp1 * ZMu_tmp2
    ZMu_tmp3 <- rowSums(ZMu_tmp3, na.rm = TRUE)

    return(ZMu_tmp3)
}

ELBO_calculation <- function(view, likelihoods, Tau, TauLn, E_ZWSq, E_ZE_W, Zeta, YTrg, YGauss, PoisRateCstnt) {
    #'
    #' Calculate the ELBO value for the current view/iterations
    #'
    #' likelihoods
    #' for poisson and bernoulli it is the bound which is used
    #' for gaussian it is expanded gaussian log likelihood
    #' it seems in MOFA they do not allow for the centering that is done for the
    #' pseudo data for non gaussian data
    #' they used the raw uncentered data for the elbo - this seems strange
    #' although i guess it acts as a lower bound for the lower bound ...
    #'
    #' @param view a character of current view name data
    #' @param likelihoods a named list of data types. The list can contain
    #' \code{gaussian}, \code{poisson} or \code{bernoulli} depending of the data
    #' type. Names are the view names.
    #' @param Tau list of Tau matrices
    #' @param TauLn list of log(Tau) matrices
    #' @param E_ZWSq ASK_DAVID
    #' @param E_ZE_W ASK_DAVID
    #' @param Zeta list of Zeta matrices
    #' @param YTrg list of data
    #' @param YGauss list of pseudo Y value matrices
    #' @param PoisRateCstnt ASK_DAVID
    #'
    #' @returns ASK_DAVID
    #'
    #' @importFrom stats plogis
    #'
    #' @examples
    #'
    #' view = "mRNA"
    #' likelihoods = c("mRNA" = "gaussian", "miRNA" = "gaussian",
    #'                 "DNAme" = "gaussian", "SNV" = "bernoulli")
    #'
    #' ELBO_L <- ELBO_calculation(view, likelihoods, Tau, TauLn, E_ZWSq, E_ZE_W, Zeta, YTrg, YGauss, PoisRateCstnt)
    #'
    #'
    #' @export

    if (likelihoods[[view]] == "poisson") {
        # b_nd is an upper bound for -log(p(y|x))
        # b_nd = k_nd/2 * (x_nd - zeta_nd)^2 + (x_nd - zeta_nd)f'(zeta_nd) + f(zeta_nd)
        # f'(a) = (1/(1+e^(-a)))(1 - y/log(1+e^a))
        # f(a) = log(1+e^a) - ylog(log(1+e^a))
        # The elbo component is -b_nd as this is a lower bound for log(p(y|x))
        ## A CONSTANT IS ADDED TO RATE CALCULATIONS HERE AS PER MOFA CODE TO AVOID ERRORS

        ELBO_L_tmpA <- 0.5 * Tau[[view]] * (E_ZWSq[[view]] - 2 * E_ZE_W[[view]] * Zeta[[view]] + Zeta[[view]]^2)
        ELBO_L_tmpB <- (E_ZE_W[[view]] - Zeta[[view]]) * plogis(Zeta[[view]]) * (1 - YTrg[[view]] / (log(1 + exp(Zeta[[view]])) + PoisRateCstnt))
        ELBO_L_tmpC <- (log(1 + exp(Zeta[[view]])) + PoisRateCstnt) - YTrg[[view]] * log((log(1 + exp(Zeta[[view]])) + PoisRateCstnt))
        ELBO_L_tmp <- -sum(ELBO_L_tmpA + ELBO_L_tmpB + ELBO_L_tmpC, na.rm = TRUE)
    } else if (likelihoods[[view]] == "bernoulli") {
        # Basedthe MOFA paper, MOFA code and eq(7) from jakoola paper
        # if g(a) = 1/(1 + e^(-a)) is the logistic (sigmoid) function
        # h_nd = (2 * y_nd - 1) * x_nd
        # lambda(a) = tanh(a/2)/(4*a)
        # b_nd = log(g(zeta_nd)) + (h_nd - zeta_nd)/2 - lambda(zeta_nd)(h_nd^2 - zeta_nd^2)
        # as tau_nd = 2 * lambda(zeta_nd) this becomes
        # b_nd = log(g(zeta_nd)) + (h_nd - zeta_nd)/2 - tau_nd/2 * (x_nd^2 - zeta_nd^2)
        # here b_nd is the lower bound for log(p(y|x))

        ELBO_L_tmpA <- log(plogis(Zeta[[view]]))
        ELBO_L_tmpB <- 0.5 * ((2 * YTrg[[view]] - 1) * E_ZE_W[[view]] - Zeta[[view]])
        ELBO_L_tmpC <- 0.5 * Tau[[view]] * (E_ZWSq[[view]] - Zeta[[view]]^2)
        ELBO_L_tmp <- sum(ELBO_L_tmpA + ELBO_L_tmpB - ELBO_L_tmpC, na.rm = TRUE)
    } else {
        # gaussian log likelihood
        # log(f(y_nd|x_nd,tau_nd)) = 1/2 * (log(tau_nd) - log(2*pi) - tau_nd * (y_nd - x_nd)^2)
        ELBO_L_tmpA <- TauLn[[view]] - log(2 * pi)
        ELBO_L_tmpB <- Tau[[view]] * (YGauss[[view]]^2 - 2 * YGauss[[view]] * E_ZE_W[[view]] + E_ZWSq[[view]])
        ELBO_L_tmp <- sum(0.5 * (ELBO_L_tmpA - ELBO_L_tmpB), na.rm = TRUE)
    }
    return(ELBO_L_tmp)
}

E_ZE_W_update <- function(view, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W) {
    #'
    #' Calculate
    #'
    #' @param view current view name
    #' @param ZMu_0 list of ZMu intercept matrices
    #' @param ZMu list of ZMu matrices
    #' @param Fctrzn_Lrn_W0 list of factorized learning set weight intercept matrices
    #' @param Fctrzn_Lrn_W list of factorized learning set weight matrices
    #'
    #' @returns ASK_DAVID
    #'
    #' @examples
    #'
    #' E_ZE_W <- E_ZE_W_update(view, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
    #'
    #'
    #' @export

    E_ZE_W <- cbind(ZMu_0, ZMu) %*% t(cbind(Fctrzn_Lrn_W0[[view]], Fctrzn_Lrn_W[[view]]))
    return(E_ZE_W)
}

E_Z_SqE_W_Sq_update <- function(view, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W) {
    #'
    #' Calculate
    #'
    #' @param view current view name
    #' @param ZMu_0 list of ZMu intercept matrices
    #' @param ZMu list of ZMu matrices
    #' @param Fctrzn_Lrn_W0 list of factorized learning set weight intercept matrices
    #' @param Fctrzn_Lrn_W list of factorized learning set weight matrices
    #'
    #' @returns ASK_DAVID
    #'
    #' @examples
    #'
    #' E_Z_SqE_W_Sq <- E_Z_SqE_W_Sq_update(view, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
    #'
    #'
    #' @export

    E_Z_SqE_W_Sq <- (cbind(ZMu_0, ZMu)^2) %*% t(cbind(Fctrzn_Lrn_W0[[view]], Fctrzn_Lrn_W[[view]])^2)
    return(E_Z_SqE_W_Sq)
}

E_ZSqE_WSq_update <- function(view, ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq) {
    #'
    #' Calculate
    #'
    #' @param view current view name
    #' @param ZMu_0 list of ZMu intercept matrices
    #' @param ZMuSq list of ZMu squared matrices
    #' @param Fctrzn_Lrn_W0 list of factorized learning set weight intercept matrices
    #' @param Fctrzn_Lrn_WSq  list of factorized learning set weight squared matrices
    #'
    #' @returns ASK_DAVID
    #'
    #' @examples
    #'
    #' E_ZSqE_WSq <- E_ZSqE_WSq_update(view, ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq)
    #'
    #'
    #' @export

    E_ZSqE_WSq <- cbind(ZMu_0^2, ZMuSq) %*% t(cbind(Fctrzn_Lrn_W0[[view]]^2, Fctrzn_Lrn_WSq[[view]]))
    return(E_ZSqE_WSq)
}

E_ZWSq_update <- function(view, E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq) {
    #'
    #' Calculate
    #'
    #' @param view current view name
    #' @param E_ZE_W ASK_DAVID
    #' @param ZMuSq list of ZMu squared matrices
    #' @param E_Z_SqE_W_Sq ASK_DAVID
    #' @param E_ZSqE_WSq ASK_DAVID
    #'
    #' @returns ASK_DAVID
    #'
    #' @examples
    #'
    #' E_ZWSq <- E_ZWSq_update(view, E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq)
    #'
    #'
    #' @export

    E_ZWSq <- (E_ZE_W[[view]]^2) - E_Z_SqE_W_Sq[[view]] + E_ZSqE_WSq[[view]]
    return(E_ZWSq)
}

VarExplFun <- function(views, YGauss, ZMu_0, Fctrzn_Lrn_W0, ZMu, Fctrzn_Lrn_W) {
    #'
    #' Calculate the variance explained by each factor for eah view
    #'
    #' @param views list of view data names
    #' @param YGauss list of pseudo Y value matrices
    #' @param ZMu_0 list of ZMu intercept matrices
    #' @param ZMu list of ZMu matrices
    #' @param Fctrzn_Lrn_W0 list of factorized learning set weight intercept matrices
    #' @param Fctrzn_Lrn_W list of factorized learning set weight matrices
    #'
    #' @returns variance explained matrix
    #'
    #' @examples
    #'
    #' VarExpl <- VarExplFun(views, YGauss, ZMu_0, Fctrzn_Lrn_W0, ZMu, Fctrzn_Lrn_W)
    #'
    #'
    #' @export

    SS_tmp <- sapply(views, function(view, YGauss, ZMu_0, Fctrzn_Lrn_W0) {
        SS_tmp <- sum((YGauss[[view]] - (matrix(ZMu_0, ncol = 1) %*% t(Fctrzn_Lrn_W0[[view]])))^2, na.rm = TRUE)
        return(SS_tmp)
    }, YGauss, ZMu_0, Fctrzn_Lrn_W0, simplify = FALSE)

    VarExpl <- sapply(views, function(view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp) {
        factorNames <- colnames(ZMu)
        var_expl_tmp <- sapply(factorNames, function(factorName, view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp) {
            RSS_tmp <- sum((YGauss[[view]] - (cbind(ZMu_0, ZMu[, factorName]) %*% t(cbind(Fctrzn_Lrn_W0[[view]], Fctrzn_Lrn_W[[view]][, factorName]))))^2, na.rm = TRUE)
            var_expl_tmp <- 1 - (RSS_tmp / SS_tmp[[view]])
            return(var_expl_tmp)
        }, view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp)
        return(var_expl_tmp)
    }, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp)

    return(VarExpl)
}

transferLearning_function <- function(TL_param, MaxIterations, MinIterations, minFactors,
                                      views, likelihoods, Fctrzn,
                                      StartDropFactor, FreqDropFactor, StartELBO, FreqELBO, DropFactorTH, ConvergenceIts, ConvergenceTH,
                                      CenterTrg, PoisRateCstnt = 0.0001, ss_start_time = NULL, outputDir = "./") {
    #'
    #' Transfer Learning with Variational Inference
    #'
    #' This function performs multi-omics matrix factorization with transfer
    #' learning. The target dataset is factorized using the latent factor values
    #' inferred from the previous factorization of a learning dataset.
    #'
    #'
    #' This function is called after target dataset is prepared (using
    #' \code{\link{TargetDataPreparation}}) and parameters initialized (using
    #' \code{\link{initTransferLearningParamaters}}).
    #'
    #'
    #' \code{TL_param} is a named list of the initialized parameters for
    #' transfer learning. It contains :
    #' 1. \code{YTrg}: a named list of matrices. Each matrix corresponds to
    #' the target dataset.
    #' 2. \code{Fctrzn_Lrn_W0}: a named list of vectors. Each vector contains
    #' the features mean weight matrix calculated for the learning dataset using
    #' MOFA.
    #' 3. \code{Fctrzn_Lrn_W}: a named list of matrices. Each matrix contains
    #' the weights matrix calculated for the learning dataset using MOFA.
    #' 4. \code{Fctrzn_Lrn_WSq}: a named list of matrices. Each matrix contains
    #' the squared weights matrix calculated for the learning dataset using
    #' MOFA.
    #' 5. \code{Tau}: a named list of matrices. Each matrix contains the Tau
    #' values matrix calculated for the learning dataset using MOFA.
    #' 6. \code{TauLn}: a named list of matrices. Each matrix contains the
    #' TauLn values matrix calculated for the learning dataset using MOFA.
    #'
    #' Names of each list should be identical (e.g.
    #' \code{c("mRNA", "miRNA", "DNAme", "SNV")}) and so each element
    #' corresponds to each omic data.
    #'
    #' To create the \code{TL_param} variable, see the
    #' \code{\link{initTransferLearningParamaters}} function.
    #'
    #' @param TL_param a named list of initialized parameters for transfer
    #' learning. It contains target dataset, weigths and scores matrices from
    #' matrix factorization of the learning dataset calculated using MOFA. See
    #' the detail section for more informations.
    #' @param MaxIterations the maximum number of iterations for the matrix
    #' factorization convergence. After this number, the factorization is
    #' stopped.
    #' @param MinIterations the minimum number of iteration for the matrix
    #' factorization convergence. Before this number, even if the function
    #' converges, the factorization is not stopped.
    #' @param minFactors the minimum number of factors we expected
    #' @param views a named vector of the target dataset. It should contains
    #' the same names used for inside the learning dataset.
    #' @param likelihoods a named vector of the target dataset types. It can
    #' contain \code{gaussian}, \code{poisson} or \code{bernoulli} depending of
    #' the data type. Names are the view names.
    #' @param Fctrzn the learning factorization model object (from \code{MOFA}).
    #' Creating using the \code{MOFA_functions.py} python script.
    #' @param StartDropFactor number after which iteration to start dropping
    #' factors
    #' @param FreqDropFactor number that corresponds to how often to drop
    #' factors
    #' @param StartELBO number after which iteration to start checking ELBO on
    #' @param FreqELBO number that correspond to how often to assess the ELBO
    #' @param DropFactorTH threshold number to drop or not factors. If factor
    #' with lowest maximum variance is below this threshold, it's dropped.
    #' @param ConvergenceIts number of consecutive checks in a row for which the
    #' change in ELBO is be
    #' @param ConvergenceTH threshold number of change in ELBO
    #' @param CenterTrg if TRUE, center the target dataset during processing, if
    #' FALSE, leave target dataset uncentered and use estimated learning dataset
    #' intercepts.
    #' @param PoisRateCstnt amount number added to the Poisson rate function
    #' to avoid error. By default is equal to 0.0001. It's used in the pseudo
    #' gaussian values calculation \code{YGauss}, see
    #' \code{\link{YGauss_calculation}} and ELBO calculation, see
    #' \code{\link{ELBO_calculation}}.
    #' @param ss_start_time time recorded before the preprocessing step starts.
    #' Generated using \code{\link{Sys.time}} function. By default is NULL.
    #' @param outputDir output directory name where to save results. By default
    #' results are saved in the current directory.
    #'
    #' @returns a named list of results. It contains
    #' 1. \code{YTrgSS} list of matrices of target dataset
    #' 2. \code{YGauss} list of matrices of pseudo gaussian target dataset
    #' 3. \code{ZMu_0} list of ZMu intercepts matrices
    #' 4. \code{ZMu} list of ZMu
    #' 5. \code{Fctrzn_Lrn_W0} list of learning features mean weight matrix
    #' 6. \code{Fctrzn_Lrn_W} list of learning weights matrix
    #' 7. \code{ELBO} numeric value of ELBO
    #' 8. \code{VarExpl} variance explained by each target dataset
    #' 9. \code{ss_start_time} time when start the analysis (i.e. before the
    #' preprocessing step)
    #' 10. \code{ss_fit_start_time} time when start the transfer learning
    #' analysis
    #' 11. \code{ss_end_time} time when finish the transfer learning.
    #'
    #' Results are also saved into \code{TL_data.rds} file.
    #'
    #' @examples
    #'
    #' TL_data <- transferLearning_function(TL_param, MaxIterations, MinIterations, minFactors,
    #'                                      views, likelihoods, Fctrzn,
    #'                                      StartDropFactor, FreqDropFactor, StartELBO,
    #'                                      FreqELBO, DropFactorTH, ConvergenceIts, ConvergenceTH,
    #'                                      CenterTrg, PoisRateCstnt = 0.0001,
    #'                                      ss_start_time = NULL, outputDir = "./")
    #'
    #'
    #' @export

    ss_fit_start_time <- Sys.time()

    ## RETREIVE PARAMETERS
    YTrgSS <- TL_param$YTrg
    Fctrzn_Lrn_W0 <- TL_param$Fctrzn_Lrn_W0
    Fctrzn_Lrn_W <- TL_param$Fctrzn_Lrn_W
    Fctrzn_Lrn_WSq <- TL_param$Fctrzn_Lrn_WSq
    Tau <- TL_param$Tau
    TauLn <- TL_param$TauLn
    smpls <- rownames(TL_param$YTrg[[1]])

    ## INIT PARAMETERS
    ELBO <- numeric()
    convergence_token <- 0

    ## FOR EACH ITERATION
    for (It in 0:MaxIterations) {
        PrintMessage <- paste0("It=", It)

        ## Drop factors explaining variance below threshold using MOFA formula:
        # SS = np.square(Y[m][gg, :]).sum()
        # Res = np.sum((Y[m][gg, :] - Ypred) ** 2.0)
        # r2[g][m, k] = 1.0 - Res / SS as per paper
        BegK <- dim(Fctrzn_Lrn_W[[1]])[2]
        if ((BegK > minFactors) & (It > 1) & (It > StartDropFactor) & (((It - StartDropFactor - 1) %% FreqDropFactor) == 0)) {
            print("Drop factors")

            VarExpl <- VarExplFun(views = views, YGauss = YGauss, ZMu_0 = ZMu_0, Fctrzn_Lrn_W0 = Fctrzn_Lrn_W0, ZMu = ZMu, Fctrzn_Lrn_W = Fctrzn_Lrn_W)

            # SS_tmp <- sapply(views, function(view, YGauss, ZMu_0, Fctrzn_Lrn_W0){
            #   SS_tmp <- sum((YGauss[[view]] - (matrix(ZMu_0,ncol = 1) %*% t(Fctrzn_Lrn_W0[[view]])))^2, na.rm=TRUE)
            #   return(SS_tmp)
            # }, YGauss, ZMu_0, Fctrzn_Lrn_W0, simplify = FALSE)

            # VarExpl = sapply(views, function(view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp){
            #   factorNames = colnames(ZMu)
            #   var_expl_tmp = sapply(factorNames, function(factorName, view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp){
            #     RSS_tmp = sum((YGauss[[view]] - (cbind(ZMu_0,ZMu[,factorName]) %*% t(cbind(Fctrzn_Lrn_W0[[view]],Fctrzn_Lrn_W[[view]][,factorName]))))^2, na.rm=TRUE)
            #     var_expl_tmp = 1-(RSS_tmp/SS_tmp[[view]])
            #     return(var_expl_tmp)
            #     }, view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp)
            #   return(var_expl_tmp)
            #   }, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp)

            var_expl_max <- apply(VarExpl, 1, max)

            ## drop factor with lowest max variance explained if below the threshold
            if (min(var_expl_max) < DropFactorTH) {
                fctrs_to_keep <- !base::rank(var_expl_max, ties.method = "first") == 1
                ZMu <- ZMu[, fctrs_to_keep]
                ZMuSq <- ZMuSq[, fctrs_to_keep]
                for (i in seq_len(length(views))) {
                    Fctrzn_Lrn_W[[i]] <- Fctrzn_Lrn_W[[i]][, fctrs_to_keep]
                    Fctrzn_Lrn_WSq[[i]] <- Fctrzn_Lrn_WSq[[i]][, fctrs_to_keep]
                }

                ## recalculate expectations based on new number of factors
                ## explanation of these formulae are further down the script
                E_ZE_W <- sapply(views, E_ZE_W_update, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
                E_Z_SqE_W_Sq <- sapply(views, E_Z_SqE_W_Sq_update, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
                E_ZSqE_WSq <- sapply(views, E_ZSqE_WSq_update, ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq)
                E_ZWSq <- sapply(views, E_ZWSq_update, E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq)
            }
        }

        if (It > 0) {
            print("Zeta, Tau and YGauss")

            ## Zeta values used for non-gaussian data
            Zeta <- sapply(views, Zeta_calculation, likelihoods, E_ZWSq, E_ZE_W)

            ## Update Tau values which only change for bernoulli data
            Tau <- sapply(views, Tau_calculation, likelihoods, Zeta, Tau)

            ## Initialise / update Pseudo Y values - called YGauss here
            YGauss <- sapply(views, YGauss_calculation, likelihoods, YTrgSS, Zeta, Tau, CenterTrg)
        }

        ## Z variances using initialised / updated tau values and W^2 values
        print("Zeta")
        ZVar <- Reduce("+", lapply(views, ZVar_calculation, Tau, Fctrzn_Lrn_WSq))
        ZVar <- 1 / (ZVar + 1)

        ## Initialize or update ZMu values
        if (It == 0) {
            print("Init Z")
            # initialise with means of learning set Z
            ZMu <- matrix(
                data = as.vector(colMeans(Fctrzn@expectations$Z$group0)),
                nrow = dim(ZVar)[1], ncol = dim(ZVar)[2],
                byrow = TRUE
            )
            rownames(ZMu) <- smpls
            colnames(ZMu) <- colnames(Fctrzn_Lrn_W[[1]])

            # vector of 1s to act as multiplier of W intercept term
            ZMu_0 <- rep(1, dim(ZVar)[1])
        } else {
            print("update Z")
            # variational updates
            for (k in seq_len(dim(ZMu)[2])) {
                ZMu_tmp <- Reduce("+", lapply(views, ZMu_calculation, k, Fctrzn_Lrn_W, Fctrzn_Lrn_W0, Tau, ZMu_0, ZMu, YGauss))
                ZMu[, k] <- ZMu_tmp * ZVar[, k] # update factor value for subsequent calculation
            }
        }

        PrintMessage <- paste0(PrintMessage, " K=", dim(ZMu)[2])

        ## Z^2 values
        print("Z squared")
        ZMuSq <- ZVar + ZMu^2

        ## Some pre calculations - results used in various parts: E_ZE_W and ZZWW
        ## E_ZE_W_nd = E[\sum_{k} z_{n,k} w_{d,k}]
        ## E_ZWSq_nd = E[(\sum_{k} z_{n,k} w_{d,k})^2] = E[\sum_k \sum_j z_{n,k} w_{d,k} z_{n,j} w_{d,j}]
        ## If A = square(ZMu%*%t(W)) - square(ZMu)%*%square(t(W))
        ## And B = ZMuSq%*%t(WSq)
        ## Then A_nd = \sum_{k} \sum_{j != k} E[z_{n,k}]E[w_{n,k}]E[z_{n,j}]E[w_{n,j}]
        ## And B_nd = \sum_{k} E[(z_{n,k})^2]E[(w_{d,k})^2]
        ## E_ZWSq_nd = (A + B)_nd = (square(ZMu%*%t(W)) - square(ZMu)%*%square(t(W)) + ZMuSq%*%t(WSq))_nd

        print("Expected")
        E_ZE_W <- sapply(views, E_ZE_W_update, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
        E_Z_SqE_W_Sq <- sapply(views, E_Z_SqE_W_Sq_update, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
        E_ZSqE_WSq <- sapply(views, E_ZSqE_WSq_update, ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq)
        E_ZWSq <- sapply(views, E_ZWSq_update, E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq)

        ## Calculate and check the ELBO

        if ((It > 0) & (It >= StartELBO) & (((It - StartELBO) %% FreqELBO) == 0)) {
            print("ELBO")

            ELBO_L <- do.call(sum, lapply(X = views, FUN = ELBO_calculation, likelihoods, Tau, TauLn, E_ZWSq, E_ZE_W, Zeta, YTrgSS, YGauss))

            ## The ELBO components for the variational and prior distributions for Z
            ELBO_P <- sum(0.5 * (-log(2 * pi) - ZMuSq))
            ELBO_Q <- sum(0.5 * (-1 - log(2 * pi) - log(ZVar)))
            ELBO <- c(ELBO, ELBO_L + ELBO_P - ELBO_Q)

            PrintMessage <- paste0(PrintMessage, " ELBO=", round(ELBO[length(ELBO)], 2))

            # Calculate and check the change in ELBO
            # I originally didn't allow negative changes in ELBO before convergence but MOFA does so I now allow it
            if (length(ELBO) >= 2) {
                ELBO_delta <- 100 * abs((ELBO[length(ELBO)] - ELBO[length(ELBO) - 1]) / ELBO[1])
                if ((ELBO_delta < ConvergenceTH)) {
                    convergence_token <- convergence_token + 1
                } else {
                    convergence_token <- 0
                }
                PrintMessage <- paste0(PrintMessage, " ELBO_delta=", round(ELBO_delta, 8), "%")
            }
        }

        print(PrintMessage)

        ## can the algorithm be stopped?
        if ((It >= 2) & (It >= MinIterations) & (convergence_token >= ConvergenceIts)) {
            print("converged")
            break
        }
    }

    ## Variance explained calculation with final factors
    VarExpl <- VarExplFun(views = views, YGauss = YGauss, ZMu_0 = ZMu_0, Fctrzn_Lrn_W0 = Fctrzn_Lrn_W0, ZMu = ZMu, Fctrzn_Lrn_W = Fctrzn_Lrn_W)

    # export the data for further analysis
    ss_end_time <- Sys.time()
    TL_data <- list(
        "YTrgSS" = YTrgSS,
        "YGauss" = YGauss,
        "ZMu_0" = ZMu_0,
        "ZMu" = ZMu,
        "Fctrzn_Lrn_W0" = Fctrzn_Lrn_W0,
        "Fctrzn_Lrn_W" = Fctrzn_Lrn_W,
        "ELBO" = ELBO,
        "VarExpl" = VarExpl,
        "ss_start_time" = ss_start_time,
        "ss_fit_start_time" = ss_fit_start_time,
        "ss_end_time" = ss_end_time
    )
    saveRDS(TL_data, file.path(outputDir, "TL_data.rds"))

    ## delete before next subset in case i've missed something in the loop
    rm(list = c(
        "YTrgSS", "YGauss", "ZMu_0", "ZMu", "Fctrzn_Lrn_W0", "Fctrzn_Lrn_W", "ELBO", "ELBO_P", "ELBO_Q", "ELBO_L",
        "E_ZE_W", "E_Z_SqE_W_Sq", "E_ZSqE_WSq", "E_ZWSq", "ZMuSq", "Fctrzn_Lrn_WSq", "Tau", "TauLn", "ZVar", "Zeta"
    ))

    invisible(gc())

    return(TL_data)
}


## --------------------------------------------------------------
