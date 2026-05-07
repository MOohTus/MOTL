## ----------------------------- TRANSFER LEARNING FUNCTIONS ----

Zeta_calculation <- function(view, likelihoods, E_ZWSq, E_ZE_W) {
    #'
    #' Calculate the Zeta matrix for the current data view
    #'
    #' For the current data view, calculate the Zeta matrix \code{Zeta}.
    #'
    #' \describe{
    #'      \item{For bernoulli data, Zeta is calculated using the expected
    #'      values of Z matrix and W squared matrix \code{E_ZWSq}. \code{E_ZWSq}
    #'      is calculated using the \code{\link{E_ZWSq_update}}
    #'      function.}{\eqn{Zeta_{nd} = \sqrt(E[(\sum_{k} z_{n,k} w_{d,k})^2])}}
    #'      \item{For other data type, Zeta is calculated using the expected
    #'      values of Z matrix and the expected values of W matrix \code{E_ZE_W}.
    #'      \code{E_ZE_W} is calculated using the \code{\link{E_ZE_W_update}}
    #'      function.}{\eqn{Zeta_{nd} = E[\sum_{k} z_{n,k} w_{d,k}]}, so
    #'      \eqn{Zeta = ZMu \%*\% t(W)}}
    #' }
    #'
    #' @param view a character of current view name data (e.g. \code{mRNA})
    #' @param likelihoods a named list of data types. The list can contain
    #' \code{gaussian}, \code{poisson} or \code{bernoulli} depending of the data
    #' type. Names are the view names.
    #' @param E_ZWSq expected values of the multiplication of the Z matrix with
    #' weight squared W matrix.
    #' @param E_ZE_W multiplication of the expected values of Z matrix with the
    #' expected values of W matrix
    #'
    #' @returns Zeta matrix for the current data view
    #'
    #' @examples
    #'
    #' data("TL_param", package = "MOTL")
    #'
    #' view <- "mRNA"
    #' ZMuSq <- TL_param$ZMuSq
    #' ZMu <- TL_param$ZMu
    #' ZMu_0 <- TL_param$ZMu_0
    #' Fctrzn_Lrn_W0 <- TL_param$Fctrzn_Lrn_W0
    #' Fctrzn_Lrn_W <- TL_param$Fctrzn_Lrn_W
    #' Fctrzn_Lrn_WSq <- TL_param$Fctrzn_Lrn_WSq
    #' likelihoods <- c("mRNA" = "gaussian", "miRNA" = "gaussian",
    #'                 "DNAme" = "gaussian", "SNV" = "bernoulli")
    #'
    #' E_ZE_W <- list()
    #' E_Z_SqE_W_Sq <- list()
    #' E_ZSqE_WSq <- list()
    #' E_ZWSq <- list()
    #'
    #' E_ZE_W$mRNA <-
    #'     E_ZE_W_update(view, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
    #' E_Z_SqE_W_Sq$mRNA <-
    #'     E_Z_SqE_W_Sq_update(view, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
    #' E_ZSqE_WSq$mRNA <-
    #'     E_ZSqE_WSq_update(view, ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq)
    #' E_ZWSq$mRNA <-
    #'     E_ZWSq_update(view, E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq)
    #'
    #' Zeta <- Zeta_calculation(view = "mRNA",
    #'                         likelihoods = likelihoods,
    #'                         E_ZWSq = E_ZWSq,
    #'                         E_ZE_W = E_ZE_W)
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
    #' Update Tau values for the current view
    #'
    #' Tau values are updated only for bernoulli data.
    #'
    #' Tau values are updated using the following equation:
    #'
    #' \eqn{Tau = (\frac{1}{2})*(\frac{1}{Zeta[[view]]})*tanh(\frac{Zeta[[view]]}{2})}
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
    #'
    #' data("TL_param", package = "MOTL")
    #'
    #' view <- "mRNA"
    #' likelihoods <- c("mRNA" = "gaussian", "miRNA" = "gaussian",
    #'                 "DNAme" = "gaussian", "SNV" = "bernoulli")
    #' Zeta <- TL_param$Zeta
    #' Tau <- TL_param$Tau
    #'
    #'
    #' Tau <- Tau_calculation(view = view,
    #'                         likelihoods = likelihoods,
    #'                         Zeta = Zeta, Tau = Tau)
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

YGauss_calculation <- function(view,
                               likelihoods,
                               YTrg,
                               Zeta,
                               Tau,
                               CenterTrg,
                               PoisRateCstnt) {
    #'
    #' Initialize or update pseudo Y values (YGauss)
    #'
    #' For gaussian data, Y values (observed data) are centered
    #' (if `CenterTrg = TRUE`) and will not change.
    #' For non gaussian, Y values are transformed and change after each
    #' update of Z matrix. The Y pseudo values are centered at each step if
    #' `CenterTrg = TRUE`.
    #' For gaussian data this is done for each iteration `It>=0`, for others
    #' it is done for each iteration, exept the first one `It>0`.
    #'
    #' @param view a character of current view name data
    #' @param likelihoods a named list of data types. The list can contain
    #' \code{gaussian}, \code{poisson} or \code{bernoulli} depending of the data
    #' type. Names are the view names.
    #' @param YTrg current data matrix
    #' @param Zeta list of Zeta matrices
    #' @param Tau list of Tau matrices
    #' @param CenterTrg if FALSE, use the estimated feature weight intercept
    #' from the \code{EstimatedIntercepts.rds} file. If TRUE, use the feature
    #' weight means.
    #' @param PoisRateCstnt small constant added when transforming Poisson data
    #' to avoid errors
    #'
    #' @returns pseudo Y values for the current view
    #'
    #' @examples
    #'
    #' data("TL_param", package = "MOTL")
    #'
    #' view <- "mRNA"
    #' likelihoods <- c("mRNA" = "gaussian", "miRNA" = "gaussian",
    #'                 "DNAme" = "gaussian", "SNV" = "bernoulli")
    #' CenterTrg <- FALSE
    #' YTrg <- TL_param$YTrg
    #' Zeta <- TL_param$Zeta
    #' PoisRateCstnt <- 0.0001
    #'
    #' YGauss <- YGauss_calculation(view = view,
    #'                                likelihoods = likelihoods,
    #'                                YTrg, Zeta, Tau, CenterTrg, PoisRateCstnt)
    #'
    #' @export

    if (likelihoods[[view]] == "poisson") {
        YGauss <- Zeta[[view]] - stats::plogis(Zeta[[view]]) *
            (1 - YTrg[[view]] / (log(1 + exp(Zeta[[view]])) + PoisRateCstnt)) / Tau[[view]]
    } else if (likelihoods[[view]] == "bernoulli") {
        YGauss <- (2 * YTrg[[view]] - 1) / (2 * Tau[[view]])
    } else {
        YGauss <- YTrg[[view]]
    }

    if (CenterTrg) {
        YGauss <-
            sweep(YGauss, 2, as.vector(colMeans(YGauss, na.rm = TRUE)), "-")
    }

    return(YGauss)
}

ZVar_calculation <- function(view, Tau, Fctrzn_Lrn_WSq) {
    #'
    #' Calculation of the Z variances for the current data
    #'
    #' Z variances is calculation using initialized or updated Tau values
    #' and the squared weight values WSq values
    #' based on the appendix of the `MOFA2` paper
    #' and [Github code](https://github.com/bioFAM/MOFA2)
    #'
    #' @param view a character of current view name data
    #' @param Tau list of Tau matrices
    #' @param Fctrzn_Lrn_WSq Factorized learning set squared weights
    #'
    #' @returns calculated Z variances matrix for the current data
    #'
    #' @examples
    #'
    #' data("TL_param", package = "MOTL")
    #'
    #' Tau <- TL_param$Tau
    #' Fctrzn_Lrn_WSq <- TL_param$Fctrzn_Lrn_WSq
    #'
    #' ZVar <- ZVar_calculation(view = "mRNA", Tau, Fctrzn_Lrn_WSq)
    #'
    #' @export
    #'
    ZVar_m <- Tau[[view]] %*% Fctrzn_Lrn_WSq[[view]]
    return(ZVar_m)
}

ZMu_calculation <-
    function(view, k, Fctrzn_Lrn_W, Fctrzn_Lrn_W0, Tau, ZMu_0, ZMu, YGauss) {
        #'
        #' Z matrix `ZMu` calculation for the current data
        #'
        #' @param view a character of current view name data
        #' @param k feature index in the current data
        #' @param Fctrzn_Lrn_W list of factorized learning set weight matrices
        #' @param Fctrzn_Lrn_W0 list of factorized learning set weight intercept
        #' matrices
        #' @param Tau list of Tau matrices
        #' @param ZMu_0 vector of coefficients for weight intercepts
        #' @param ZMu matrix of Z values
        #' @param YGauss list of pseudo Y value matrices
        #'
        #' @returns ZMu values for the current view (Z matrix)
        #'
        #' @examples
        #'
        #' data("TL_param", package = "MOTL")
        #'
        #' k <- 10
        #' view <- "mRNA"
        #' Fctrzn_Lrn_W <- TL_param$Fctrzn_Lrn_W
        #' Fctrzn_Lrn_W0 <- TL_param$Fctrzn_Lrn_W0
        #' Tau <- TL_param$Tau
        #' ZMu_0 <- TL_param$ZMu_0
        #' ZMu <- TL_param$ZMu
        #' YGauss <- TL_param$YTrg
        #' ZMu <- TL_param$ZMu
        #'
        #' ZMu <- ZMu_calculation(view,
        #'                         k,
        #'                         Fctrzn_Lrn_W,
        #'                         Fctrzn_Lrn_W0,
        #'                         Tau,
        #'                         ZMu_0,
        #'                         ZMu,
        #'                         YGauss)
        #'
        #' @export

        ZMu_tmp1 <- matrix(Fctrzn_Lrn_W[[view]][, k],
                           nrow = dim(Tau[[view]])[1],
                           ncol = dim(Tau[[view]])[2],
                           byrow = TRUE)
        ZMu_tmp1 <- Tau[[view]] * ZMu_tmp1
        ZMu_tmp2 <- cbind(ZMu_0, ZMu[, -k]) %*%
            t(cbind(Fctrzn_Lrn_W0[[view]], Fctrzn_Lrn_W[[view]][, -k]))
        ZMu_tmp2 <- YGauss[[view]] - ZMu_tmp2
        ZMu_tmp3 <- ZMu_tmp1 * ZMu_tmp2
        ZMu_tmp3 <- rowSums(ZMu_tmp3, na.rm = TRUE)

        return(ZMu_tmp3)
    }

ELBO_calculation <-
    function(view, likelihoods, Tau, TauLn, E_ZWSq, E_ZE_W, Zeta, YTrg, YGauss,
             PoisRateCstnt) {
        #'
        #' Calculate the ELBO value for the current view/iterations
        #'
        #' for poisson and bernoulli it is the bound which is used
        #' for gaussian it is expanded gaussian log likelihood
        #'
        #' @param view a character of current view name data
        #' @param likelihoods a named list of data types. The list can contain
        #' \code{gaussian}, \code{poisson} or \code{bernoulli} depending of the
        #' data type. Names are the view names.
        #' @param Tau list of Tau matrices
        #' @param TauLn list of log(Tau) matrices
        #' @param E_ZWSq expected values of the multiplication of the Z matrix
        #' with weight squared W matrix. See \code{\link{E_ZWSq_update}}
        #' function.
        #' @param E_ZE_W multiplication of the expected values of Z matrix
        #' with the expected values of W matrix. Seed
        #' \code{\link{E_ZE_W_update}} function.
        #' @param Zeta list of Zeta matrices
        #' @param YTrg list of data
        #' @param YGauss list of pseudo Y value matrices
        #' @param PoisRateCstnt small constant added for Poisson data to
        #' avoid errors
        #'
        #' @returns the ELBO value for the current view/iteration
        #'
        #' @examples
        #'
        #' \donttest{
        #'
        #' data("TL_param", package = "MOTL")
        #'
        #' view <- "mRNA"
        #' likelihoods <- c("mRNA" = "gaussian", "miRNA" = "gaussian",
        #'                 "DNAme" = "gaussian", "SNV" = "bernoulli")
        #' Tau <- TL_param$Tau
        #' TauLn <- TL_param$TauLn
        #' Zeta <- TL_param$Zeta
        #' YTrg <- TL_param$YTrg
        #' ZMu <- TL_param$ZMu
        #' CenterTrg <- FALSE
        #' PoisRateCstnt <- 0.0001
        #' YGauss <- TL_param$YTrg
        #'
        #' YGauss <- YGauss_calculation(view = view,
        #'                                likelihoods = likelihoods,
        #'                                YTrg, Zeta, Tau, CenterTrg, PoisRateCstnt)
        #'
        #' ZMuSq <- TL_param$ZMuSq
        #' ZMu_0 <- TL_param$ZMu_0
        #' Fctrzn_Lrn_W0 <- TL_param$Fctrzn_Lrn_W0
        #' Fctrzn_Lrn_W <- TL_param$Fctrzn_Lrn_W
        #' Fctrzn_Lrn_WSq <- TL_param$Fctrzn_Lrn_WSq
        #'
        #' E_ZE_W <- list()
        #' E_Z_SqE_W_Sq <- list()
        #' E_ZSqE_WSq <- list()
        #' E_ZWSq <- list()
        #'
        #' E_ZE_W$mRNA <-
        #'     E_ZE_W_update(view, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
        #' E_Z_SqE_W_Sq$mRNA <-
        #'     E_Z_SqE_W_Sq_update(view, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
        #' E_ZSqE_WSq$mRNA <-
        #'     E_ZSqE_WSq_update(view, ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq)
        #' E_ZWSq$mRNA <-
        #'     E_ZWSq_update(view, E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq)
        #'
        #' ELBO_L <-
        #'     ELBO_calculation(view, likelihoods, Tau, TauLn, E_ZWSq, E_ZE_W,
        #'                         Zeta, YTrg, YGauss, PoisRateCstnt)
        #' }
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
            ELBO_L_tmpB <- (E_ZE_W[[view]] - Zeta[[view]]) * stats::plogis(Zeta[[view]]) * (1 - YTrg[[view]] / (log(1 + exp(Zeta[[view]])) + PoisRateCstnt))
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

            ELBO_L_tmpA <- log(stats::plogis(Zeta[[view]]))
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
    #' Calculate `E_ZE_W`
    #'
    #' `E_ZE_W` is the multiplication of the expected values of Z matrix with the
    #' expected values of W matrix \eqn{E[Z]E[W]}.
    #'
    #' @param view current view name
    #' @param ZMu_0 vector of coefficients for weight intercepts
    #' @param ZMu matrix of Z values
    #' @param Fctrzn_Lrn_W0 list of factorized learning set weight intercept
    #' matrices
    #' @param Fctrzn_Lrn_W list of factorized learning set weight matrices
    #'
    #' @returns `E_ZE_W` for current view
    #'
    #' @examples
    #'
    #' data("TL_param", package = "MOTL")
    #' view <- "mRNA"
    #' ZMu <- TL_param$ZMu
    #' ZMu_0 <- TL_param$ZMu_0
    #' Fctrzn_Lrn_W0 <- TL_param$Fctrzn_Lrn_W0
    #' Fctrzn_Lrn_W <- TL_param$Fctrzn_Lrn_W
    #'
    #' E_ZE_W <- E_ZE_W_update(view, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
    #'
    #'
    #' @export

    E_ZE_W <-
        cbind(ZMu_0, ZMu) %*% t(cbind(Fctrzn_Lrn_W0[[view]], Fctrzn_Lrn_W[[view]]))
    return(E_ZE_W)
}

E_Z_SqE_W_Sq_update <- function(view, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W) {
    #'
    #' Calculate `E_Z_SqE_W_Sq`
    #'
    #' `E_Z_SqE_W_Sq` is the multiplication of the squared expected values of Z
    #' matrix with the squared expected values of W matrix
    #' \eqn{((E[Z])^2)*((E[W])^2)}.
    #'
    #' @param view current view name
    #' @param ZMu_0 vector of coefficients for weight intercepts
    #' @param ZMu matrix of Z values
    #' @param Fctrzn_Lrn_W0 list of factorized learning set weight intercept
    #' matrices
    #' @param Fctrzn_Lrn_W list of factorized learning set weight matrices
    #'
    #' @returns `E_Z_SqE_W_Sq` for current view
    #'
    #' @examples
    #'
    #' data("TL_param", package = "MOTL")
    #' view <- "mRNA"
    #' ZMu <- TL_param$ZMu
    #' ZMu_0 <- TL_param$ZMu_0
    #' Fctrzn_Lrn_W0 <- TL_param$Fctrzn_Lrn_W0
    #' Fctrzn_Lrn_W <- TL_param$Fctrzn_Lrn_W
    #'
    #' E_Z_SqE_W_Sq <-
    #'     E_Z_SqE_W_Sq_update(view, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
    #'
    #' @export

    E_Z_SqE_W_Sq <-
        (cbind(ZMu_0, ZMu)^2) %*% t(cbind(Fctrzn_Lrn_W0[[view]], Fctrzn_Lrn_W[[view]])^2)

    return(E_Z_SqE_W_Sq)
}

E_ZSqE_WSq_update <-
    function(view, ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq) {
        #'
        #' Calculate `E_ZSqE_WSq`
        #'
        #' `E_ZSqE_WSq` is the multiplication of the expected values of the
        #' squared Z matrix with the expected values of the squared W matrix
        #' \eqn{E[Z^2]*E[W^2]}
        #'
        #' @param view current view name
        #' @param ZMu_0 vector of coefficients for weight intercepts
        #' @param ZMuSq matrix of squared Z values
        #' @param Fctrzn_Lrn_W0 list of factorized learning set weight intercept
        #' matrices
        #' @param Fctrzn_Lrn_WSq  list of factorized learning set weight squared
        #' matrices
        #'
        #' @returns E_ZSqE_WSq for current view
        #'
        #' @examples
        #'
        #' data("TL_param", package = "MOTL")
        #' view <- "mRNA"
        #' ZMuSq <- TL_param$ZMuSq
        #' ZMu_0 <- TL_param$ZMu_0
        #' Fctrzn_Lrn_W0 <- TL_param$Fctrzn_Lrn_W0
        #' Fctrzn_Lrn_WSq <- TL_param$Fctrzn_Lrn_WSq
        #'
        #' E_ZSqE_WSq <-
        #'     E_ZSqE_WSq_update(view, ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq)
        #'
        #' @export

        E_ZSqE_WSq <-
            cbind(ZMu_0^2, ZMuSq) %*% t(cbind(Fctrzn_Lrn_W0[[view]]^2, Fctrzn_Lrn_WSq[[view]]))

        return(E_ZSqE_WSq)
    }

E_ZWSq_update <- function(view, E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq) {
    #'
    #' Calculate `E_ZWSq`
    #'
    #' `E_ZWSq` is the expected values of the multiplication of the Z matrix with
    #' weight squared W matrix \eqn{E[(ZW)^2]}.
    #'
    #' @param view current view name
    #' @param E_ZE_W matrix of \eqn{E[Z]E[W]} values for current view
    #' @param ZMuSq matrix of squared ZMu values for current view
    #' @param E_Z_SqE_W_Sq matrix of \eqn{((E[Z])^2)((E[W])^2)} values for
    #' current view
    #' @param E_ZSqE_WSq matrix of \eqn{E[Z^2]E[W^2]} values for current view
    #'
    #' @returns `E_ZWSq` values for current view
    #'
    #' @examples
    #'
    #' data("TL_param", package = "MOTL")
    #' view <- "mRNA"
    #' ZMuSq <- TL_param$ZMuSq
    #' ZMu <- TL_param$ZMu
    #' ZMu_0 <- TL_param$ZMu_0
    #' Fctrzn_Lrn_W0 <- TL_param$Fctrzn_Lrn_W0
    #' Fctrzn_Lrn_W <- TL_param$Fctrzn_Lrn_W
    #' Fctrzn_Lrn_WSq <- TL_param$Fctrzn_Lrn_WSq
    #'
    #' E_ZE_W <- list()
    #' E_Z_SqE_W_Sq <- list()
    #' E_ZSqE_WSq <- list()
    #' E_ZWSq <- list()
    #'
    #'
    #' E_ZE_W$mRNA <- E_ZE_W_update(view, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
    #' E_Z_SqE_W_Sq$mRNA <-
    #'     E_Z_SqE_W_Sq_update(view, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
    #' E_ZSqE_WSq$mRNA <-
    #'     E_ZSqE_WSq_update(view, ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq)
    #'
    #' E_ZWSq$mRNA <- E_ZWSq_update(view, E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq)
    #'
    #' @export

    E_ZWSq <- (E_ZE_W[[view]]^2) - E_Z_SqE_W_Sq[[view]] + E_ZSqE_WSq[[view]]
    return(E_ZWSq)
}

VarExplFun <- function(views, YGauss, ZMu_0, Fctrzn_Lrn_W0, ZMu, Fctrzn_Lrn_W) {
    #'
    #' Calculate the variance explained by each factor for each view
    #'
    #' @param views list of view names
    #' @param YGauss list of pseudo Y value matrices
    #' @param ZMu_0 vector of coefficients for weight intercepts
    #' @param ZMu matrix of Z values
    #' @param Fctrzn_Lrn_W0 list of factorized learning set weight intercept
    #' matrices
    #' @param Fctrzn_Lrn_W list of factorized learning set weight matrices
    #'
    #' @returns variance explained matrix
    #'
    #' @examples
    #'
    #' \donttest{
    #' VarExpl <-
    #'     VarExplFun(views, YGauss, ZMu_0, Fctrzn_Lrn_W0, ZMu, Fctrzn_Lrn_W)
    #'}
    #'
    #' @export

    names(views) <- views

    SS_tmp <- lapply(views, function(view, YGauss, ZMu_0, Fctrzn_Lrn_W0) {
        SS_tmp <- sum((YGauss[[view]] - (matrix(ZMu_0, ncol = 1) %*% t(Fctrzn_Lrn_W0[[view]])))^2, na.rm = TRUE)
        return(SS_tmp)
    }, YGauss, ZMu_0, Fctrzn_Lrn_W0)

    VarExpl <- do.call(cbind, lapply(views, function(view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp) {
        factorNames <- colnames(ZMu)
        names(factorNames) <- factorNames
        var_expl_tmp <- unlist(lapply(factorNames, function(factorName, view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp) {
            RSS_tmp <- sum((YGauss[[view]] - (cbind(ZMu_0, ZMu[, factorName]) %*% t(cbind(Fctrzn_Lrn_W0[[view]], Fctrzn_Lrn_W[[view]][, factorName]))))^2, na.rm = TRUE)
            var_expl_tmp <- 1 - (RSS_tmp / SS_tmp[[view]])
            return(var_expl_tmp)
        }, view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp))
        return(var_expl_tmp)
    }, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp))

    return(VarExpl)
}

transferLearning_function <- function(TL_param, MaxIterations, MinIterations,
                                      minFactors, views, likelihoods, Fctrzn,
                                      StartDropFactor, FreqDropFactor,
                                      StartELBO, FreqELBO, DropFactorTH,
                                      ConvergenceIts, ConvergenceTH,
                                      CenterTrg, PoisRateCstnt = 0.0001,
                                      ss_start_time = NULL, outputDir = "./") {
    #'
    #' Transfer Learning with Variational Inference
    #'
    #' This function performs multi-omics matrix factorization with transfer
    #' learning. The target dataset is factorized using the latent factor values
    #' inferred from the previous factorization of a learning dataset.
    #'
    #' This function is called after target dataset is prepared (using
    #' \code{\link{TargetDataPreparation}}) and parameters initialized (using
    #' \code{\link{initTransferLearningParamaters}}).
    #'
    #' \code{TL_param} is a named list of the initialized parameters and
    #' data objects for transfer learning. It contains :
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
    #' @param TL_param a named list of initialized parameters and data objects
    #' for transfer learning. It contains target dataset, weigths and scores
    #' matrices from matrix factorization of the learning dataset calculated
    #' using MOFA. See the detail section for more informations.
    #' @param MaxIterations the maximum number of iterations for the matrix
    #' factorization convergence. After this number, the factorization is
    #' stopped.
    #' @param MinIterations the minimum number of iterations for the matrix
    #' factorization convergence. Before this number, even if the function
    #' converges, the factorization is not stopped.
    #' @param minFactors the minimum number of factors to retain
    #' @param views a named vector of the target dataset. It should contains
    #' the same names used for inside the learning dataset.
    #' @param likelihoods a named vector of the target dataset types. It can
    #' contain \code{gaussian}, \code{poisson} or \code{bernoulli} depending of
    #' the data type. Names are the view names.
    #' @param Fctrzn the learning factorization model object (from \code{MOFA}).
    #' Creating using the \code{MOFA_functions.py} python script.
    #' @param StartDropFactor number after which iteration to start dropping
    #' factors
    #' @param FreqDropFactor number that corresponds to how often to check
    #' whether to drop factors
    #' @param StartELBO number after which iteration to start checking ELBO on
    #' @param FreqELBO number that correspond to how often to assess the ELBO
    #' @param DropFactorTH threshold number to drop or not factors. If factor
    #' with lowest maximum variance explained is below this threshold, it's
    #' dropped.
    #' @param ConvergenceIts number of consecutive iterations that change in
    #' ELBO is below threshold before convergence
    #' @param ConvergenceTH threshold number for change in ELBO for checking
    #' convergence
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
    #' \donttest{
    #'
    #' data("TL_param", package = "MOTL")
    #'
    #' ss_start_time <- Sys.time()
    #' minFactors <- 13
    #' StartDropFactor <- 1
    #' FreqDropFactor <- 1
    #' StartELBO <- 1
    #' FreqELBO <- 5
    #' DropFactorTH <- 0.01
    #' MaxIterations <- 10
    #' MinIterations <- 2
    #' ConvergenceIts <- 2
    #' ConvergenceTH <- 0.0005
    #' PoisRateCstnt <- 0.0001
    #'
    #' TL_data <- transferLearning_function(TL_param, MaxIterations,
    #'                                     MinIterations, minFactors,
    #'                                     views, likelihoods, Fctrzn,
    #'                                     StartDropFactor, FreqDropFactor,
    #'                                     StartELBO, FreqELBO,
    #'                                     DropFactorTH, ConvergenceIts,
    #'                                     ConvergenceTH,
    #'                                     CenterTrg, PoisRateCstnt = 0.0001,
    #'                                     ss_start_time = NULL,
    #'                                     outputDir = "./")
    #' }
    #'
    #'
    #' @export

    ss_fit_start_time <- Sys.time()

    ## RETRIEVE PARAMETERS
    YTrgSS <- TL_param$YTrg
    Fctrzn_Lrn_W0 <- TL_param$Fctrzn_Lrn_W0
    Fctrzn_Lrn_W <- TL_param$Fctrzn_Lrn_W
    Fctrzn_Lrn_WSq <- TL_param$Fctrzn_Lrn_WSq
    Tau <- TL_param$Tau
    TauLn <- TL_param$TauLn
    smpls <- rownames(TL_param$YTrg[[1]])
    names(views) <- views

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
            message("Drop factors")

            VarExpl <- VarExplFun(views = views, YGauss = YGauss, ZMu_0 = ZMu_0, Fctrzn_Lrn_W0 = Fctrzn_Lrn_W0, ZMu = ZMu, Fctrzn_Lrn_W = Fctrzn_Lrn_W)

            # SS_tmp <- sapply(views, function(view, YGauss, ZMu_0, Fctrzn_Lrn_W0){
            #   SS_tmp <- sum((YGauss[[view]] - (matrix(ZMu_0,ncol = 1) %*% t(Fctrzn_Lrn_W0[[view]])))^2, na.rm=TRUE)
            #   return(SS_tmp)
            # }, YGauss, ZMu_0, Fctrzn_Lrn_W0, simplify = FALSE)

            # VarExpl = sapply(views, function(view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp){
            #     factorNames = colnames(ZMu)
            #     var_expl_tmp = sapply(factorNames, function(factorName, view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp){
            #       RSS_tmp = sum((YGauss[[view]] - (cbind(ZMu_0,ZMu[,factorName]) %*% t(cbind(Fctrzn_Lrn_W0[[view]],Fctrzn_Lrn_W[[view]][,factorName]))))^2, na.rm=TRUE)
            #       var_expl_tmp = 1-(RSS_tmp/SS_tmp[[view]])
            #       return(var_expl_tmp)
            #       }, view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp)
            #     return(var_expl_tmp)
            #     }, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp)

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
                E_ZE_W <- lapply(views, E_ZE_W_update, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
                E_Z_SqE_W_Sq <- lapply(views, E_Z_SqE_W_Sq_update, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
                E_ZSqE_WSq <- lapply(views, E_ZSqE_WSq_update, ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq)
                E_ZWSq <- lapply(views, E_ZWSq_update, E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq)
            }
        }

        if (It > 0) {
            message("Zeta, Tau, YGauss")

            ## Zeta values used for non-gaussian data
            Zeta <- lapply(views, Zeta_calculation, likelihoods, E_ZWSq, E_ZE_W)

            ## Update Tau values which only change for bernoulli data
            Tau <- lapply(views, Tau_calculation, likelihoods, Zeta, Tau)

            ## Initialise / update Pseudo Y values - called YGauss here
            YGauss <- lapply(views, YGauss_calculation, likelihoods, YTrgSS, Zeta, Tau, CenterTrg)
        }

        ## Z variances using initialised / updated tau values and W^2 values
        message("Zeta")
        ZVar <- Reduce("+", lapply(views, ZVar_calculation, Tau, Fctrzn_Lrn_WSq))
        ZVar <- 1 / (ZVar + 1)

        ## Initialize or update ZMu values
        if (It == 0) {
            message("Init Z")
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
            message("update Z")
            # variational updates
            for (k in seq_len(dim(ZMu)[2])) {
                ZMu_tmp <- Reduce("+", lapply(views, ZMu_calculation, k, Fctrzn_Lrn_W, Fctrzn_Lrn_W0, Tau, ZMu_0, ZMu, YGauss))
                ZMu[, k] <- ZMu_tmp * ZVar[, k] # update factor value for subsequent calculation
            }
        }

        PrintMessage <- paste0(PrintMessage, " K=", dim(ZMu)[2])

        ## Z^2 values
        message("Z squared")
        ZMuSq <- ZVar + ZMu^2

        ## Some pre calculations - results used in various parts: E_ZE_W and ZZWW
        ## E_ZE_W_nd = E[\sum_{k} z_{n,k} w_{d,k}]
        ## E_ZWSq_nd = E[(\sum_{k} z_{n,k} w_{d,k})^2] = E[\sum_k \sum_j z_{n,k} w_{d,k} z_{n,j} w_{d,j}]
        ## If A = square(ZMu%*%t(W)) - square(ZMu)%*%square(t(W))
        ## And B = ZMuSq%*%t(WSq)
        ## Then A_nd = \sum_{k} \sum_{j != k} E[z_{n,k}]E[w_{n,k}]E[z_{n,j}]E[w_{n,j}]
        ## And B_nd = \sum_{k} E[(z_{n,k})^2]E[(w_{d,k})^2]
        ## E_ZWSq_nd = (A + B)_nd = (square(ZMu%*%t(W)) - square(ZMu)%*%square(t(W)) + ZMuSq%*%t(WSq))_nd

        message("Expected")
        E_ZE_W <- lapply(views, E_ZE_W_update, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
        E_Z_SqE_W_Sq <- lapply(views, E_Z_SqE_W_Sq_update, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
        E_ZSqE_WSq <- lapply(views, E_ZSqE_WSq_update, ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq)
        E_ZWSq <- lapply(views, E_ZWSq_update, E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq)

        ## Calculate and check the ELBO

        if ((It > 0) & (It >= StartELBO) & (((It - StartELBO) %% FreqELBO) == 0)) {
            message("ELBO")

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

        message(PrintMessage)

        ## can the algorithm be stopped?
        if ((It >= 2) & (It >= MinIterations) & (convergence_token >= ConvergenceIts)) {
            message("converged")
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
