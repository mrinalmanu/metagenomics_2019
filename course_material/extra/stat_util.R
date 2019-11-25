library(philr)
library(robustlmm)
library(plyr)
library(docstring)  # enables Python-like docstring documentation
# library(functional)

# add a signed log1p transform to make balance boxplots more readable
signed.log1p.transform <- scales::trans_new(
    'signed.log',
    transform=function(x) sign(x) * log1p(abs(x)),
    inverse=function(x) sign(x) * expm1(abs(x))
)


select.features <- function(glmnet.coeffs) {
    #' Select features based on coefficients from a logistic glmnet model

    # selecting featurs
    selected <- as.matrix(glmnet.coeffs)
    selected <- rownames(selected)[which(selected != 0)]
    # removing the intercept
    selected[2:length(selected)]
}

select.features.multinomial <- function(glmnet.coeffs) {
    #' Select features based on coefficients from a multinomial glmnet model

    featues <- lapply(glmnet.coeffs, select.features)
    unique(unlist(featues))
}

fit.mixed.models <- function(data, responses, predictors, ...) {
    #' Fit a series of mixed models
    #' @param predictors a vector of strings representing the right-hand side of the model formula;
    #' when length(predictors) == 1, the same predictor design is used for all response variables

    # formulae <- sapply(paste(responses, paste(factors, collapse=" + "), sep=' ~ '), as.formula)
    formulae <- sapply(paste(responses, predictors, sep=' ~ '), as.formula)
    sapply(formulae, function(formula_) rlmer(formula_, data=data, ...))
    # sapply(formulae, function(formula_) lmer(formula_, data=data, ...))
}

validate.residuals <- function(models, grouping.factor, alpha) {
    #' Test residuals for heteroskedasticity and normality
    #' @param models a list of fitted models
    #' @param grouping.factor a single vector of grouping factor values
    #' @param alpha a null-hypothesis rejection threshold

    # heteroskedasticity
    equation <- paste('residuals')
    levene <- lapply(models, residuals) %>%
        lapply(function(x) leveneTest(x, grouping.factor)$`Pr(>F)`[[1]]) %>%
        unlist
    # normality
    shapiro <- lapply(models, residuals) %>%
        lapply(shapiro.test) %>%
        lapply(function(x) x$p.value) %>%
        unlist

    (levene > alpha) & (shapiro > alpha)
}

qqplot.lmer <- function(model) {
    qqnorm(model %>% resid)
    qqline(model %>% resid)
}

extract.balance.names <- function(tree, taxonomy, balances, clade, threshold) {
    #' Extract names given a phylogenetic tree, a taxonomy table, a vector of balance names
    #' and a clade type (up or down)

    sapply(balances, function(x) name.balance(tree, taxonomy, x, thresh=threshold, return.votes=clade))
}

pack.balance.names <- function(up.names, down.names, ranks) {
    #' Pack names extracted by extract.balance.names into a single long df

    # match lengths
    stopifnot(length(up.names) == length(down.names))
    # match features
    stopifnot(dimnames(up.names)[[2]] == dimnames(down.names)[[2]])

    balance.dimnames <- dimnames(up.names)[[2]]

    flatten.rank.votes <- function(balance, rank) {
        votes <- as.data.frame(balance[[rank]])
        votes$rank <- rank
        votes
    }

    flatten.balance.votes <- function(balance) {
        ldply(lapply(ranks, function(rank) flatten.rank.votes(balance, rank)), data.frame)
    }

    bind.directions <- function(dimension) {
        balance.name <- up.names[[dimension*2 - 1]]
        dim.name <- balance.dimnames[dimension]
        up.votes <- flatten.balance.votes(up.names[[dimension*2]])
        up.votes$clade <- '+'
        down.votes <- flatten.balance.votes(down.names[[dimension*2]])
        down.votes$clade <- '-'
        votes <- rbind(up.votes, down.votes)
        votes$feature <- dim.name
        votes$balance <- balance.name
        names(votes)[names(votes) == 'Freq'] <- 'nleaves'
        votes
    }
    name.positions <- seq(1, length(balance.dimnames))
    ldply(lapply(name.positions, bind.directions), data.frame)
 }

