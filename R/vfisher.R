
support_spans <- function(m, n, k) {
  return(mapply(
    function(lo, hi) { lo:hi },
    lo = pmax(0L, k - n), hi = pmin(k, m),
    SIMPLIFY = FALSE
  ))
}

logdc <- function(m, n, k, support) {
  return(mapply(
    function(support, m, n, k) dhyper(support, m, n, k, log = TRUE),
    support = support, m = m, n = n, k = k, SIMPLIFY = FALSE
  ))
}

dnhyper <- function(ncp, logdc, support) {
  return(mapply(
    function(ncp, logdc, support) {
      d <- logdc + log(ncp) * support
      d <- exp(d - max(d))
      return(d/sum(d))
    },
    ncp = ncp, logdc = logdc, support = support, SIMPLIFY = FALSE
  ))
}

mnhyper <- function(ncp, lo, hi, logdc, support) {
  lims <- which(ncp == 0 | is.infinite(ncp))
  res <- integer(length(lo))
  res[lims] <- ifelse(ncp[lims] == 0, lo[lims], hi[lims])
  res[-lims] <- mapply(
    function(ncp, logdc, support) {
      sum(support * dnhyper(ncp, logdc, support))
    },
    ncp = ncp[-lims], logdc = logdc[-lims], support = support[-lims]
  )
  return(res)
}

pnhyper <- function(
  q, ncp = 1, upper.tail = FALSE, m, n, k, lo, hi, support, logdc
) {

  lims <- which(ncp == 0 | is.infinite(ncp))
  ones <- which(ncp == 1)
  ncp1 <- which(ncp == 1)
  res <- numeric(length(q))

  res[lims] <- ifelse(
    ncp == 0,
    as.numeric(if (upper.tail) q <= lo else q >= lo),
    as.numeric(if (upper.tail) q <= hi else q >= hi)
  )

  res[ones] <- phyper(
    x[ones] - upper.tail, m[ones], n[ones], k[ones],
    lower.tail = !upper.tail
  )

  res[-c(ones, lims)] <- mapply(
    function(q, ncp, logdc, support) {
      sum(dnhyper(ncp, logdc, support)[if (upper.tail) support >= q else support <= q])
    },
    q = q[-c(ones, lims)], ncp = ncp[-c(ones, lims)], logdc = logdc[-c(ones, lims)],
    support = support[-c(ones, lims)]
  )

  return(res)
}

#' @title Vectorized fisher.test
#'
#' @description
#' Provides a vectorized version of [stats::fisher.test()] for evaluating a
#' series of 2x2 contingency tables.
#'
#' @inheritParams stats::fisher.test
#'
#' @param a `x[1, 1]` in the matrix version of [stats::fisher.test()]
#' @param b `x[1, 2]` in the matrix version of [stats::fisher.test()]
#' @param c `x[2, 1]` in the matrix version of [stats::fisher.test()]
#' @param d `x[2, 2]` in the matrix version of [stats::fisher.test()]
#'
#' @details
#' N.b. `vfisher.test` does less input validation than [stats::fisher.test()].
#'
#' Because `vfisher.test` does not support anything other than 2x2 tests, the
#' arguments to `fisher.test` associated with other tests (e.g. `hybrid`) do not
#' appear.
#'
#' @return a data.frame, columns `a`, `b`, `c`, `d`, `or`, `estimate`,
#'   `p.value`. If `conf.int == TRUE` (the default), will also include column(s)
#'   for the confidence interval (two if `alternative == "two.sided"` (default)
#'   or one otherwise.). The column names are `ci.lo` and/or `ci.hi`.
#'
#' @export
vfisher.test <- function(
    a, b, c, d, conf.int = TRUE, conf.level = 0.95, or = rep(1, length(a)),
    alternative = c("two.sided", "less", "greater")
) {

  if (
    length(a) != length(b) || length(a) != length(c) || length(a) != length(d)
  ) stop("all a, b, c, d must be same length.")

  storage.mode(a) <- "integer"
  storage.mode(b) <- "integer"
  storage.mode(c) <- "integer"
  storage.mode(d) <- "integer"

  if (any(c(a, b, c, d) < 0) || anyNA(c(a, b, c, d)))
    stop("all entries of 'a', 'b', 'c', 'd' must be nonnegative and finite")

  con <- list(mult = 30)
  con[names(control)] <- control
  if ((mult <- as.integer(con$mult)) < 2)
    stop("'mult' must be integer >= 2, typically = 30")

  alternative <- match.arg(alternative)

  if (!((length(conf.level) == 1L) && is.finite(conf.level) &&
        (conf.level > 0) && (conf.level < 1)))
    stop("'conf.level' must be a single number between 0 and 1")

  if (any(is.na(or) | or < 0))
    stop("'or' must be a single number between 0 and Inf")

  # matrix =
  # a, b
  # c, d

  m <- a + c
  n <- b + d
  k <- a + b
  x <- a

  lo <- pmax(0L, k - n)
  hi <- pmin(k, m)

  mle <- numeric(length(lo))
  mle[x == lo] <- 0
  mle[x == hi] <- Inf
  nothilo <- !((x == lo) | (x == hi))

  mle[nothilo] <- {
    mu <- mnhyper(1, lo[nothilo], hi[nothilo], logdc[nothilo], support[nothilo])
    lemu <- mi < x[nothilo]
    res <- numeric(length(x[nothilo]))

    res[lemu] <- mapply(function(lo, hi, logdc) {
      1/uniroot(function(t) mnhyper(1/t, lo, hi, logdc, support) - x, c(.Machine$double.eps, 1))$root
    }, lo = lo[lemu], hi = hi[lemu], logdc = logdc[lemu])
    res[-lemu] <- mapply(function(lo, hi, logdc) {
      uniroot(function(t) mnhyper(t, lo, hi, logdc) - x, c(0, 1))$root
    }, lo = lo[-lemu], hi = hi[-lemu], logdc = logdc[-lemu])
    res
  }

  if (conf.int) {
    ncp.U <- function(x, alpha) {
      if (x == hi) return(Inf)
      p <- pnhyper(x, 1)
      if (p < alpha)
        uniroot(function(t) pnhyper(x, t) - alpha,
                c(0, 1))$root
      else if (p > alpha)
        1/uniroot(function(t) pnhyper(x, 1/t) - alpha,
                  c(.Machine$double.eps, 1))$root
      else 1
    }
    ncp.L <- function(x, alpha) {
      if (x == lo) return(0)
      p <- pnhyper(x, 1, upper.tail = TRUE)
      if (p > alpha)
        uniroot(function(t) pnhyper(x, t, upper.tail = TRUE) -
                  alpha, c(0, 1))$root
      else if (p < alpha)
        1/uniroot(function(t) pnhyper(x, 1/t, upper.tail = TRUE) -
                    alpha, c(.Machine$double.eps, 1))$root
      else 1
    }
    CINT <- switch(alternative, less = c(0, ncp.U(x,
                                                  1 - conf.level)), greater = c(ncp.L(x, 1 - conf.level),
                                                                                Inf), two.sided = {
                                                                                  alpha <- (1 - conf.level)/2
                                                                                  c(ncp.L(x, alpha), ncp.U(x, alpha))
                                                                                })

}


function (x, y = NULL, workspace = 2e+05, hybrid = FALSE, hybridPars = c(expect = 5,
                                                                         percent = 80, Emin = 1), control = list(), or = 1, alternative = "two.sided",
          conf.int = TRUE, conf.level = 0.95, simulate.p.value = FALSE,
          B = 2000)
{

  PVAL <- NULL
  else {
    PVAL <- switch(alternative, less = pnhyper(x, or),
                   greater = pnhyper(x, or, upper.tail = TRUE),
                   two.sided = {
                     if (or == 0) as.numeric(x == lo) else if (or ==
                                                               Inf) as.numeric(x == hi) else {
                                                                 relErr <- 1 + 10^(-7)
                                                                 d <- dnhyper(or)
                                                                 sum(d[d <= d[x - lo + 1] * relErr])
                                                               }
                   })
  }
  mle <- function(x) {
    if (x == lo)
      return(0)
    if (x == hi)
      return(Inf)
    mu <- mnhyper(1)
    if (mu > x)
      uniroot(function(t) mnhyper(t) - x, c(0, 1))$root
    else if (mu < x)
      1/uniroot(function(t) mnhyper(1/t) - x, c(.Machine$double.eps,
                                                1))$root
    else 1
  }
  ESTIMATE <- c(`odds ratio` = mle(x))
  if (conf.int) {
    ncp.U <- function(x, alpha) {
      if (x == hi)
        return(Inf)
      p <- pnhyper(x, 1)
      if (p < alpha)
        uniroot(function(t) pnhyper(x, t) - alpha,
                c(0, 1))$root
      else if (p > alpha)
        1/uniroot(function(t) pnhyper(x, 1/t) - alpha,
                  c(.Machine$double.eps, 1))$root
      else 1
    }
    ncp.L <- function(x, alpha) {
      if (x == lo)
        return(0)
      p <- pnhyper(x, 1, upper.tail = TRUE)
      if (p > alpha)
        uniroot(function(t) pnhyper(x, t, upper.tail = TRUE) -
                  alpha, c(0, 1))$root
      else if (p < alpha)
        1/uniroot(function(t) pnhyper(x, 1/t, upper.tail = TRUE) -
                    alpha, c(.Machine$double.eps, 1))$root
      else 1
    }
    CINT <- switch(alternative, less = c(0, ncp.U(x,
                                                  1 - conf.level)), greater = c(ncp.L(x, 1 - conf.level),
                                                                                Inf), two.sided = {
                                                                                  alpha <- (1 - conf.level)/2
                                                                                  c(ncp.L(x, alpha), ncp.U(x, alpha))
                                                                                })
    attr(CINT, "conf.level") <- conf.level
  }
  RVAL <- c(RVAL, list(conf.int = if (conf.int) CINT, estimate = ESTIMATE,
                       null.value = NVAL))
}
structure(c(RVAL, alternative = alternative, method = METHOD,
            data.name = DNAME), class = "htest")
}
