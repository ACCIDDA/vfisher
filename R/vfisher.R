
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
  if (ncp == 0) {
    return(lo)
  } else if (ncp == Inf) {
    return(hi)
  } else {
    return(sum(support * dnhyper(ncp, logdc, support)))
  }
}

pnhyper <- function(
  q, ncp = 1, upper.tail = FALSE, m, n, k
) {
  if (ncp == 1) {
    phyper(x - upper.tail, m, n, k, lower.tail = !upper.tail)
  } else if (ncp == 0) {
    return(as.numeric(if (upper.tail) q <= lo else q >= lo))
  } else if (ncp == Inf) {
    return(as.numeric(if (upper.tail) q <= hi else q >= hi))
  } else {
    return(sum(dnhyper(ncp)[if (upper.tail) support >= q else support <= q]))
  }
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

  core <- function(m, n, k, x, lo, hi) {

  }

  mapply(function(
    lo, hi
  ) {
    support <- lo:hi
    logdc <- dhyper(support, m, n, k, log = TRUE)
    dnhyper <- function(ncp) {
      d <- logdc + log(ncp) * support
      d <- exp(d - max(d))
      d/sum(d)
    }
    mnhyper <- function(ncp) {
      if (ncp == 0) return(lo)
      if (ncp == Inf) return(hi)
      sum(support * dnhyper(ncp))
    }
    pnhyper <- function(q, ncp = 1, upper.tail = FALSE) {
      if (ncp == 1) {
        phyper(x - upper.tail, m, n, k, lower.tail = !upper.tail)
      }
      if (ncp == 0) {
        return(as.numeric(if (upper.tail) q <= lo else q >= lo))
      }
      if (ncp == Inf) {
        return(as.numeric(if (upper.tail) q <= hi else q >= hi))
      }
      sum(dnhyper(ncp)[if (upper.tail) support >= q else support <= q])
    }
    p.value <- switch(alternative,
                      less = pnhyper(x, or), greater = pnhyper(x, or, upper.tail = TRUE),
                      two.sided = {
                        if (or == 0) as.numeric(x == lo) else if (or == Inf) as.numeric(x == hi) else {
                          relErr <- 1 + 10^(-7)
                          d <- dnhyper(or)
                          sum(d[d <= d[x - lo + 1] * relErr])
                        }
                      }
    )
    mle <- if (x == lo) { 0 } else if (x == hi) { Inf } else {
      mu <- mnhyper(1)
      if (mu > x)
        uniroot(function(t) mnhyper(t) - x, c(0, 1))$root
      else if (mu < x)
        1/uniroot(function(t) mnhyper(1/t) - x, c(.Machine$double.eps, 1))$root
      else 1
    }
    list(estimate = mle, p.value = p.value)
  }, lo = lo, hi = hi)


  mapply(
    \(lo, hi) {
      support <- lo:hi
      logdc <- dhyper(support, m, n, k, log = TRUE)
      dnhyper <- function(ncp) {
        d <- logdc + log(ncp) * support
        d <- exp(d - max(d))
        d/sum(d)
      }
      mnhyper <- function(ncp) {
        if (ncp == 0)
          return(lo)
        if (ncp == Inf)
          return(hi)
        sum(support * dnhyper(ncp))
      }
      pnhyper <- function(q, ncp = 1, upper.tail = FALSE) {
        if (ncp == 1) {
          return(if (upper.tail) phyper(x - 1, m, n, k,
                                        lower.tail = FALSE) else phyper(x, m, n, k))
        }
        if (ncp == 0) {
          return(as.numeric(if (upper.tail) q <= lo else q >=
                              lo))
        }
        if (ncp == Inf) {
          return(as.numeric(if (upper.tail) q <= hi else q >=
                              hi))
        }
        sum(dnhyper(ncp)[if (upper.tail) support >= q else support <=
                           q])
      }


      list(estimate = mle, p.value = p.value)
    }, lo, hi)
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
