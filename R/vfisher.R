
logdc_calc <- function(m, n, k, support) {
  return(mapply(
    function(support, m, n, k) dhyper(support, m, n, k, log = TRUE),
    support = support, m = m, n = n, k = k, SIMPLIFY = FALSE
  ))
}

dnhyper <- function(ncp, logdc, support) {
    d <- logdc + log(ncp) * support
    d <- exp(d - max(d))
    return(d/sum(d))
}

mnhyper <- function(ncp, lo, hi, logdc, support) {
  if (ncp == 0) {
    return(lo)
  } else if (is.infinite(ncp)) {
    return(hi)
  } else {
    return(sum(support * dnhyper(ncp, logdc, support)))
  }
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

# ncp.U == ncp_ci(..., lower = FALSE)
# ncp.L == ncp_ci(..., lower = TRUE)
ncp_ci <- function(x, alpha, m, n, k, lo, hi, support, logdc, lower = FALSE) {
  if (x == hi) {
    return(Inf)
  } else {
    p <- pnhyper(x, 1, upper.tail = lower, m, n, k, lo, hi, support, logdc)
    ple <- p < alpha
    if (p == alpha) {
      return(1)
    } else if (ple != lower) {
      return(uniroot(function(t) pnhyper(x, t, upper.tail = lower, m, n, k, lo, hi, support, logdc) - alpha, c(0, 1))$root)
    } else { # ple == lower
      return(1/uniroot(function(t) pnhyper(x, 1/t, upper.tail = lower, m, n, k, lo, hi, support, logdc) - alpha, c(.Machine$double.eps, 1))$root)
    }
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
#' @return data.table, columns `a`, `b`, `c`, `d`, `or`, `p.value`, `estimate`.
#'
#' If `conf.int == TRUE` (the default), will also include columns for the
#' confidence interval, `ci.lo` and `ci.hi`. If `alternative == "less"`, `ci.lo`
#' will be `0`, and similarly for `alternative == "greater"`, `ci.hi` will be
#' `Inf`. Otherwise (the default), the CI will be centered and both low and high
#' ends will take on values between 0 and Inf.
#'
#' @export
#' @import data.table
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

  if (any(!is.numeric(or) | is.na(or) | or < 0))
    stop("'or' must be a non-NA number between 0 and Inf")

  if (length(or) != length(a)) {
    warning("`length(or) != length(a)`; using `rep(or, length.out = length(a))` to extend.")
    or <- rep(or, length.out = length(a))
  }

  # matrix =
  # a, b
  # c, d

  result_dt <- data.table(a = a, b = b, c = c, d = d, or = or)
  result_dt[,
    m := a + c
  ][,
    n := b + d
  ][,
    k := a + b
  ][,
    lo := pmax(0L, k - n)
  ][,
    hi := pmin(k, m)
  ]

  result_dt[, rowid := 1:.N]

  result_dt[,
    support := .(list(lo:hi)), by = rowid
  ][,
    logdc :=  .(list(dhyper(support[[1]], m, n, k, log = TRUE))), by = rowid
  ]

  # x == a

  result_dt[a == lo, estimate := 0]
  result_dt[a == hi, estimate := Inf]
  result_dt[!(a == lo | a == hi), mnhyper1 := mnhyper(1, lo, hi, logdc[[1]], support[[1]]), by = rowid]
  result_dt[mnhyper1 == a, estimate := 1]
  result_dt[mnhyper1 < a, estimate := 1/uniroot(function(t) mnhyper(1/t, lo, hi, logdc[[1]], support[[1]]) - a, c(.Machine$double.eps, 1))$root, by = rowid]
  result_dt[mnhyper1 > a, estimate := uniroot(function(t) mnhyper(t, lo, hi, logdc[[1]], support[[1]]) - a, c(0, 1))$root, by = rowid ]

  sdcols <- c("a", "b", "c", "d", "or", "p.value", "estimate")

  if (conf.int) {

    alternative <- match.arg(alternative)

    if (!((length(conf.level) == 1L) && is.finite(conf.level) &&
          (conf.level > 0) && (conf.level < 1)))
      stop("'conf.level' must be a single number between 0 and 1")

    sdcols <- c(sdcols, c("ci.lo", "ci.hi"))

    setattr(result_dt, "conf.level", conf.level)

    if (alternative == "less") {
      result_dt[, ci.lo := 0]
      result_dt[, ci.hi := ncp_ci(a, 1 - conf.level, m, n, k, lo, hi, support[[1]], logdc[[1]], lower = FALSE), by = rowid]
    } else if (alternative == "greater") {
      result_dt[, ci.lo := ncp_ci(a, 1 - conf.level, m, n, k, lo, hi, support[[1]], logdc[[1]], lower = TRUE), by = rowid]
      result_dt[, ci.hi := Inf]
    } else {
      alpha <- (1 - conf.level)/2
      result_dt[, ci.lo := ncp_ci(a, alpha, m, n, k, lo, hi, support[[1]], logdc[[1]], lower = TRUE), by = rowid]
      result_dt[, ci.hi := ncp_ci(a, alpha, m, n, k, lo, hi, support[[1]], logdc[[1]], lower = FALSE), by = rowid]
    }

  }

  result_dt[, .SD, .SDcols = sdcols]

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
