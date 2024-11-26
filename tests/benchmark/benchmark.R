
library(vfisher)
library(microbenchmark)

# create an OR + n1, n2, "sample" around it, run both ways

n1 <- rnbinom(1000, 30, 0.7)
n2 <- rnbinom(1000, 30, 0.7)
or <- rgamma(1000, shape = 2.5, scale = 2)
p1 <- rbeta(1000, 3, 6)
p2 = ((p1/(1-p1))/or)/(1 + (p1/(1-p1))/or)

a <- rbinom(1000, n1, p1)
b <- n1 - a
c <- rbinom(1000, n2, p2)
d <- n2 - c

matrix(c(a[1], c[1], b[1], d[1]), nrow = 2)

t(mapply(\(a, b, c, d) {
  res <- fisher.test(
    x = matrix(c(a, c, b, d), nrow = 2)
  )
  c(res$estimate, res$conf.int)
}, a, b, c, d))

vfisher.test(a, b, c, d)
