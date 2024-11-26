
test_that("vfisher gives same results as fisher.test example 1", {
  TeaTasting <- matrix(
    c(3, 1, 1, 3), nrow = 2, dimnames = list(
      Guess = c("Milk", "Tea"),
      Truth = c("Milk", "Tea")
  ))
  fref <- fisher.test(TeaTasting, alternative = "greater")
  vres <- vfisher.test(
    TeaTasting[1,1], TeaTasting[1, 2], TeaTasting[2, 1], TeaTasting[2, 2], alternative = "greater"
  )
  expect_equal(fref$p.value, vres[, p.value])
  expect_equal(unname(fref$estimate), vres[, estimate])
  expect_equal(fref$conf.int[1], vres[, ci.lo])
  expect_equal(fref$conf.int[2], vres[, ci.hi])
})

test_that("vfisher gives same results as fisher.test example 2", {
  Convictions <- matrix(c(2, 10, 15, 3), nrow = 2,
                        dimnames =
                          list(c("Dizygotic", "Monozygotic"),
                               c("Convicted", "Not convicted")))
  fref <- fisher.test(Convictions, alternative = "less")
  vres <- vfisher.test(
    Convictions[1,1], Convictions[1, 2], Convictions[2, 1], Convictions[2, 2], alternative = "less"
  )
  expect_equal(fref$p.value, vres[, p.value])
  expect_equal(unname(fref$estimate), vres[, estimate])
  expect_equal(fref$conf.int[1], vres[, ci.lo])
  expect_equal(fref$conf.int[2], vres[, ci.hi])
})

test_that("vfisher gives same results as fisher.test example 1 & 2", {
  TeaTasting <- matrix(
    c(3, 1, 1, 3), nrow = 2, dimnames = list(
      Guess = c("Milk", "Tea"),
      Truth = c("Milk", "Tea")
    ))
  Convictions <- matrix(c(2, 10, 15, 3), nrow = 2,
                        dimnames =
                          list(c("Dizygotic", "Monozygotic"),
                               c("Convicted", "Not convicted")))
  fref1 <- fisher.test(TeaTasting)
  fref2 <- fisher.test(Convictions)

  vres <- vfisher.test(
    c(TeaTasting[1,1], Convictions[1,1]),
    c(TeaTasting[1,2], Convictions[1,2]),
    c(TeaTasting[2,1], Convictions[2,1]),
    c(TeaTasting[2,2], Convictions[2,2])
  )

  expect_equal(fref1$p.value, vres[1, p.value])
  expect_equal(fref2$p.value, vres[2, p.value])
  expect_equal(unname(fref1$estimate), vres[1, estimate])
  expect_equal(unname(fref2$estimate), vres[2, estimate])
  expect_equal(as.numeric(fref1$conf.int), vres[1, c(ci.lo, ci.hi)])
  expect_equal(as.numeric(fref2$conf.int), vres[2, c(ci.lo, ci.hi)])
})
