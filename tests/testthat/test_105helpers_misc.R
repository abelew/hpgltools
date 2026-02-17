start <- as.POSIXlt(Sys.time())
context("105helpers_misc.R")

## This function actually makes untenable assumptions about where the repository lives.
testing <- get_git_commit(gitdir = "")
expected <- c("glue", "character")
actual <- class(testing)
test_that("get_git_commit() gave me a commit id?", {
  expect_equal(expected, actual)
})

## So, a bunch of functions exported in helpers are difficult to test because
## they are pretty specific to their little domains.  I probably therefore will
## skip testing some (many) of them.
end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 105helpers_misc.R in ", elapsed,  " seconds.")
