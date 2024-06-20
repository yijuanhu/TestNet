context("Testing `TestNet` function")
library(TestNet)
library(testthat)

data(sim.otu.tab)

# test
test_that("`TestNet` function provides expected results", {
    system.time(res.TestNet <- TestNet(otu.tab=sim.otu.tab))
    p.omni.lower.tri <- res.TestNet$p.omni[lower.tri(res.TestNet$p.omni)]
    p.omni.lower.tri <- signif(p.omni.lower.tri, 3)
    expect_equivalent(p.omni.lower.tri[1:3], c(0.38, 0.65, 0.80))
})
