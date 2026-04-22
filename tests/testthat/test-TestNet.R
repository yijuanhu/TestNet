context("Testing `TestNet` function")
library(TestNet)
library(testthat)

data(sim.otu.tab)
sim.otu.tab.sub <- sim.otu.tab[1:20,]

# test
test_that("`TestNet` function provides expected results", {
    system.time(res.TestNet <- TestNet(otu.tab=sim.otu.tab.sub))
    p.omni.lower.tri <- res.TestNet$p.omni[lower.tri(res.TestNet$p.omni)]
    p.omni.lower.tri <- signif(p.omni.lower.tri, 3)
    expect_equivalent(p.omni.lower.tri[1:3], c(0.205, 0.400, 0.330))
})
