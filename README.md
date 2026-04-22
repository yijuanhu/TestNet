# TestNet

This package implements the testing method, TestNet, for inferring microbial networks. It differs from existing microbial network analyses in that it provides calibrated results by controlling the false discovery rate. TestNet accounts for the features of compositionality, sparsity, and overdispersion in taxa count data. It also accommodates both independent and clustered samples, offers separate linear and nonlinear tests for each pair of taxa, and includes an omnibus test that bypasses the need to pre-specify the type of relationship for each pair of taxa.


To install the package:

devtools::install_github("yijuanhu/TestNet")

