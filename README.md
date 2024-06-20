# TestNet

The TestNet package implements the testing method, called TestNet, for inference of microbial networks. TestNet provides calibrated results by controlling the false discovery rate. It accounts for the features of compositionality, sparsity, and overdispersion in microbiome read count data. It accommodates both independent and clustered samples, offers separate linear and nonlinear tests, and includes an omnibus test that eliminates the need to pre-specify the type of relationship.


To install the package:

devtools::install_github("yijuanhu/TestNet", build_vignettes=TRUE)

