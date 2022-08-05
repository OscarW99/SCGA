# Used this to generate API token before downloading package: https://combiz.github.io/scFlow/ .. section == 'How to contribure'


devtools::install_github("combiz/scFlow", auth_token = "HIDDEN")

# Debugging the scFlow downlaod
# https://stackoverflow.com/questions/59303030/devtools-error-in-read-dcfpath-desc-line-starting-this-corresponds-to 



# Github open issue:

# > devtools::install_github("combiz/scFlow")
# Downloading GitHub repo combiz/scFlow@HEAD
# Error: Failed to install 'scFlow' from GitHub:
#   Line starting 'Roxyg ...' is malformed!


# > sessionInfo()
# R version 4.1.1 (2021-08-10)
# Platform: x86_64-conda-linux-gnu (64-bit)
# Running under: Red Hat Enterprise Linux

# Matrix products: default
# BLAS/LAPACK: /ahg/regevdata/projects/ICA_Lung/Oscar/conda/conda_env/lib/libopenblasp-r0.3.18.so

# locale:
# [1] C

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base

# loaded via a namespace (and not attached):
#  [1] ps_1.7.0          prettyunits_1.1.1 rprojroot_2.0.3   crayon_1.5.1
#  [5] withr_2.5.0       brio_1.1.3        R6_2.5.1          lifecycle_1.0.1
#  [9] magrittr_2.0.3    rlang_1.0.2       cachem_1.0.6      cli_3.3.0
# [13] curl_4.3.2        remotes_2.4.2     fs_1.5.2          testthat_3.1.4
# [17] callr_3.7.0       ellipsis_0.3.2    desc_1.4.1        devtools_2.4.3
# [21] tools_4.1.1       glue_1.6.2        purrr_0.3.4       pkgload_1.2.4
# [25] fastmap_1.1.0     compiler_4.1.1    processx_3.6.0    pkgbuild_1.3.1
# [29] sessioninfo_1.2.2 memoise_2.0.1     usethis_2.1.6
