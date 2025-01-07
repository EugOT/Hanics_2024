#!/usr/bin/env Rscript
# render.R

workflowr::wflow_build(here::here("analysis/01A-eda.Rmd"), verbose = T)
