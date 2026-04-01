library(usethis)

use_roxygen_md()
use_testthat(3)
use_mit_license("Yilin Chiu")
use_readme_rmd()
use_news_md()
use_git()



usethis::use_git_config(
  user.name = "Yilin Chiu",
  user.email = "yilin1107@mail.ndmutsgh.edu.tw"
)

system("git config --global --get user.name")
system("git config --global --get user.email")



system("git status")
system('git add .')
system('git commit -m "Initial package skeleton for L1000PAI"')

usethis::use_git()


devtools::document()
devtools::load_all()
devtools::check()

devtools::document()
devtools::load_all()
devtools::test()
devtools::check()



if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("clusterProfiler")


usethis::use_package("clusterProfiler", type = "Imports")







devtools::install()



install.packages(c("usethis", "devtools", "roxygen2", "testthat", "knitr", "rmarkdown"))

library(usethis)

create_package("D:/Rpackages/L1000PAIdata")


devtools::document()
devtools::build_vignettes()
devtools::check()



gmt_crispr <- get_l1000_library("crispr")
head(gmt_crispr)

gmt_chem <- get_l1000_library("chem")
head(gmt_chem)

gmt_kg <- get_l1000_library("kinase_gpcr")
head(gmt_kg)



devtools::check()




remove.packages("L1000PAI")


getwd()
list.files()


devtools::install()
library(L1000PAI)
devtools::document()
devtools::build_readme()
devtools::build_vignettes()
devtools::check()


resolve_score_direction("chem", "auto")


devtools::document()
devtools::load_all()
devtools::test()
devtools::check()







