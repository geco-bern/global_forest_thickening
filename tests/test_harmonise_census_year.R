suppressPackageStartupMessages(library(dplyr))
source("R/harmonise_census_year.R")
x <- tibble(
  plotID=c("a","a","a","a","a","b","b","c","c","d","d","d"),
  census=c(1,1,1,1,2,1,1,1,1,1,1,1),
  year=c(2001,2000,2001,NA,2010,2005,2000,NA,NA,1999,NA,1999),
  row_id=1:12
)
y <- harmonise_census_year(x)
stopifnot(
  identical(y$row_id,x$row_id),
  identical(y$year,c(2001,2001,2001,2001,2010,2000,2000,NA,NA,1999,1999,1999)),
  identical(harmonise_census_year(x[12:1,]) |> arrange(row_id),y),
  identical(harmonise_census_year(y),y)
)
# Support source column names without renaming the input.
raw <- rename(x,Quadrat=plotID,Census=census,Year=year)
stopifnot(identical(harmonise_census_year(raw,"Quadrat","Census","Year")$Year,y$year))
cat("Modal census years: modes, ties, missing values, groups and order passed.\n")
