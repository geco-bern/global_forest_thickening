suppressPackageStartupMessages({library(testthat); library(dplyr)})
source('R/forestplots_management.R')
meta <- tibble(
 plotID = LETTERS[1:10], forest_type = c(rep('Lowland Tropical Forest',8),'Savanna',''),
 anthropogenic_effects=c(' Old-Growth ', 'Reduced impact logging: 10+ year post-logging',
 'Lightly logged > 30 years ago','3-400 year old fire/logging','Logged- Logging experiment 1994',
 'Reduced Impact Logging: Logged  around 2016','Not logged','5 years since distubance','Intact','Intact'),
 Monitoring.Start.Date=rep('2012',10))
x <- classify_forestplots_management(meta)
expect_equal(nrow(x),8L)
expect_equal(x$management_cat,c(2L,1L,1L,1L,1L,3L,3L,1L))
expect_equal(x$management_since_census1_yrs,c(NA,10,30,300,18,NA,NA,5))
expect_equal(nrow(classify_forestplots_management(meta,keep_forest=FALSE)),10L)
# Actual census years override monitoring-start metadata for dated events.
trees <- tidyr::expand_grid(plotID=LETTERS[1:10],year=c(2010,2015,2020))
y <- join_forestplots_metadata(trees,x)
expect_equal(nrow(y),24L)
expect_equal(unique(y$management_since_census1_yrs[y$plotID=='E']),16)
expect_true(all(y$management[y$plotID=='F']==1L))
expect_true(all(y$management_cat[y$plotID=='F']==3L))
expect_true(all(is.na(y$management_since_census1_yrs[y$plotID=='F'])))
expect_true(all(y$management[y$plotID=='A']==0L))
# Metadata survives plot/subplot IDs and generic aggregation.
y$plotID <- paste0(y$plotID,'_cell1')
stands <- distinct(y,plotID,year) |> mutate(management_since_census1_yrs=NaN)
z <- forestplots_stand_management(stands,y)
expect_equal(z$management_cat,y$management_cat)
expect_equal(z$management_since_census1_yrs,y$management_since_census1_yrs)
expect_equal(nrow(z),nrow(y))
# Unknown future descriptions stay category 3; a URL is not a census year.
unknown <- meta[1,];unknown$anthropogenic_effects <- 'New undocumented description'
expect_equal(classify_forestplots_management(unknown)$management_cat,3L)
unknown$anthropogenic_effects <- 'Logged 2015'; unknown$Monitoring.Start.Date <- 'https://example/2018'
expect_true(is.na(classify_forestplots_management(unknown)$management_since_census1_yrs))
expect_equal(nrow(classify_forestplots_management(meta[FALSE,])),0L)
cat('History lookup, forest selection, dated events, baseline anchoring and aggregation checks passed.\n')
