library(dplyr)
aws.s3::s3load(
  object = file.path("sierra-2020-bm-test", "models", "percent_canopy_data.Rda"),
  bucket = "silviaterra-sequoia"
)

canopy_cover_table <- training_data %>%
  select(
    bapa, tpa, cc_prop = percent_cc
  )

usethis::use_data(canopy_cover_table)
