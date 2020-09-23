lefse_out <- lefse(
  oxygen,
  normalization = 1e6,
  bootstrap_n = 5,
  summarize = FALSE,
  class = "oxygen_availability",
  subclass = "body_site"
)
