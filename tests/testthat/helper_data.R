data("oxygen")
data("pediatric_ibd")

lefse_out <- lefse(
  oxygen,
  normalization = 1e6,
  bootstrap_n = 5,
  summarize = "lefse",
  class = "oxygen_availability",
  subclass = "body_site"
)
