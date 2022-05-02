# test post hoc test result

    Code
      print(res_test, digits = 5)
    Output
                                              comparisons diff_mean     pvalue
      Enterotype 2-Enterotype 1 Enterotype 2-Enterotype 1  -0.28139 4.7701e-08
      Enterotype 3-Enterotype 1 Enterotype 3-Enterotype 1  -0.26045 1.6364e-09
      Enterotype 3-Enterotype 2 Enterotype 3-Enterotype 2   0.02094 7.8899e-01
                                 ci_lower  ci_upper
      Enterotype 2-Enterotype 1 -0.371347 -0.191443
      Enterotype 3-Enterotype 1 -0.331229 -0.189681
      Enterotype 3-Enterotype 2 -0.057576  0.099457

# test visualization of posthoc test, data of signicance level annotation

    Code
      annotation_single
    Output
                xmin         xmax y_position annotation
      1 Enterotype 2 Enterotype 1    0.54668        ***
      2 Enterotype 3 Enterotype 1    0.61035        ***
      3 Enterotype 3 Enterotype 2    0.30558        NS.

---

    Code
      head(annotation_all)
    Output
                xmin         xmax y_position annotation           feature
      1 Enterotype 2 Enterotype 1    0.13116        NS. p__Actinobacteria
      2 Enterotype 3 Enterotype 1     0.2694        NS. p__Actinobacteria
      3 Enterotype 3 Enterotype 2    0.29804        NS. p__Actinobacteria
      4 Enterotype 2 Enterotype 1    0.57527        NS.  p__Bacteroidetes
      5 Enterotype 3 Enterotype 1    0.63693        ***  p__Bacteroidetes
      6 Enterotype 3 Enterotype 2    0.56712         **  p__Bacteroidetes

