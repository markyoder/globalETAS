### Global (not-)ETAS

This repository contains Python for an earthquake seismicity model, based on aftershock scaling. The model was published in
[Yoder et al. (2015)](https://doi.org/10.1007/s00024-014-0785-z). It introduces a near-field constraint based on the finite size
and temporeal duration of an individual earthquake, which facilitates estimates of immediately post-mainshock seismicity rates. 
This can significantly increase compute time performance by avoiding fitting algorithms.

A map of expected seismicity is producec by summing the contributions of all events at all sites on the map. Because all events can be used to 
directly contribute to rate estimates, we can achieve higher resolutions than most contemporary implementations of ETAS and 
ETAS-like methods. Additionally, because there is no data fitting required, this method is relatively insensitive to catalog 
incompleteness and near-field catalog anomalies, immediately following and near-by large earthquakes. 
