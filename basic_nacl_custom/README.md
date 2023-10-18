# Tutorial 5.5: Creating "binless" resampling schemes: Na+/Cl- association simulations
This folder contains the same basic NaCl tutorial from Tutorial 5.5.

The custom distance scheme is enabled in this version such that segments are organized pair-wise based on similarity and merged based on that hierarchy, as indicated in sort.py. You will have to use the `custom_order_resampler` branch from https://github.com/jeremyleung/westpa to use sort.py. The lack of "80" in the name indicate that the gamma_ln collision frequency is set to 5 ps^{-1}, a value used when there is explicit water modeled. 




## Authors

* **Jeremy Leung** - *Primary work* - [jeremyleung521](https://github.com/jeremyleung521)
* **Anthony Bogetti** - *Primary work* - [atbogetti](https://github.com/atbogetti)
