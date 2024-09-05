# T2D Discontinuation

This repository contains the code for investigating prediction of discontinuation in T2D therapies.

## Files included:

0.  `set_up_data_and_functions.R`: This file contains the function used to select the right patients for analysis.
1.  `characteristics_tables.R`: This file contains the code required to create characteristics tables used in the manuscript.
2.  `propensity_score_lm.R`: This file contains the code required to fit a multivariate / multi-outcome propensisty score model (not used for this analysis).
3.  `univariate_analysis.R`: This file contains the code required to conduct univariate analysis between clinical features and discontinuation rates.
4.  `multivariate_analysis.R`: This file contains the code required to fit a multivariate outcome model using linear regression (not used for this analysis).
5.  `heterogeneity_analysis.R`: This file contains the code for validating average treatment effects overall (not used for this analysis).
6.  `heterogeneity_absolute_relative_risk.R`: This file contains the code for validating average treatment effects on a drug vs drug basis (not used for this analysis).
7.  `bart_multivariate_analysis.R`: This file contains the code for fitting BART models to the data at 3-month/6-month/12-month outcomes.
