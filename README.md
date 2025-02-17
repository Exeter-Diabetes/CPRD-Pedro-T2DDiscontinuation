# T2D Discontinuation

This repository contains the code for investigating prediction of discontinuation in T2D therapies.

## Files included:

-   `00.set_up_data_and_functions.R`: This file contains the function used to select the right patients for analysis.
-   `01.characteristics_tables.R`: This file contains the code required to create characteristics tables used in the manuscript.
-   `02.propensity_score_lm.R`: This file contains the code required to fit a multivariate / multi-outcome propensisty score model (not used for the final analysis).
-   `03.univariate_analysis.R`: This file contains the code required to conduct univariate analysis between clinical features and discontinuation rates.
-   `04.multivariate_analysis.R`: This file contains the code required to fit a multivariate outcome model using linear regression (not used for the final analysis).
-   `05.heterogeneity_analysis.R`: This file contains the code for validating average treatment effects overall (not used for the final analysis).
-   `06.bart_multivariate_analysis.R`: This file contains the code for fitting BART models to the data at 3-month/6-month/12-month outcomes.
-   `07.01.heterogeneity_bart_model_3m.R`: This file contains the code for analysing heterogeneity in 3-months discontinuation.
-   `07.02.heterogeneity_bart_model_6m.R`: This file contains the code for analysing heterogeneity in 6-months discontinuation.
-   `07.03.heterogeneity_bart_model_12m.R`: This file contains the code for analysing heterogeneity in 12-months discontinuation.
