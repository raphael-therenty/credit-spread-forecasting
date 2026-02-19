# Credit Conditions and Economic Activity in the US: The Role of the Corporate Spread (2010–2019)

## Project Overview
This research examines the predictive power of the corporate credit spread on US economic activity during the post-subprime crisis period (**2010–2019**). 

The study utilizes monthly time-series data, focusing on:
* **Industrial Production (IP):** Proxy for real economic activity.
* **BAA-10Y Spread:** Proxy for credit risk and financing conditions.
* **Fed Balance Sheet:** Accounting for unconventional monetary policy.

## Methodology
The project follows a rigorous econometric workflow to ensure robust results:

1.  **Stationarity Testing:** Augmented Dickey-Fuller (ADF) and Phillips-Perron tests to determine integration orders.
2.  **Univariate Modeling:** Estimation of an **ARIMA** model for baseline forecasting.
3.  **Multivariate Dynamics:** Implementation of a **Vector Autoregression (VAR)** model to capture interdependencies.
4.  **Long-term Equilibrium:** * **Engle-Granger** & **Johansen** Cointegration tests.
    * Estimation of an **Error Correction Model (ECM)**.



## Key Findings
* **Cointegration:** A significant long-term relationship exists between credit spreads and industrial production.
* **Economic Impact:** While short-term fluctuations in the spread have limited immediate effects, a **structural degradation** of credit conditions acts as a major, persistent drag on real activity.
* **Error Correction:** The adjustment mechanism back to equilibrium is relatively slow, with a coefficient of $\alpha \approx -0.071$. This implies that roughly **7.1%** of the deviation from the long-term trend is corrected each month.

## Project Structure
* `data/`: Raw and processed CSV files.
* `scripts/`: R scripts for data cleaning, VAR/VECM estimation, and diagnostic tests.
* `results/`: The final research paper/memo in PDF format.

## Requirements
To replicate this analysis, you will need **R** and the following packages:
```R
install.packages(c("urca", "vars", "tseries", "forecast", "tidyverse", "tsDyn"))
