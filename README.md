The R Script does a number of items to analyze various market segments across multiple time frames.
 * Tracks 7 different ETFs, following segments such as large cap stocks, REITs, corporate bonds, etc.
 * Breaks the data into 5 different time frames: Daily, Weekly, Monthly, Quarterly, and Annual
 * Sets 3 different return targets (annualized for smaller time frames)
     - 2%
     - 4%
     - 6%
 * Uses a solver to optimize the portfolios with the following parameters:
     - Full investment requirement
     - Long-only portfolio
     - 10 year data period
 * Then completes separate solutions to minimize the following risk measures:
     - Variance
     - Mean Absolute Deviation
     - Conditional Value at Risk

In total, 45 different portfolio allocations are created.
  
The obtained portfolio weights are then backtested against the previous data and produce the following statistics for analysis:
  * Expected Return
  * Standard Deviation
  * Mean Absolute Deviation
  * Value at Risk
  * Conditional Value at Risk
