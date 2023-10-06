# Time-limited Balanced Truncation within Incremental Four-Dimensional Variational Data Assimilation

This MATLAB repository contains code for the numerical results of the following paper:

1.  König, J. & Freitag, M. A. "[Time-limited balanced truncation within incremental four-dimensional variational data assimilation](https://onlinelibrary.wiley.com/doi/10.1002/pamm.202300019)". Proceedings in Applied Mathematics and Mechanics, e202300019 (2023).

## Summary
The work in [1] applies time-limited balanced truncation within the inner loop of incremental four-dimensional variational data assimilation. Results and code for (time-limited) balanced truncation for Bayesian inference from [2] (to be found at
[https://github.com/joskoUP/TLBTforDA](https://github.com/joskoUP/TLBTforDA)) and [3] (to be found at
[https://github.com/elizqian/balancing-bayesian-inference](https://github.com/elizqian/balancing-bayesian-inference)) allow to simplify the problem and use a Quasi-Newton minimisation with a constant Hessian.

## Examples
To run this code, you need the MATLAB Control System Toolbox.

To generate the plots comparing the full incremental 4D-Var method to time-limited balanced truncation and $\alpha$-bounded balanced truncation, run the **Comparison_lorenzN_4dvar.m** script.<br />
There will be a pop-up menu. Choose the observation points and prior covariance as follows:
* Figure 1: Observations _in all N variables_, prior covariance matrix _Gaussian exponential_.
* Figure 2: Observations _every s points_, in command window: `Take observations every ... points.` - enter _9_, prior covariance matrix _scaled identity_.

## References
2. König, J., Freitag, M. A. "[Time-limited Balanced Truncation for Data Assimilation Problems](https://link.springer.com/article/10.1007/s10915-023-02358-4)." Journal of Scientific Computing 97.47 (2023).
3. Qian, E., Tabeart, J. M., Beattie, C., Gugercin, S., Jiang, J., Kramer, P. R., & Narayan, A.
"[Model reduction for linear dynamical systems via balancing for Bayesian inference](https://link.springer.com/article/10.1007/s10915-022-01798-8)." Journal of Scientific Computing 91.29 (2022).

### Contact
Please feel free to contact [Josie König](https://www.math.uni-potsdam.de/professuren/datenassimilation/personen/josie-koenig) with any questions about this repository or the associated paper.
