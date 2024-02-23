# OneStep - Version 0.9.3

## NEW

-   add a optim control parameter in control list of `onestep()`. supported control are trace, fnscale, maxit, abstol, reltol, alpha, beta, gamma, REPORT, warn.1d.NelderMead, type, lmm, factr, pgtol, tmax, temp
-   use `log=TRUE` argument when computing log-likelihood
-   add student `t` distribution
-   add location scaled student `lst` distribution
-   add `nbstep` component in `onestep()` output

# OneStep - Version 0.9.2

## NEW

-   allow initial guess supplied by the user.
-   add defensive programming.
-   update manual page.
-   take into account weights in `benchonestep()`.
-   new function in `benchonestep.replicate()`.

## BUG FIX

-   bug fixed in `benchonestep()` for one-parameter distribution.

# OneStep - Version 0.9.1

## BUG FIX

-   bug fixed in `onestep_closedformula()` for the Weibull distribution.

# OneStep - Version 0.9.0

-   Initial CRAN release.
