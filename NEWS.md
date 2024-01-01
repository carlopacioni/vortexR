# vortexR 2.0.0
* Add CI checks via GitHub Actions.
* Add test coverage reporting at codecov.io via GitHub Actions. 
  Requires a repository secret via GitHub repo > Settings > Secrets and Variables > Actions > New repository secret
  at the [upstream Gh repo settings](https://github.com/carlopacioni/vortexR/settings/secrets/actions) using the 
  Token from [codecov](https://app.codecov.io/gh/carlopacioni/vortexR/settings).

# vortexR 1.1.8
collate_* functions gain dec_sep argument to control decimal separator charcter in order to solve reading numerical data on European machines
Minor edits to help files

# vortexR 1.1.7
Ensured copatibility with R 4.*
fit_regression gains a new argument: `links`, which allows for a better control
on how the Beta regression is fitted.
Improved testing for fit_regression (this is for development only, not visible to users)
Updated references

# vortexR 1.1.6
Ensured copatibility with irr 0.84.1

# vortexR 1.1.5
Ensured copatibility with R 3.6.*

# vortexR 1.1.4
Implemented SSMD_matrix. Calculates a matrix of SSMD and associated pvalues

# vortexR 1.0.4
Fix compatibility issue with new version of GGally

# vortexR 1.0.3
CRAN release
