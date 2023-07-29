# Integration step

- subroutine `adjust_eps()` in `hmc.F` adjusts integration step `eps` after every 30 calls to get acceptance rate of about 80%
- graph `eps.pdf` illustrates how integration step changes every 30 measurement steps (1000 measurements at `beta` 2, 4 and 6 and `kappa = 0.25`)
- graph `eps-kappa-beta.pdf` shows how integration step depends on `kappa` and `beta` when it is adjusted to acceptance rate of about 80%

*Hip, 2023-07-29*