# Acceptance check

- subroutine `adjust_eps()` from `hmc.F` adjusts `eps` after every 30 calls with the aim to get about 80% acceptance rate
- to check if it is working properly the runs were performed with `eps` input parameter equal to 1 (this is very large, so the acceptance rate is expected to be zero at the start of the iteration)
- graph `acceptance.pdf` shows that already after about 150 iterations `eps` is properly adjusted to get about 80% acceptance rate and it seems that this is `beta` independent
- **conclusion:** input parameter `eps` is rather useles because it is quickly adjusted - usually already during the thermalization

*Hip, 2023-07-29*
