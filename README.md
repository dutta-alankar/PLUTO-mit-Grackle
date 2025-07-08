# `PLUTO` astrophysics MHD code integrated with external `Grackle` radiative cooling library
PLUTO version: 4.4 patch 3

## Test plot with equilibrium cooling at solar metallicity and FG2011 photo-ionizing UV background radiation
Details of the test:
- A few gas cells are initialized at $T=2\times 10^6 \ \rm K$ in one-dimension (3D also possible).
- The gas density is set to $10^{-2}\ \rm cm^{-3}$ at rest.
- The gas cools according to $\dot{e}=-n_H^2 \Lambda(T)$, where $e$ is the internal energy density.
- Gas temperature and $\dot{e}$ (averaged over all cells) are calculated and plotted over a time of about 1 Gyr.
- Set velocity in the positive x-direction by 1 code unit with periodic boundary to test the flux of ions.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://github.com/user-attachments/assets/aa3e5a9c-0f19-4d26-8392-ed94a73d47dd">
  <img alt="" src="https://github.com/user-attachments/assets/f6c64fdd-1a2f-48f6-ae8f-3c52e91a7c8d">
</picture>

