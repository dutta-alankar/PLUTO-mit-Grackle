# `PLUTO` astrophysics MHD code integrated with external `Grackle` radiative cooling library
PLUTO version: 4.4 patch 3

## Test plot with equilibrium cooling at solar metallicity and FG2011 photo-ionizing UV background radiation
Details of the test:
- A few gas cells are initialized at $T=2\times 10^6 \ \rm K$ in one-dimension.
- The gas density is set to $10^{-3}\ \rm cm^{-3}$ at rest.
- The gas cools according to $\dot{e}=-n_H^2 \Lambda(T)$, where $e$ is the internal energy density.
- Gas temperature and $\dot{e}$ (averaged over all cells) are calculated and plotted over a time of about 1 Gyr.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://github.com/user-attachments/assets/739840b2-1a82-4b7d-8b1b-eaf31af46442">
  <img alt="" src="https://github.com/dutta-alankar/PLUTO-mit-Grackle/blob/31ad1a410f2d141100b155f66631a9f60c571d2e/cooling-function.png">
</picture>
