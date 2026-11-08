# Licence and attribution

The files in this folder are by David B. Davidson, Curtin University, from:

D.B. Davidson, "MATLAB files for reconstructing antenna patterns from
spherical mode data", Curtin University research data, July 2024.
https://doi.org/10.25917/N05J-TH75

They accompany D.B. Davidson and A.T. Sutinjo, "Efficient storage of
embedded element patterns for low frequency radio telescopes" (prepared
for *Radio Science*).

The dataset is licensed under the Creative Commons Attribution 4.0
International licence (CC BY 4.0):
https://creativecommons.org/licenses/by/4.0/

Changes made in this repository:

- `ReadSphModes_freqs.m`: extended from the original single-frequency
  reader to read several frequencies from one `.sph` file.
- `FFfromSph.m`, `ReadSphModes_freqs.m`, `SphModes_demo.m`: comments
  tidied (H1 line, usage example, a note that the demo's data files are
  not included). The code is otherwise unchanged.
- The `.sph` example files of the dataset are not included.
