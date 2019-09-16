Demo1 runs an example of a spherical lens. Two dataset are available, lensConv loads the data of a lens under a thermal gradient on a three-point mount, lensConv4 loads the data of a lens under a thermal gradient on a three-point mount.

Three different illumination conditions are proposed, the user can change the source location and number of rays, it is recommended to keep the illumination aperture since surfaces identification is particularly sensitive to that parameter. The user should also ensure that rays will not reach the lens edge.

Demo2 runs an example of a square optical window mounted on a three point supported, in which one introduce a punctual load. WindowExpBK7 loads the dataset for a window made of BK7 glass, and WindowExpFS loads the data for a window of fused silica both subjected to the same load. Stress optical coefficient are set for a bk7 window, for a FS window k11 = -0.7 e-12 and k11 = -4.2 e-12. Other glasses stress optical coefficients are in “Handbook of optical material. Weber M. CRC Press, 2003“.
If verbosity is 1, the script plots the algorithm outputs.

Data is saved in the output folder, so it can be loaded for further analysis.
