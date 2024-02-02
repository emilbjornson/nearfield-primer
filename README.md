[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=emilbjornson/nearfield-primer)

A Primer on Near-Field Beamforming for Arrays and Reconfigurable Intelligent Surfaces
==================

This is a code package is related to the following scientific article:

Emil Björnson, Özlem Tuğfe Demir, and Luca Sanguinetti “[A Primer on Near-Field Beamforming for Arrays and Reconfigurable Intelligent Surfaces](https://arxiv.org/pdf/2110.06661.pdf
),” Asilomar Conference on Signals, Systems, and Computers, Virtual conference, October-November 2021.

The package contains a simulation environment, based on Matlab, that reproduces some of the numerical results and figures in the article. We encourage you to also perform reproducible research!


## Abstract of Article

Wireless communication systems have almost exclusively operated in the far-field of antennas and antenna arrays,
which is conventionally characterized by having propagation
distances beyond the Fraunhofer distance. This is natural since
the Fraunhofer distance is normally only a few wavelengths.
With the advent of active arrays and passive reconfigurable
intelligent surfaces (RIS) that are physically large, it is plausible
that the transmitter or receiver is located in between the
Fraunhofer distance of the individual array/surface elements and
the Fraunhofer distance of the entire array. An RIS then can be
configured to reflect the incident waveform towards a point in
the radiative near-field of the surface, resulting in a beam with
finite depth, or as a conventional angular beam with infinity
focus, which only results in amplification in the far-field. To
understand when these different options are viable, an accurate
characterization of the near-field behaviors is necessary. In this
paper, we revisit the motivation and approximations behind the
Fraunhofer distance and show that it is not the right metric for
determining when near-field focusing is possible. We obtain the
distance range where finite-depth beamforming is possible and
the distance where the beamforming gain tapers off.

## Content of Code Package

The article contains 6 simulation figures, numbered 3-5 and 7-9. Figures 3, 4, 5, 7, 8, and 9 are generated respectively by the Matlab scripts simulateFigure3.m, simulateFigure4.m, simulateFigure5.m, simulateFigure7.m, simulateFigure8.m, and simulateFigure9.m.

See each file for further documentation.

## Acknowledgements

This work was supported by the FFL18-0277 grant from the SSF.

## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
