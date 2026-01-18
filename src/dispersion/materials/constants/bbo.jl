const ALPHA_BBO_SOURCE = """
# Source

  These coefficients are taken from Applied Optics, 41(13), 2474.

# Reference

  @587.6
    n_o = 1.533
    n_e = 1.673

# Another Source

   The following coefficients are taken from
   http://www.newlightphotonics.com/Birefringent-Crystals/alpha-BBO-Crystals
 
   * coeffs_o = (2.67579, 0.02099, 0.00470, 0.00528)
   * coeffs_e = (2.31197, 0.01184, 0.016070, 0.00400)
"""

"""
Sellmeier coefficients for α-BBO, biaxal crystal

$ALPHA_BBO_SOURCE
"""
const AlphaBBO = (
    coeffs_o = (2.7471, 0.01878, 0.01822, 0.01354),
    coeffs_e = (2.3174, 0.01224, 0.01667, 0.001516),
    range = (0.19, 3.5),
)


const BETA_BBO_SOURCE = """
# Source
  
https://refractiveindex.info/?shelf=main&book=BaB2O4&page=Tamosauskas-e
https://refractiveindex.info/?shelf=main&book=BaB2O4&page=Tamosauskas-o

1. G. Tamošauskas, G. Beresnevičius, D. Gadonas, A. Dubietis. Transmittance and phase matching of BBO crystal in the 3−5 μm range and its application for the characterization of mid-infrared laser pulses. Opt. Mater. Express 8, 1410-1418 (2018)
2. G. Tamošauskas. β-barium borate (BBO) absorption in the 0.188-6.22 μm range. arXiv:2111.01212 [physics.optics] (2021)


# Other Souce (Not used here):

  K. Kato, IEEE J. Quantum Electron QE-22, 1013 (1986)

  Model equation is different from the conventional Sellmeier equation.

 * coeffs_o = (2.7359, 0.01878, 0.01822, 0.01354)
 * coeffs_e = (2.3753, 0.01224, 0.01667, 0.01516)

  D. Eimerl et al., J. Appl. 62, 1968 (1987).

  Model equation is different from the conventional Sellmeier equation.

 * coeffs_o = (2.7405, 0.0184, 0.0179, 0.0155)
 * coeffs_e = (2.3730, 0.0128, 0.0156, 0.0044)
"""

"""
Sellmeier coefficients for β-BBO

$BETA_BBO_SOURCE
"""
const BetaBBO = (
    e = (A = 0.0, B = (1.151075, 0.21803, 0.656), C = (0.007142, 0.02259, 263.0)),
    o = (A = 0.0, B = (0.90291, 0.83155, 0.76536), C = (0.003926, 0.018786, 60.01)),
    range = (0.188, 5.2),
)

