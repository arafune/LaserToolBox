const MgF2_SOURCE = """
# Source

  Appl. Opt. 23, 1980-1985 (1984).

## Reference value
  @587.6nm:
    n_e = 1.390
    n_o = 1.378
"""

const CALCITE_SOURCE = """
# Source

  Opt. Commun. 163, 95-102 (1999).

## Reference value
  @1.064:
    n_e= 1.480
    n_o= 1.642

* Negative birefringence
"""


const QUARTZ_SOURCE = """
# Source

    Opt. Commun. 163, 95-102 (1999).
    https://refractiveindex.info/?shelf=main&book=SiO2&page=Ghosh-e
 
"""

"""
Sellmeier coefficients for MgF₂ (biaxial crystal).

$MgF2_SOURCE
"""
const MgF2 = (
    e = (
        A = 0.0,
        B = (0.41344023, 0.50497499, 2.4904862),
        C = (0.03684262^2, 0.09076162^2, 23.771995^2),
    ),
    o = (
        A = 0.0,
        B = (0.48755108, 0.39875031, 2.3120353),
        C = (0.04338408^2, 0.09461442^2, 23.793604^2),
    ),
    range = (0.2, 7.0),
)

"""
Sellmeier coefficients for Calcite (CaCO3).

$CALCITE_SOURCE
"""
const Calcite = (
    e = (A = 0.35859695, B = (0.82427380, 0.14429128), C = (1.06689543 * 1e-2, 120.0)),
    o = (A = 0.73358749, B = (0.96464345, 1.82831454), C = (1.94325203 * 1e-2, 120.0)),
    range = (0.204, 2.172),
)


"""
Sellmeier coefficients for crystal quartz.

$QUARTZ_SOURCE
"""
const Quartz = (
    e = (A = 0.28851804, B = (1.09509924, 1.15662475), C = (1.02101864 * 1e-2, 100.0)),
    o = (A = 0.28604141, B = (1.07044083, 1.10202242), C = (1.00585997 * 1e-2, 100.0)),
    range = (0.198, 2.0531),
)

