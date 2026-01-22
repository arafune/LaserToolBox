const BK7_SOURCE = """
# Source

  https://refractiveindex.info/?shelf=glass&book=BK7&page=SCHOTT

The values of the coefficients are taken from the SCHOTT catalog,
as listed on the Thorlabs web.
 """

const CaF2_SOURCE = """
# Source

  J. Phys. Chem. Ref. Data 9 161 (1980).

## Reference

  0.15 - 12 μm
"""

const FUSED_SILICA_SOURCE = """
# Source

  https://refractiveindex.info/?shelf=glass&book=fused_silica&page=Malitson

"""

const SF10_SOURCE = """
# Source

   https://refractiveindex.info/?shelf=popular_glass&book=SF10&page=SCHOTT

 Reference:
   0.15 - 12 µm
 """

const SF11_SOURCE = """
# Source

    https://refractiveindex.info/?shelf=specs&book=SCHOTT-optical&page=N-SF11

## Reference

  0.37 - 2.5 µm

# Another Source

     https://refractiveindex.info/?shelf=specs&book=SCHOTT-optical&page=SF11
     B = (1.73848403, 0.311168974, 1.17490871)
     C = (0.0136068604, 0.0615960463, 121.922711)
  """


const CORNING7980_SOURCE = """
# Source
  https://www.corning.com/media/worldwide/csm/documents/HPFS_Product_Brochure_All_Grades_2015_07_21.pdf

frequently refer as UVFS

"""

const CORNING7979_SOURCE = """
# Source
  https://www.corning.com/media/worldwide/csm/documents/HPFS_Product_Brochure_All_Grades_2015_07_21.pdf

frequently refer as FS for IR for IR

"""

# ---------------------------
"""
Sellmeier coefficients for BK7.

$BK7_SOURCE
"""
const BK7 = (
    A = 0.0,
    B = (1.03961212, 0.231792344, 1.01046945),
    C = (0.00600069867, 0.0200179144, 103.560653),
    range = (0.3, 2.5),
)


"""
Sellmeier coefficients for CaF2.

$CaF2_SOURCE
"""
const CaF2 = (
    A = 0.33973,
    B = (0.69913, 0.11994, 4.35181),
    C = (0.09374^2, 21.18^2, 38.46^2),
    range = (0.15, 12.0),
)

"""
Sellmeier coefficients for fused silica (UVFS)

$FUSED_SILICA_SOURCE
"""
const FusedSilica = (
    A = 0.0,
    B = (0.6961663, 0.4079426, 0.8974794),
    C = (0.06840432^2, 0.11624142^2, 9.8961612^2),
    range = (0.21, 6.7),
)

"""
Sellmeier coefficients for SF10.

$SF10_SOURCE
"""
const SF10 = (
    A = 0.0,
    B = (1.62153902, 0.256287842, 1.64447552),
    C = (0.0122241457, 0.0595736775, 147.468793),
    range = (0.38, 2.5),
)

"""
Sellmeier coefficients for SF11.

$SF11_SOURCE
"""
const SF11 = (
    A = 0.0,
    B = (1.73759695, 0.313747346, 1.89878101),
    C = (0.013188707, 0.0623068142, 155.23629),
    range = (0.37, 2.5),
)



"""
Sellmeier coefficients for Corning 7979 Fused Silica (UVFS)

$CORNING7980_SOURCE
"""
const Corning7980 = (
    A = 0.0,
    B = (0.68374049400, 0.42032361300, 0.58502748000),
    C = (0.00460352869, 0.01339688560, 64.49327320000),
    range = (0.185, 1.1),
)

"""
Sellmeier coefficients for Corning 7979 Fused Silica (FS for IR) 

$CORNING7979_SOURCE
"""
const Corning7979 = (
    A = 0.0,
    B = (3.550277875E-02, 7.306029048E-01, 3.475321572E-01, 9.216052441E-01),
    C = (-5.783959035E-03, 5.600103210E-03, 1.389808930E-02, 1.006578079E+02),
    range = (0.185, 2.3),
)

