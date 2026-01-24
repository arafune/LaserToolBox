using Test
using LaserToolBox

const ri = LaserToolBox.n

@test beta_n(ri.sf11, 0.8, order = 2) ≈ 187.50 atol=1e-1
