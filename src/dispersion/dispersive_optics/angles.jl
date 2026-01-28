brewster_angle_deg(n::Real) = rad2deg(atan(n))
brewster_angle_deg(n) = rad2deg.(atan.(n))

ideal_apex_deg(n::Real) = rad2deg(2asin(sin(deg2rad(brewster_angle_deg(n)))/n))
ideal_apex_deg(n) = @. rad2deg(2asin(sin(deg2rad(brewster_angle_deg(n)))/n))
