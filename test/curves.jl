using Dates

effective_date = Date(2020, 01, 15)

days = [369, 462, 553, 735, 1100, 1465]
df = [0.97678403, 0.97153083, 0.96655364, 0.95695995, 0.93829085, 0.91954718]
test_date = effective_date + Day(900)

const dc_act_360 = QuantLib.Time.Actual360()

# test log-linear interpolation
test_curve_loglinear = QuantLib.InterpolatedDiscountCurve(effective_date .+ Day.(days), df, dc_act_360, QuantLib.Math.LogLinear())
disc_factor = QuantLib.discount(test_curve_loglinear, test_date)
x1 = (test_date - test_curve_loglinear.dates[4]).value
dx = (test_curve_loglinear.dates[5]-test_curve_loglinear.dates[4]).value
dy = log(test_curve_loglinear.data[5]) - log(test_curve_loglinear.data[4])
df_recalc = exp(log(test_curve_loglinear.data[4]) + x1*dy/dx)
@test disc_factor == df_recalc

test_curve_linear = QuantLib.InterpolatedDiscountCurve(effective_date .+ Day.(days), df, dc_act_360, QuantLib.Math.LinearInterpolation())
disc_factor = QuantLib.discount(test_curve_linear, test_date)
x1 = (test_date - test_curve_linear.dates[4]).value
dx = (test_curve_linear.dates[5]-test_curve_linear.dates[4]).value
dy = (test_curve_linear.data[5]) - (test_curve_linear.data[4])
df_recalc = ((test_curve_linear.data[4]) + x1*dy/dx)
@test disc_factor == df_recalc
