using Optim

function lmfit{I <: Integer}(f::Function, p0::Vector{Float64}, maxIter::I)
	# this is a convenience function for the curve_fit() methods
	# which assume f(p) is the cost functionj i.e. the residual of a
	# model where
	#   model(xpts, params...) = ydata + error (noise)

	# this minimizes f(p) using a least squares sum of squared error:
	#   sse = sum(f(p)^2)
	# This is currently embedded in Optim.levelberg_marquardt()
	# which calls sumabs2
	#
	# returns p, f(p), g(p) where
	#   p    : best fit parameters
	#   f(p) : function evaluated at best fit p, (weighted) residuals
	#   g(p) : estimated Jacobian at p (Jacobian with respect to p)

	# construct Jacobian function, which uses finite difference method
	g = Calculus.jacobian(f)

	results = Optim.levenberg_marquardt(f, g, p0; maxIter=maxIter, tolG=1e-8, lambda=1.0)
	# p = results.minimum
	# resid = f(p)
	# dof = length(resid) - length(p)

  return results
end
