function [gammaFRA, varianceIntegral] = fraConvexityAdjustment(a, sigma, gamma, t_start, t_end)
%% FRACONVEXITYADJUSTMENT  MHW FRA convexity adjustment, Appendix C.
%
%  Appendix C defines:
%    gamma_FRA = exp(-int_0^{t_i} sigma_tilde_i(t) * eta_i(t) dt).
%
%  In the one-factor MHW model:
%    eta_i(t) = gamma * sigma_tilde_i(t),
%  hence:
%    gamma_FRA = exp(-gamma * int_0^{t_i} sigma_tilde_i(t)^2 dt).

delta = t_end - t_start;
varianceIntegral = sigma^2 * hwB(a, delta)^2 * hwVarianceFactor(a, t_start);
gammaFRA = exp(-gamma * varianceIntegral);
end
