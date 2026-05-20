%% =========================================================================
%  swaption/ folder - MHW calibration on diagonal ATM CS swaptions
%  =========================================================================
%
%  REFERENCE: Baviera (2019) - "Back-of-the-envelope swaptions in a very
%  parsimonious multicurve interest rate model" - IJTAF 22(5), 1950027
%
%  VALUE DATE: 10 September 2015
%
%  =========================================================================
%  FILE STRUCTURE:
%
%  Pricing/calibration:
%  bachelierCS.m              Bachelier CS swaption formula, Appendix B
%  bachelierCS_implVol.m      Normal volatility inversion via bachelierCS
%  mhwPrice.m                 MHW ATM CS receiver swaption price, Prop. 3.2
%  calibrateMHW.m             MHW calibration, eq. (4.1), on 1y9y,...,9y1y
%  convexityAdjFRA.m          Optional FRA adjustment step, Appendix C
%
%  External callable utilities:
%  getDF.m                    OIS discount factor interpolation
%  getPD.m                    Euribor 6m pseudo-discount interpolation
%  cashAnnuity.m              Cash annuity C_{alpha,omega}(S), eq. (2.17)
%  swapRateATM.m              Multicurve ATM swap rate
%  mhwZeta.m                  zeta_alpha, eq. (3.3)
%  mhwV.m                     Coefficient v_{alpha,i}, Lemma 3.1
%  hwB.m                      HW coefficient B(a,tau)
%  hwVarianceFactor.m         HW variance factor
%  fraConvexityAdjustment.m   gamma_FRA, Appendix C
%  swaptionUtils.m            Utility index; no duplicated formulas
%
%  =========================================================================
%  PIPELINE (Section 4 of [B2019]):
%
%  Step 1 - Bootstrap curves (in bootstrap/)
%    OIS   = bootstrapOIS(mktData.OIS, params)
%    EUR6M = bootstrapEuribor6m(mktData.EUR6M, mktData.FRA, OIS, params)
%
%  Step 2 - MHW calibration (this folder)
%    [MHWparams, calibRes] = calibrateMHW(mktData.swaptVol, OIS, EUR6M, params)
%
%  Step 3 - Convexity adjustments (this folder)
%    EUR6M_adj = convexityAdjFRA(EUR6M, mktData.FRA, OIS, MHWparams, params)
%
%  =========================================================================
%  KEY MATHEMATICS (from [B2019]):
%
%  MHW model (eq. 3.1):
%    sigma(t,T) = (1-gamma) * sigma_tilde(t,T)   [OIS discount volatility]
%    eta(t,T)   = gamma     * sigma_tilde(t,T)   [spread volatility]
%    sigma_tilde(t,T) = sigma * (1-exp(-a(T-t))) / a   [pseudo-discount volatility]
%
%  Parameters:
%    a, sigma > 0, gamma in [0,1]
%    Calibrated values in the paper: a=12.94%, sigma=1.26%, gamma=0.07%
%
%  Key coefficients (Lemma 3.1, eq. 3.3-3.7):
%    zeta_alpha^2 = sigma^2(1-exp(-2a*t_alpha))/(2a)
%    v_{alpha,i} = zeta_alpha(1-exp(-a(t_i-t_alpha)))/a
%    xi_{alpha,i} = (1-gamma)v_{alpha,i}          [for B_{alpha,i}]
%    nu_{alpha,i} = v_{alpha,i}-gamma*v_{alpha,i+1} [for beta_i B_{alpha,i}]
%
%  Objective function (eq. 4.1):
%    Err^2(p) = sum_{i=1}^{M=9} [R^{C,MHW}_i(p) - R^{C,mkt}_i]^2
%
%  =========================================================================
