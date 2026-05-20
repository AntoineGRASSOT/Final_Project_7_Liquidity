function files = swaptionUtils()
%% SWAPTIONUTILS  Index of external utility functions in the MHW swaption library.
%
%  This file does not contain duplicated implementations. The operative
%  functions are standalone MATLAB files and can be called directly:
%
%    getDF.m                    - OIS discount factor interpolation
%    getPD.m                    - Euribor 6m pseudo-discount interpolation
%    cashAnnuity.m              - CS cash annuity, eq. (2.17)
%    swapRateATM.m              - multicurve ATM swap rate
%    mhwZeta.m                  - zeta_alpha, eq. (3.3)
%    mhwV.m                     - coefficient v_{alpha,i}, Lemma 3.1
%    hwB.m                      - coefficient B(a,tau)
%    hwVarianceFactor.m         - HW variance integral
%    fraConvexityAdjustment.m   - gamma_FRA, Appendix C
%
%  Usage:
%    files = swaptionUtils();

files = {'getDF', 'getPD', 'cashAnnuity', 'swapRateATM', ...
         'mhwZeta', 'mhwV', 'hwB', 'hwVarianceFactor', ...
         'fraConvexityAdjustment'};
end
