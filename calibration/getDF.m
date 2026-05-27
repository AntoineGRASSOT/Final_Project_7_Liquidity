function DF = getDF(OIS, T_query)
%% GETDF  Interpolate OIS discount factors.
%
%  Linear interpolation on zero rates z(t) = -ln(DF)/t  (ACT/365).
%  Flat extrapolation outside the pillar range; DF = 1 for t <= 0.
%
%  Inputs:
%    OIS     - from bootstrap:
%                OIS.T   year fractions from settlement (ACT/365)  
%                OIS.DF  OIS discount factors at pillars           
%    T_query - year fractions (ACT/365) at which to evaluate DF  
%
%  Output:
%    DF      - interpolated discount factors, same shape as T_query

%% --- pillars ----------------------------------------------------------
t = OIS.T(:);
z = -log(OIS.DF(:)) ./ t;          % zero rates at pillars

%% -----------------------------
sz      = size(T_query);
T_query = T_query(:);               % flatten to column

%% --- vectorised interpolation with flat extrapolation -----------------
%  Clamping T_query to [t(1), t(end)] before interp1 achieves flat
%  extrapolation on both sides in a single call (no NaN, no branches).
T_clamp = max(t(1), min(t(end), T_query));
z_query = interp1(t, z, T_clamp, 'linear');

DF = exp(-z_query .* T_query);

%% --- override: t <= 0  =>  DF = 1  ------------------------------------
DF(T_query <= 0) = 1.0;

%% --- restore original shape -------------------------------------------
DF = reshape(DF, sz);
end
