function PD = getPD(EUR6M, T_query)
%% GETPD  Interpolate Euribor 6m pseudo-discount factors.
%
%  Linear interpolation on pseudo-zero rates q(t) = -ln(PD)/t  (ACT/365).
%  Flat extrapolation outside the pillar range; PD = 1 for t <= 0.
%
%  Inputs:
%    EUR6M   - from bootstrap:
%                EUR6M.T   year fractions from settlement (ACT/365)  
%                EUR6M.PD  pseudo-discount factors at pillars         
%    T_query - year fractions (ACT/365) at which to evaluate PD    
%
%  Output:
%    PD      - interpolated pseudo-discount factors, same shape as T_query

%% --- pillars ----------------------------------------------------------
t = EUR6M.T(:);
q = -log(EUR6M.PD(:)) ./ t;        % pseudo-zero rates at pillars

%% -----------------------------
sz      = size(T_query);
T_query = T_query(:);               % flatten to column

%% --- interpolation with flat extrapolation -----------------
%  Clamping T_query to [t(1), t(end)] before interp1 achieves flat
%  extrapolation on both sides.
T_clamp = max(t(1), min(t(end), T_query));
q_query = interp1(t, q, T_clamp, 'linear');

PD = exp(-q_query .* T_query);

%% --- override: t <= 0  =>  PD = 1  ------------------------------------
PD(T_query <= 0) = 1.0;

%% --- restore original shape -------------------------------------------
PD = reshape(PD, sz);
end
