function EUR6M_adj = convexityAdjFRA(EUR6M_raw, FRAdata, OIS, MHWparams, params)
%% CONVEXITYADJFRA  FRA convexity adjustments, Appendix C in Baviera (2019).
%
%  Step 3 of the calibration cascade: after calibrating the MHW parameters,
%  the Euribor 6m curve buckets affected by FRA convexity adjustments may be
%  rebuilt. In the paper's numerical case these adjustments are negligible,
%  but the implementation follows relations C.3-C.4.
%
%  Appendix C:
%    F_FRA(t0) = gamma_FRA * L(t0;t_i,t_{i+1}) + (gamma_FRA - 1)/delta
%    gamma_FRA = exp(-int_0^{t_i} sigma_tilde_i(t) * eta_i(t) dt)
%
%  In the one-factor MHW model:
%    eta_i(t) = gamma * sigma_tilde_i(t)
%  hence the FRA adjustment depends on a, sigma and gamma.
%
%  Inputs:
%    EUR6M_raw - pseudo-discount curve before the adjustment
%    FRAdata   - structure with startT, endT, mid and optionally startMon
%    OIS       - unused here; kept for pipeline compatibility
%    MHWparams - structure with a, sigma, gamma
%    params    - unused here; kept for pipeline compatibility
%
%  Output:
%    EUR6M_adj - pseudo-discount curve updated where the correction is
%                numerically significant

%#ok<*INUSD>
a     = MHWparams.a;
sigma = MHWparams.sigma;
gamma = MHWparams.gamma;

fprintf('>>> FRA convexity adjustments (Appendix C, Baviera 2019):\n');

EUR6M_adj = EUR6M_raw;
nFRA = length(FRAdata.mid);

for k = 1:nFRA
    t_start = FRAdata.startT(k);
    t_end   = FRAdata.endT(k);
    delta   = t_end - t_start;
    F_FRA   = FRAdata.mid(k);

    [gammaFRA, varianceIntegral] = fraConvexityAdjustment(a, sigma, gamma, t_start, t_end);

    % Exact solution of C.3 with respect to the forward Libor/pseudo-discount rate.
    F_fwd = (F_FRA - (gammaFRA - 1) / delta) / gammaFRA;
    correction_bps = (F_fwd - F_FRA) * 1e4;

    if isfield(FRAdata, 'startMon')
        labelStart = round(FRAdata.startMon(k));
    else
        labelStart = round(12 * t_start);
    end

    fprintf('   FRA %dx%d: gammaFRA=%.12f, int=%.4e, deltaF=%.6f bps  ', ...
            labelStart, labelStart + 6, gammaFRA, varianceIntegral, correction_bps);

    if abs(correction_bps) < 0.01
        fprintf('[NEGLIGIBLE]\n');
        continue;
    end

    fprintf('[APPLIED]\n');

    if ~(isfield(EUR6M_adj, 'T_knots') && isfield(EUR6M_adj, 'PD_knots'))
        warning('convexityAdjFRA:missingKnots', ...
                'Curve has no T_knots/PD_knots fields: FRA bucket update skipped.');
        continue;
    end

    [~, idx_ts] = min(abs(EUR6M_adj.T_knots - t_start));
    [~, idx_te] = min(abs(EUR6M_adj.T_knots - t_end));

    if idx_ts == idx_te
        warning('convexityAdjFRA:sameKnot', ...
                'FRA start/end collapse onto the same knot; update skipped.');
        continue;
    end

    PD_ts = EUR6M_adj.PD_knots(idx_ts);
    delta_eff = EUR6M_adj.T_knots(idx_te) - EUR6M_adj.T_knots(idx_ts);
    PD_te_new = PD_ts / (1 + delta_eff * F_fwd);

    delta_PD = abs(PD_te_new - EUR6M_adj.PD_knots(idx_te));
    if delta_PD > 1e-10
        EUR6M_adj.PD_knots(idx_te) = PD_te_new;
        if isfield(EUR6M_adj, 'ZR_knots')
            EUR6M_adj.ZR_knots(idx_te) = -log(PD_te_new) / EUR6M_adj.T_knots(idx_te);
        end
        fprintf('      PD knot %.4fy adjusted by %.3e\n', EUR6M_adj.T_knots(idx_te), delta_PD);
    end
end

if isfield(EUR6M_adj, 'T_knots') && isfield(EUR6M_adj, 'ZR_knots') && ...
        isfield(EUR6M_adj, 'T') && ~isempty(EUR6M_adj.ZR_knots)
    EUR6M_adj.ZR = interp1(EUR6M_adj.T_knots, EUR6M_adj.ZR_knots, ...
                           EUR6M_adj.T, 'pchip', 'extrap');
    EUR6M_adj.PD = exp(-EUR6M_adj.ZR .* EUR6M_adj.T);
end

fprintf('   Convexity adjustments completed. For 10 Sep 2015 they are negligible.\n\n');
end
