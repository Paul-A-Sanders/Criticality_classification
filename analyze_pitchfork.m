function [theta, R, rmse] = analyze_pitchfork(q_time, tau)
% analyze_pitchfork  Fit the reduced pitchfork drift from time-series data
%
%   dq/dtau = mu*q + s*q^3
%
% The returned coefficient vector is
%   theta = [mu_hat; s_hat].

% Ensure column vectors
q_time = q_time(:);
tau = tau(:);

N = numel(q_time);
if numel(tau) ~= N
    error('q_time and tau must have the same length.');
end
if N < 2
    error('Need at least two time samples.');
end

% Time differences (allow nonuniform dtau)
dtau = diff(tau);

% Time derivative (forward difference)
q_tau = diff(q_time) ./ dtau;

% Predictors evaluated at the left time point
q_left = q_time(1:end-1);
q3_left = q_left.^3;

% Design matrix [q, q^3]
B = [q_left, q3_left];

% Stable least-squares solve
theta = nan(2,1);
tol_rcond = 1e-12;

try
    [Q,Rqr] = qr(B,0);
    if min(abs(diag(Rqr))) > eps*max(abs(diag(Rqr)))*100
        theta = Rqr \ (Q' * q_tau);
    else
        rc = rcond(B' * B);
        if rc > tol_rcond
            theta = (B' * B) \ (B' * q_tau);
        else
            theta = pinv(B) * q_tau;
            warning('analyze_pitchfork:IllConditioned', ...
                'Design matrix nearly singular (rcond=%g). Used pinv fallback.', rc);
        end
    end
catch
    theta = pinv(B) * q_tau;
    warning('analyze_pitchfork:FallbackPINV', ...
        'QR solve failed; used pinv fallback.');
end

% Residuals and RMSE
prediction = B * theta;
R = q_tau - prediction;
rmse = sqrt(mean(R.^2));

end
