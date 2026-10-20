function [theta, R, rmse] = analyze_pitchfork(q_time, tau)
% ANALYZE_PITCHFORK  Estimate the pitchfork coefficients from a time series.
%
%   Fits the reduced normal form
%
%       dq/dtau = mu*q + s*q^3
%
%   by ordinary least squares: the finite-difference time derivative of the
%   data is regressed on the candidate terms [q, q^3].
%
%   Syntax
%       [theta, R, rmse] = analyze_pitchfork(q_time, tau)
%
%   Inputs
%       q_time   N-element vector, the observed trajectory q(tau_j).
%       tau      N-element vector of sample times (need not be uniform).
%
%   Outputs
%       theta    2-by-1 vector [mu_hat; s_hat] of estimated coefficients.
%       R        (N-1)-by-1 residuals of the regression.
%       rmse     Root-mean-square of R.
%
%   Method
%     The regression uses forward differences, matching the Euler-Maruyama
%     scheme in simulate_pitchfork:
%         (q_{j+1} - q_j)/(tau_{j+1} - tau_j)  ~  mu*q_j + s*q_j^3
%     The least-squares problem is solved with a cascade of increasingly
%     robust solvers: QR decomposition, then the normal equations, then the
%     pseudo-inverse (with a warning) if the design matrix is nearly singular.
%
%   See also SIMULATE_PITCHFORK.

%% Input checks and shaping
q_time = q_time(:);
tau    = tau(:);

N = numel(q_time);
if numel(tau) ~= N
    error('q_time and tau must have the same length.');
end
if N < 2
    error('Need at least two time samples.');
end

%% Regression problem: q_tau = B*theta
dtau  = diff(tau);                 % Time steps (may be non-uniform)
q_tau = diff(q_time) ./ dtau;      % Forward-difference derivative

q_left  = q_time(1:end-1);         % Predictors are evaluated at the left end of each step
q3_left = q_left.^3;

B = [q_left, q3_left];             % Design matrix [q, q^3]

%% Least-squares solve
tol_rcond = 1e-12;                 % Minimum acceptable rcond for the normal equations

try
    [Q,Rqr] = qr(B,0);             % Economy-size QR
    if min(abs(diag(Rqr))) > eps*max(abs(diag(Rqr)))*100
        % Well-conditioned: solve via QR
        theta = Rqr \ (Q' * q_tau);
    else
        rc = rcond(B' * B);
        if rc > tol_rcond
            % Moderately conditioned: normal equations
            theta = (B' * B) \ (B' * q_tau);
        else
            % Nearly singular: minimum-norm solution via pseudo-inverse
            theta = pinv(B) * q_tau;
            warning('analyze_pitchfork:IllConditioned', ...
                'Design matrix nearly singular (rcond=%g). Used pinv fallback.', rc);
        end
    end
catch
    % Any failure above (e.g. non-finite data): fall back to pinv
    theta = pinv(B) * q_tau;
    warning('analyze_pitchfork:FallbackPINV', ...
        'QR solve failed; used pinv fallback.');
end

%% Residuals and RMSE
prediction = B * theta;
R = q_tau - prediction;
rmse = sqrt(mean(R.^2));

end
