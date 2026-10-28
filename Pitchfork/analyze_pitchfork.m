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
%     The least-squares problem is solved with MATLAB's backslash operator,
%     as in analyze_SH_scaled. For a rectangular design matrix this uses a
%     QR factorisation and warns if the matrix is rank deficient.
%
%   See also SIMULATE_PITCHFORK, ANALYZE_SH_SCALED.

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

%% Least-squares fit
theta = B \ q_tau;                 % theta = [mu_hat; s_hat]

%% Residuals and RMSE
prediction = B * theta;
R = q_tau - prediction;
rmse = sqrt(mean(R.^2));

end
