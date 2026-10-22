function [tau, q] = simulate_pitchfork(mu, s, tmax, tsteps, ...
    noiseType, NoiseScale, InitialScale)
% SIMULATE_PITCHFORK  Euler-Maruyama simulation of the stochastic pitchfork normal form.
%
%   Model (in scaled time tau):
%
%       dq = (mu*q + s*q^3) dtau + NoiseScale*dW
%
%   Syntax
%       [tau, q] = simulate_pitchfork(mu, s, tmax, tsteps, noiseType, NoiseScale, InitialScale)
%
%   Inputs
%       mu           Linear coefficient (mu < 0: q = 0 is linearly stable).
%       s            Cubic coefficient / sign: +1 subcritical, -1 supercritical.
%       tmax         Final scaled time tau_max.
%       tsteps       Number of stored time points, including tau = 0 and tau = tmax.
%       noiseType    'White'   Gaussian increments,
%                    'Uniform' unit-variance uniform increments,
%                    'None'    deterministic run.
%                    Case-insensitive; char vector or string scalar.
%       NoiseScale   Noise amplitude.
%       InitialScale Standard deviation of the random initial condition,
%                    q(0) = InitialScale*randn.
%
%   Outputs
%       tau          tsteps-by-1 vector of times, linspace(0, tmax, tsteps)'.
%       q            tsteps-by-1 vector with the simulated trajectory.
%
%   Notes
%     * The time step is dtau = tmax/(tsteps-1). To obtain a desired step
%       dt_target exactly, the calling script should use
%           tsteps = round(tmax/dt_target) + 1.
%     * Random numbers come from MATLAB's global stream, so seed it in the
%       caller (e.g. rng(seed,'combRecursive')) for reproducible runs.
%     * randn() is always called for the initial condition, even when
%       InitialScale = 0, so the noise sequence for a given seed does not
%       depend on InitialScale.
%
%   See also ANALYZE_PITCHFORK.

if nargin < 7
    error('All 7 arguments must be provided.');
end
if tsteps < 2
    error('tsteps must be at least 2.');
end

noiseMode = lower(char(noiseType));
if ~ismember(noiseMode, {'white','uniform','none'})
    error('Unknown noiseType: %s', noiseType);
end

%% Time grid
dtau = tmax/(tsteps-1);
tau  = linspace(0,tmax,tsteps)';

%% Initial condition
q = zeros(tsteps,1);
q(1) = InitialScale*randn();

%% Euler-Maruyama time stepping
% q_{i+1} = q_i + dtau*(mu*q_i + s*q_i^3) + NoiseScale*sqrt(dtau)*xi_i
sigma_dtau = NoiseScale*sqrt(dtau);    % Prefactor of the Gaussian increment

for ti = 1:tsteps-1
    qi = q(ti);

    % Deterministic drift
    drift = mu*qi + s*qi^3;

    % Stochastic increment
    switch noiseMode
        case 'white'
            dW = sigma_dtau*randn();
        case 'uniform'
            % Uniform on [-sqrt(3), sqrt(3)] has unit variance.
            dW = NoiseScale*sqrt(3)*(2*rand()-1)*sqrt(dtau);
        case 'none'
            dW = 0;
    end

    q(ti+1) = qi + dtau*drift + dW;
end

end
