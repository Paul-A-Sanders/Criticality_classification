function [tau, X, q] = simulate_scaled_SH_Fourier(mu, s, gamma, tmax, tsteps, L, xsteps, ...
    noiseType, correlation_length, NoiseScale, InitialScale)
% SIMULATE_SCALED_SH_FOURIER  Stochastic Swift-Hohenberg equation, FFT spectral solver.
%
%   Simulates the reduced stochastic Swift-Hohenberg equation
%
%       dq = [mu*q - gamma*(1+d_xx)^2*q + s*q^3] dtau + NoiseScale*dW
%
%   on a one-dimensional periodic domain, using Euler-Maruyama in time and
%   a Fourier spectral discretisation in space. The linear part is
%   advanced explicitly in Fourier space, the cubic term is evaluated in
%   physical space.
%
%   Syntax
%       [tau, X, q] = simulate_scaled_SH_Fourier(mu, s, gamma, tmax, tsteps, L, xsteps, ...
%                         noiseType, correlation_length, NoiseScale, InitialScale)
%
%   Inputs
%       mu, gamma           Linear coefficients (see equation above).
%       s                   Cubic coefficient / sign (+1 or -1).
%       tmax                Final scaled time tau_max.
%       tsteps              Number of stored time points (including tau = 0
%                           and tau = tmax); must be >= 2.
%       L                   Length of the periodic domain.
%       xsteps              Number of unique grid points Nx on [0, L); an
%                           integer >= 4.
%       noiseType           'White'      space-time white noise,
%                           'Correlated' white noise smoothed with a Gaussian
%                                        filter (needs correlation_length > 0),
%                           'Uniform'    unit-variance uniform noise,
%                           'None'       deterministic run.
%                           Case-insensitive; char vector or string scalar.
%       correlation_length  Width of the Gaussian smoothing filter; only
%                           used for 'Correlated'.
%       NoiseScale          Noise amplitude.
%       InitialScale        Standard deviation of the random initial
%                           condition, q(x,0) = InitialScale*randn.
%
%   Outputs
%       tau   tsteps-by-1 vector of times.
%       X     1-by-Nx vector of positions, X_j = j*dX, j = 0..Nx-1, dX = L/Nx.
%             X = L is not stored because it is identical to X = 0.
%       q     tsteps-by-Nx array; row i is the field at time tau(i).
%
%   Notes
%     * Stability: explicit Euler requires dtau*|lambda_k| < 2 for every
%       damped mode. The function raises an error if this is violated and a
%       warning when dtau*max|lambda_k| >= 1.5.
%     * Random numbers come from MATLAB's global stream; seed it in the
%       caller (e.g. rng(seed,'combRecursive')) for reproducible runs.
%
%   See also ANALYZE_SH_SCALED.

if tsteps < 2
    error('tsteps must be at least 2.');
end
if xsteps < 4 || xsteps ~= round(xsteps)
    error('xsteps must be an integer of at least 4.');
end

%% Check the noise type before doing any work
noiseMode = lower(char(noiseType));
switch noiseMode
    case {'white','uniform','none'}
        % Nothing to prepare.
    case 'correlated'
        if correlation_length <= 0
            error('noiseType ''Correlated'' requires correlation_length > 0.');
        end
    otherwise
        error('Unknown noiseType: %s', noiseType);
end

%% Time and spatial grids
dtau = tmax/(tsteps-1);
tau  = linspace(0,tmax,tsteps)';

Nx = xsteps;
dX = L/Nx;
X  = (0:Nx-1)*dX;

%% Fourier wavenumbers
% Define the physical wavenumbers in centred order, then rearrange them
% into MATLAB's native (unshifted) FFT storage order. This keeps fft(q)
% unshifted during time stepping while making the convention explicit.
if mod(Nx,2) == 0
    k_centered = (2*pi/L)*(-Nx/2:Nx/2-1)';
else
    k_centered = (2*pi/L)*(-(Nx-1)/2:(Nx-1)/2)';
end
k_fft = ifftshift(k_centered);

%% Linear Swift-Hohenberg operator
% d_xx -> -k^2, hence  lambda(k) = mu - gamma*(1-k^2)^2.
Lk = mu - gamma*(1-k_fft.^2).^2;

% Explicit-Euler stability check for the damped linear modes: for
% lambda < 0, forward Euler requires dtau*|lambda| < 2.
dampedRates = -Lk(Lk < 0);
if ~isempty(dampedRates)
    stabilityNumber = dtau*max(dampedRates);
    if stabilityNumber >= 2
        error(['Explicit Euler time step is linearly unstable: ', ...
            'dtau*max(-lambda_k) = %.6g >= 2.'], stabilityNumber);
    elseif stabilityNumber >= 1.5
        warning(['Explicit Euler time step is close to the linear stability limit: ', ...
            'dtau*max(-lambda_k) = %.6g.'], stabilityNumber);
    end
end

%% Initial condition
q = zeros(tsteps,Nx);
q(1,:) = InitialScale*randn(1,Nx);

%% Noise
% sqrt(dtau/dX) is the grid scaling for unit-amplitude space-time white noise.
epsilon = NoiseScale*sqrt(dtau/dX);

% Gaussian smoothing filter for correlated noise (applied via FFT).
if strcmp(noiseMode,'correlated')
    dist = min((0:Nx-1)'*dX, L-(0:Nx-1)'*dX);     % Periodic distance to x = 0
    filt = exp(-(dist.^2)/(correlation_length^2));
    filt = filt/sqrt(sum(filt.^2));               % L2 normalisation keeps the variance comparable
    filt_hat = fft(filt);
end

%% Euler-Maruyama time stepping
for ti = 1:tsteps-1
    q_phys = q(ti,:)';

    % Nonlinear term (physical space)
    NL = s*q_phys.^3;

    % Deterministic step (Fourier space)
    q_hat  = fft(q_phys);
    NL_hat = fft(NL);

    q_hat_next = q_hat + dtau*(Lk.*q_hat + NL_hat);
    q_next = real(ifft(q_hat_next));

    % Stochastic forcing
    switch noiseMode
        case 'white'
            Noise = epsilon*randn(Nx,1);

        case 'correlated'
            eta   = epsilon*randn(Nx,1);
            Noise = real(ifft(fft(eta).*filt_hat));

        case 'uniform'
            % Unit-variance uniform noise before multiplication by epsilon.
            Noise = epsilon*sqrt(3)*(2*rand(Nx,1)-1);

        case 'none'
            Noise = zeros(Nx,1);
    end

    q(ti+1,:) = (q_next + Noise)';
end

end
