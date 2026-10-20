function [theta, k_max, max_eigenvalue] = analyze_SH_scaled(q, tau, X)
% ANALYZE_SH_SCALED  Estimate the Swift-Hohenberg coefficients from space-time data.
%
%   Fits the reduced spatial model
%
%       q_tau = a0*q + a2*q_XX + a4*q_XXXX + s*q^3
%
%   by least squares over all spatial locations and time steps, and
%   derives the dominant mode of the fitted linear operator.
%
%   Syntax
%       [theta, k_max, max_eigenvalue] = analyze_SH_scaled(q, tau, X)
%
%   Inputs
%       q     Nt-by-Nx array of the field, one row per time, one column per
%             spatial point (periodic domain, uniform spacing).
%       tau   Nt-element vector of times (may be non-uniform).
%       X     Nx-element vector of spatial positions (uniform spacing).
%
%   Outputs
%       theta           [a0_hat; a2_hat; a4_hat; s_hat] (4-by-1).
%       k_max           Wavenumber k >= 0 at which the estimated growth rate
%                       lambda(k) is largest.
%       max_eigenvalue  lambda(k_max), the largest estimated growth rate.
%
%   Dispersion relation
%     In Fourier space the fitted linear operator is
%
%         lambda(k) = a0 - a2*k^2 + a4*k^4.
%
%     For the scaled Swift-Hohenberg equation
%         dq = [mu*q - gamma*(1+d_xx)^2*q + s*q^3] dtau + dW
%     the true coefficients are a0 = mu - gamma, a2 = -2*gamma, a4 = -gamma,
%     so lambda(k) = mu - gamma*(1-k^2)^2 with its maximum lambda = mu at k = 1.
%
%   Notes
%     * Spatial derivatives are evaluated spectrally (FFT) on the periodic
%       domain; the time derivative is a forward difference.
%     * The data are centred by subtracting the global mean before fitting.
%     * k_max is located with fminbnd on [0, max|k|], which returns a local
%       maximiser. This is the global maximum when a4_hat < 0 (as for the
%       true model); otherwise lambda(k) is largest at an end point of the
%       interval and the result should be interpreted with care.
%     * The dispersion-function handle is only needed here and is not returned.
%
%   See also SIMULATE_SCALED_SH_FOURIER.

%% Data and grids
tau = tau(:);
X   = X(:);

[Nt,Nx] = size(q);

dtau_vec = diff(tau);        % Time steps (Nt-1 values)

dX = mean(diff(X));          % Uniform spatial spacing
L  = Nx*dX;                  % Periodic domain length

%% Fourier wavenumbers
% Define the physical wavenumbers in centred order, then rearrange them
% into MATLAB's native (unshifted) FFT storage order.
if mod(Nx,2) == 0
    k_centered = (2*pi/L)*(-Nx/2:Nx/2-1)';
else
    k_centered = (2*pi/L)*(-(Nx-1)/2:(Nx-1)/2)';
end
k_fft = ifftshift(k_centered);

%% Centre the data
Q = q' - mean(q(:));         % Nx-by-Nt, one column per time

%% Regression matrices
% A holds the time derivative (response); the four slices of B hold the
% regressors q, q_XX, q_XXXX and q^3, each evaluated at the left end of
% every time step.
A = zeros(Nx,Nt-1);
B = zeros(Nx,Nt-1,4);

for j = 1:Nt-1

    % Time derivative (forward difference)
    A(:,j) = (Q(:,j+1)-Q(:,j))/dtau_vec(j);

    % q
    B(:,j,1) = Q(:,j);

    % Spatial derivatives (spectral): d_xx -> -k^2, d_xxxx -> +k^4
    Qhat = fft(Q(:,j));
    B(:,j,2) = real(ifft(-k_fft.^2 .* Qhat));
    B(:,j,3) = real(ifft( k_fft.^4 .* Qhat));

    % Cubic term
    B(:,j,4) = Q(:,j).^3;

end

%% Least-squares fit over all space-time points
A = A(:);
B = reshape(B,[],4);

theta = B\A;

%% Estimated dispersion relation and its maximum
a0_hat = theta(1);
a2_hat = theta(2);
a4_hat = theta(3);

dispersion_hat = @(k) a0_hat - a2_hat*k.^2 + a4_hat*k.^4;

% Maximise lambda(k) over 0 <= k <= largest resolved wavenumber.
k_max = fminbnd(@(k)-dispersion_hat(k), 0, max(abs(k_fft)));

max_eigenvalue = dispersion_hat(k_max);

end
