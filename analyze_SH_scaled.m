function [theta, k_max, max_eigenvalue] = analyze_SH_scaled(q, tau, X)
% analyze_SH_scaled
%
% Fits the reduced spatial model
%
%   q_tau = a0*q + a2*q_XX + a4*q_XXXX + s*q^3
%
% using all spatial locations and time steps in the dataset.
% The returned coefficient vector is
%   theta = [a0_hat; a2_hat; a4_hat; s_hat].
% k_max and max_eigenvalue are retained because they are used in the
% dispersion-summary plots. The dispersion-function handle itself is
% only needed locally and is therefore not returned.
%
% Spatial derivatives are evaluated spectrally using the FFT on the
% periodic spatial domain.

%% Data and grids
tau = tau(:);
X = X(:);

[Nt,Nx] = size(q);

dtau_vec = diff(tau);

dX = mean(diff(X));
L = Nx*dX;

%% Fourier wavenumbers
% Define physical wavenumbers in centered order and rearrange them to
% MATLAB's native unshifted FFT storage order.
if mod(Nx,2) == 0
    k_centered = (2*pi/L)*(-Nx/2:Nx/2-1)';
else
    k_centered = (2*pi/L)*(-(Nx-1)/2:(Nx-1)/2)';
end
k_fft = ifftshift(k_centered);

%% Centre data
Q = q' - mean(q(:));     % Nx x Nt

%% Construct regression matrices
A = zeros(Nx,Nt-1);
B = zeros(Nx,Nt-1,4);

for j = 1:Nt-1

    % Time derivative
    A(:,j) = (Q(:,j+1)-Q(:,j))/dtau_vec(j);

    % q
    B(:,j,1) = Q(:,j);

    % Spatial derivatives
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

%% Estimated dispersion relation
a0_hat = theta(1);
a2_hat = theta(2);
a4_hat = theta(3);

dispersion_hat = @(k) ...
    a0_hat - a2_hat*k.^2 + a4_hat*k.^4;

[k_max,~] = fminbnd(@(k)-dispersion_hat(k), ...
    0,max(abs(k_fft)));

max_eigenvalue = dispersion_hat(k_max);

end
