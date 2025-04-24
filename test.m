% sanity check

f = @(n, x) hoFunction(n, x);

fplot(@(x) f(0, x), [-0.2, 0.2])

f2 = @(x) abs(f(0, x)).^2;
integral(f2, -1, 1)

function psi = hoFunction(n, x, omega, m, hbar)
% ho_wavefunction  Returns the nth quantum harmonic oscillator wavefunction
%   psi = ho_wavefunction(n, x, omega, m, hbar)
%   n     - quantum number (non-negative integer)
%   x     - position array
%   m     - mass
%   omega - angular frequency
%   hbar  - reduced Planck constant

    if nargin < 5
        hbar = 1.054571817e-34;
    end
    if nargin < 4
        m = 2.2069391e-25;
    end
    if nargin < 3
        omega = 75e3*2*pi;
    end

    x = 1e-6 * x; % make x unit in um

    % Rescale factor
    alpha = sqrt(m * omega / hbar);
    xi = alpha * x;

    % Hermite polynomial
    Hn = hermiteH(n, xi);  % Symbolic toolbox function

    % Normalization factor
    norm_factor = 1 / sqrt(2^n * factorial(n)) * (alpha / pi)^(1/4);

    % Wavefunction
    psi = norm_factor * Hn .* exp(-xi.^2 / 2);
end
