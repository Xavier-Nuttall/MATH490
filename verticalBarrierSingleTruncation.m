%% Parameters
function [A,B,C,D] = verticalBarrierSingleTruncation(maxDepth, barrierDepth, period, colocationPointCount, truncation, modeTruncation, positiveIncidentAmplitude, negativeIncidentAmplitude)
if nargin < 9
    negativeIncidentAmplitude = 0;
end
if nargin < 8
    positiveIncidentAmplitude = 1;
end
if nargin < 7
    modeTruncation = 10;
end
if nargin < 6
    truncation = 1000;
end
if nargin < 5
    colocationPointCount = 50;
end

A = zeros(1, modeTruncation).';
D = zeros(1, modeTruncation).';
A(1) = positiveIncidentAmplitude;
D(1) = negativeIncidentAmplitude;

gravity = 9.81;
frequency = 2*pi/period;
alpha = frequency^2/gravity;
%% Wave numbers

waveNumbers = dispersion_free_surface(alpha, truncation-1, maxDepth) * 1i;
waveNumbers(1) = -waveNumbers(1);
reducedWaveNumbers = waveNumbers(1:modeTruncation);

colocationPoints = linspace(-barrierDepth, 0, colocationPointCount).';

deltaZ = colocationPoints(2) - colocationPoints(1);
weights = eye(colocationPointCount) * deltaZ;

weights(1,1) = weights(1,1)/2;
weights(end,end) = weights(end,end)/2;

%% Calculating u and kernel
kernel = getKernel(waveNumbers, colocationPoints, maxDepth, barrierDepth);


functionValue = f(A(1), colocationPoints, waveNumbers(1), maxDepth);

u = (kernel * weights) \ functionValue;

%% Calculate qcoefficients for B and C

% Calculate coefficents for phi-
B = A - diag(1./(1i * reducedWaveNumbers .* phi_norm_square(reducedWaveNumbers, maxDepth, barrierDepth))) * phi(colocationPoints, reducedWaveNumbers, maxDepth).' * weights * u;

% Calculate coefficients for phi+
C = D + diag(1./(1i * reducedWaveNumbers .* phi_norm_square(reducedWaveNumbers, maxDepth, barrierDepth))) * phi(colocationPoints, reducedWaveNumbers, maxDepth).' * weights * u;

end
%% Functions
function output = phi_norm_square(waveNumbers, maxDepth, barrierDepth) %#ok<INUSD>

output = (...
    maxDepth + sinh(2 * waveNumbers * maxDepth) ...
    ./ (2 * waveNumbers)...
    ) /2;
end

function output = getKernel(waveNumbers, colocationPoints, maxDepth, barrierDepth)
Phi = phi(colocationPoints, waveNumbers, maxDepth);
K = diag(1./(phi_norm_square(waveNumbers, maxDepth, barrierDepth) .* 1i .* waveNumbers));

output = Phi * K * Phi.';
end

function output = phi(z,waveNumbers, maxDepth)
output = cosh((z + maxDepth) * waveNumbers);
end

function output = f(A, z, waveNumbers, maxDepth)
output = A * phi(z, waveNumbers, maxDepth);
end

