%% Parameters
truncation = 1000;


colocationPointCount = 200;
maxDepth = 20;
barrierDepth = 5;


A = zeros(1, truncation).';
D = zeros(1, truncation).';
A(1) = 1;


gravity = 9.81;
period = 5;
frequency = 2*pi/period;
alpha = frequency^2/gravity;
%%

waveNumbers = dispersion_free_surface(alpha, truncation-1, maxDepth) * 1i;
waveNumbers(1) = -waveNumbers(1);

colocationPoints = linspace(-barrierDepth, 0, colocationPointCount).';

deltaZ = colocationPoints(2) - colocationPoints(1);
weights = eye(colocationPointCount) * deltaZ;

weights(1,1) = weights(1,1)/2;
weights(end,end) = weights(end,end)/2;



%% For testing the Kernel function
% kernelSlow should be the same as kernel
% tic
% kernelSlow = zeros(colocationPointCount);
% for i = 1:colocationPointCount
%     for j = 1:colocationPointCount
%         for k = 1:truncation
%             kernelSlow(j,i) = kernelSlow(j,i) + ...
%                 phi(colocationPoints(i), waveNumbers(k), maxDepth) ...
%                 * phi(colocationPoints(j), waveNumbers(k), maxDepth) ...
%                 / (phi_norm_square(waveNumbers(k), maxDepth, barrierDepth) * 1i * waveNumbers(k));
%         end
%     end
% end
% toc
%%

%% Calculating u and kernel
tic
kernel = getKernel(waveNumbers, colocationPoints, maxDepth, barrierDepth);
toc


functionValue = f(A(1), colocationPoints, waveNumbers(1), maxDepth);

u = (kernel * weights) \ functionValue;

%%

%% Calculate qcoefficients for B and C

% Calculate coefficents for phi-
B = A - diag(1./(1i * waveNumbers .* phi_norm_square(waveNumbers, maxDepth, barrierDepth))) * phi(colocationPoints, waveNumbers, maxDepth).' * weights * u; 

% Calculate coefficients for phi+
C = D + diag(1./(1i * waveNumbers .* phi_norm_square(waveNumbers, maxDepth, barrierDepth))) * phi(colocationPoints, waveNumbers, maxDepth).' * weights * u;
%%

%% Printouts
fprintf("--------------------------------\n")
fprintf("Truncaion %d:\n\n", truncation)
fprintf("Energy Conservation for truncation %d: |%d|^2 + |%d|^2 - |%d|^2 =  %d\n",truncation, B(1), C(1), A(1), abs(B(1)^2) + abs(C(1)^2) - abs(A(1))^2)

fprintf("Difference phi+ and phi- along x=0 with truncation %d: %d\n", truncation, sum(A) - sum(B) - sum(C))

values = zeros(truncation);
bZeros = zeros(truncation);
cZeros = zeros(truncation);

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

%%