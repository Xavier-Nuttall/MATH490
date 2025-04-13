%% Parameters
truncation = 5;
colocationPointCount = 3;
maxDepth = 40;
barrierDepth = 40;


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

weights = eye(colocationPointCount) * (barrierDepth)/(colocationPointCount-1);
colocationPoints = linspace(-barrierDepth, 0, colocationPointCount).';

weights(1,1) = weights(1,1)/2;
weights(end,end) = weights(end,end)/2;



%% For testing the Kernel function
% kernelSlow should be the same as kernel
kernelSlow = zeros(colocationPointCount);
for i = 1:colocationPointCount
    for j = 1:colocationPointCount
        for k = 1:truncation
            kernelSlow(j,i) = kernelSlow(j,i) + ...
                phi(colocationPoints(i), waveNumbers(k), maxDepth) ...
                * phi(colocationPoints(j), waveNumbers(k), maxDepth) ...
                / (phi_norm_square(waveNumbers(k), maxDepth, barrierDepth) * 1i * waveNumbers(k));
        end
    end
end
%%

%% For testing the K function
% KSlow should be the same as K
kSlowFirst = 0;
for k = 1:truncation
    kSlowFirst = kSlowFirst + ...
        phi(-20, waveNumbers(k), maxDepth) ...
        * phi(-20, waveNumbers(k), maxDepth) ...
        / (phi_norm_square(waveNumbers(k), maxDepth, barrierDepth) * 1i * waveNumbers(k));
end
kFirst = K(-20,-20, waveNumbers, maxDepth, barrierDepth);
%%

%% Calculating u and kernel
kernel = zeros(colocationPointCount);

kernel = getKernel(waveNumbers, colocationPoints, maxDepth, colocationPointCount, barrierDepth, kernel);

functionValue = f(A(1), colocationPoints, waveNumbers(1), maxDepth);

u = (weights * kernel) \ functionValue;

diagU = diag(u);
%%

%% Weights x kernel should be identity
kernelDifference = weights * kernel - eye(colocationPointCount);
%%

%% Calculate qcoefficients for B and C

% Calculate coefficents for phi-
B = A - sum(( ...
    weights * diagU * phi(colocationPoints, waveNumbers, maxDepth) ...
    ...
    )).';

% Calculate coefficients for phi+
C = D + sum(( ...
    weights * diagU * phi(colocationPoints, waveNumbers, maxDepth) ...
    ...
    )).';
%%

%% Printouts
fprintf("--------------------------------\n")
fprintf("Truncaion %d:\n\n", truncation)
fprintf("Energy Conservation for truncation %d: |%d|^2 + |%d|^2 - |%d|^2 =  %d\n",truncation, B(1), C(1), A(1), abs(B(1)^2) + abs(C(1)^2) - abs(A(1))^2)

fprintf("Difference phi+ and phi- along x=0 with truncation %d: %d\n", truncation, sum(A) - sum(B) - sum(C))

% Save to values
% values(truncation/25) = abs(B(1)) + abs(C(1)) - abs(A(1));
% bZeros(truncation/25) = B(1);
% cZeros(truncation/25) = C(1);
%%

values = zeros(truncation);
bZeros = zeros(truncation);
cZeros = zeros(truncation);

uCorrect = A(1) * waveNumbers(1) * 1i *  cosh(waveNumbers(1) *(colocationPoints + maxDepth))/cosh(waveNumbers(1) * maxDepth);

difference = uCorrect - u;


%% Functions
function output = phi_norm_square(waveNumbers, maxDepth, barrierDepth) %#ok<INUSD>
    N = cosh(waveNumbers * maxDepth);
    

    % output = 1/2 * (barrierDepth + ...
    %     (...
    %         sinh(2 * waveNumbers * maxDepth) - sinh(2 * waveNumbers * maxDepth - 2* waveNumbers * barrierDepth)...
    %     )./ (2 * waveNumbers)) ./ N.^2;


    output = (...
        maxDepth + sinh(2 * waveNumbers * maxDepth) ...
            ./ (2 * waveNumbers)... 
        ) /2 ...
        ./(N.^2);
end

function output = getKernel2(waveNumbers, colocationPoints, maxDepth, barrierDepth)
    Phi = phi(colocationPoints, waveNumbers, maxDepth);
    K = phi_norm_square(waveNumbers, maxDepth, barrierDepth) .* 1i .* waveNumbers;

    output = Phi * K * Phi.';
end

function output = phi(z,waveNumbers, maxDepth)
    % normalisation constant
    N = cosh(waveNumbers * maxDepth);

    output = cosh((z + maxDepth) * waveNumbers) ... 
    ./ N;
end

function output = f(A, z, waveNumbers, maxDepth)
    output = A * phi(z, waveNumbers, maxDepth);
end


% Working? K function
function output = K(z, xi, waveNumbers, maxDepth, barrierDepth)
    output = sum( ...
        phi(z,waveNumbers, maxDepth) .* phi(xi,waveNumbers, maxDepth) ...
        ./ ...
        (phi_norm_square(waveNumbers, maxDepth, barrierDepth) * 1i .* waveNumbers) ...
        );
end

function output = getKernel(waveNumbers, colocationPoints, maxDepth, colocationPointCount, barrierDepth, previous)
    output = previous;
    for i = 1:colocationPointCount
    
        for j = 1:colocationPointCount
            output(i,j) = output(i,j) +  K(colocationPoints(i), colocationPoints(j), waveNumbers, maxDepth, barrierDepth);
        end
    end
end
%%