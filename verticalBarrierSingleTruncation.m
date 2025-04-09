truncation = 25;
colocationPointCount = 250;
maxDepth = 40;
barrierDepth = 40;


A = zeros(1, truncation)';
D = zeros(1, truncation)';
A(1) = 1;


gravity = 9.81;
period = 5;
frequency = 2*pi/period;
alpha = frequency^2/gravity;

% alpha = 1
% epsilon = 1e-6;


maveNumbers = dispersion_free_surface(alpha, truncation-1, maxDepth) * 1i;

maveNumbers(1) = -maveNumbers(1);


weights = eye(colocationPointCount) * 2 * (barrierDepth)/(colocationPointCount-1);
colocationPoints = linspace(-barrierDepth, 0, colocationPointCount);

weights(1,1) = weights(1,1)/2;
weights(end,end) = weights(end,end)/2;

kernel = zeros(colocationPointCount);

kernel = getKernel(maveNumbers, colocationPoints, maxDepth, colocationPointCount, barrierDepth, kernel);

u = (weights * kernel) \ f(A(1), colocationPoints, maveNumbers(1), maxDepth)';
diag_u = diag(u);


% Calculate coefficents for phi-
B = A - sum(( ...
    weights * diag(u) * phi(colocationPoints, maveNumbers', maxDepth)' ...
    ...
    ))';

% Calculate coefficients for phi+
C = D + sum(( ...
    weights * diag(u) * phi(colocationPoints, maveNumbers', maxDepth)' ...
    ...
    ))';


fprintf("--------------------------------\n")
fprintf("Truncaion %d:\n\n", truncation)
fprintf("Energy Conservation for truncation %d: |%d|^2 + |%d|^2 - |%d|^2 =  %d\n",truncation, B(1), C(1), A(1), abs(B(1)^2) + abs(C(1)^2) - abs(A(1))^2)

fprintf("Difference phi+ and phi- along x=0 with truncation %d: %d\n", truncation, sum(A) - sum(B) - sum(C))

% Save to values
% values(truncation/25) = abs(B(1)) + abs(C(1)) - abs(A(1));
% bZeros(truncation/25) = B(1);
% cZeros(truncation/25) = C(1);

values = zeros(size(25:25:2000));
bZeros = zeros(size(25:25:2000));
cZeros = zeros(size(25:25:2000));

uCorrect = phi(colocationPoints, maveNumbers', maxDepth) * A(1);

difference = abs(uCorrect - u);



function output = phi_norm_square(maveNumbers, maxDepth, barrierDepth) %#ok<INUSD>
    output = (...
        sinh(maveNumbers * maxDepth) ...
            ./ (2 * maveNumbers)... 
        );
end



function output = phi(z,maveNumbers, maxDepth)
    output = cosh(maveNumbers*(z + maxDepth));
end

function output = f(A, z, maveNumbers, maxDepth)
    output = A * phi(z, maveNumbers, maxDepth);
end

function output = K(z, xi, maveNumbers, maxDepth, barrierDepth)
    output = sum( ...
        phi(z,maveNumbers, maxDepth) .* phi(xi,maveNumbers, maxDepth) ...
        ./ ...
        (phi_norm_square(maveNumbers, maxDepth, barrierDepth) * 1i .* maveNumbers) ...
        )/maxDepth;
end

function output = getKernel(maveNumbers, colocationPoints, maxDepth, colocationPointCount, barrierDepth, previous)
    output = previous;
    for i = 1:colocationPointCount
    
        for j = 1:colocationPointCount
            output(i,j) = output(i,j) +  K(colocationPoints(i), colocationPoints(j), maveNumbers, maxDepth, barrierDepth);
        end
    end
end