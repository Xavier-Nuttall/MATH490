%% Parameters
truncation = 10000;
modeTruncation = 200;


colocationPointCount = 200;
maxDepth = 20;
barrierDepth = 4;


A = zeros(1, modeTruncation).';
D = zeros(1, modeTruncation).';
A(1) = 1;

gravity = 9.81;

%% Values
colocationPoints = linspace(-barrierDepth, 0, colocationPointCount).';

deltaZ = colocationPoints(2) - colocationPoints(1);
weights = eye(colocationPointCount) * deltaZ;

weights(1,1) = weights(1,1)/2;
weights(end,end) = weights(end,end)/2;


%% Plotting params
N = 20;
offset = 5;
phiPlusLegend = strings(size(1:N));
phiMinusLegend = strings(size(1:N));

% Create a color spectrum
colours = winter(N);

phiPlusFigure = figure();
hold on
title("$\phi^+$ for various periods", 'Interpreter', 'latex')


phiMinusFigure = figure();
hold on
title("$\phi^-$ for various periods", 'Interpreter', 'latex')

%% Iterate over different periods
for period = 1:N
    frequency = 2*pi/(period + offset);
    alpha = frequency^2/gravity;

    % Calculate new wave numbers
    waveNumbers = dispersion_free_surface(alpha, truncation-1, maxDepth) * 1i;
    waveNumbers(1) = -waveNumbers(1);
    reducedWaveNumbers = waveNumbers(1:modeTruncation);


    %% Calculating u and kernel
    tic
    kernel = getKernel(waveNumbers, colocationPoints, maxDepth, barrierDepth);
    toc
    
    
    functionValue = f(A(1), colocationPoints, waveNumbers(1), maxDepth);
    
    u = (kernel * weights) \ functionValue;
    
    
    
    %% Calculate qcoefficients for B and C
    
    % Calculate coefficents for phi-
    B = A - diag(1./(1i * reducedWaveNumbers .* phi_norm_square(reducedWaveNumbers, maxDepth, barrierDepth))) * phi(colocationPoints, reducedWaveNumbers, maxDepth).' * weights * u; 
    
    % Calculate coefficients for phi+
    C = D + diag(1./(1i * reducedWaveNumbers .* phi_norm_square(reducedWaveNumbers, maxDepth, barrierDepth))) * phi(colocationPoints, reducedWaveNumbers, maxDepth).' * weights * u;
    
    
    %% Plotting
    fullHeightPoints = linspace(-maxDepth, 0, round(fix(colocationPointCount/barrierDepth) * maxDepth)).';
    
    phiPlus = phi(fullHeightPoints, reducedWaveNumbers, maxDepth) * (A + B) ;
    
    phiMinus = phi(fullHeightPoints, reducedWaveNumbers, maxDepth) * (C + D);

    figure(phiPlusFigure)
    plot(phiPlus,fullHeightPoints, 'r.', 'MarkerSize', 1, 'Color', colours(period,:))

    phiPlusLegend(period) = sprintf('$\\phi^+_{%d}$', period+offset);
    
    figure(phiMinusFigure)
    plot(phiMinus,fullHeightPoints, 'b.', 'MarkerSize', 1, 'Color', colours(period,:))
    phiMinusLegend(period) = sprintf('$\\phi^-_{%d}$', period+offset);


end
figure(phiMinusFigure)
legend(phiMinusLegend, 'Interpreter', 'latex')

figure(phiPlusFigure)
legend(phiPlusLegend, 'Interpreter', 'latex')

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