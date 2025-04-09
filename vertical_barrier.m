truncation = 25;
colocation_point_count = 250;
max_depth = 40;
barrier_depth = 40;


A = zeros(1, truncation)';
D = zeros(1, truncation)';
A(1) = 1;


gravity = 9.81;
period = 5;
frequency = 2*pi/period;
alpha = frequency^2/gravity;

% alpha = 1
% epsilon = 1e-6;


wave_numbers = dispersion_free_surface(alpha, truncation-1, max_depth) * 1i;

wave_numbers(1) = -wave_numbers(1);


weights = eye(colocation_point_count) * 2 * (barrier_depth)/(colocation_point_count-1);
colocation_points = linspace(-barrier_depth, 0, colocation_point_count);

weights(1,1) = weights(1,1)/2;
weights(end,end) = weights(end,end)/2;

kernel = zeros(colocation_point_count);

values = zeros(size(25:25:2000));
B_zeros = zeros(size(25:25:2000));
C_zeros = zeros(size(25:25:2000));

previous = zeros(colocation_point_count);
for truncation = 25:25:2000
    wave_numbers = dispersion_free_surface(alpha,truncation-1,max_depth) * 1i;
    wave_numbers(1) = -wave_numbers(1);

    A = zeros(1,truncation)';
    A(1) = 1;
    new_wave_numbers = wave_numbers(end-24:end);
    D = zeros(1,truncation)';
    

    kernel = getKernel(new_wave_numbers, colocation_points, max_depth, colocation_point_count, barrier_depth, previous);
    previous = kernel;
    % Solve for u at colocation points
    % fprintf("\n")
    % rcond(kernel)
    
    % kernel(1,:)
    % weighted_kernel = weights * kernel;
    % weighted_kernel(1,:)
    % kernel;
    % size(weights)
    u = (weights * kernel) \ f(A(1), colocation_points, wave_numbers(1), max_depth)';
    diag_u = diag(u);
    % weights * diag_u;
    % size(weights)
    % size(diag_u)
    % size(phi(colocation_points, wave_numbers', max_depth))
    % colo_phi = phi(colocation_points, wave_numbers', max_depth)

    % Calculate coefficents for phi-
    B = A - sum(( ...
        weights * diag(u) * phi(colocation_points, wave_numbers', max_depth)' ...
        ...
        ))';
    
    % Calculate coefficients for phi+
    C = D + sum(( ...
        weights * diag(u) * phi(colocation_points, wave_numbers', max_depth)' ...
        ...
        ))';
    
    
    fprintf("--------------------------------\n")
    fprintf("Truncaion %d:\n\n", truncation)
    fprintf("Energy Conservation for truncation %d: |%d|^2 + |%d|^2 - |%d|^2 =  %d\n",truncation, B(1), C(1), A(1), abs(B(1)^2) + abs(C(1)^2) - abs(A(1))^2)
    
    fprintf("Difference phi+ and phi- along x=0 with truncation %d: %d\n", truncation, sum(A) - sum(B) - sum(C))

    % Save to values
    values(truncation/25) = abs(B(1)) + abs(C(1)) - abs(A(1));
    B_zeros(truncation/25) = B(1);
    C_zeros(truncation/25) = C(1);
end


function output = phi_norm_square(wave_numbers, max_depth, barrier_depth)
    output = (...
        sinh(wave_numbers * max_depth) ...
            ./ (2 * wave_numbers)... 
        );
end



function output = phi(z,wave_numbers, max_depth)
    output = cosh(wave_numbers*(z + max_depth));
end

function output = f(A, z, wave_numbers, max_depth)
    output = A * phi(z, wave_numbers, max_depth);
end

function output = K(z, xi, wave_numbers, max_depth, barrier_depth)
    output = sum( ...
        phi(z,wave_numbers, max_depth) .* phi(xi,wave_numbers, max_depth) ...
        ./ ...
        (phi_norm_square(wave_numbers, max_depth, barrier_depth) * 1i .* wave_numbers) ...
        )/max_depth;
end

function output = getKernel(wave_numbers, colocation_points, max_depth, colocation_point_count, barrier_depth, previous)
    output = previous;
    for i = 1:colocation_point_count
    
        for j = 1:colocation_point_count
            output(i,j) = output(i,j) +  K(colocation_points(i), colocation_points(j), wave_numbers, max_depth, barrier_depth);
        end
    end
end