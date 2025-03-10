close all

% Set number of iterations to be performed
nk = 500

% Set parameters alpha and beta
alpha = 2;
beta  = 3;

% Set the number of meshpoints so the interior has N x N such points
N = 50;

% Compute the distance between mesh points, in each direction
h = 1/(N+1);

% Compute the optimal relaxation parameter omega
w = 2 / (1 + sin(pi * h));

% Convergence tolerance - using 0.000001
converge = 1e-6;

% Set convergence flags
SOR_converge = false;
GS_converge = false;

% We will have arrays that capture the boundary as well as the interior
% meshpoints.  As a result, we need those arrays to be of size (N+2) x
% (N+2) so that those indexed 2:N+1, 2:N+1 represent the interior.  

% Compute the x-values at each point i,j, including the boundary
x = h * [ 0:N+1 ];   % Notice this creates a row vector

% Compute the y-values at each point i,j, including the boundary
y = h * [ 0:N+1 ];   % Notice this creates a row vector

% Create an array that captures the load at each point i,j
for i=1:N+2
    for j=1:N+2
        F( i,j ) = ...
            ( alpha^2 + beta^2 ) * pi^2 * sin( alpha * pi * x( i ) ) * sin( beta * pi * y( j ) );
    end
end

% Set the initial values at the mesh points.  Use SOR to indicate
% the array for the SOR iteration.  Use GS to indicate the array
% is for the Gauss-Seidel iteration.
SOR = zeros( N+2, N+2 );
GS = zeros( N+2, N+2 );

% Perform nk iterations
for k = 1:nk
    k           % print current iteration index

    % Store old SOR and GS values for convergence tolerance computation
    SOR_old = SOR
    GS_old = GS
    
    % update all the interior points (Jacobi iteration)
    for i=2:N+1
        for j=2:N+1
             SOR( i,j ) = ((1-w) * SOR( i,j )) + (w * (( SOR( i, j-1 ) + SOR( i-1, j ) + SOR( i+1, j ) + SOR( i, j+1 ) + h^2 * F( i, j ) ) / 4));
        end
    end 
    
    subplot( 3, 1, 1 );  % plot in top graph
    mesh( x, y, SOR );
    axis( [ 0 1 0 1 -1.5 1.5 ]);

    % update all the interior points (Gauss-Seidel iteration)
    for i=2:N+1
        for j=2:N+1
            GS( i,j ) = ( GS( i, j-1 ) + GS( i-1, j ) + GS( i+1, j ) + GS( i, j+1 ) + h^2 * F( i, j ) ) / 4;
        end
    end
    
    subplot( 3, 1, 2 );  % plot in bottom graph
    mesh( x, y, GS );
    axis( [ 0 1 0 1 -1.5 1.5 ]);
    
    subplot( 3, 1, 3);
    mesh( x, y, SOR - GS );
    axis( [ 0 1 0 1 -0.05 0.05 ])

    % Compute when the convergence has reached 0 (or close to 0 in this case for an estiamtion) ---
    max_diff_SOR = max(max(abs(SOR - SOR_old)));
    max_diff_GS = max(max(abs(GS - GS_old)));
    
    if ~SOR_converge && max_diff_SOR < converge
        fprintf('SOR converged at iteration %d\n', k);
        sor_converged = true;
    end
    
    if ~GS_converge && max_diff_GS < converge
        fprintf('GS converged at iteration %d\n', k);
        gs_converged = true;
    end
    
    % If both have converged, stop early
    if SOR_converge && GS_converge
        fprintf('Both SOR and GS converged at iteration %d\n', k);
        break;
    end
    
    % wait to continue to the next iteration
    next = input( 'press RETURN to continue' );
end
