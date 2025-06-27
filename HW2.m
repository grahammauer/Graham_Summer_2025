%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Nystrom Discretization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u, nodes] = nystrom_discretization(N, quadrature_method)

%   Solves the 2nd-kind Fredholm equation given a number of nodes, kernel function, right-hand-side, and quadrature method. 
%
%   Inputs
%   ------
%       N : int, number of quadrature nodes to evaluate the function at
%       quadrature_method : str, either 'comp_trap' or 'gaussian'
%
%   Outputs
%   -------
%       u : vector (length N), numerical solution to the IE
    
    % Compute weights and nodes for quadrature method
    switch quadrature_method
        case 'comp_trap'
            nodes = linspace(0,1,N);
            weights = 2 * ones(1, N);
            weights(1) = 1;
            weights(end) = 1;
            weights = weights / (2*(N-1));

        case 'gaussian'
            [nodes, weights] = legpts(N, [0, 1]);
            nodes = nodes';

        otherwise
            error('Non-existent quadrature scheme: %s', quadrature_method);
    end

    % Construct the Kernel matrix (I know it is slow)

    K_matrix = zeros(N);
    for i = 1:N
        for j = 1:N
            K_matrix(i,j) = exp(nodes(i) * nodes(j)) * weights(j);
            %% Why use the same nodes for the s and t directions? As long as it matches with how RHS is discretized, does it matter?
        end
    end

    K_matrix = K_matrix + eye(N);

    % Construct the RHS
    RHS = exp(nodes) + 1 ./ (nodes + 1) .* (exp(nodes + 1) - 1);

    % Matrix solve for u
    u = K_matrix \ RHS';

end

%% 2a

n_values = [4,8,16,32,64,128,256];
n_values = linspace(4,256,253);
max_error_array = zeros(2,length(n_values));

for i = 1:length(n_values)
    [trap_values, trap_nodes] = nystrom_discretization(n_values(i), 'comp_trap');
    [gauss_values, gauss_nodes] = nystrom_discretization(n_values(i), 'gaussian');

    trap_values = trap_values(:);
    trap_nodes = trap_nodes(:);

    gauss_values = gauss_values(:);
    gauss_nodes = gauss_nodes(:);

    max_error_array(1, i) = log2(max(abs(trap_values - exp(trap_nodes))));
    max_error_array(2, i) = log2(max(abs(gauss_values - exp(gauss_nodes))));
end

% figure;
% plot(log2(n_values), max_error_array(1, :), '-o', 'DisplayName', 'Composite Trapezoid');
% hold on;
% plot(log2(n_values), max_error_array(2, :), '-s', 'DisplayName', 'Gaussian');
% xlabel('Log2 n');
% ylabel('Log2 Maximum Error');
% title('Log2 Max Error');
% legend('Location', 'southwest');
% grid on;
% saveas(gcf, 'HW2_2a.png')

%% 2d
[u, x] = nystrom_discretization(5, 'gaussian');
poly = polyfit(x, u, length(x) - 1);
xx = linspace(0,1,10000);
poly_val = polyval(poly, xx);
% figure;
% plot(xx, log10(abs(poly_val - exp(xx))), '-', 'DisplayName', 'Interpolant Error');
% hold on;
% scatter(x, log10(abs(polyval(poly, x) - exp(x))), 300, '.', 'DisplayName', 'Nystrom Error');
% xlabel('x')
% ylabel('u')
% title('Difference between interpolated and analytical solution')
% legend('Location', 'northwest');
% grid on;
% saveas(gcf, 'HW2_2d.png')

% figure;
% plot(log10(n_values), log10(2.^max_error_array(1, :)), '-o', 'DisplayName', 'Composite Trapezoid');
% hold on;
% plot(log10(n_values), log10(2.^max_error_array(2, :)), '-s', 'DisplayName', 'Gaussian');
% xlabel('Log10 n');
% ylabel('Log10 Maximum Error');
% title('Log10 Max Error');
% legend('Location', 'southwest');
% grid on;
% saveas(gcf, 'HW2_2d.2.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Directional Derivative of Fundamental Solution %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dd = Direct_Derivative_Fundamental_Laplace(x, y, n_y)

%   Returns the directional derivative at N x-points for one given y and
%   n_y value.
%
%   Inputs
%   ------
%       x : vector, (2,N)
%       y : vector, (2,1)
%       n_y : vector (2,1)
%
%   Outputs
%   -------
%       dd : vector (1, N)

    dd = (n_y(1).*(x(1,:)-y(1)) + n_y(2).*(x(2,:)-y(2))) ./ (2 * pi * ((x(1,:)-y(1)).^2 + (x(2,:)-y(2)).^2));

end

[x_1, x_2] = meshgrid(-1:0.0101:1, -1:0.0101:1);

x = [x_1(:)'; x_2(:)'];

y = [0.5; -0.2];
n_y = [1/2; sqrt(3)/2];

dd = Direct_Derivative_Fundamental_Laplace(x, y, n_y);

dd_grid = reshape(dd, size(x_1));

% figure;
% contourf(x_1, x_2, dd_grid, linspace(-3, 3, 21));
% colorbar;
% axis equal;
% title('Contour Plot of Directional Derivative');
% xlabel('$x_1$', 'Interpreter', 'latex');
% ylabel('$x_2$', 'Interpreter', 'latex');
% saveas(gcf, 'HW2_3.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4) Setting up Closed Curves %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = boundary(R, Rp, theta)

%   Returns a $2\pi$-periodic parameterization of a closed curve given
%   $R(\theta)$ and $R'(\theta)$.
%
%   Inputs
%   ------
%       R : func, radius as a function of theta
%       Rp : func, derivative of radius as a function of theta
%       theta : vector (1,n) of points to evaluate z at
%
%   Outputs
%   -------
%       z : vector (z_1, z_2) representing the boundary for the closed
%       curve.

    r = R(theta);
    z = [r .* cos(theta);
         r .* sin(theta)];
end

function zp = boundary_derivative(R, Rp, theta)

%   Returns a $2\pi$-periodic parameterization of the derivative of a 
%   closed curve given $R(\theta)$ and $R'(\theta)$.
%
%   Inputs
%   ------
%       R : func, radius as a function of theta
%       Rp : func, derivative of radius as a function of theta
%       theta : vector (1,n) of points to evaluate z at
%
%   Outputs
%   -------
%       zp : vector (zp_1, zp_2) representing the derivative of the 
%           boundary for the closed curve.
    
    r = R(theta);
    rp = Rp(theta);
    zp = [rp .* cos(theta); rp .* sin(theta)] + [r .* -sin(theta); r .* cos(theta)];

end

function n = normal(R, Rp, theta)

%   Returns the unit normal vector to some boundary parameterized by 
%   $R(\theta)$ and $R'(\theta)$.
%
%   Inputs
%   ------
%       R : func, radius as a function of theta
%       Rp : func, derivative of radius as a function of theta
%       theta : vector (1,n) of points to evaluate z at
%
%   Outputs
%   -------
%       n : vector (n_1, n_2) representing the derivative of the boundary for the closed
%       curve.

    zp = boundary_derivative(R, Rp, theta);
    n = [zp(2, :); -zp(1, :)] ./ vecnorm(zp, 2, 1);

end

function k = boundary_curvature(R, Rp, Rpp, theta)

%   Returns a $2\pi$-periodic parameterization of the curvature of a 
%   closed curve given $R''(\theta)$.
%
%   Inputs
%   ------
%       Rpp : func, second derivative of radius as a function of theta
%       theta : vector (1,n) of points to evaluate z at
%
%   Outputs
%   -------
%       k : vector (1, n) representing the curvature of the 
%           boundary for the closed curve.

    r = R(theta);
    rp = Rp(theta);
    rpp = Rpp(theta);

    k = (r.^2 + 2 * rp.^2 - r .* rpp) ./ ((r.^2 + rp.^2).^(3/2));

end

theta = linspace(0, 2*pi, 30+1);
theta(end) = [];
xx = linspace(0, 2*pi, 1000);

R = @(theta) 1 + 0.3*cos(3*theta);
Rp = @(theta) -0.9*sin(3*theta);

z = boundary(R, Rp, theta);
zp = boundary_derivative(R, Rp, theta);
z_bdd = boundary(R, Rp, xx);
n = normal(R, Rp, theta);

% figure;
% plot(z_bdd(1,:), z_bdd(2,:), 'DisplayName', 'Boundary')
% hold on;
% quiver(z(1,:), z(2,:), n(1,:), n(2,:), 'DisplayName', 'Normal Vector');
% quiver(z(1,:), z(2,:), zp(1,:), zp(2,:), 'DisplayName', 'Speed Function');
% axis equal;
% xlabel('$x_1$', 'Interpreter', 'latex');
% ylabel('$x_2$', 'Interpreter', 'latex');
% title('Boundary with normal vectors');
% grid on;
% legend;
% saveas(gcf, 'HW2_4.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5) Checking Gauss' Law %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the boundary
theta = linspace(0, 2*pi, 30+1);
theta(end) = [];
xx = linspace(0, 2*pi, 1000);

R = @(theta) 1 + 0.3*cos(3*theta);
Rp = @(theta) -0.9*sin(3*theta);

z = boundary(R, Rp, theta);
zp = boundary_derivative(R, Rp, theta);
z_bdd = boundary(R, Rp, xx);
n = normal(R, Rp, theta);

% Create [-2,2]^2

[x_1, x_2] = meshgrid(-2:0.015:2, -2:0.015:2);
x = [x_1(:)'; x_2(:)'];

% Evaluate the surface integral at points in [0,2]^2

dd = zeros(1, size(x, 2));
for i = 1:length(theta)

    y = z(:,i);
    n_y = n(:,i);

    dd = dd + Direct_Derivative_Fundamental_Laplace(x, y, n_y) * norm(zp(:,i));

end
dd = 2 * pi * dd / length(theta);
dd_grid = reshape(dd, size(x_1));

% figure;
% contourf(x_1, x_2, dd_grid, linspace(-1.5, 0.5, 11));
% hold on;
% plot(z_bdd(1,:), z_bdd(2,:), 'r');
% scatter(z(1,:), z(2,:), 30, 'r', 'filled');
% colorbar;
% axis equal;
% grid on;
% title("Gauss' Law Solution")
% saveas(gcf, 'HW2_5b.2.png')

theta_flat = atan2(x_2(:), x_1(:));
r_flat = sqrt(x_1(:).^2 + x_2(:).^2);

R_vals = R(theta_flat);

inside_mask_grid = reshape(r_flat < R_vals, size(x_1));

comparison_mask = double(inside_mask_grid);

abs_error = log10(abs(dd_grid + comparison_mask));

% figure;
% contourf(x_1, x_2, abs_error, 15);
% hold on;
% plot(z_bdd(1,:), z_bdd(2,:), 'r');
% scatter(z(1,:), z(2,:), 30, 'r', 'filled');
% xlabel('x_{1}');
% ylabel('x_{2}');
% colorbar;
% grid on;
% axis equal;
% title("log10 Gauss' Law Total Error");
% saveas(gcf, 'HW2_5b.png')

% Determine convergence of specific point (0.2, 0.1) for n = [30,60,120,240,480]

% n_vals = [15, 30,60,120,240,480,960];
n_vals = 4:120;
u = zeros(size(n_vals));
x = [0.2, 0.1]';
for j = 1:length(n_vals)

    n_val = n_vals(j);
    
    % Define the boundary
    theta = linspace(0, 2*pi, n_val+1);
    theta(end) = [];

    % R = @(theta) 1;
    % Rp = @(theta) 0;
    
    R = @(theta) 1 + 0.3*cos(3*theta);
    Rp = @(theta) -0.9*sin(3*theta);
    
    z = boundary(R, Rp, theta);
    zp = boundary_derivative(R, Rp, theta);
    n = normal(R, Rp, theta);

    for i = 1:n_val

        y = z(:,i);
        n_y = n(:,i);

        u(j) = u(j) + Direct_Derivative_Fundamental_Laplace(x, y, n_y) * norm(zp(:,i));

    end

    u(j) = 2 * pi * u(j) / n_val;

end

% figure;
% loglog(n_vals, abs(u + 1), 'ro-', 'DisplayName', 'Error at (0.2, 0.1)');
% hold on;
% loglog(n_vals, 0.1 * n_vals.^(-1), 'b--', 'DisplayName', 'N^{-1}');
% legend;
% xlabel('N');
% ylabel('log10|u + 1|');
% grid on;
% title('Convergence of u(0.2, 0.1)');
% saveas(gcf, 'HW2_5c.png');

%%%%%%%%%%%%%%%%%%%%%%%%%
% 6) Laplace's Equation %
%%%%%%%%%%%%%%%%%%%%%%%%%