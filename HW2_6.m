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

    k = (r.^2 + 2 * rp.^2 - r .* rpp) ./ ((r.^2 + rp.^2).^(3/2)); % I ended up just using the definition for a polar defined function

end

function tau = solve_laplace_density(n, R, Rp, Rpp, f)

%   Returns the density of the dipole layer for the Laplace equation solved over a boundary 
%   parameterized by $R(\theta)$, $R'(\theta)$, and $R''(\theta)$ with 
%   $u(x) = f$ for $x\in\Gamma$.
%
%   Inputs
%   ------
%       n : int, number of boundary nodes 
%       R : func, radius as a function of theta
%       Rp : func, first derivative of radius as a function of theta
%       Rpp : func, second derivative of radius as a function of theta
%       f : func, solution on the boundary, evaluable at point (x1, x2)
%
%   Outputs
%   -------
%       tau : vector (1, n) representing the density of the dipole layer on
%       the boundary

    theta = linspace(0,2*pi,n+1);
    theta(end) = [];
    
    % Define boundary quantities
    z = boundary(R, Rp, theta);
    zp = boundary_derivative(R, Rp, theta);
    speed = vecnorm(zp, 2, 1);
    n_y = normal(R, Rp, theta);
    kappa = boundary_curvature(R, Rp, Rpp, theta);
    
    % Construct Nystrom matrix
    K = zeros(n,n);
    for i = 1:n
    
        for j = 1:n
            
            % Fill diagonal
            if i == j
    
                K(i,j) = kappa(i) / (4 * pi);
            
            % Fill off-diagonal
            else
    
                r1 = z(1,i)-z(1,j);
                r2 = z(2,i)-z(2,j);
    
                K(i,j) = 1 / (2*pi) * (n_y(1,i)*r1 + n_y(2,i)*r2) / (r1^2 + r2^2);
    
            end
    
        end
    
    end
    
    % Take a look at the kernel
    figure;
    imagesc(K);
    colorbar;
    title('Kernel Function');
    
    % Construct operator
    D = eye(n) - 2 * K * diag(speed .* (2 * pi / n * ones(1,n)));
    
    % Initialize boundary conditions
    boundary_data = f(z(1,:), z(2,:));
    
    % Construct dipole density
    tau = D \ (-2 * boundary_data)';

    % Take a look at tau
    figure;
    plot(linspace(0,2*pi,n), tau);
    colorbar;
    title('\tau Dipole Density');
    xlabel('\theta');
    ylabel('\tau');

end

% Define functions and initial test case for n=30
n = 30;

R = @(theta) 1 + 0.3*cos(3*theta);
Rp = @(theta) -0.9*sin(3*theta);
Rpp = @(theta) -2.7*cos(3*theta);

f = @(x1, x2) cos(x1) .* exp(x2);

% Define useful boundary information
theta = linspace(0,2*pi,n+1);
theta(end) = []; % Ensure that there is not a repeated point
z = boundary(R, Rp, theta);
zp = boundary_derivative(R, Rp, theta);
speed = vecnorm(zp, 2, 1);
n_y = normal(R, Rp, theta);
kappa = boundary_curvature(R, Rp, Rpp, theta);
boundary_data = f(z(1,:), z(2,:));

%% Construct solution
tau = solve_laplace_density(n, R, Rp, Rpp, f);

% Create [-2,2]^2
[x_1, x_2] = meshgrid(-2:0.015:2, -2:0.015:2);
x = [x_1(:)'; x_2(:)'];

% Evaluate the surface integral at points in [-2,2]^2 and compute the solution in the domain
dd = zeros(1, size(x, 2));
for i = 1:n
    
    % For each node on the boundary, add its contribution to each x value
    dd = dd + Direct_Derivative_Fundamental_Laplace(x, z(:,i), n_y(:,i)) * tau(i) * speed(i) * 2 * pi / n; % Exactly the same as in the Gauss' law test, except $\tau$ is now variable

end

dd_grid = reshape(dd, size(x_1));

% Create an interior mask
theta_flat = atan2(x_2(:), x_1(:));
r_flat = sqrt(x_1(:).^2 + x_2(:).^2);

R_vals = R(theta_flat);

inside_mask_grid = reshape(r_flat < R_vals, size(x_1));

masked_dd_grid = dd_grid;
masked_dd_grid(~inside_mask_grid) = NaN;

% For creating a denser boundary for plotting
xx = linspace(0, 2*pi, 1000);
z_bdd = boundary(R, Rp, xx);

%% Four potentially useful plots for visualization
figure;
contourf(x_1, x_2, masked_dd_grid, linspace(-5,5,21));
hold on;
plot(z_bdd(1,:), z_bdd(2,:), 'r');
scatter(z(1,:), z(2,:), 30, 'r', 'filled');
xlabel('x_{1}');
ylabel('x_{2}');
colorbar;
grid on;
axis equal;
title('Interior Laplace Dirichlet');
saveas(gcf, 'HW2_6c.png')

figure;
contourf(x_1, x_2, dd_grid, linspace(-5,5,21));
hold on;
plot(z_bdd(1,:), z_bdd(2,:), 'r');
scatter(z(1,:), z(2,:), 30, 'r', 'filled');
xlabel('x_{1}');
ylabel('x_{2}');
colorbar;
grid on;
axis equal;
title('Interior Laplace Dirichlet Unmasked');

figure;
surf(x_1, x_2, masked_dd_grid, 'EdgeColor', 'none');
view(3);
colorbar;
xlabel('x_{1}');
ylabel('x_{2}');
zlabel('u(x)');
title('Interior Laplace Dirichlet Surface');
grid on;
hold on;
z_boundary = zeros(1, size(z_bdd, 2));
plot3(z_bdd(1,:), z_bdd(2,:), z_boundary, 'r-', 'LineWidth', 1.5);
scatter3(z(1,:), z(2,:), boundary_data, 30, 'r', 'filled');

% Compute the error over the domain
solution = cos(x_1) .* exp(x_2);
error_grid = dd_grid - solution;
error_grid(~inside_mask_grid) = NaN;

figure;
contourf(x_1, x_2, log10(abs(error_grid)), 21);
hold on;
plot(z_bdd(1,:), z_bdd(2,:), 'r');
scatter(z(1,:), z(2,:), 30, 'r', 'filled');
xlabel('x_{1}');
ylabel('x_{2}');
colorbar;
grid on;
axis equal;
title('Interior Laplace Dirichlet Log10 Error');

%% Pointwise spectral convergence test, waiting until better results in the n=30 case
% 
% n_vals = 10:50;
% u = zeros(size(n_vals));
% x = [0.2, 0.1]';
% for l = 1:length(n_vals)
% 
%     n = n_vals(l);
% 
%     % Define boundary quantities
%     theta = linspace(0,2*pi,n+1);
%     theta(end) = [];
% 
%     z = boundary(R, Rp, theta);
%     zp = boundary_derivative(R, Rp, theta);
%     speed = vecnorm(zp, 2, 1);
%     n_y = normal(R, Rp, theta);
%     kappa = boundary_curvature(R, Rp, Rpp, theta);
% 
%     % Construct density
%     tau = solve_laplace_density(n, R, Rp, Rpp, f);
% 
%     x = [0.2; 0.1];
% 
%     dd = zeros(1, size(x, 2));
%     for i = 1:n
% 
%         dd = dd + Direct_Derivative_Fundamental_Laplace(x, z(:,i), n_y(:,i)) * tau(i) * speed(i) * 2 * pi / n;
% 
%     end
% 
%     u(l) = dd;
% 
% end
% 
% figure;
% loglog(n_vals, abs(u - cos(0.2)*exp(0.1)), 'ro-', 'DisplayName', 'Error at (0.2, 0.1)');
% legend;
% xlabel('N');
% ylabel('log10$|u - \\cos(0.2)\\exp(0.1)|$');
% grid on;
% title('Convergence of u(0.2, 0.1)');
% saveas(gcf, 'HW2_6c3.png');
