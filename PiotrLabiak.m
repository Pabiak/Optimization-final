R = 5;

fun = @(a) -1 * a * sqrt(R^2 - a^2); % Objective function
lb = 0; % Lower bound of the variable 'a'
ub = R; % Upper bound of the variable 'a'

function [x, fval] = calculate_minbnd(fun, lb, ub)
    tolerance = 1e-8; % Tolerance for convergence
    max_iter = 1000; % Maximum number of iterations

    % Golden ratio
    phi = (1 + sqrt(5)) / 2;

    % Initial points
    x1 = lb;
    x4 = ub;
    x2 = x4 - (x4 - x1) / phi;
    x3 = x1 + (x4 - x1) / phi;

    % Evaluate function values at initial points
    f1 = fun(x1);
    f2 = fun(x2);
    f3 = fun(x3);
    f4 = fun(x4);

    iter = 0;
    while abs(x4 - x1) > tolerance && iter < max_iter
        if f2 < f3
            x4 = x3;
            x3 = x2;
            x2 = x4 - (x4 - x1) / phi;
            f3 = f2;
            f2 = fun(x2);
        else
            x1 = x2;
            x2 = x3;
            x3 = x1 + (x4 - x1) / phi;
            f2 = f3;
            f3 = fun(x3);
        end
        iter = iter + 1;
    end

    % Choose the final solution
    if f2 < f3
        x = x2;
        fval = f2;
    else
        x = x3;
        fval = f3;
    end
end

[a, fval] = calculate_minbnd(fun, lb, ub);

% Solving b
b = sqrt(R^2 - a^2);

% Calculating the area
area = 4 * a * b;
disp('The largest area: ');
disp(area);

a = R * sqrt(2);
b = R * sqrt(2);

% plotting
theta = linspace(0, 2*pi, 100);
x_circle = R * cos(theta);
y_circle = R * sin(theta);

figure;
hold on;
axis equal;
xlabel('x');
ylabel('y');
title('Diagram');

plot(x_circle, y_circle, 'b');
plot([0, R], [0, 0], 'r--', 'LineWidth', 1);

rectangle('Position', [-a/2, -b/2, a, b], 'EdgeColor', 'r', 'FaceColor', 'none');

text(R/2, 0.2, sprintf('R = %.2f', R), 'Color', 'r', 'HorizontalAlignment', 'center');
text(-a/2, -b/2-0.2, sprintf('a = %.2f', a), 'Color', 'r', 'HorizontalAlignment', 'left');
text(-a/2-0.2, b/2, sprintf('b = %.2f', b), 'Color', 'r', 'Rotation', 90, 'VerticalAlignment', 'top');
