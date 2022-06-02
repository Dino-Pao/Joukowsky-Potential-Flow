% Script created and edited by DANILO PAOLINO, from the university of
% Naples Federico II. The script generate an airfoil using conformal
% mapping proposed by Joukowsy, ad evaluate the velocity field around it.
% The theory behind it was taken from the wikipedia page:
% https://en.wikipedia.org/wiki/Joukowsky_transform

clc
clear all
close all

% Properties of uniform flow 
Vinfty = 10;                % Flow velocity: [m/s]
alpha = deg2rad(+7);        % Angle of attack: [rad]
pinfty = 1.01325e5;         % Static pressure: [Pa]
rho = 1.18;                 % Density: [kg/m^3]

Ptot = pinfty + 0.5 * rho * Vinfty^2;       % Total pressure: [Pa]

% Next lines will define the real part, the imaginary part and the complex
% number which identify the center of the Joukowsky circle. This
% coordinates cannot equal zero, or the airfoil will collapse to a segment.

mu_x = -0.08;               % x–coordinate of the center: [m]
mu_y = +0.08;               % y–coordinate of the center: [m]
mu = mu_x + 1i*mu_y;

% The radius of the circle, on the complex plane, cannot be selected
% arbitrarily, but it has to satisfy the following condition:

R = sqrt( (1 - mu_x)^2 + mu_y^2);       % Radius: [m]

% The same concept applies for the circulation around the airfoil: in order
% to satisfy the Kutta condition (stagnation point coincident with the
% trailing edge) Gamma is evaluated with the following formula:

Gamma = 4*pi*Vinfty*R*sin(alpha + asin(mu_y/R));    % Circulation: [m^2/s]

% The first thing to do is to determine the mesh to use, to compute the
% airfoil geometry, as well as the flow field around it. It will be created
% varing the angle and the radius:
%
% theta:    angle, from 0 to 2*pi;
% r:        radius, from a minimum of R to 5*R.

theta = linspace(0, 2*pi,  200);        %[rad]
r = linspace(R, 5*R, 200);              %[m]

% Coordinates of the circle, on the complex plane, can be evaluated:
xi = mu_x + R*cos(theta);               %[m]
eta = mu_y + R*sin(theta);              %[m]
zeta = xi + 1i * eta;

% Applying the conformal mapping technique real coordinates of the profile,
% as two row arrays of x and y, can be evaluated:
x_profilo = real(zeta + 1./zeta);       % x array of airfoil: [m]
y_profilo = imag(zeta + 1./zeta);       % y array of airfoil: [m]


figure(1)
subplot(121)
plot(real(zeta), imag(zeta), 'k', 'LineWidth',1)
grid on
hold on
axis equal
xlabel('Real part Re$(z)$', 'Interpreter','latex')
ylabel('Imaginary part Im$(z)$', 'Interpreter','latex')
plot(mu_x, mu_y, 'k+', 'LineWidth',1)


subplot(122)
plot(x_profilo, y_profilo, 'k', 'LineWidth',1)
xlabel('Axis $x$', 'Interpreter','latex')
ylabel('Axis $y$', 'Interpreter','latex')
axis equal
grid on

%% This section will evaluate the flow field around the airfoil.

% The first thing to do is create the grid from the intersection of theta
% and r values, and traslate the result alongside the x–axis and the y–axis
% because of the center of the circonference:
[Radius, Theta] = meshgrid(r, theta);

A = mu_x + Radius.*cos(Theta);
B = mu_y + Radius.*sin(Theta);

Z = A + 1i*B;

% This is the velocity field in the complex plane
Wtilde = ...
    Vinfty * exp(-1i * alpha) + ...
    1i * Gamma ./ (2 * pi * (Z - mu) ) - ...
    Vinfty * R^2 * exp(1i * alpha) ./ ((Z - mu).^2) ;

figure(3)
subplot(121)
contourf(A, B, real(Wtilde), -20:20, 'ShowText','on')
hold on
grid on
axis equal
plot(xi, eta, 'k')
title('Complex plane: $u$ velocity', 'Interpreter','latex')
xlabel('Real part')
ylabel('Imaginary part')
view(2)
subplot(122)
contourf(A, B, imag(Wtilde), -20:20, 'ShowText','on')
hold on
grid on
axis equal
plot(xi, eta, 'k')
title('Complex plane: $v$ velocity', 'Interpreter','latex')
xlabel('Real part')
ylabel('Imaginary part')
view(2)

% The results gained from the previous parts are now mapped to the real
% {x,y}-plane.
W = Wtilde ./ (1 - Z.^(-2));

X = real(Z + 1./Z);
Y = imag(Z + 1./Z);

% This part is used to move and scale the airfoil, in order to set its
% minimum x value to 0, and its maximum to 1, like NACA airfoils.

min_index = find(x_profilo == min(x_profilo));
max_index = find(x_profilo == max(x_profilo));

scale_factor = x_profilo(max_index) - x_profilo(min_index);
X = X - x_profilo(min_index);
X = X/scale_factor;
Y = Y/scale_factor;

x_profilo_scaled = ( x_profilo - x_profilo(min_index) )/scale_factor;

y_profilo_scaled = y_profilo /scale_factor;

ux = real(W);
uy = - imag(W);

en_cin = rho * (ux.^2 + uy.^2) / 2;
p = Ptot - en_cin;

%% Plotting the results

figure(4)
contourf(X, Y, en_cin, 100, 'EdgeColor', 'none')
hold on
grid on
axis equal
plot(x_profilo_scaled, y_profilo_scaled, 'r')
axis([-.5 1.5 -.5 .5])
view(2)
legend('Kinetic energy contour')
xlabel('$x$ axis', 'Interpreter','latex')
ylabel('$y$ axis', 'Interpreter','latex')

% The following figure shows the pressure contours close to the airfoil
figure(5)
contourf(X, Y, p, 100, 'EdgeColor', 'k')
grid on
axis equal
hold on
plot(x_profilo_scaled, y_profilo_scaled, 'r')
axis([-.5 1.5 -.5 .5])
view(2)
legend('Static pressure contour')
xlabel('$x$ axis', 'Interpreter','latex')
ylabel('$y$ axis', 'Interpreter','latex')


% The following figure plots the vector map around the airfoil.
figure(6)
scale = 0.5;
quiver(X, Y, ux, uy, scale)
hold on
axis equal
plot(x_profilo_scaled, y_profilo_scaled, 'r')
axis([-0.5 1.5 -.5 0.5])
view(2)
xlabel('$x$ axis', 'Interpreter','latex')
ylabel('$y$ axis', 'Interpreter','latex')

%% XFoil setup
% The last part of the script will save the airfoil geometry to a .txt file
% to import in a session of XFoil, in order to validate the results.

x_profilo_XFoil = ...
    [x_profilo(max_index:end), x_profilo(2:max_index)] - x_profilo(min_index);
x_profilo_XFoil = x_profilo_XFoil/scale_factor;

y_profilo_XFoil = ...
    [y_profilo(max_index:end), y_profilo(2:max_index)];
y_profilo_XFoil = y_profilo_XFoil/scale_factor;

figure(9)
plot(x_profilo_XFoil, y_profilo_XFoil, 'k')
hold on
grid on
axis equal
scatter(x_profilo_XFoil(1), y_profilo_XFoil(1), 'r')
scatter(x_profilo_XFoil(end), y_profilo_XFoil(end), 'b')
legend('Airfoil', 'TE First point', 'TE Last Point')

airfoil_XFoil = [x_profilo_XFoil', y_profilo_XFoil'];

writematrix(airfoil_XFoil, 'airfoil_sim0000.txt','Delimiter',' ')

% These are the values of the pressure coefficient, get from an XFoil run.
% As we can see they validate the script.

XFoil_data = readmatrix("airfoil_sim0000.cp", 'FileType','text');
XFoil_x_data = XFoil_data(:,1);
XFoil_y_data = XFoil_data(:,2);
XFoil_cp_data = XFoil_data(:,3);

% Pressure coefficient plot, evaluated on the border of the profile
figure(8)
cp = (p(:,1) - pinfty) / (0.5 * rho * Vinfty^2);
plot(X(:, 1), cp)
hold on
plot(XFoil_x_data, XFoil_cp_data, 'o')
grid on
xlabel('Axis $x$ [m]', 'Interpreter','latex')
ylabel('Pressure coefficient $c_p = \frac{p - p_\infty}{\frac{1}{2}\rho c_\infty^2}$', 'Interpreter','latex')
