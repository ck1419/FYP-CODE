%% STARTUP

clc
clear
close all


%% INITIAL VARIABLES

%Settings for Newton-Rhapson
iterations = 25;
tolerance = 0.01;   %Used to check results
variable_count = 12;

%Operating Points
Pconu = 0;
Pconl = 0;
Pgrid = 500 * 1e6;
Qgrid = 100 * 1e6;
Vgrid_RE = 400 * 1e3;
Vgrid_IM = 100 * 1e3;
Vhvdc = 200 * 1e3;
Xarm = 2;
Xl = 2;
R = 2;


%% PRE-ITERATION CALCULATIONS

%Combines imaginary and real components
Vgrid = Vgrid_RE + (Vgrid_IM * 1i);
Sgrid = Pgrid + (Qgrid * 1i);


%% NEWTON-RHAPSON CALCULATION

%Initial matrix to Solve newton-Raphson with
x = zeros(variable_count, iterations);
x(:,1) = ones(variable_count,1) * 100;

%Loop to execute Newton-Raphson
for n = 2:iterations
    f12_value = f12(x(:,n-1), R, Xl, Xarm, Vhvdc, Vgrid, Pconu, Pconl, Pgrid, 1, 1, 1);
    f12_delta_value = f12_delta(x(:,n-1), R, Xl, Xarm, Vhvdc, Vgrid, Pconu, Pconl, Pgrid, 1, 1, 1);
    x(:,n) = x(:,n-1) - (f12_delta_value^-1 * f12_value);
end

