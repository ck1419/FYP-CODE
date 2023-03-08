%% STARTUP

clc
clear
close all


%% INITIAL VARIABLES

%Settings for Newton-Rhapson
iterations = 25;
tolerance = 0.01;   %Used to check results
variable_count = 6;

%Operating Points
Pcon = 0;
Pgrid = 500 * 1e6;
Qgrid = 100 * 1e6;
Vgrid_RE = 400 * 1e3;
Vgrid_IM = 100 * 1e3;
Vhvdc = 200 * 1e3;
Xarm_PU = 0.15;
R_PU = 0.01;

%Converter Limits
rated_voltage = 600e3;
rated_current = 100e3;
rated_power = 2000e6;


%% PRE-ITERATION CALCULATIONS

%Combines imaginary and real components
Vgrid = Vgrid_RE + (Vgrid_IM * 1i);
Sgrid = Pgrid + (Qgrid * 1i);

%PU Conversion
Z_PUBase = abs(Vgrid)^2 / abs(Sgrid);
Xarm = Xarm_PU * Z_PUBase;
R = R_PU * Z_PUBase;


%% NEWTON-RHAPSON CALCULATION

%Initial matrix to Solve newton-Raphson with
x = zeros(variable_count, iterations);
x(:,1) = ones(variable_count,1) * 100;

%Loop to execute Newton-Raphson
for n = 2:iterations
    f11_value = f11(x(:,n-1), Pcon, Xarm, R, Vgrid_RE, Vgrid_IM, Vhvdc, Pgrid, Qgrid);
    f11_delta_value = f11_delta(x(:,n-1), Pcon, Xarm, R, Vgrid_RE, Vgrid_IM, Vhvdc, Pgrid, Qgrid);
    x(:,n) = x(:,n-1) - (f11_delta_value^-1 * f11_value);
end


%% OTHER CALCULATIONS FOR FINAL ITERATION

%Finds combined Vac/Iac values
Vac = x(1, iterations) + (x(2, iterations)*1i);
Iac = x(3, iterations) + (x(4, iterations)*1i);
Vdc = x(5, iterations);
Idc = x(6, iterations);

%Finds the phase of Vac/Iac
phase_vac = rad2deg( angle(Vac) );
phase_iac = rad2deg( angle(Iac) );
phase_dif = phase_vac - phase_iac;

%Finds the reactive power in the converter
Scon = (Vac * conj(Iac)) + (Vdc * Idc);
Qcon = Scon - Pcon;


%% RESULTS

fprintf('CALCULATION CHECK: \n')

%Check for f(x) = 0
pass_eqn = f11_eqn_check(x(:,iterations), Pcon, Xarm, R, Vgrid_RE, Vgrid_IM, Vhvdc, Pgrid, Qgrid, tolerance);
%Check for no imaginary components of each variable
pass_im = check_im(x(:,iterations));
%Check for power flow direction
pass_flow = check_powerflow(Vdc, Idc);

fprintf('\nCONVERTER PHYSICAL LIMITS CHECK: \n')

%Voltage limit check
pass_voltage_lim = check_voltage_limit(Vac,Vdc,rated_voltage);
%Current limit check
pass_current_lim = check_current_limit(Iac,Idc,rated_current);
%Power limit check
pass_power_lim = check_power_limit(Vgrid,Iac,rated_power);

%Displays results
f11_results_display(Vac, Iac, Vdc, Idc, phase_vac, phase_iac, phase_dif, Qcon)


%% PLOTS

msg_Pcon = ['Pcon = ' num2str(Pcon, '%.2e')];
msg_Sgrid = ['Sgrid = ' num2str(Sgrid, '%.2e')];
msg_Vgrid = ['Vgrid = ' num2str(Vgrid, '%.2e')];
msg_Vhvdc = ['Vhvdc = ' num2str(Vhvdc, '%.2e')];
msg_XarmPU = ['Xarm PU = ' num2str(Xarm_PU)];
msg_RPU = ['R PU = ' num2str(R_PU)];

msg = {msg_Pcon msg_Sgrid msg_Vgrid msg_Vhvdc msg_XarmPU msg_RPU};

%AC Phasor Plot
% plot_AC(Vac, Iac, msg)