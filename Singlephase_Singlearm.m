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
Pgrid = 100 * 1e6;
Qgrid = 300 * 1e6;
Vgrid_RE = 400 * 1e3;
Vgrid_IM = 100 * 1e3;
Vhvdc = 600 * 1e3;
Xarm_PU = 0.2;
R_PU = 0.05;
Rarm_PU = 0.1;

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
Rarm = Rarm_PU * Z_PUBase;

%For runs where Rarm and Xarm needs to remain the same
Xarm = 100;
Rarm = 50;
R = 25;

%% NEWTON-RHAPSON CALCULATION

%Initial matrix to Solve newton-Raphson with
x = zeros(variable_count, iterations);
x(:,1) = ones(variable_count,1) * 100;

%Loop to execute Newton-Raphson
for n = 2:iterations
    f11_value = f11(x(:,n-1), Pcon, Xarm, R, Rarm, Vgrid_RE, Vgrid_IM, Vhvdc, Pgrid, Qgrid);
    f11_delta_value = f11_delta(x(:,n-1), Xarm, R, Rarm, Vgrid_RE, Vgrid_IM);
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
plot_AC(Vac, Iac, 'Single Phase Single Arm', msg)