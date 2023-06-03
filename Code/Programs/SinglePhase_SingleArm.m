%% STARTUP

clc
clear
close all
addpath("../Functions/")


%% CHANGEABLE VARIABLES

%Settings for Newton-Rhapson
max_iterations = 25;
tolerance = 0.01;

%Operating Points
Pcon = 0;
Pgrid = 640 * 1e6;
Qgrid = 000 * 1e6;
Vgrid_RE = -525 * 1e3;
Vgrid_IM = -250 * 1e3;
Vhvdc = 800 * 1e3;
Xarm_PU = 0.04;
R_PU = 0.02;
Rarm_PU = 0.01;


%% PRE-ITERATION CALCULATIONS

%Combines imaginary and real components
Vgrid = Vgrid_RE + (Vgrid_IM * 1i);
Sgrid = Pgrid + (Qgrid * 1i);

%PU Conversion
Z_PUBase = abs(Vgrid)^2 / abs(Sgrid);
Xarm = Xarm_PU * Z_PUBase;
R = R_PU * Z_PUBase;
Rarm = Rarm_PU * Z_PUBase;


%% NEWTON-RHAPSON CALCULATION

%Initial matrix to Solve newton-Raphson with
in = ones(6,1) * 1000;
final = SinglePhase_SingleArm_Calc(in, max_iterations, tolerance, Pcon, Xarm, R, Rarm, Vgrid_RE, Vgrid_IM, Vhvdc, Pgrid, Qgrid);


%% OTHER CALCULATIONS FOR FINAL ITERATION

%Finds combined Vac/Iac values
Vac = final(1) + (final(2)*1i);
Iac = final(3) + (final(4)*1i);
Vdc = final(5);
Idc = final(6);

%Finds the phase of Vac/Iac
phase_vac = rad2deg( angle(Vac) );
phase_iac = rad2deg( angle(Iac) );
phase_dif = phase_vac - phase_iac;

%Finds the reactive power in the converter
Scon = (Vac * conj(Iac)) + (Vdc * Idc);
Qcon = Scon - Pcon;


%% RESULTS

%Displays results
fprintf('\nFINAL ITERATION RESULTS: \n')
disp(['Vac = ' num2str(real(Vac), '%3.3e') disp_sign(Vac) num2str(abs(imag(Vac)), '%3.3e') 'i'])
disp(['Iac = ' num2str(real(Iac), '%3.3e') disp_sign(Iac) num2str(abs(imag(Iac)), '%3.3e') 'i'])
disp(['Vdc = ' num2str(Vdc, '%3.3e')])
disp(['Idc = ' num2str(Idc, '%3.3e')])
fprintf('\nCALCULATED RESULTS: \n')
disp(['Vac Phase = ' num2str(phase_vac) '°'])
disp(['Iac Phase = ' num2str(phase_iac) '°'])
disp(['Phase Difference = ' num2str(phase_dif) '°'])
disp(['Converter Reactive Power = ' num2str(imag(Qcon), '%3.3e')])
fprintf('\n')


%% PLOTS

%Creates the textbox message
msg_Pcon = ['Pcon = ' num2str(Pcon, '%.2e')];
msg_Sgrid = ['Sgrid = ' num2str(Sgrid, '%.2e')];
msg_Vgrid = ['Vgrid = ' num2str(Vgrid, '%.2e')];
msg_Vhvdc = ['Vhvdc = ' num2str(Vhvdc, '%.2e')];
msg_XarmPU = ['Xarm PU = ' num2str(Xarm_PU)];
msg_RPU = ['R PU = ' num2str(R_PU)];
msg_RarmPU = ['Rarm PU = ' num2str(Rarm_PU)];
msg = {msg_Pcon msg_Sgrid msg_Vgrid msg_Vhvdc msg_XarmPU msg_RarmPU msg_RPU};

%AC Phasor Plot
plot_AC(Vac, Iac, 'Single Phase Single Arm', [.131 .131 .795 .795], msg)