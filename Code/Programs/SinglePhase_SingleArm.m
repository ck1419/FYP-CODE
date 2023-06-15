%% STARTUP

clc
clear
close all
addpath("../Functions/")


%% CHANGEABLE VARIABLES

%Settings for Newton-Rhapson
max_iterations = 50;
tolerance = 0.01;

%Operating Points
Pcon = 0;
Pgrid = 640/3 * 1e6;
Qgrid = 000 * 1e6;
Vgrid_RE = 220 * 1e3;
Vgrid_IM = 0 * 1e3;
Vhvdc = 800 * 1e3;
Xarm_PU = 0.168;
R_PU = 0.005;
Rarm_PU = 0.005;


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

%Finds combined Vac/Iac values
Vac = final(1) + (final(2)*1i);
Iac = final(3) + (final(4)*1i);
Vdc = final(5);
Idc = final(6);



%% OTHER CALCULATIONS FOR FINAL ITERATION

%Finds the phase of Vac/Iac
phase_vac = rad2deg( angle(Vac) );
phase_iac = rad2deg( angle(Iac) );
phase_dif = phase_vac - phase_iac;

%Finds the reactive power in the converter
Qcon = imag(Vac * conj(Iac));


%% RESULTS

%Displays results
fprintf('\nFINAL ITERATION RESULTS: \n')
disp(['VAC = ' num2str(real(Vac)/1e3) disp_sign(Vac) num2str(abs(imag(Vac))/1e3) 'i kV'])
disp(['IAC = ' num2str(real(Iac)) disp_sign(Iac) num2str(abs(imag(Iac))) 'i A'])
disp(['VDC = ' num2str(Vdc/1e3) ' kV'])
disp(['IDC = ' num2str(Idc) ' A'])
fprintf('\nCALCULATED RESULTS: \n')
disp(['QCON = ' num2str(Qcon/1e6') ' MVAR'])
disp(['VAC Phase = ' num2str(phase_vac) '°'])
disp(['IAC Phase = ' num2str(phase_iac) '°'])
fprintf('\n')


%% PLOTS

%Creates the textbox message
msg_Pcon = ['Pcon = ' num2str(Pcon/1e6) ' MW'];
msg_Sgrid = ['Sgrid = ' num2str(Sgrid/1e6) ' MVA'];
msg_Vgrid = ['Vgrid = ' num2str(Vgrid/1e3) ' kV'];
msg_Vhvdc = ['Vhvdc = ' num2str(Vhvdc/1e3) ' kV'];
msg_XarmPU = ['Xarm = ' num2str(Xarm_PU) ' PU'];
msg_RPU = ['R = ' num2str(R_PU) ' PU'];
msg_RarmPU = ['Rarm = ' num2str(Rarm_PU) ' PU'];
msg = {msg_Pcon msg_Sgrid msg_Vgrid msg_Vhvdc msg_XarmPU msg_RarmPU msg_RPU};

%AC Phasor Plot
plot_AC(Vac, Iac, 'Single Phase Single Arm', [.131 .131 .795 .795], msg)