%% STARTUP

clc
clear
close all
addpath("../Functions/")


%% CHANGEABLE VARIABLES

%Settings for Newton-Rhapson
max_iteration = 50;
tolerance = 0.05;

%Operating Points
Pconu = 0;
Pconl = 0;
Pgrid = 500 * 1e6;
Qgrid = 0 * 1e6;
Vgrid_RE = 400 * 1e3;
Vgrid_IM = 00 * 1e3;
Vhvdc = 600 * 1e3;
Xarm_PU = 0.2;
Xl_PU = 0.02;
R_PU = 0.01;
Rl_PU = 0.01;
idcdif_ref = 0;
imiacsum_ref = 0;


%% PRE-ITERATION CALCULATIONS

%Combines imaginary and real components
Vgrid = Vgrid_RE + (Vgrid_IM * 1i);
Sgrid = Pgrid + (Qgrid * 1i);

%PU Conversion
Z_PU = abs(Vgrid)^2 / abs(Sgrid);
Xl = Xl_PU * Z_PU;
Xarm = Xarm_PU * Z_PU;
R = R_PU * Z_PU;
Rl = Rl_PU * Z_PU;


%% NEWTON-RHAPSON CALCULATION

in = ones(12,1) * 1000;
final = SinglePhase_TwoArm_Calc(in, max_iteration, tolerance, R, Rl, Xl, Xarm, Vhvdc, Vgrid_RE, Vgrid_IM, Pconu, Pconl, Pgrid, Qgrid, idcdif_ref, imiacsum_ref);


%% SEPARATE VARIABLES

vdcsum = final(1);
vdcdif = final(2); 
idcdif = final(3);  
idcsum = final(4);  
revacsum = final(5);   
imvacsum = final(6); 
revacdif = final(7);       
imvacdif = final(8);
reiacsum = final(9);      
imiacsum = final(10);    
reiacdif = final(11);  
imiacdif = final(12);    


%% OTHER CALCULATIONS FOR FINAL ITERATION

%Parses variables
vacsum = revacsum + (imvacsum * 1i);
vacdif = revacdif + (imvacdif * 1i);
iacdif = reiacdif + (imiacdif * 1i);
iacsum = reiacsum + (imiacsum * 1i);

%Finds the reactive power in the converter
Qconu = imag( ((vdcsum/2) * (idcsum - idcdif/2)) - (vacdif*conj(iacdif)/2) - (vacsum*conj(iacsum)) );
Qconl = imag( ((vdcsum/2) * (idcsum - idcdif/2)) - (vacdif*conj(iacdif)/2) + (vacsum*conj(iacsum)) );

%Finds the phases
phase_vacdif = rad2deg( angle(vacdif) );
phase_iacdif = rad2deg( angle(iacdif) );
phase_vacsum = rad2deg( angle(vacsum) );
phase_iacsum = rad2deg( angle(iacsum) );


%% DISPLAY OUTPUTS

%Displays results
fprintf('\nFINAL ITERATION RESULTS: \n')
disp(['VDC SUM = ' num2str(vdcsum/1e3) ' kV'])
disp(['VAC DIF = ' num2str(real(vacdif/1e3)) disp_sign(vacdif) num2str(abs(imag(vacdif/1e3))) 'i kV'])
disp(['IDC SUM = ' num2str(idcsum) ' A'])
disp(['IAC DIF = ' num2str(real(iacdif)) disp_sign(iacdif) num2str(abs(imag(iacdif))) 'i A'])
disp(['VAC SUM = ' num2str(real(vacsum/1e3)) disp_sign(vacsum) num2str(abs(imag(vacsum))/1e3') 'i kV'])
disp(['VDC DIF = ' num2str(vdcdif/1e3) ' kV'])
disp(['IAC SUM = ' num2str(real(iacsum)) disp_sign(iacsum) num2str(abs(imag(iacsum))) 'i A'])
disp(['IDC DIF = ' num2str(idcdif) ' A'])
fprintf('\nCALCULATED VALUES: \n')
disp(['QCON U = ' num2str(Qconu/1e6') ' MVAR'])
disp(['QCON L = ' num2str(Qconl/1e6) ' MVAR'])
disp(['VAC DIF Phase = ' num2str(phase_vacdif) '°'])
disp(['IAC DIF Phase = ' num2str(phase_iacdif) '°'])
disp(['VAC SUM Phase = ' num2str(phase_vacsum) '°'])
disp(['IAC SUM Phase = ' num2str(phase_iacsum) '°'])

%Creates the textbox message
msg_Pconu = ['Pconu = ' num2str(Pconu/1e6) ' MW'];
msg_Pconl = ['Pconl = ' num2str(Pconl/1e6) ' MW'];
msg_Sgrid = ['Sgrid = ' num2str(Sgrid/1e6) ' MVA'];
msg_Vgrid = ['Vgrid = ' num2str(Vgrid/1e3) ' kV'];
msg_Vhvdc = ['Vhvdc = ' num2str(Vhvdc/1e3) ' kV'];
msg_XarmPU = ['Xarm = ' num2str(Xarm_PU) ' PU'];
msg_RPU = ['R = ' num2str(R_PU) ' PU'];
msg_RlPU = ['Rl = ' num2str(Rl_PU) ' PU'];
msg_XlPU = ['Xl = ' num2str(Xl_PU) ' PU'];
msg_Idc = ['Idcdif ref = ' num2str(idcdif_ref) ' A'];
msg_Iac = ['Im(Iacsum) ref = ' num2str(imiacsum_ref) ' A'];
msg = {msg_Pconu msg_Pconl msg_Sgrid msg_Vgrid msg_Vhvdc msg_XarmPU msg_XlPU msg_RlPU msg_RPU msg_Idc msg_Iac};

plot_AC(vacdif, iacdif, 'Single Phase Two Arm Differential Values', [.2685 .13 .795 .795], msg)


%% DEBUG DATA

Vacu = -vacdif + vacsum/2;
Vdcu = -vdcdif + vdcsum/2;
Iacu = iacdif/2 + iacsum;
Idcu = idcdif/2 + idcsum;
Vacl = vacdif + vacsum/2;
Vdcl = vdcdif + vdcsum/2;
Iacl = -iacdif/2 + iacsum;
Idcl = -idcdif/2 + idcsum;