%% STARTUP

clc
clear
close all
addpath("../Functions/")


%% CHANGEABLE VARIABLES

%Settings for Newton-Rhapson
max_iteration = 50;
tolerance = 0.05;   %Tolerance percentage in difference between iterations for final answer

%Global Operating Points
Pgrid = 500 * 1e6;
Qgrid = 0 * 1e6;
Vgrid_RE = 400 * 1e3;
Vgrid_IM = 0 * 1e3;
Vhvdc = 600 * 1e3;
Xarm_PU = 0.015;
Xl_PU = 0.02;
R_PU = 0.01;
Rl_PU = 0.01;

%Phase A Operating Points
Pconu_a = 0;
Pconl_a = 0;
idcdif_ref_a = 0;
imiacsum_ref_a = 0;

%Phase B Operating Points
Pconu_b = 0;
Pconl_b = 0;
idcdif_ref_b = 0;
imiacsum_ref_b = 0;
b_phase = 120;  %Relative to A in degrees

%Phase C Operating Points
Pconu_c = 0;
Pconl_c = 0;
idcdif_ref_c = 0;
imiacsum_ref_c = 0;
c_phase = -120; %Relative to A in degrees


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

%Combines variables into groups
Pcon_array = [Pconu_a, Pconl_a, Pconu_b, Pconl_b, Pconu_c, Pconl_c];
Ref_array = [idcdif_ref_a, imiacsum_ref_a, idcdif_ref_b, imiacsum_ref_b, idcdif_ref_c, imiacsum_ref_c];


%% NEWTON-RHAPSON CALCULATION

in = ones(36,1) * 1000;
final = ThreePhase_SixArm_Calc(in, max_iteration, tolerance, R, Rl, Xl, Xarm, Vhvdc, Vgrid_RE, Vgrid_IM, Pcon_array, Pgrid, Qgrid, Ref_array, b_phase, c_phase);


%% SEPARATE VARIABLES

vdcsum_a = final(1);
vdcdif_a = final(2); 
idcdif_a = final(3);  
idcsum_a = final(4);  
revacsum_a = final(5);   
imvacsum_a = final(6); 
revacdif_a = final(7);       
imvacdif_a = final(8);
reiacsum_a = final(9);      
imiacsum_a = final(10);    
reiacdif_a = final(11);  
imiacdif_a = final(12);   
vdcsum_b = final(13);
vdcdif_b = final(14); 
idcdif_b = final(15);  
idcsum_b = final(16);  
revacsum_b = final(17);   
imvacsum_b = final(18); 
revacdif_b = final(19);       
imvacdif_b = final(20);
reiacsum_b = final(21);      
imiacsum_b = final(22);    
reiacdif_b = final(23);  
imiacdif_b = final(24);  
vdcsum_c = final(25);
vdcdif_c = final(26); 
idcdif_c = final(27);  
idcsum_c = final(28);  
revacsum_c = final(29);   
imvacsum_c = final(30); 
revacdif_c = final(31);       
imvacdif_c = final(32);
reiacsum_c = final(33);      
imiacsum_c = final(34);    
reiacdif_c = final(35);  
imiacdif_c = final(36);  


%% CALCULATE VARIABLES

%Phase A
vacsum_a = revacsum_a + (imvacsum_a * 1i);
vacdif_a = revacdif_a + (imvacdif_a * 1i);
iacdif_a = reiacdif_a + (imiacdif_a * 1i);
iacsum_a = reiacsum_a + (imiacsum_a * 1i);

%Phase B
vgrid_b = Vgrid * exp(1i*deg2rad(b_phase));
vacsum_b = revacsum_b + (imvacsum_b * 1i);
vacdif_b = revacdif_b + (imvacdif_b * 1i);
iacdif_b = reiacdif_b + (imiacdif_b * 1i);
iacsum_b = reiacsum_b + (imiacsum_b * 1i);

%Phase C
vgrid_c = Vgrid * exp(1i*deg2rad(c_phase));
vacsum_c = revacsum_c + (imvacsum_c * 1i);
vacdif_c = revacdif_c + (imvacdif_c * 1i);
iacdif_c = reiacdif_c + (imiacdif_c * 1i);
iacsum_c = reiacsum_c + (imiacsum_c * 1i);


%% DISPLAY OUTPUTS

%Displays results
fprintf('\nPHASE A: \n')
disp(['VGRID = ' num2str(real(Vgrid), '%3.3e') disp_sign(Vgrid) num2str(abs(imag(Vgrid)), '%3.3e') 'i'])
disp(['VDC SUM = ' num2str(vdcsum_a, '%3.3e')])
disp(['VAC DIF = ' num2str(real(vacdif_a), '%3.3e') disp_sign(vacdif_a) num2str(abs(imag(vacdif_a)), '%3.3e') 'i'])
disp(['IDC SUM = ' num2str(idcsum_a, '%3.3e')])
disp(['IAC DIF = ' num2str(real(iacdif_a), '%3.3e') disp_sign(iacdif_a) num2str(abs(imag(iacdif_a)), '%3.3e') 'i'])
disp(['VAC SUM = ' num2str(real(vacsum_a), '%3.3e') disp_sign(vacsum_a) num2str(abs(imag(vacsum_a)), '%3.3e') 'i'])
disp(['VDC DIF = ' num2str(vdcdif_a, '%3.3e')])
disp(['IAC SUM = ' num2str(real(iacsum_a), '%3.3e') disp_sign(iacsum_a) num2str(abs(imag(iacsum_a)), '%3.3e') 'i'])
disp(['IDC DIF = ' num2str(idcdif_a, '%3.3e')])

fprintf('\nPHASE B: \n')
disp(['VGRID = ' num2str(real(vgrid_b), '%3.3e') disp_sign(vgrid_b) num2str(abs(imag(vgrid_b)), '%3.3e') 'i'])
disp(['VDC SUM = ' num2str(vdcsum_b, '%3.3e')])
disp(['VAC DIF = ' num2str(real(vacdif_b), '%3.3e') disp_sign(vacdif_b) num2str(abs(imag(vacdif_b)), '%3.3e') 'i'])
disp(['IDC SUM = ' num2str(idcsum_b, '%3.3e')])
disp(['IAC DIF = ' num2str(real(iacdif_b), '%3.3e') disp_sign(iacdif_b) num2str(abs(imag(iacdif_b)), '%3.3e') 'i'])
disp(['VAC SUM = ' num2str(real(vacsum_b), '%3.3e') disp_sign(vacsum_b) num2str(abs(imag(vacsum_b)), '%3.3e') 'i'])
disp(['VDC DIF = ' num2str(vdcdif_b, '%3.3e')])
disp(['IAC SUM = ' num2str(real(iacsum_b), '%3.3e') disp_sign(iacsum_b) num2str(abs(imag(iacsum_b)), '%3.3e') 'i'])
disp(['IDC DIF = ' num2str(idcdif_b, '%3.3e')])

fprintf('\nPHASE C: \n')
disp(['VGRID = ' num2str(real(vgrid_c), '%3.3e') disp_sign(vgrid_c) num2str(abs(imag(vgrid_c)), '%3.3e') 'i'])
disp(['VDC SUM = ' num2str(vdcsum_c, '%3.3e')])
disp(['VAC DIF = ' num2str(real(vacdif_c), '%3.3e') disp_sign(vacdif_c) num2str(abs(imag(vacdif_c)), '%3.3e') 'i'])
disp(['IDC SUM = ' num2str(idcsum_c, '%3.3e')])
disp(['IAC DIF = ' num2str(real(iacdif_c), '%3.3e') disp_sign(iacdif_c) num2str(abs(imag(iacdif_c)), '%3.3e') 'i'])
disp(['VAC SUM = ' num2str(real(vacsum_c), '%3.3e') disp_sign(vacsum_c) num2str(abs(imag(vacsum_c)), '%3.3e') 'i'])
disp(['VDC DIF = ' num2str(vdcdif_c, '%3.3e')])
disp(['IAC SUM = ' num2str(real(iacsum_c), '%3.3e') disp_sign(iacsum_c) num2str(abs(imag(iacsum_c)), '%3.3e') 'i'])
disp(['IDC DIF = ' num2str(idcdif_c, '%3.3e')])

