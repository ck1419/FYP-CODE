%% STARTUP

clc
clear
close all


%% INITIAL VARIABLES

%Settings for Newton-Rhapson
iterations = 50;
tolerance = 0.05;   %Tolerance percentage in difference between iterations for final answer
variable_count = 12;

%Operating Points
Pconu = 0;
Pconl = 0;
Pgrid = 500 * 1e6;
Qgrid = 100 * 1e6;
Vgrid_RE = 400 * 1e3;
Vgrid_IM = 100 * 1e3;
Vhvdc = 600 * 1e3;
Xarm_PU = 0.015;
Xl_PU = 0.02;
R_PU = 0.01;
Rl_PU = 0.01;

idcdif_ref = 0;
imiacsum_ref = 0;


%% PRE-ITERATION CALCULATIONS

%Allows system to converge faster without affecting results
idcdif_ref_temp = idcdif_ref;
if idcdif_ref == 0
    idcdif_ref_temp = 1e-9;
end

imiacsum_ref_temp = imiacsum_ref;
if imiacsum_ref == 0
    imiacsum_ref_temp = 1e-9;
end

Vgrid_IM_temp = Vgrid_IM;
if Vgrid_IM == 0
    Vgrid_IM_temp = 1e-9;
end
    
Qgrid_temp = Qgrid;
if Qgrid == 0
    Qgrid_temp = 1e-9;
end

%Combines imaginary and real components
Vgrid = Vgrid_RE + (Vgrid_IM_temp * 1i);
Sgrid = Pgrid + (Qgrid_temp * 1i);
Igrid = abs(Sgrid)/abs(Vgrid);

%PU Conversion
Z_PU = abs(Vgrid)^2 / abs(Sgrid);
Xl = Xl_PU * Z_PU;
Xarm = Xarm_PU * Z_PU;
R = R_PU * Z_PU;
Rl = Rl_PU * Z_PU;


%% NEWTON-RHAPSON CALCULATION

%Initial matrix to Solve newton-Raphson with
x = zeros(variable_count, iterations);
x(:,1) = ones(variable_count,1) * 100;

%Loop to execute Newton-Raphson
for n = 2:iterations
    f12_value = f12(x(:,n-1), R, Rl, Xl, Xarm, Vhvdc, Vgrid, Pconu, Pconl, Sgrid, idcdif_ref_temp, imiacsum_ref_temp);
    f12_delta_value = f12_delta(x(:,n-1), R, Rl, Xl, Xarm, Vgrid);
    x(:,n) = x(:,n-1) - (f12_delta_value^-1 * f12_value);
    if all((x(:,n)./x(:,n-1)) <= 1+tolerance) && all((x(:,n)./x(:,n-1)) >= 1-tolerance)
        iterated = n;
        final = x(:,n);
        break
    end
    if n == iterations
        final = x(:,n);
        fprintf('WARNING: MAX ITERATIONS REACHED \n')
    end
end


%% CLEAN UP VARIABLES

for n = 1:variable_count
    if abs(final(n)) <= 1
        final(n) = 0;
    end
end


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


%% CALCULATE VARIABLES

vacsum = revacsum + (imvacsum * 1i);
vacdif = revacdif + (imvacdif * 1i);
iacdif = reiacdif + (imiacdif * 1i);
iacsum = reiacsum + (imiacsum * 1i);

Qconu = imag( ((vdcsum/2) * (idcsum - idcdif/2)) - (vacdif*conj(iacdif)/2) - (vacsum*conj(iacsum)) );
Qconl = imag( ((vdcsum/2) * (idcsum - idcdif/2)) - (vacdif*conj(iacdif)/2) + (vacsum*conj(iacsum)) );


%% DISPLAY OUTPUTS

fprintf(['ITERATIONS = ' num2str(iterated) '\n\n'])

fprintf('\nFINAL ITERATION RESULTS: \n')
disp(['VDC SUM = ' num2str(vdcsum, '%3.3e')])
disp(['VAC DIF = ' num2str(real(vacdif), '%3.3e') disp_sign(vacdif) num2str(abs(imag(vacdif)), '%3.3e') 'i'])
disp(['IDC SUM = ' num2str(idcsum, '%3.3e')])
disp(['IAC DIF = ' num2str(real(iacdif), '%3.3e') disp_sign(iacdif) num2str(abs(imag(iacdif)), '%3.3e') 'i'])

disp(['VAC SUM = ' num2str(real(vacsum), '%3.3e') disp_sign(vacsum) num2str(abs(imag(vacsum)), '%3.3e') 'i'])
disp(['VDC DIF = ' num2str(vdcdif, '%3.3e')])
disp(['IAC SUM = ' num2str(real(iacsum), '%3.3e') disp_sign(iacsum) num2str(abs(imag(iacsum)), '%3.3e') 'i'])
disp(['IDC DIF = ' num2str(idcdif, '%3.3e')])

fprintf('\nCALCULATED VALUES: \n')
disp(['QCON U = ' num2str(Qconu, '%3.3e') 'i'])
disp(['QCON L = ' num2str(Qconl, '%3.3e') 'i'])

msg_Pconu = ['Pconu = ' num2str(Pconu, '%.2e')];
msg_Pconl = ['Pconu = ' num2str(Pconl, '%.2e')];
msg_Sgrid = ['Sgrid = ' num2str(Sgrid, '%.2e')];
msg_Vgrid = ['Vgrid = ' num2str(Vgrid, '%.2e')];
msg_Vhvdc = ['Vhvdc = ' num2str(Vhvdc, '%.2e')];
msg_XarmPU = ['Xarm PU = ' num2str(Xarm_PU)];
msg_RPU = ['R PU = ' num2str(R_PU)];
msg_RlPU = ['Rl PU = ' num2str(Rl_PU)];
msg_XlPU = ['Xl PU = ' num2str(Xl_PU)];
msg_Idc = ['Idcdif ref = ' num2str(idcdif_ref)];
msg_Iac = ['Iacsum ref = ' num2str(imiacsum_ref)];

msg = {msg_Pconu msg_Pconl msg_Sgrid msg_Vgrid msg_Vhvdc msg_XarmPU msg_XlPU msg_RlPU msg_RPU msg_Idc msg_Iac};

plot_AC(vacdif, iacdif, 'Single Phase Two Arm Differential Values', [.565 .2895 .565 .286], msg)

abs(iacdif/2)*sqrt(2) + abs(idcsum)