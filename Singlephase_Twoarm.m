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
Qgrid = 500 * 1e6;
Vgrid_RE = 400 * 1e3;
Vgrid_IM = 100 * 1e3;
Vhvdc = 600 * 1e3;
Xarm_PU = 0.2;
Xl_PU = 0.4;
R_PU = 0.05;
Rl_PU = 0.1;

idcdif_ref = 1e-3;
reiacsum_ref = 1e-3;


%% PRE-ITERATION CALCULATIONS

%Combines imaginary and real components
Vgrid = Vgrid_RE + (Vgrid_IM * 1i);
Sgrid = Pgrid + (Qgrid * 1i);
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
    f12_value = f12(x(:,n-1), R, Rl, Xl, Xarm, Vhvdc, Vgrid, Pconu, Pconl, Sgrid, idcdif_ref, reiacsum_ref);
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

%This allows us to relax the tolerances before the values properly converges to 0
for n = 1:variable_count
    if n <= 2 && abs(final(n)) <= Vgrid*tolerance/2
        final(n) = 0;
    elseif n > 2 && n <= 4 && abs(final(n)) <= Igrid*tolerance/2
        final(n) = 0;
    elseif n > 4 && n <= 8 && abs(final(n)) <= Vgrid*tolerance/2
        final(n) = 0;
    elseif n > 8 && abs(final(n)) <= Igrid*tolerance/2
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

plot_AC(vacdif, iacdif, 'AC Grid Values', 'temp')

abs(iacdif/2)*sqrt(2) + abs(idcsum)