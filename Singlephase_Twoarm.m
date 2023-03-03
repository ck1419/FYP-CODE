%% STARTUP

clc
clear
close all


%% INITIAL VARIABLES

%Settings for Newton-Rhapson
iterations = 25;
tolerance = 1;   %Used to check results
variable_count = 12;

%Operating Points
Pconu = 0;
Pconl = 0;
Pgrid = 500 * 1e6;
Qgrid = 100 * 1e6;
Vgrid_RE = 400 * 1e3;
Vgrid_IM = 100 * 1e3;
Vhvdc = 200 * 1e3;
Xarm_PU = 0.2;
Xl_PU = 0.15;
R_PU = 0.05;

idcac_ref = 1e-3;
reiacdc_ref = 1e-3;


%% PRE-ITERATION CALCULATIONS

%Combines imaginary and real components
Vgrid = Vgrid_RE + (Vgrid_IM * 1i);
Sgrid = Pgrid + (Qgrid * 1i);

%PU Conversion
Z_PU = abs(Vgrid)^2 / abs(Sgrid);
Xl = Xl_PU * Z_PU;
Xarm = Xarm_PU * Z_PU;
R = R_PU * Z_PU;


%% NEWTON-RHAPSON CALCULATION

%Initial matrix to Solve newton-Raphson with
x = zeros(variable_count, iterations);
x(:,1) = ones(variable_count,1) * 100;

%Loop to execute Newton-Raphson
for n = 2:iterations
    f12_value = f12(x(:,n-1), R, Xl, Xarm, Vhvdc, Vgrid, Pconu, Pconl, Sgrid, idcac_ref, reiacdc_ref);
    f12_delta_value = f12_delta(x(:,n-1), R, Xl, Xarm, Vhvdc, Vgrid, Pconu, Pconl, Sgrid, idcac_ref, reiacdc_ref);
    x(:,n) = x(:,n-1) - (f12_delta_value^-1 * f12_value);
    if (abs(f12_value)) <= tolerance
        final = x(:,n);
        iterated = n;
        break
    end
    if n == iterations
        final = x(:,n);
        fprintf('WARNING: MAX ITERATIONS REACHED \n')
    end
end


%% SEPARATE VARIABLES

%Rounds any magnitude less than 0.5 to the closest integer
%Helps to clean up values that hasnt converged properly yet
for n = 1:variable_count
    if abs(final(n)) <= 0.5
        final(n) = round(final(n));
    end
end

vdcsum = final(1);
vdcdif = final(2);         % 0
revacsum = final(3);       % 0
imvacsum = final(4);       % 0
revacdif = final(5);       
imvacdif = final(6);
reiacac = final(7);        %Iac of AC grid
imiacac = final(8);       
reiacdc = final(9);        %Iac of DC grid - 0
imiacdc = final(10);       % 0
idcac = final(11);         %Idc of AC grid - 0
idcdc = final(12);         %Idc of DC grid

vacsum = revacsum + (imvacsum * 1i);
vacdif = revacdif + (imvacdif * 1i);
iacac = reiacac + (imiacac * 1i);
iacdc = reiacdc + (imiacdc * 1i);


%% CALCULATE VARIABLES

Qconu = imag( ((vdcsum/2) * (idcdc + idcac/2)) - (vacdif*conj(iacac)/2) - (vacsum*conj(iacdc)) );
Qconl = imag( ((vdcsum/2) * (idcdc + idcac/2)) - (vacdif*conj(iacac)/2) + (vacsum*conj(iacdc)) );


%% DISPLAY OUTPUTS

fprintf(['ITERATIONS = ' num2str(iterated) '\n\n'])

disp(['VDC SUM = ' num2str(vdcsum, '%3.3e')])
disp(['VAC DIF = ' num2str(real(vacdif), '%3.3e') disp_sign(vacdif) num2str(abs(imag(vacdif)), '%3.3e') 'i'])
disp(['IAC AC = ' num2str(real(iacac), '%3.3e') disp_sign(iacac) num2str(abs(imag(iacac)), '%3.3e') 'i'])
disp(['IDC DC = ' num2str(idcdc, '%3.3e')])

disp(['VAC SUM = ' num2str(real(vacsum), '%3.3e') disp_sign(vacsum) num2str(abs(imag(vacsum)), '%3.3e') 'i'])
disp(['VDC DIF = ' num2str(vdcdif, '%3.3e')])
disp(['IAC DC = ' num2str(real(iacdc), '%3.3e') disp_sign(iacdc) num2str(abs(imag(iacdc)), '%3.3e') 'i'])
disp(['IDC AC = ' num2str(idcac, '%3.3e')])

fprintf('\nCALCULATED VALUES: \n')
disp(['QCON U = ' num2str(Qconu, '%3.3e') 'i'])
disp(['QCON L = ' num2str(Qconl, '%3.3e') 'i'])

plot_AC(vacdif, iacac, 'AC Grid Values', 'temp')