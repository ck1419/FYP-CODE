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
Vgrid_RE = 400 * 1e3;
Vgrid_IM = 150 * 1e3;
Vhvdc = 200 * 1e3;
Xarm_PU = 0.2;
Xl_PU = 0.4;
R_PU = 0.05;
Rl_PU = 0.1;

voltage_lim = 1000*1e3;
current_lim = 1.5*1e3;

idcac_ref = 1*1e-3;
reiacdc_ref = 1*1e-3;


%% SWEEP SETTINGS

angle_size = 10;

exponent_mat = linspace(0.3,1.35,30);
magnitude_coefficient = (10 .^ exponent_mat - 0.9)/10;


%% NEWTON-RHAPSON SWEEP

%Initial matrix to Solve newton-Raphson with
x = zeros(variable_count, iterations);
x(:,1) = ones(variable_count,1) * 100;

failed_voltage_angle = [];
failed_voltage_magnitude = [];
failed_current_angle = [];
failed_current_magnitude = [];
failed_max = [];
data_collection = zeros(variable_count, (360/angle_size)-1);

disp('ITERATION / ANGLE / MAGNITUDE MULTIPLIER')
%Loop to change operating conditions
for angle_loop = 0:(360/angle_size)-1
    angle = angle_loop * angle_size;
    magnitude = 500 * 1e6;

    for change = magnitude_coefficient
        Pgrid = magnitude * cosd(angle) * change;
        Qgrid = magnitude * sind(angle) * change;

        %Combines imaginary and real components
        Vgrid = Vgrid_RE + (Vgrid_IM * 1i);
        Sgrid = Pgrid + (Qgrid * 1i);
        
        %PU Conversion
        Z_PU = abs(Vgrid)^2 / abs(Sgrid);
        Xl = Xl_PU * Z_PU;
        Xarm = Xarm_PU * Z_PU;
        R = R_PU * Z_PU;
        Rl = Rl_PU * Z_PU;

        %Loop to execute Newton-Raphson
        for n = 2:iterations
            f12_value = f12(x(:,n-1), R, Rl, Xl, Xarm, Vhvdc, Vgrid, Pconu, Pconl, Sgrid, idcac_ref, reiacdc_ref);
            f12_delta_value = f12_delta(x(:,n-1), R, Rl, Xl, Xarm, Vhvdc, Vgrid, Pconu, Pconl, Sgrid, idcac_ref, reiacdc_ref);
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

        %Cleans and assign variables
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

        %Check for limits
        if check_voltage_limit(vacdif, vdcsum, voltage_lim) == 0 %FAILED CHECK
            failed_voltage_angle = [failed_voltage_angle, angle];
            failed_voltage_magnitude = [failed_voltage_magnitude, magnitude*change];
            disp([num2str(angle_loop) ', ' num2str(angle) ', ' num2str(change) ': VOLTAGE LIMIT'])
            data_collection(:,angle_loop+1) = final;
            break
        elseif check_current_limit(iacac/2, idcdc, current_lim) == 0 %FAILED CHECK
            failed_current_angle = [failed_current_angle, angle];
            failed_current_magnitude = [failed_current_magnitude, magnitude*change];
            disp([num2str(angle_loop) ', ' num2str(angle) ', ' num2str(change) ': CURRENT LIMIT'])
            data_collection(:,angle_loop+1) = final;
            break
        elseif change == magnitude_coefficient(end)
            failed_max = [failed_max, angle];
            disp([num2str(angle_loop) ', ' num2str(angle) ': MAXED'])
            data_collection(:,angle_loop+1) = final;
        end
    end
end


%% PLOT SWEEP RESULTS

failed_current_p = failed_current_magnitude .* cosd(failed_current_angle);
failed_current_q = failed_current_magnitude .* sind(failed_current_angle);

failed_voltage_p = failed_voltage_magnitude .* cosd(failed_voltage_angle);
failed_voltage_q = failed_voltage_magnitude .* sind(failed_voltage_angle);

failed_max_p = 500*1e6*2 * cosd(failed_max);
failed_max_q = 500*1e6*2 * sind(failed_max);

max_magnitude = 500 * 1e6 * 3;

figure
hold on
grid on
axis equal
plot(failed_current_p, failed_current_q)
plot(failed_voltage_p, failed_voltage_q)
plot(failed_max_p, failed_max_q, 'x')
xlabel('Pgrid')
ylabel('Qgrid')
% xlim([-max_magnitude, max_magnitude])
% ylim([-max_magnitude, max_magnitude])
legend('Current Limit', 'Voltage Limit', 'Max Iterations')


%% DATA FOR DEBUG

vdcsum = data_collection(1,:);
vdcdif = data_collection(2,:);         % 0
revacsum = data_collection(3,:);       % 0
imvacsum = data_collection(4,:);       % 0
revacdif = data_collection(5,:);       
imvacdif = data_collection(6,:);
reiacac = data_collection(7,:);        %Iac of AC grid
imiacac = data_collection(8,:);       
reiacdc = data_collection(9,:);        %Iac of DC grid - 0
imiacdc = data_collection(10,:);       % 0
idcac = data_collection(11,:);         %Idc of AC grid - 0
idcdc = data_collection(12,:);         %Idc of DC grid
vacsum = revacsum + (imvacsum * 1i);
vacdif = revacdif + (imvacdif * 1i);
iacac = reiacac + (imiacac * 1i);
iacdc = reiacdc + (imiacdc * 1i);

debug_voltage = abs(vacdif)*sqrt(2) + abs(vdcsum);
debug_current = abs(iacac/2)*sqrt(2) + abs(idcdc);