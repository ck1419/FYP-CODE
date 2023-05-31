%% STARTUP

clc
clear
close all


%% INITIAL VARIABLES

%Settings for Newton-Rhapson
iterations = 25;
tolerance = 1;   %Used to check results
variable_count = 6;

%Operating Points
Pcon = 0;
Vgrid_RE = 400 * 1e3;
Vgrid_IM = 0 * 1e3;
Vhvdc = 600 * 1e3;

%Converter Limits
voltage_lim = 1175e3;
current_lim = 2500;


%% SWEEP SETTINGS

angle_size = 0.5;

exponent_mat = linspace(0.1,1.5,1000);
magnitude_coefficient = (10 .^ exponent_mat - 0.9)/10;

%For runs where Rarm and Xarm needs to remain the same
Xarm = 8;
Rarm = 4;
R = 4;


%% NEWTON-RHAPSON SWEEP

figure
hold on
grid on
axis equal

for nominal_change = 1:3

    if nominal_change == 2
        Vgrid_RE = Vgrid_RE * 1.05;
    elseif nominal_change == 3
        Vgrid_RE = (Vgrid_RE / 1.05) * 0.95;
    end

    %Initial matrix to Solve newton-Raphson with
    x = zeros(variable_count, iterations);
    x(:,1) = ones(variable_count,1) * 100;
    
    failed_voltage_angle = 0;
    failed_voltage_magnitude = 0;
    failed_current_angle = 0;
    failed_current_magnitude = 0;
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
    
            %Loop to execute Newton-Raphson
            for n = 2:iterations
                f11_value = f11(x(:,n-1), Pcon, Xarm, R, Rarm, Vgrid_RE, Vgrid_IM, Vhvdc, Pgrid, Qgrid);
                f11_delta_value = f11_delta(x(:,n-1), Xarm, R, Rarm, Vgrid_RE, Vgrid_IM);
                x(:,n) = x(:,n-1) - (f11_delta_value^-1 * f11_value);
                if (abs(f11_value)) <= tolerance
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
    
            %Finds combined Vac/Iac values
            Vac = final(1) + (final(2)*1i);
            Iac = final(3) + (final(4)*1i);
            Vdc = final(5);
            Idc = final(6);
    
            %Check for limits
            if check_voltage_limit(Vac, Vdc, voltage_lim) == 0 %FAILED CHECK
                failed_voltage_angle = [failed_voltage_angle, angle];
                failed_voltage_magnitude = [failed_voltage_magnitude, magnitude*change];
                disp([num2str(angle_loop) ', ' num2str(angle) ', ' num2str(change) ': VOLTAGE LIMIT'])
                data_collection(:,angle_loop+1) = final;
                break
            elseif check_current_limit(Iac, Idc, current_lim) == 0 %FAILED CHECK
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
    
    if nominal_change == 1
        failed_current_p = failed_current_magnitude .* cosd(failed_current_angle);
        failed_current_q = failed_current_magnitude .* sind(failed_current_angle);
        failed_voltage_p = failed_voltage_magnitude .* cosd(failed_voltage_angle);
        failed_voltage_q = failed_voltage_magnitude .* sind(failed_voltage_angle);
        failed_max_p = 500*1e6*2 * cosd(failed_max);
        failed_max_q = 500*1e6*2 * sind(failed_max);
    elseif nominal_change == 2
        failed_current_p_increase = failed_current_magnitude .* cosd(failed_current_angle);
        failed_current_q_increase = failed_current_magnitude .* sind(failed_current_angle);
        failed_voltage_p_increase = failed_voltage_magnitude .* cosd(failed_voltage_angle);
        failed_voltage_q_increase = failed_voltage_magnitude .* sind(failed_voltage_angle);
        failed_max_p_increase = 500*1e6*2 * cosd(failed_max);
        failed_max_q_increase = 500*1e6*2 * sind(failed_max);
    elseif nominal_change == 3
        failed_current_p_decrease = failed_current_magnitude .* cosd(failed_current_angle);
        failed_current_q_decrease = failed_current_magnitude .* sind(failed_current_angle);
        failed_voltage_p_decrease = failed_voltage_magnitude .* cosd(failed_voltage_angle);
        failed_voltage_q_decrease = failed_voltage_magnitude .* sind(failed_voltage_angle);
        failed_max_p_decrease = 500*1e6*2 * cosd(failed_max);
        failed_max_q_decrease = 500*1e6*2 * sind(failed_max);
    end
end

    
%% PLOT SWEEP RESULTS

Pcon = 0;
Vgrid_RE = 400 * 1e3;
Vgrid_IM = 0 * 1e3;
Vhvdc = 600 * 1e3;
Vgrid = Vgrid_RE + (Vgrid_IM * 1i);

plot(failed_current_p_decrease, failed_current_q_decrease, '.', 'color', "#008000", 'markersize', 5)
plot(failed_current_p, failed_current_q, '.', 'color', "#FF0000", 'markersize', 5)
plot(failed_current_p_increase, failed_current_q_increase, '.', 'color', "#0000FF", 'markersize', 5)

plot(failed_voltage_p_decrease, failed_voltage_q_decrease, '.', 'color', "#7CFC00", 'markersize', 5)
plot(failed_voltage_p, failed_voltage_q, '.', 'color', "#FF00FF", 'markersize', 5)
plot(failed_voltage_p_increase, failed_voltage_q_increase, '.', 'color', "#0096FF", 'markersize', 5)

xlabel('Pgrid')
ylabel('Qgrid')
title('Single-Phase Single-Arm (V_{GRID} Changed)')
legend('5% Decrease', 'Nominal Value', '5% Increase')

msg_Vgrid = ['Vgrid = ' num2str(Vgrid, '%.2e')];
msg_Vhvdc = ['Vhvdc = ' num2str(Vhvdc, '%.2e')];
msg_RPU = ['R = ' num2str(R)];
msg_RarmPU = ['Rarm = ' num2str(Rarm)];
msg_XarmPU = ['Xarm = ' num2str(Xarm)];
msg_Vlim = ['Voltage Limit = ' num2str(voltage_lim, '%.2e')];
msg_Ilim = ['Current Limit = ' num2str(current_lim, '%.2e')];

msg = {msg_Vgrid msg_Vhvdc msg_RPU msg_RarmPU msg_XarmPU msg_Vlim msg_Ilim};
annotation('textbox', [.131 .131 .795 .795],'String',msg,'FitBoxToText','on');


%% DATA FOR DEBUG

%Finds combined Vac/Iac values
Vac = data_collection(1,:) + (data_collection(2,:)*1i);
Iac = data_collection(3,:) + (data_collection(4,:)*1i);
Vdc = data_collection(5,:);
Idc = data_collection(6,:);

debug_voltage = abs(Vac)*sqrt(2) + abs(Vdc);
debug_current = abs(Iac)*sqrt(2) + abs(Idc);