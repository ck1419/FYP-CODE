%% STARTUP

clc
clear
close all
addpath("../Functions/")


%% CHANGEABLE VARIABLES

%Settings for Newton-Rhapson
max_iteration = 25;
tolerance = 0.05;   %Used to check results

%Operating Points
Pconu = 0;
Pconl = 0;
Vgrid_RE = 400 * 1e3;
Vgrid_IM = 0 * 1e3;
Vhvdc = 600 * 1e3;
idcdif_ref = 0;
imiacsum_ref = 0;
Xarm = 5;
Xl = 8;
R = 2;
Rl = 4;

%Converter Limits
voltage_lim = 1180*1e3;
current_lim = 1500;

%Sweep Settings
angle_size = 1;
magnitude_steps = 500;


%% NEWTON-RHAPSON SWEEP

exponent_mat = linspace(0.1,1.5,magnitude_steps);
magnitude_coefficient = (10 .^ exponent_mat - 0.9)/10;
    
for nominal_change = 1:3

    if nominal_change == 2
        Vgrid_RE = Vgrid_RE * (1+change_percentage);
    elseif nominal_change == 3
        Vgrid_RE = (Vgrid_RE / (1+change_percentage)) * (1-change_percentage);
    end

    %Initial matrix to Solve newton-Raphson with
    in = ones(12,1) * 1000;
    failed_voltage_angle = 0;
    failed_voltage_magnitude = 0;
    failed_current_angle = 0;
    failed_current_magnitude = 0;
    failed_max = [];
    data_collection = zeros(12, (360/angle_size)-1);
    
    disp('ITERATION / ANGLE / MAGNITUDE MULTIPLIER')
    %Loop to change operating conditions
    for angle_loop = 0:(360/angle_size)-1
        angle = angle_loop * angle_size;
        magnitude = 500 * 1e6;
    
        for change = magnitude_coefficient
            Pgrid = magnitude * cosd(angle) * change;
            Qgrid = magnitude * sind(angle) * change;
            
            final = SinglePhase_TwoArm_Calc(in, max_iteration, tolerance, R, Rl, Xl, Xarm, Vhvdc, Vgrid_RE, Vgrid_IM, Pconu, Pconl, Pgrid, Qgrid, idcdif_ref, imiacsum_ref);
    
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
            vacsum = revacsum + (imvacsum * 1i);
            vacdif = revacdif + (imvacdif * 1i);
            iacdif = reiacdif + (imiacdif * 1i);
            iacsum = reiacsum + (imiacsum * 1i);
    
            %Check for limits
            if check_limit(vacdif, vdcsum, voltage_lim) == 0 %FAILED CHECK
                failed_voltage_angle = [failed_voltage_angle, angle];
                failed_voltage_magnitude = [failed_voltage_magnitude, magnitude*change];
                disp([num2str(angle_loop) ', ' num2str(angle) ', ' num2str(change) ': VOLTAGE LIMIT'])
                data_collection(:,angle_loop+1) = final;
                break
            elseif check_limit(iacdif/2, idcsum, current_lim) == 0 %FAILED CHECK
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

Pconu = 0;
Pconl = 0;
Vgrid_RE = 400 * 1e3;
Vgrid_IM = 0 * 1e3;
Vhvdc = 600 * 1e3;

figure
hold on
grid on
axis equal

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

msg_Pconu = ['Pconu = ' num2str(Pconu, '%.2e')];
msg_Pconl = ['Pconu = ' num2str(Pconl, '%.2e')];
msg_Vgrid = ['Vgrid = ' num2str(Vgrid, '%.2e')];
msg_Vhvdc = ['Vhvdc = ' num2str(Vhvdc, '%.2e')];
msg_XarmPU = ['Xarm = ' num2str(Xarm)];
msg_RPU = ['R = ' num2str(R)];
msg_RlPU = ['Rl = ' num2str(Rl)];
msg_XlPU = ['Xl = ' num2str(Xl)];
msg_Idc = ['Idcdif ref = ' num2str(idcdif_ref)];
msg_Iac = ['Im(Iacsum) ref = ' num2str(imiacsum_ref)];

msg = {msg_Pconu msg_Pconl msg_Vgrid msg_Vhvdc msg_XarmPU msg_XlPU msg_RlPU msg_RPU msg_Idc msg_Iac};
annotation('textbox', [.57 .25 .565 .286],'String',msg,'FitBoxToText','on');


%% DATA FOR DEBUG

vdcsum = data_collection(1,:);
vdcdif = data_collection(2,:);         % 0
idcdif = data_collection(3,:);       % 0
idcsum = data_collection(4,:);       % 0
revacsum = data_collection(5,:);       
imvacsum = data_collection(6,:);
revacdif = data_collection(7,:);        %Iac of AC grid
imvacdif = data_collection(8,:);       
imiacsum = data_collection(9,:);        %Iac of DC grid - 0
reiacsum = data_collection(10,:);       % 0
reiacdif = data_collection(11,:);         %Idc of AC grid - 0
imiacdif = data_collection(12,:);         %Idc of DC grid
vacsum = revacsum + (imvacsum * 1i);
vacdif = revacdif + (imvacdif * 1i);
iacdif = reiacdif + (imiacdif * 1i);
iacsum = reiacsum + (imiacsum * 1i);

debug_voltage = abs(vacdif)*sqrt(2) + abs(vdcsum);
debug_current = abs(iacdif/2)*sqrt(2) + abs(idcsum);