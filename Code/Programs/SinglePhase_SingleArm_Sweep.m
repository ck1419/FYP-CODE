%% STARTUP

clc
clear
close all
addpath("../Functions/")


%% CHANGEABLE VARIABLES

%Settings for Newton-Rhapson
iterations = 25;
tolerance = 0.05;   %Used to check results

%Operating Points
Pcon = 0;
Vgrid_RE = 400 * 1e3;
Vgrid_IM = 0 * 1e3;
Vhvdc = 600 * 1e3;
Xarm = 15.25;
Rarm = 0.4;
R = 2;

%Converter Limits
voltage_lim = 1250e3;
current_lim = 6050;

%Sweep Settings
angle_size = 0.5;
magnitude_steps = 1250;
min_magnitude = 0;
max_magnitude = 2000*1e6;
change_percentage = 0.045;
varying = 1; %Vgrid = 0; Vhvdc = 1;
halfbridge = 1; %fullbridge = 0; halfbridge = 1;


%% NEWTON-RHAPSON SWEEP

%Creates the multipliers for the sweep
magnitude_mat = linspace(min_magnitude, max_magnitude, magnitude_steps);

%First loop for change in operating condition
for nominal_change = 1:3
    
    %Changes operating condition each loop
    Vgrid_RE_temp = Vgrid_RE;
    Vgrid_IM_temp = Vgrid_IM;
    Vhvdc_temp = Vhvdc;
    if varying == 1
        if nominal_change == 2
            Vhvdc_temp = Vhvdc * (1+change_percentage);
        elseif nominal_change == 3
            Vhvdc_temp = Vhvdc * (1-change_percentage);
        end
    else
        if nominal_change == 2
            Vgrid_RE_temp = Vgrid_RE * (1+change_percentage);
            Vgrid_IM_temp = Vgrid_IM * (1+change_percentage);
        elseif nominal_change == 3
            Vgrid_RE_temp = Vgrid_RE * (1-change_percentage);
            Vgrid_IM_temp = Vgrid_IM * (1-change_percentage);
        end
    end


    %Initial matrix to Solve newton-Raphson with
    in = ones(6,1) * 1000;
    failed_voltage_angle = [];
    failed_voltage_magnitude = [];
    failed_current_angle = [];
    failed_current_magnitude = [];
    failed_max = [];
    data_collection = zeros(6, (360/angle_size)-1);
    
    %Loop to change Sgrid angle
    disp('CALCULATION / ITERATION / ANGLE / MAGNITUDE')
    for angle_loop = 0:(360/angle_size)-1
        angle = angle_loop * angle_size;
    
        %Loop to change Sgrid magnitude
        for magnitude = magnitude_mat
            Pgrid = magnitude * cosd(angle);
            Qgrid = magnitude * sind(angle);
            
            %Newton-Raphson calculation
            final = SinglePhase_SingleArm_Calc(in, iterations, tolerance, Pcon, Xarm, R, Rarm, Vgrid_RE_temp, Vgrid_IM_temp, Vhvdc_temp, Pgrid, Qgrid);
    
            %Extracts variables from Newton-Raphson output
            Vac = final(1) + (final(2)*1i);
            Iac = final(3) + (final(4)*1i);
            Vdc = final(5);
            Idc = final(6);
    
            %Check for limits
            if check_limit(Vac, Vdc, voltage_lim, halfbridge) %FAILED CHECK
                failed_voltage_angle = [failed_voltage_angle, angle];
                failed_voltage_magnitude = [failed_voltage_magnitude, magnitude];
                disp([num2str(nominal_change) ', ' num2str(angle_loop) ', ' num2str(angle) ', ' num2str(magnitude) ': VOLTAGE LIMIT'])
                data_collection(:,angle_loop+1) = final;
                in = final;
                break
            elseif check_limit(Iac, Idc, current_lim, 0) %FAILED CHECK
                failed_current_angle = [failed_current_angle, angle];
                failed_current_magnitude = [failed_current_magnitude, magnitude];
                disp([num2str(nominal_change) ', ' num2str(angle_loop) ', ' num2str(angle) ', ' num2str(magnitude) ': CURRENT LIMIT'])
                data_collection(:,angle_loop+1) = final;
                in = final;
                break
            elseif magnitude == magnitude_mat(end)
                failed_max = [failed_max, angle];
                disp([num2str(nominal_change) ', ' num2str(angle_loop) ', ' num2str(angle) ': MAXED'])
                data_collection(:,angle_loop+1) = final;
                in = final;
            end
        end
    end
    
    %Stores variable for each operating condition
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

%Calculates variables for output plot
Vgrid = Vgrid_RE + (Vgrid_IM * 1i);

%Plots operating region
figure
hold on
grid on
axis equal
plot(failed_current_p_decrease/1e6, failed_current_q_decrease/1e6, '.', 'color', "#008000", 'markersize', 5)
plot(failed_current_p/1e6, failed_current_q/1e6, '.', 'color', "#FF0000", 'markersize', 5)
plot(failed_current_p_increase/1e6, failed_current_q_increase/1e6, '.', 'color', "#0000FF", 'markersize', 5)
plot(failed_voltage_p_decrease/1e6, failed_voltage_q_decrease/1e6, '.', 'color', "#7CFC00", 'markersize', 5)
plot(failed_voltage_p/1e6, failed_voltage_q/1e6, '.', 'color', "#FF00FF", 'markersize', 5)
plot(failed_voltage_p_increase/1e6, failed_voltage_q_increase/1e6, '.', 'color', "#0096FF", 'markersize', 5)

%Creates legend
h = zeros(3, 1);
h(1) = plot(NaN,NaN, '.', 'color', "#008000", 'markersize', 15);
h(2) = plot(NaN,NaN, '.', 'color', "#FF0000", 'markersize', 15);
h(3) = plot(NaN,NaN, '.', 'color', "#0000FF", 'markersize', 15);
legend(h, [num2str(change_percentage*100) '% Decrease'], 'Nominal Value', [num2str(change_percentage*100) '% Increase'])

%Apply axis labels
xlabel('Pgrid [MW]')
ylabel('Qgrid [MVAR]')
if varying == 0
    title('Single-Phase Single-Arm (V_{GRID} Changed)')
else
    title('Single-Phase Single-Arm (V_{HVDC} Changed)')
end

%Adds textbox with nominal operating conditions
msg_Pcon = ['Pconu = ' num2str(Pcon/1e6) ' MW'];
msg_Vgrid = ['Vgrid = ' num2str(Vgrid/1e3) ' kV'];
msg_Vhvdc = ['Vhvdc = ' num2str(Vhvdc/1e3) ' kV'];
msg_R = ['R = ' num2str(R) ' Ω'];
msg_Rarm = ['Rarm = ' num2str(Rarm) ' Ω'];
msg_Xarm = ['Xarm = ' num2str(Xarm) ' Ω'];
msg_Vlim = ['Voltage Limit = ' num2str(voltage_lim/1e3) ' kV'];
msg_Ilim = ['Current Limit = ' num2str(current_lim) ' A'];
msg_half = ['Halfbridge = ' num2str(halfbridge)];
msg = {msg_Pcon msg_Vgrid msg_Vhvdc msg_R msg_Rarm msg_Xarm msg_Vlim msg_Ilim msg_half};
annotation('textbox', [.131 .131 .795 .795],'String',msg,'FitBoxToText','on');



%% DATA FOR DEBUG

%Finds combined Vac/Iac values
Vac = data_collection(1,:) + (data_collection(2,:)*1i);
Iac = data_collection(3,:) + (data_collection(4,:)*1i);
Vdc = data_collection(5,:);
Idc = data_collection(6,:);

debug_voltage = abs(Vac)*sqrt(2) + abs(Vdc);
debug_current = abs(Iac)*sqrt(2) + abs(Idc);