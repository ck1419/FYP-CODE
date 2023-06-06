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
voltage_lim = 900*1e3;
current_lim = 2500;

%Sweep Settings
angle_size = 0.5;
magnitude_steps = 1000;
min_exponent = 0.9;
max_exponent = 2.5;
change_percentage = 0.05;
varying = 0; %Vgrid = 0; Vhvdc = 1;
halfbridge = 0; %0 = fullbridge; 1 = halfbridge


%% NEWTON-RHAPSON SWEEP

%Creates the multipliers for the sweep
exponent_mat = linspace(min_exponent,max_exponent,magnitude_steps);
magnitude_coefficient = (10 .^ exponent_mat - 0.9)/10;
    
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
    in = ones(12,1) * 1000;
    failed_voltage_angle = 0;
    failed_voltage_magnitude = 0;
    failed_current_angle = 0;
    failed_current_magnitude = 0;
    failed_max = [];
    data_collection = zeros(12, (360/angle_size)-1);
    
    %Loop to change Sgrid angle
    disp('ITERATION / ANGLE / MAGNITUDE MULTIPLIER')
    for angle_loop = 0:(360/angle_size)-1
        angle = angle_loop * angle_size;
        magnitude = 500 * 1e6;
        
        %Loop to change Sgrid magnitude
        for change = magnitude_coefficient
            Pgrid = magnitude * cosd(angle) * change;
            Qgrid = magnitude * sind(angle) * change;
            
            %Newton-Raphson calculation
            final = SinglePhase_TwoArm_Calc(in, max_iteration, tolerance, R, Rl, Xl, Xarm, Vhvdc_temp, Vgrid_RE_temp, Vgrid_IM_temp, Pconu, Pconl, Pgrid, Qgrid, idcdif_ref, imiacsum_ref);
    
            %Extracts variables from Newton-Raphson output
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
    
            %Transforms variables
            Vacu = -vacdif + vacsum/2;
            Vdcu = -vdcdif + vdcsum/2;
            Iacu = iacdif/2 + iacsum;
            Idcu = idcdif/2 + idcsum;
            Vacl = vacdif + vacsum/2;
            Vdcl = vdcdif + vdcsum/2;
            Iacl = -iacdif/2 + iacsum;
            Idcl = -idcdif/2 + idcsum;

            %Check for limits
            if check_limit(Vacu, Vdcu, voltage_lim, halfbridge) || check_limit(Vacl, Vdcl, voltage_lim, halfbridge) %FAILED CHECK
                failed_voltage_angle = [failed_voltage_angle, angle];
                failed_voltage_magnitude = [failed_voltage_magnitude, magnitude*change];
                disp([num2str(angle_loop) ', ' num2str(angle) ', ' num2str(change) ': VOLTAGE LIMIT'])
                data_collection(:,angle_loop+1) = final;
                break
            elseif check_limit(Iacu, Idcu, current_lim, halfbridge) || check_limit(Iacl, Idcl, current_lim, halfbridge) %FAILED CHECK
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
plot(failed_current_p_decrease, failed_current_q_decrease, '.', 'color', "#008000", 'markersize', 5)
plot(failed_current_p, failed_current_q, '.', 'color', "#FF0000", 'markersize', 5)
plot(failed_current_p_increase, failed_current_q_increase, '.', 'color', "#0000FF", 'markersize', 5)
plot(failed_voltage_p_decrease, failed_voltage_q_decrease, '.', 'color', "#7CFC00", 'markersize', 5)
plot(failed_voltage_p, failed_voltage_q, '.', 'color', "#FF00FF", 'markersize', 5)
plot(failed_voltage_p_increase, failed_voltage_q_increase, '.', 'color', "#0096FF", 'markersize', 5)

%Apply axis labels and legend
xlabel('Pgrid')
ylabel('Qgrid')
legend([num2str(change_percentage*100) '% Decrease'], 'Nominal Value', [num2str(change_percentage*100) '% Increase'])
if varying == 0
    title('Single-Phase Two-Arm (V_{GRID} Changed)')
else
    title('Single-Phase Two-Arm (V_{HVDC} Changed)')
end

%Adds textbox with nominal operating conditions
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
msg_Vlim = ['Voltage Limit = ' num2str(voltage_lim, '%.2e')];
msg_Ilim = ['Current Limit = ' num2str(current_lim, '%.2e')];
msg = {msg_Pconu msg_Pconl msg_Vgrid msg_Vhvdc msg_XarmPU msg_XlPU msg_RlPU msg_RPU msg_Idc msg_Iac msg_Vlim msg_Ilim};
annotation('textbox', [.131 .131 .795 .795],'String',msg,'FitBoxToText','on');




%% DATA FOR DEBUG

vdcsum = data_collection(1,:);
vdcdif = data_collection(2,:);   
idcdif = data_collection(3,:);  
idcsum = data_collection(4,:);  
revacsum = data_collection(5,:);       
imvacsum = data_collection(6,:);
revacdif = data_collection(7,:);  
imvacdif = data_collection(8,:);       
imiacsum = data_collection(9,:);  
reiacsum = data_collection(10,:); 
reiacdif = data_collection(11,:); 
imiacdif = data_collection(12,:); 
vacsum = revacsum + (imvacsum * 1i);
vacdif = revacdif + (imvacdif * 1i);
iacdif = reiacdif + (imiacdif * 1i);
iacsum = reiacsum + (imiacsum * 1i);

Vacu = -vacdif + vacsum/2;
Vdcu = -vdcdif + vdcsum/2;
Iacu = iacdif/2 + iacsum;
Idcu = idcdif/2 + idcsum;
Vacl = vacdif + vacsum/2;
Vdcl = vdcdif + vdcsum/2;
Iacl = -iacdif/2 + iacsum;
Idcl = -idcdif/2 + idcsum;

debug_voltage = abs(Vacu)*sqrt(2) + abs(Vdcu);
debug_current = abs(iacdif/2)*sqrt(2) + abs(idcsum);