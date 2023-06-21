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
Vgrid_RE = 220 * 1e3;
Vgrid_IM = 0 * 1e3;
Vhvdc = 800 * 1e3;
idcdif_ref = 0;
imiacsum_ref = 0;
Xarm = 38.11;
Xl = 19.05;
R = 1.15;
Rl = 0.55;

%Converter Limits
voltage_lim = 1200e3;
current_lim = 1800;

%Sweep Settings
angle_size = 1;
magnitude_steps = 600;
min_magnitude = 0;
max_magnitude = 1200*1e6;
change_percentage = 0.1;
varying = 0; %Vgrid = 0; Vhvdc = 1;
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
    in = ones(12,1) * 1000;
    failed_voltage_angle = [];
    failed_voltage_magnitude = [];
    failed_current_angle = [];
    failed_current_magnitude = [];
    failed_max = [];
    data_collection = zeros(12, (360/angle_size)-1);
    
    %Loop to change Sgrid angle
    disp('CALCULATION / ITERATION / ANGLE / MAGNITUDE')
    for angle_loop = 0:(360/angle_size)-1
        angle = angle_loop * angle_size;
        magnitude = 500 * 1e6;
        
        %Loop to change Sgrid magnitude
        for magnitude = magnitude_mat
            Pgrid = magnitude * cosd(angle);
            Qgrid = magnitude * sind(angle);
            
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
                failed_voltage_magnitude = [failed_voltage_magnitude, magnitude];
                disp([num2str(nominal_change) ', ' num2str(angle_loop) ', ' num2str(angle) ', ' num2str(magnitude/1e6) ': VOLTAGE LIMIT'])
                data_collection(:,angle_loop+1) = final;
                in = final;
                break
            elseif check_limit(Iacu, Idcu, current_lim, 0) || check_limit(Iacl, Idcl, current_lim, 0) %FAILED CHECK
                failed_current_angle = [failed_current_angle, angle];
                failed_current_magnitude = [failed_current_magnitude, magnitude];
                disp([num2str(nominal_change) ', ' num2str(angle_loop) ', ' num2str(angle) ', ' num2str(magnitude/1e6) ': CURRENT LIMIT'])
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
if halfbridge == 1
    halfbridge_text = ' Half-Bridge ';
else
    halfbridge_text = ' Full-Bridge ';
end
if varying == 1
    varying_text = '(V_{HVDC} Changed)';
else
    varying_text = '(V_{GRID} Changed)';
end
title(['Single-Phase Two-Arm' halfbridge_text varying_text])

%Adds textbox with nominal operating conditions
msg_Pconu = ['Pconu = ' num2str(Pconu/1e6) ' MW'];
msg_Pconl = ['Pconl = ' num2str(Pconl/1e6') ' MW'];
msg_Vgrid = ['Vgrid = ' num2str(Vgrid/1e3) ' kV'];
msg_Vhvdc = ['Vhvdc = ' num2str(Vhvdc/1e3) ' kV'];
msg_Xarm = ['Xarm = ' num2str(Xarm) ' Ω'];
msg_R = ['R = ' num2str(R) ' Ω'];
msg_Rl = ['Rl = ' num2str(Rl) ' Ω'];
msg_Xl = ['Xl = ' num2str(Xl) ' Ω'];
msg_Idc = ['Idcdif ref = ' num2str(idcdif_ref) ' A'];
msg_Iac = ['Im(Iacsum) ref = ' num2str(imiacsum_ref) ' A'];
msg_Vlim = ['Voltage Limit = ' num2str(voltage_lim/1e3) ' kV'];
msg_Ilim = ['Current Limit = ' num2str(current_lim) ' A'];
msg = {msg_Pconu msg_Pconl msg_Vgrid msg_Vhvdc msg_Xarm msg_Xl msg_Rl msg_R msg_Idc msg_Iac msg_Vlim msg_Ilim};
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