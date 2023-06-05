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
Xarm = 5;
Xl = 8;
R = 2;
Rl = 4;

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

%Converter Limits
voltage_lim = 900*1e3;
current_lim = 2500;

%Sweep Settings
angle_size = 0.5;
magnitude_steps = 100;
min_exponent = 0.9;
max_exponent = 2.5;
change_percentage = 0.05;
varying = 1; %Vgrid = 0; Vhvdc = 1;


%% PRE-SWEEP CALCULATIONS

%Combines variables into groups
Pcon_array = [Pconu_a, Pconl_a, Pconu_b, Pconl_b, Pconu_c, Pconl_c];
Ref_array = [idcdif_ref_a, imiacsum_ref_a, idcdif_ref_b, imiacsum_ref_b, idcdif_ref_c, imiacsum_ref_c];


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
    in = ones(36,1) * 1000;
    failed_voltage_angle = 0;
    failed_voltage_magnitude = 0;
    failed_current_angle = 0;
    failed_current_magnitude = 0;
    failed_max = [];
    data_collection = zeros(36, (360/angle_size)-1);
    
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
            final = ThreePhase_SixArm_Calc(in, max_iteration, tolerance, R, Rl, Xl, Xarm, Vhvdc_temp, Vgrid_RE_temp, Vgrid_IM_temp, Pcon_array, Pgrid, Qgrid, Ref_array, b_phase, c_phase);

            %Variable extraction
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
            vacsum_a = revacsum_a + (imvacsum_a * 1i);
            vacdif_a = revacdif_a + (imvacdif_a * 1i);
            iacdif_a = reiacdif_a + (imiacdif_a * 1i);
            iacsum_a = reiacsum_a + (imiacsum_a * 1i);
            vacsum_b = revacsum_b + (imvacsum_b * 1i);
            vacdif_b = revacdif_b + (imvacdif_b * 1i);
            iacdif_b = reiacdif_b + (imiacdif_b * 1i);
            iacsum_b = reiacsum_b + (imiacsum_b * 1i);
            vacsum_c = revacsum_c + (imvacsum_c * 1i);
            vacdif_c = revacdif_c + (imvacdif_c * 1i);
            iacdif_c = reiacdif_c + (imiacdif_c * 1i);
            iacsum_c = reiacsum_c + (imiacsum_c * 1i);
    
            %Transforms variables
            Vacu_a = -vacdif_a + vacsum_a/2;
            Vdcu_a = -vdcdif_a + vdcsum_a/2;
            Iacu_a = iacdif_a/2 + iacsum_a;
            Idcu_a = idcdif_a/2 + idcsum_a;
            Vacl_a = vacdif_a + vacsum_a/2;
            Vdcl_a = vdcdif_a + vdcsum_a/2;
            Iacl_a = -iacdif_a/2 + iacsum_a;
            Idcl_a = -idcdif_a/2 + idcsum_a;
            Vacu_b = -vacdif_b + vacsum_b/2;
            Vdcu_b = -vdcdif_b + vdcsum_b/2;
            Iacu_b = iacdif_b/2 + iacsum_b;
            Idcu_b = idcdif_b/2 + idcsum_b;
            Vacl_b = vacdif_b + vacsum_b/2;
            Vdcl_b = vdcdif_b + vdcsum_b/2;
            Iacl_b = -iacdif_b/2 + iacsum_b;
            Idcl_b = -idcdif_b/2 + idcsum_b;
            Vacu_c = -vacdif_c + vacsum_c/2;
            Vdcu_c = -vdcdif_c + vdcsum_c/2;
            Iacu_c = iacdif_c/2 + iacsum_c;
            Idcu_c = idcdif_c/2 + idcsum_c;
            Vacl_c = vacdif_c + vacsum_c/2;
            Vdcl_c = vdcdif_c + vdcsum_c/2;
            Iacl_c = -iacdif_c/2 + iacsum_c;
            Idcl_c = -idcdif_c/2 + idcsum_c;


            %% FINISHED UP TO THIS POINT

            %Check for limits
            if check_limit(Vacu_a, Vdcu_a, voltage_lim) == 0 || check_limit(Vacl_a, Vdcl_a, voltage_lim) == 0 || check_limit(Vacu_b, Vdcu_b, voltage_lim) == 0 || check_limit(Vacl_b, Vdcl_b, voltage_lim) == 0 || check_limit(Vacu_c, Vdcu_c, voltage_lim) == 0 || check_limit(Vacl_c, Vdcl_c, voltage_lim) == 0 %FAILED CHECK
                failed_voltage_angle = [failed_voltage_angle, angle];
                failed_voltage_magnitude = [failed_voltage_magnitude, magnitude*change];
                disp([num2str(angle_loop) ', ' num2str(angle) ', ' num2str(change) ': VOLTAGE LIMIT'])
                data_collection(:,angle_loop+1) = final;
                break
            elseif check_limit(Iacu_a, Idcu_a, current_lim) == 0 || check_limit(Iacl_a, Idcl_a, current_lim) == 0 || check_limit(Iacu_b, Idcu_b, current_lim) == 0 || check_limit(Iacl_b, Idcl_b, current_lim) == 0 || check_limit(Iacu_c, Idcu_c, current_lim) == 0 || check_limit(Iacl_c, Idcl_c, current_lim) == 0 %FAILED CHECK
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

% %Adds textbox with nominal operating conditions
% msg_Pconu = ['Pconu = ' num2str(Pconu, '%.2e')];
% msg_Pconl = ['Pconu = ' num2str(Pconl, '%.2e')];
% msg_Vgrid = ['Vgrid = ' num2str(Vgrid, '%.2e')];
% msg_Vhvdc = ['Vhvdc = ' num2str(Vhvdc, '%.2e')];
% msg_XarmPU = ['Xarm = ' num2str(Xarm)];
% msg_RPU = ['R = ' num2str(R)];
% msg_RlPU = ['Rl = ' num2str(Rl)];
% msg_XlPU = ['Xl = ' num2str(Xl)];
% msg_Idc = ['Idcdif ref = ' num2str(idcdif_ref)];
% msg_Iac = ['Im(Iacsum) ref = ' num2str(imiacsum_ref)];
% msg_Vlim = ['Voltage Limit = ' num2str(voltage_lim, '%.2e')];
% msg_Ilim = ['Current Limit = ' num2str(current_lim, '%.2e')];
% msg = {msg_Pconu msg_Pconl msg_Vgrid msg_Vhvdc msg_XarmPU msg_XlPU msg_RlPU msg_RPU msg_Idc msg_Iac msg_Vlim msg_Ilim};
% annotation('textbox', [.131 .131 .795 .795],'String',msg,'FitBoxToText','on');




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
debug_current = abs(Iacu)*sqrt(2) + abs(Idcu);