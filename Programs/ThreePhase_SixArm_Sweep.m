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
Vgrid_RE = 220 * 1e3;           %Phase A reference voltage
Vgrid_IM = 0 * 1e3;
Vhvdc = 800 * 1e3;
Xarm = 12.7;
Xl = 6.35;
R = 0.38;
Rl = 0.19;

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
b_phase = -120;  %Relative to A in degrees

%Phase C Operating Points
Pconu_c = 0;
Pconl_c = 0;
idcdif_ref_c = 0;
imiacsum_ref_c = 0;
c_phase = 120; %Relative to A in degrees

%Converter Limits
voltage_lim = 1200e3;
current_lim = 1800;

%Sweep Settings
angle_size = 2;
magnitude_steps = 100;
min_magnitude = 0*1e6;
max_magnitude = 3000*1e6;
change_percentage = 0.25;
varying = 0; %Vgrid = 0; Vhvdc = 1;
halfbridge = 1; %fullbridge = 0; halfbridge = 1;


%% PRE-SWEEP CALCULATIONS

%Combines variables into groups
Pcon_array = [Pconu_a, Pconl_a, Pconu_b, Pconl_b, Pconu_c, Pconl_c];
Ref_array = [idcdif_ref_a, imiacsum_ref_a, idcdif_ref_b, imiacsum_ref_b, idcdif_ref_c, imiacsum_ref_c];


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
    in = ones(36,1) * 1000;
    failed_voltage_angle = [];
    failed_voltage_magnitude = [];
    failed_current_angle = [];
    failed_current_magnitude = [];
    failed_max = [];
    data_collection = zeros(36, (360/angle_size)-1);
    
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
            if check_limit(Vacu_a, Vdcu_a, voltage_lim, halfbridge) || check_limit(Vacl_a, Vdcl_a, voltage_lim, halfbridge) || check_limit(Vacu_b, Vdcu_b, voltage_lim, halfbridge) || check_limit(Vacl_b, Vdcl_b, voltage_lim, halfbridge) || check_limit(Vacu_c, Vdcu_c, voltage_lim, halfbridge) || check_limit(Vacl_c, Vdcl_c, voltage_lim, halfbridge) %FAILED CHECK
                failed_voltage_angle = [failed_voltage_angle, angle];
                failed_voltage_magnitude = [failed_voltage_magnitude, magnitude];
                disp([num2str(nominal_change) ', ' num2str(angle_loop) ', ' num2str(angle, '%.1f') ', ' num2str(magnitude/1e6, '%.2f') ': VOLTAGE LIMIT'])
                data_collection(:,angle_loop+1) = final;
                in = final;
                break
            elseif check_limit(Iacu_a, Idcu_a, current_lim, 0) || check_limit(Iacl_a, Idcl_a, current_lim, 0) || check_limit(Iacu_b, Idcu_b, current_lim, 0) || check_limit(Iacl_b, Idcl_b, current_lim, 0) || check_limit(Iacu_c, Idcu_c, current_lim, 0) || check_limit(Iacl_c, Idcl_c, current_lim, 0) %FAILED CHECK
                failed_current_angle = [failed_current_angle, angle];
                failed_current_magnitude = [failed_current_magnitude, magnitude];
                disp([num2str(nominal_change) ', ' num2str(angle_loop) ', ' num2str(angle, '%.1f') ', ' num2str(magnitude/1e6, '%.2f') ': CURRENT LIMIT'])
                data_collection(:,angle_loop+1) = final;
                in = final;
                break
            elseif magnitude == magnitude_mat(end)
                failed_max = [failed_max, angle];
                disp([num2str(nominal_change) ', ' num2str(angle_loop) ', ' num2str(angle, '%.1f') ': MAXED'])
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
title(['Three-Phase Six-Arm' halfbridge_text varying_text])


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