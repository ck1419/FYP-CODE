%Checks for the converter's current limitations
%0 = Fail
function pass = check_voltage_limit(Vac, Vdc, rated_voltage)
    if abs(Vac)*sqrt(2) + abs(Vdc) >= rated_voltage
        pass = 0;
%     elseif (abs(Vdc) - abs(Vac)) < 0
%         pass = 0;
    else
        pass = 1;
    end
end