%Checks for the converter's current limitations
function pass = check_current_limit(Iac, Idc, rated_voltage)
    if abs(Iac)*sqrt(2) + abs(Idc) <= rated_voltage
        pass = 1;
    else
        pass = 0;
    end
end