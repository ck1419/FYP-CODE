%Checks for the converter's current limitations
function pass = check_voltage_limit(Vac, Vdc, rated_voltage)
    if abs(Vac)*sqrt(2) + Vdc <= rated_voltage
        pass = 1;
    else
        pass = 0;
    end
end