%Checks for the converter's current limitations
function pass = check_current_limit(Iac, Idc, rated_voltage)
    if abs(Iac)*sqrt(2) + Idc <= rated_voltage
        pass = 1;
        fprintf('PASSED: converter rated current limit\n')
    else
        pass = 0;
        fprintf('FAILED: converter rated current limit\n')
    end
end