%Checks for the direction of the power flow
function pass = check_powerflow(Vdc, Idc)
    if Vdc * Idc <= 0
        pass = 1;
        fprintf('PASSED: power flow direction\n')
    else
        pass = 0;
        fprintf('FAILED: power flow direction\n')
    end
end

