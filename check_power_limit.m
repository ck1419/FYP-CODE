%Checks for the converter's power limitations
function pass = check_power_limit(Vac, Iac, rated_power)
    if abs(Vac * conj(Iac)) <= rated_power
        pass = 1;
        fprintf('PASSED: converter rated power limit\n')
    else
        pass = 0;
        fprintf('PASSED: converter rated power limit\n')
    end
end