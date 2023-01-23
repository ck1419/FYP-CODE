function pass = check_voltage_limit(Vac, Vdc, rated_voltage)
    if abs(Vac)*sqrt(2) + Vdc <= rated_voltage
        pass = 1;
        fprintf('PASSED: converter rated voltage limit\n')
    else
        pass = 0;
        fprintf('FAILED: converter rated voltage limit\n')
    end
end