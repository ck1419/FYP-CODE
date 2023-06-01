%Checks for the converter's limitations
%0 = Fail
function pass = check_limit(ac, dc, rated_voltage)
    if abs(ac)*sqrt(2) + abs(dc) >= rated_voltage
        pass = 0;
    else
        pass = 1;
    end
end