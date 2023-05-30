%Checks for the converter's current limitations
%0 = Fail
function pass = check_current_limit(Iac, Idc, rated_current)
    if abs(Iac)*sqrt(2) + abs(Idc) >= rated_current
        pass = 0;
%     elseif (abs(Idc) - abs(Iac)) < 0
%         pass = 0;
    else
        pass = 1;
    end
end