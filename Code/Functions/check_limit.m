%Checks for the converter's limitations
%0 = Fail
function pass = check_limit(ac, dc, rated, halfbridge)
    if (abs(ac)*sqrt(2) + abs(dc) >= rated) || ((abs(dc) - abs(ac)*sqrt(2) < 0) && halfbridge)
        pass = 1;
    else
        pass = 0;
    end
end