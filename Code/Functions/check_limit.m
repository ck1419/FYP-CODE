%Checks for the converter's limitations
function pass = check_limit(ac, dc, rated, halfbridge)
    if (abs(ac)*sqrt(2) + abs(dc) >= rated) || ((abs(dc) - abs(ac)*sqrt(2) < 0) && halfbridge)
        pass = true;
    else
        pass = false;
    end
end