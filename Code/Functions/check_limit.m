%Checks for the converter's limitations
function pass = check_limit(ac, dc, rated, is_halfbridge)
    if (abs(ac)*sqrt(2) + abs(dc) >= rated) || (abs(ac)*sqrt(2) - abs(dc) < 0 && is_halfbridge)
        pass = 1;
    else
        pass = 0;
    end
end