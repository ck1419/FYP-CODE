%Checks whether the equations for newton raphson are fulfilled (f = 0)
function pass = f11_eqn_check(x, Pcon, Xarm, R, Vgrid_RE, Vgrid_IM, Vhvdc, Pgrid, Qgrid, tolerance)
    %Checks the final iteration results with the tolerance given
    check = f11(x, Pcon, Xarm, R, Vgrid_RE, Vgrid_IM, Vhvdc, Pgrid, Qgrid);
    pass = check < tolerance;
    
    %Counts the amount of passes
    counter = 0;
    for i = 1:length(pass)
        counter = counter + pass;
    end

    %Displays results of check
    if counter == length(pass)
        fprintf('PASSED: f(x) = 0 \n')
    else
        fprintf('FAILED: f(x) = 0 \n')
    end
end