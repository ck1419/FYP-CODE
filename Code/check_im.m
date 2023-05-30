%Checks to make sure the variables only contains a magnitude
function pass = check_im(x)
    %Checks for no imaginary components
    pass = imag(x) == 0;

    %Counts the amount of passes
    counter = 0;
    for i = 1:length(pass)
        counter = counter + pass;
    end

    %Displays results of check
    if counter == length(pass)
        fprintf('PASSED: no imaginary component \n')
    else
        fprintf('FAILED: no imaginary component \n')
    end
end


