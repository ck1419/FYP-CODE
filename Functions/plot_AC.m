%Simple AC Phasor Plot
function plot_AC(Vac, Iac, heading, pos, msg)
    %Initializes plot and forces equal axis
    figure
    hold on
    grid on
    axis equal

    %Plots points as arrows with a legend
    quiver(0, 0, real(Vac)/1000, imag(Vac)/1000)
    quiver(0, 0, real(Iac), imag(Iac))
    xlabel('Real')
    ylabel('Imaginary')
    legend('Vac (kV)', 'Iac (A)')
    title(heading)

    %Adds text box with custom message at a custom position
    annotation('textbox', pos,'String',msg,'FitBoxToText','on');
end


