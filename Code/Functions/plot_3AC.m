%3-Phase AC Phasor Plot
function plot_3AC(Vac_a, Iac_a, Vac_b, Iac_b, Vac_c, Iac_c, heading)
    %Initializes plot and forces equal axis
    figure
    hold on
    grid on
    axis equal

    %Plots points as arrows with a legend
    quiver(0, 0, real(Vac_a)/100, imag(Vac_a)/100)
    quiver(0, 0, real(Iac_a), imag(Iac_a))
    quiver(0, 0, real(Vac_b)/100, imag(Vac_b)/100)
    quiver(0, 0, real(Iac_b), imag(Iac_b))
    quiver(0, 0, real(Vac_c)/100, imag(Vac_c)/100)
    quiver(0, 0, real(Iac_c), imag(Iac_c))
    xlabel('Real')
    ylabel('Imaginary')
    legend('[A] Vac/100 (V)', '[A] Iac (A)', '[B] Vac/100 (V)', '[B] Iac (A)', '[C] Vac/100 (V)', '[C] Iac (A)')
    title(heading)

%     %Adds text box with custom message at a custom position
%     annotation('textbox', pos,'String',msg,'FitBoxToText','on');
end