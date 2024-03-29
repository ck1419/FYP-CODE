%3-Phase AC Phasor Plot
function plot_3AC(Vac_a, Iac_a, Vac_b, Iac_b, Vac_c, Iac_c, heading)
    %Initializes plot and forces equal axis
    figure
    hold on
    grid on
    axis equal

    %Plots points as arrows with a legend
    quiver(0, 0, real(Vac_a)/1000, imag(Vac_a)/1000)
    quiver(0, 0, real(Iac_a), imag(Iac_a))
    quiver(0, 0, real(Vac_b)/1000, imag(Vac_b)/1000)
    quiver(0, 0, real(Iac_b), imag(Iac_b))
    quiver(0, 0, real(Vac_c)/1000, imag(Vac_c)/1000)
    quiver(0, 0, real(Iac_c), imag(Iac_c))
    xlabel('Real')
    ylabel('Imaginary')
    legend('[A] Vac (kV)', '[A] Iac (A)', '[B] Vac (kV)', '[B] Iac (A)', '[C] Vac (kV)', '[C] Iac (A)')
    title(heading)
end