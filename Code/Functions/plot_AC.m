%Simple AC Phasor Plot
function plot_AC(Vac, Iac, heading, pos, msg)
    figure
    hold on
    grid on
    axis equal

    quiver(0, 0, real(Vac)/1e2, imag(Vac)/1e2)
    quiver(0, 0, real(Iac), imag(Iac))

    annotation('textbox', pos,'String',msg,'FitBoxToText','on');
    xlabel('Real')
    ylabel('Imaginary')
    legend('Vac/100 (V)', 'Iac (A)')
    title(heading)
end


