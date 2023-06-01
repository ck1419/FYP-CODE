%Function to display results for single-phase single-arm
function f11_results_display(Vac, Iac, Vdc, Idc, phase_vac, phase_iac, phase_dif, Qcon)
    if sign(imag(Iac)) == -1
        sign_iac = ' - ';
    else
        sign_iac = ' + ';
    end

    if sign(imag(Vac)) == -1
        sign_vac = ' - ';
    else
        sign_vac = ' + ';
    end

    fprintf('\nFINAL ITERATION RESULTS: \n')
    disp(['Vac = ' num2str(real(Vac), '%3.3e') sign_vac num2str(abs(imag(Vac)), '%3.3e') 'i'])
    disp(['Iac = ' num2str(real(Iac), '%3.3e') sign_iac num2str(abs(imag(Iac)), '%3.3e') 'i'])
    disp(['Vdc = ' num2str(Vdc, '%3.3e')])
    disp(['Idc = ' num2str(Idc, '%3.3e')])

    fprintf('\nCALCULATED RESULTS: \n')
    disp(['Vac Phase = ' num2str(phase_vac) '°'])
    disp(['Iac Phase = ' num2str(phase_iac) '°'])
    disp(['Phase Difference = ' num2str(phase_dif) '°'])
    disp(['Converter Reactive Power = ' num2str(imag(Qcon), '%3.3e')])
    fprintf('\n')
end
