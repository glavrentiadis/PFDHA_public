function [ f,famp,fphase,c1s ] = CalcFourierFromTrigPolyn( f_trigpolyn,coeffs_trigpolyn )
%CalcFourierFromTrigPolyn calculates the Fourier amplitude and phase spectra from the trigonometric polynomial coefficients

%input assertions
assert(size(f_trigpolyn,2) == 1,'Error. f_trigpolyn must be a column array')
assert(all(size(f_trigpolyn) == size(coeffs_trigpolyn)),'Error. f_trigpolyn and coeffs_trigpolyn must have the same size')

%initialize arrays
f = unique(f_trigpolyn);
famp = nan(size(f));
fphase = nan(size(f));

for j = 1:length(f)
    lgi_freq = f_trigpolyn == f(j);
    switch sum(lgi_freq) 
        case 1
            famp(j) = abs(coeffs_trigpolyn(lgi_freq));
            if sign(coeffs_trigpolyn(lgi_freq)) == 1
                fphase(j) = 0;
            else
                fphase(j) = pi;
            end
        case 2
            amp_sin_cos = coeffs_trigpolyn(lgi_freq);
            famp(j) = 0.5*sqrt(sum(amp_sin_cos.^2));
            fphase(j) = atan2(amp_sin_cos(1),amp_sin_cos(2));
        otherwise
            error('Error. Invalid entries in f array')
    end
end

c1s = famp.*exp(1i*fphase);

end

