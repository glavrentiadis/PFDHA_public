function [ fphase ] = IntPhaseDer(k_array,fphase_deriv,dx,integ_order)
%IntPhaseDer creates a single sided phase spectra based on the phase and phase derivatives
%Version 2: Allows for higher order of integration
% Input Arguments:
%   k_array = wavenumber array
%   fphase_deriv = phase derivative array
%   dx = distance interval
%   integ_order = order of integration

assert(size(k_array,2)==1,'Error. k_array must be a vector')
assert(size(k_array,1) == size(fphase_deriv,1),'Error. k_array and fphase_deriv must have the same length')
assert((max(diff(k_array))-min(diff(k_array))<1e-9),'Error. wavenumber array must have constant spacing')
assert(abs(k_array(1))<1e-5,'Error. k_array array must start with the zero wavenumber')

dk = mean(diff(k_array));

if ~exist('integ_order','var')
    integ_order = 1;
end
integ_order = min(integ_order,length(k_array));

%first wavenumber
fphase = nan(length(k_array),1);
fphase(1) = 0;

%other wavenumber
switch integ_order
    case 1
        for j = 2:length(k_array)
			fphase(j) = fphase(j-1)+dk*fphase_deriv(j);
        end
    case 2
        %2nd wavenumber phase
        fphase(2) = fphase(2-1)+dk*fphase_deriv(2);
        %3rd till last wavenumber phase
        for j = 3:length(k_array)
			fphase(j) = 2/3*(-1/2*fphase(j-2)+2*fphase(j-1)+dk*fphase_deriv(j));
		end
    case 3
        %2nd wavenumber phase
        fphase(2) = fphase(2-1)+dk*fphase_deriv(2);
        %3rd wavenumber phase
        fphase(3) = 2/3*(-1/2*fphase(3-2)+2*fphase(3-1)+dk*fphase_deriv(3));
        %4rd till last wavenumber phase
        for j = 4:length(k_array)
            fphase(j) = 6/11*(1/3*fphase(j-3)-3/2*fphase(j-2)+3*fphase(j-1)+dk*fphase_deriv(j));
        end
    case 4
        %2nd wavenumber phase
        fphase(2) = fphase(2-1)+dk*fphase_deriv(2);
        fphase(3) = 2/3*(-1/2*fphase(3-2)+2*fphase(3-1)+dk*fphase_deriv(3));
        fphase(4) = 6/11*(1/3*fphase(4-3)-3/2*fphase(4-2)+3*fphase(4-1)+dk*fphase_deriv(4));
        %5th till last wavenumber phase
        for j = 5:length(k_array)
            fphase(j) = 12/25*(-1/4*fphase(j-4)+4/3*fphase(j-3)-3*fphase(j-2)+4*fphase(j-1)+dk*fphase_deriv(j));
        end
    otherwise
        error('Invalid integ_order')
end
fphase = mod(fphase,2*pi);


if abs(k_array(end)-1/(2*dx))<1e-9
	fphase(end) = 0; %Nyquist wavenumber
elseif (k_array(end)>-1/(2*dx))<1e-4
	assert(k_array(end)<=1/(2*dx),'Error. Invalid wavenumber')
end

end

