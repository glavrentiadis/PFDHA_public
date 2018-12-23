function [ Fphase1s ] = IntPhaseDer(k_array,fphase_deriv,dx,integ_order)
%IntPhaseDer creates a single sided phase spectra based on the phase and phase derivatives
%Version 2: Allows for higher order of integration
%   k_array = wavenumber array
%   fphase = phase array
%	fphase_deriv = phase derivatives array
%	dx = distance interval
%   integ_order = order of integration

assert(size(k_array,2)==1,'Error. k_array must be a vector array')
assert(size(k_array,1) == size(fphase_deriv,1),'Error. k_array and fphase_deriv must have the same length')
assert((max(diff(k_array))-min(diff(k_array))<1e-9),'Error. wavenumber array must have a constant spacing')
assert(abs(k_array(1))<1e-5,'Error. k_array array must start with the zero wavenumber')

dk = mean(diff(k_array));

if ~exist('integ_order','var')
    integ_order = 1;
end
integ_order = min(integ_order,length(k_array));

%first wavenumber
Fphase1s = nan(length(k_array),1);
Fphase1s(1) = 0;

%other wavenumber
switch integ_order
    case 1
        for j = 2:length(k_array)
			Fphase1s(j) = Fphase1s(j-1)+dk*fphase_deriv(j);
        end
    case 2
        %2nd wavenumber phase
        Fphase1s(2) = Fphase1s(2-1)+dk*fphase_deriv(2);
        %3rd till last wavenumber phase
        for j = 3:length(k_array)
			Fphase1s(j) = 2/3*(-1/2*Fphase1s(j-2)+2*Fphase1s(j-1)+dk*fphase_deriv(j));
		end
    case 3
        %2nd wavenumber phase
        Fphase1s(2) = Fphase1s(2-1)+dk*fphase_deriv(2);
        %3rd wavenumber phase
        Fphase1s(3) = 2/3*(-1/2*Fphase1s(3-2)+2*Fphase1s(3-1)+dk*fphase_deriv(3));
        %4rd till last wavenumber phase
        for j = 4:length(k_array)
            Fphase1s(j) = 6/11*(1/3*Fphase1s(j-3)-3/2*Fphase1s(j-2)+3*Fphase1s(j-1)+dk*fphase_deriv(j));
        end
    case 4
        %2nd wavenumber phase
        Fphase1s(2) = Fphase1s(2-1)+dk*fphase_deriv(2);
        Fphase1s(3) = 2/3*(-1/2*Fphase1s(3-2)+2*Fphase1s(3-1)+dk*fphase_deriv(3));
        Fphase1s(4) = 6/11*(1/3*Fphase1s(4-3)-3/2*Fphase1s(4-2)+3*Fphase1s(4-1)+dk*fphase_deriv(4));
        %5th till last wavenumber phase
        for j = 5:length(k_array)
            Fphase1s(j) = 12/25*(-1/4*Fphase1s(j-4)+4/3*Fphase1s(j-3)-3*Fphase1s(j-2)+4*Fphase1s(j-1)+dk*fphase_deriv(j));
        end
    otherwise
        error('Invalid integ_order')
end
Fphase1s = mod(Fphase1s,2*pi);


if abs(k_array(end)-1/(2*dx))<1e-9
	Fphase1s(end) = 0; %Nyquist wavenumber
elseif (k_array(end)>-1/(2*dx))<1e-4
	assert(k_array(end)<=1/(2*dx),'Error. Invalid wavenumber')
end

end

