function [ B_amp ] = ButterAmp(freq_array,Bamp_k0,KC,Npoles)
%ButterAmp creates a Butterworth amplitude spectrum based on Bamp_k0, KC
%and Npoles for the frequencies in freq_array
% Input arguments:
%	freq_array:	frequency/wavenumber array
%	Bamp_k0:	amplitude at the zero wavenumber
%	KC:			corner frequency/wavenumber
%	Npoles:		number of poles
% Output arguments:
%	B_amp:		Butterworth amplitude spectrum

B_amp = (1./sqrt(1+(freq_array./KC).^(2*Npoles))')';

B_amp = Bamp_k0.*B_amp;

end

