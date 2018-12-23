function [f2s,C2s,f_opt] = Create2SidedSepct_RealSig(freq1s,famp1s,fphase1s,dx,flag_Nyq_freq )
%Create2SidedSepct returns the complex Fourier coefficients of two-sided spectrum based on
%amplitude and phase single-sided spectra of a real signal
% Input arguments
%	freq1s:		frequency array
%	famp1s:		amplitude array
%	fphase1s:	Phase array
%	dx: 		distance interval
%	flag_Nyq_freq: if true check if Nyquist frequency exist in freq and create C2s appropriately 

assert(all(freq1s>=0),'This function works with the positive 1sided spectrum');
assert(all(diff(freq1s)>0),'freq1s vector must be in an ascending order')
assert(freq1s(1)<1e-5,'First entry in freq1s, famp1s and fphase1s must correspond to the zero frequency/wavenumber');


if flag_Nyq_freq
	if abs(freq1s(end)-1/(2*dx))<1e-9 %last frequency corresponds to Nyquist frequency/wavenumber
		f_opt = 1;
	else
		f_opt = 2;
	end
else
	f_opt = 2;
end

%complex Fourier coefficients
C_1s = famp1s(2:end).*exp(1i*fphase1s(2:end)); %excluding the coefficient for f=0
C0 =famp1s(1)*exp(1i*fphase1s(1)); %complex coefficient for f=0

%confirm that the array of complex coefficient is consistent with wavenumber array
assert(abs(imag(C0))<1e-8,'Error. Invalid Fourier coefficient for zero frequency/wavenumber')
if f_opt == 1;
	assert(abs(imag(C_1s(end)))<1e-8,'Error. Invalid Fourier coefficient for Nyquist wavenumber')
end

if f_opt == 1 %even number of points unless flag_Nyq_freq = false
	%complex spectrum
	C2s = [conj(flipud(C_1s));C0;C_1s(1:end-1)]; %order compatible with fftshift, f_Nyq is in negative frequencies
	%frequency array
	f2s = [-flipud(freq1s);freq1s(2:end-1)]; 
elseif f_opt == 2 %odd number of points or flag_Nyq_freq = false
	%complex spectrum
	C2s = [conj(flipud(C_1s));C0;C_1s];
	%frequency array
	f2s = [-flipud(freq1s);freq1s(2:end)];
end


end

