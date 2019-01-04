function [ fphase_dev ] = PhaseDerivBoore03(t,y,f2s,fun_fourier)
%PhaseDerivBoore03 computes the phase derivatives with the method outlined in Boore 2003
%   Reference: Boore, David M. "Phase derivatives and simulation of strong ground motions." Bulletin of the Seismological Society of America 93.3 (2003): 1132-1143.
% Input Arguments:
%   t:              time array
%   y:              signal array
%   f2s:            double sided frequency array
%   fun_fourier:	function handle to compute Fourier, format: c2s = fun_fourier(t,y)
% Output Arguments:
%   fphase_dev:     phase derivative vector


%number of frequencies
n = length(f2s);

%offset t to zero
t = t - min(t);

%product of time and signal
ty = t.*y;

%compute Fourier transforms
% [~,~,~,~,~,~,~,Y] = fun_fourier(t,y);
% [~,~,~,~,~,~,~,TY] = fun_fourier(t,ty); 
Y = fun_fourier(t,y);
TY = fun_fourier(t,ty); 


%numerator and denominator: equation 9, Boore 03
fphase_deriv_num = -2*pi()*(real(Y).*real(TY)+imag(Y).*imag(TY)); %numerator, add minus to be compatible with matlab-fft (e^(-2*pi*k*f)) 3/11/2018
fphase_deriv_den = abs(Y).^2; %denominator

%phase deriavtive without smoothing
fphase_dev_no_smth = fphase_deriv_num./fphase_deriv_den;
%smoothing window size
win_mov_avg = max(round(0.005*length(t)),3);

%phase deriavtive with smoothing
fphase_dev = TriangMovMean(fphase_dev_no_smth,win_mov_avg);


end

