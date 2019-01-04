function [f,famp,fphase,fphder,c1s,f_Nyq,f2s,c2s] = Rft(t_array,s_array,dt_mean,flag_phderiv,oversamp_ratio,wt_reg,wt_pwr)
%Rft performs a regularized least squares Fourier transform of uniformly spaced signals 
% Input Arguments:
%   t_array:        time/distance array
%   s_array:        signal array
%   dt_mean:        (optional) time/distance interval for the development of frequency/wavenumber array
%   flag_phderiv:   if true return phase derivatives
%   oversamp_ratio: frequency/wavenumber oversampling ratio
%   wt_reg:         regularization weight (parameter alpha)
%   wt_pwr:         relative weight low vs high frequencies/wavenumbers (parameter alpha)
% Output Arguments:
%   f:              1-sided frequency/wavenumber vector
%   famp:           1-sided amplitude vector
%   fphase:         1-sided phase vector
%   fphder:         1-sided phase derivative vector
%   c1s:            1-sided Fourier complex coefficients
%   f_Nyq:          Nyquist frequency/wavenumber
%   f2s:             2-sided frequency/wavenumber vector [-f_nyq,f_nyq), compatible with fftshift 
%   c2s:            2-sided Fourier complex coefficients


%assign default values if not specified
%time/distance interval
if isnan(dt_mean)
    dt_mean = mean(diff(t_array));
end

%check input validity
assert(all(size(t_array)==size(s_array)),'Error. signal and time/distance arrays must have the same size')
assert(abs(oversamp_ratio-ceil(oversamp_ratio))<1e-4,'Error. oversamp_ratio must be an integer')

%offset time or distance array to start from zero
t_array = t_array - min(t_array); 

%frequency sampling
dt = dt_mean/oversamp_ratio;
fs = 1/dt;

%Number of points
n_pt = length(t_array)*oversamp_ratio;

%frequency vector
f = (fs*(0:floor(n_pt/2))/n_pt)';
f_Nyq = fs/2;
f_Nyq(2) = f_Nyq/oversamp_ratio;

%Regularized Fourier Transform
%LS-regression matrix
[f4reg,trig_polyn4reg] = TrigPolyn(t_array,f,dt);

%Regularization weights
wt = abs(f4reg'); %frequency proportional weights
wt(wt==0) = min(wt(wt~=0)); %for zero weights assign them to smallest non-zero
wt = wt.^wt_pwr;
wt = wt_reg*f(2)*wt./min(wt); %normalization 
%perform regression
trig_polyn_amp = (trig_polyn4reg'*trig_polyn4reg+diag(wt)'*diag(wt))\(trig_polyn4reg'*s_array);

%compute Fourier amplitude and phase
[f_out,famp,fphase,c1s] = CalcFourierFromTrigPolyn(f4reg',trig_polyn_amp);

%compute two sided complex spectrum (-f_Nyq, f_Nyq-df)
f2s = fs/n_pt*(-floor(n_pt/2):floor((n_pt-1)/2))';
c2s = nan(n_pt,1);
c2s(1+(0:floor(n_pt/2)),1) = conj(flipud(c1s)); % -f_Nyq to 0
c2s(end-(floor((n_pt-1)/2)-1):end,1) = c1s(2:(floor((n_pt-1)/2)+1)); % 0 to f_Nyq-df

%compute phase derivatives
if flag_phderiv == true(1)
    func4Fourier = @(t,x) Rft(t,x,dt_mean,false(1),oversamp_ratio,wt_reg,wt_pwr);
    func4FourierC2s = @(t,x) CompCoefs(func4Fourier,t,x);

    %compute phase derivatives 2sided spectrum
    fphase_deriv2s_boore = PhaseDerivBoore03(t_array,s_array,f2s,func4FourierC2s);    
    %phase derivative single-sided
    fphder = flipud(fphase_deriv2s_boore(1:floor(n_pt/2)+1)); %03/05/18, flip the array because it corresponds to the negative freq, don't use the positive because f_Nyq is missing.
else
    fphder = nan(size(f));
end


end


function [c2s,f2s] = CompCoefs(fun_FT,t,x)

[~,~,~,~,~,~,f2s,c2s] = fun_FT(t,x);

end