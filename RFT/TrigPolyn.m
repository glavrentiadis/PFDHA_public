function [ freq_array_output,trig_polyn_mat ] = TrigPolyn( t,f,dt_mean )
%TrigPolyn returns the trigonometric polynomials at times in t for frequencies in f
% Input Arguments:
%   t:  time/distance vector
%   f:  frequency/wavenumber vector

%Input assertions
assert(size(t,2)==1,'Error. t must be a column vector')
assert(size(f,2) == 1,'Error. f must be a column vector')
assert(max(diff(f)) - min(diff(f)) <= 1e-4,'Error. Frequencies in f must be sampled at a constant interval')
assert(abs(f(1)) <1e-10,'Error. f must start with the zero frequency')

%frequency interval
df = f(2);

%initialize trigonometric polynomial matrix
trig_polyn_mat = nan(length(t),length(f)*2);
lgi_only_cos = false(size(f)); %if true compute only cosine
if any(f==0) %if zero frequency exist remove zero frequency sine
    trig_polyn_mat = trig_polyn_mat(:,1:end-1);
    lgi_only_cos(f==0) = true(1);
end
%if Nyquist frquency exist remove Nyquist sine wave
if abs(f(end)-1/(2*dt_mean))<1e-9
    trig_polyn_mat = trig_polyn_mat(:,1:end-1);
    lgi_only_cos((abs(f-1/(2*dt_mean))<1e-9)) = true(1);
end

jj = 0; %index for tirg_polyn_mat 
for j = 1:length(f)
    if lgi_only_cos(j)
        jj = jj +1; %index of right most non-nan column
    	trig_polyn_mat(:,jj) = [cos(2*pi*f(j)*t)];
        freq_array_output(1,jj) = f(j);
    else
        jj = jj +2; %index of right most non-nan column
        trig_polyn_mat(:,jj+[-1,0]) = [-sin(2*pi*f(j)*t),cos(2*pi*f(j)*t)];
        freq_array_output(1,jj+[-1,0]) = f(j);
    end
end

trig_polyn_mat = df * trig_polyn_mat;

end

