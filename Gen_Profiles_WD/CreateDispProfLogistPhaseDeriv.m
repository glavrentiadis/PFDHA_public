function [dist_array,disp_gen_all,iter2pass,k2s,C2s_all] = CreateDispProfLogistPhaseDeriv(dist_array,Bk0,KC,Np,st_fa,mu,s,nprof,th2acpt)
%CreateDispProfLogistPhaseDeriv generates slip profiles in the wavenumber
%domain, phase derivatives are assumed to follow a logistic distribution
% Input Arguments:
%   dist_array:	  along strike distance array
%   Bk0:		  amplitude of Butterworth filter
%   KC:			  corner wavenumber of Butterworth filter
%   Np:			  number of poles of Butterworth filter
%	st_fa:		  standard deviation of Fourier amplitude variates 
%   mu:			  center of logistic distribution
%   s:			  scale of logistic distribution
%   nprof:		  number of profiles to be generated
%   th2acpt:      threshold criteria for acceptance
%         th2acpt(1): maximum ratio of slip at the edge to average slip to be accepted
%         th2acpt(1): maximum ratio of negative slip to be accepted
% Output Arguments
%   dist_array:   distance array (same as input)
%   disp_gen_all: all generated slip profiles
%   iter2pass:    number of iterations to generate each profile
%   k2s:          2-sided wavenumber array of generated profiles
%   C2s_all:      2-sided complex Fourier coefficient array of gen. profiles

%amplitude ratio for significant wavenumbers
ratio_sig = 0.02; 
%number of iterations before interruption
limit_iters = 5000;
limit_iters_ns = 1000;

%distance array
assert(abs(dist_array(1)) <1e-9,'Error. dist_array must start at zero')
assert(max(diff(dist_array)) - min(diff(dist_array)) < 1e-6,'Error. distance array must be sampled at a constant interval.')
dx = mean(diff(dist_array)); %distance interval
n_points = length(dist_array)-1; %number of points
srl = (n_points)*dx; %surface rupture length

%wavenumber and Butterworth arrays
k_array = 1/(dx*n_points)*(0:floor(n_points/2))';
dk = 1/(dx*n_points);
bamp_array = ButterAmp(k_array,Bk0,KC,Np);
nk = length(k_array);
ik_sig = bamp_array > ratio_sig * Bk0; %significant wavenumbers
nk_sig = sum(ik_sig);

%initialize matrices to store generated profiles
disp_gen_all = nan(n_points+1,nprof);
C2s_all = nan(n_points,nprof);
iter2pass = nan(1,nprof);

%generate slip profiles
for j = 1:nprof
    flag_premature_stop = false(1);
    flag_repeat_both = true(1);
    n_counts = 0;
	%break wavenumber amplitude array into significant and non-significant parts
    k_s = k_array(ik_sig);
	% k_ns = k_array(~ik_sig);
    ba_s = bamp_array(ik_sig);
    ba_ns = bamp_array(~ik_sig);
	
	while and(~flag_premature_stop,flag_repeat_both)
        flag_repeat_sp = true(1); %flag to repeat significant part
        flag_repeat_nsp = true(1); %flag to repeat non-significant part
        while flag_repeat_sp
			%create amplitude, phase spectra for significant wavenumbers
			fa_s = ba_s.*10.^(normrnd(0,st_fa,[nk_sig,1])); fa_s(1) = ba_s(1);
			fpd_s = random('logistic',mu,s,[nk_sig,1]);
			%integrate phase derivatives
			fp_s = IntPhaseDer(k_s,fpd_s,dx,4);
            %create two sided complex spectrum
            [k2s_s,C2s_s] = Create2SidedSepct_RealSig(k_s,fa_s,fp_s,dx,true(1));
            if abs(sum(C2s_s))/Bk0 > ratio_sig/2 %if displacement at edge exceeds Bk0*ratio_sig/2 perform another iteration
                continue
            end
            
            %generate slip profile significant wavenumbers
			dist_s = linspace(0,srl,length(k2s_s)+1)';
            dx_s = dist_s(2);
            disp_gen_s = ifft(ifftshift(C2s_s))*length(k2s_s)*dk;
			%area under slip in negative direction
			area_negsp = -1*sum(disp_gen_s(disp_gen_s<0))*dx_s;
            %area under slip profile
			area_sp = sum(abs(disp_gen_s))*dx_s;
            
            %conditions to exit
            n_counts = n_counts+1;
            if area_negsp/area_sp  <= th2acpt(2)
                flag_repeat_sp = false(1);
                iter2pass(j) = n_counts;
                fprintf('\tDisp. Prof %i of %i (Low wavenumbers) \n',j,nprof)
            elseif n_counts >= limit_iters
                flag_premature_stop = true(1);
                flag_repeat_sp = false(1);
                flag_repeat_nsp = false(1);
                save('interupted_dispprof_generation.mat')
                fprintf('\tInterupted. Disp. Prof %i of %i\n',j,nprof)
            end
        end
		
		%iterate at the non-significant frequencies
        n_counts_s = n_counts;
        n_counts_ns = 0;
        n_counts_nz_edge = 0;
        
        while flag_repeat_nsp
			%create amplitude, phase spectra for non-significant wavenumbers
			fa_ns = ba_ns.*10.^(normrnd(0,st_fa,[nk-nk_sig,1]));
			fpd_ns = random('logistic',mu,s,[nk-nk_sig,1]);
			%create amplitude, phase spectra for non-significant all wavenumbers
			fa = [fa_s;fa_ns];
			fpd = [fpd_s;fpd_ns];
			%integrate phase derivatives
			fp = IntPhaseDer(k_array,fpd,dx,4);
		    %create two sided complex spectrum
            [k2s,C2s] = Create2SidedSepct_RealSig(k_array,fa,fp,dx,true(1));
            if abs(sum(C2s_s))/Bk0  > th2acpt(1) %if normalized displacement at edge exceeds threshold iterate perform another iteration
				n_counts_nz_edge = n_counts_nz_edge+1;
                if n_counts_nz_edge>=1000
                    flag_repeat_nsp = false(1);
                end
                continue
            end
            
            %generate slip profile significant wavenumbers
            disp_gen = ifft(ifftshift(C2s))*length(k2s)*dk;
            %disp_gen_s = invNUDFT_sine_cosine(f2s_s,C2s_s,dist_array);
            %area under slip in negative direction
			area_negsp = -1*sum(disp_gen(disp_gen<0))*dx;
            %area under slip profile
			area_sp = sum(abs(disp_gen))*dx;
            
            %conditions to exit
            n_counts_ns = n_counts_ns+1;
            n_counts = n_counts+1;
            if area_negsp/area_sp < th2acpt(2)
                flag_repeat_both = false(1);
                flag_repeat_nsp = false(1);
                iter2pass(j) = n_counts;
                fprintf('\tDisp. Prof %i of %i (High wavenumbers) \n',j,nprof)
            elseif n_counts >= limit_iters
                flag_premature_stop = true(1);
                flag_repeat_nsp = false(1);
                save('interupted_dispprof_generation.mat')
                fprintf('\tInterupted. Disp. Prof %i of %i\n',j,nprof)
            elseif n_counts_ns >= limit_iters_ns
                flag_repeat_nsp = false(1);
                fprintf('\tRepeating Disp. Prof %i of %i\n',j,nprof)
            end
        end		
		
	end
    
    if ~flag_premature_stop
        disp_gen_all(:,j) = [disp_gen;disp_gen(1)];
        C2s_all(:,j) = C2s;
    else
        disp_gen_all(:,j) = nan(size(dist_array));
        break
    end
end

end		


