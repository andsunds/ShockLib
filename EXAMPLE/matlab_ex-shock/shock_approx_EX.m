%% Shock approximation script
% Example script on how to use the method find_approx to find the
% approximation coefficients from discrete data of psi, E and rho.
%
% (c) Andréas Sundström, 2018
clc;clf;clear all

%Adding the path to the library containing the methods of interest.
addpath(strcat(getenv('SHOCKLIB'), '/matlabLib'));
% Path to the data
DATA_path = strcat(getenv('SHOCKLIB'), '/EXAMPLE/in-data/');
data_name = 'psi-E-rho_M1-300_tau200_MB.tsv';
data      = load(strcat(DATA_path,data_name));
% Human-readable variable names
X=data(:,1); psi=data(:,2); E=data(:,3); rho=data(:,4);



FFT_sample_points = 1024;   % Numer of sample points to unse in the FFT
exp_cutoff        = 5.4;     % x-value above which psi(x) is asympotically exponential (decay)
poly_deg          = 14;     % degree of polynomial to use in the upstream approximation

% Calculating the approximatino coefficients
% If the fit is bad (extra noticable in E and rho), try adjusting
% exp_cutoff and/or poly_deg
[Cm,K, Cp]=find_approx(X,psi,E,rho, FFT_sample_points, exp_cutoff, poly_deg);
% Cm - is the fourier coefficient for the downstream
% K  - is the base frequency of the Fourier series
% Cp - is the polynomial coefficient of the numerator polynomial in the
%      upstream approximation


% Plotting the approximated functions
x=linspace(min(X),max(X),5000); % dense sample in x

% The number of Fourier (cosine) terms to use in the approximation,
% anything above 10 should work fine.
N_Fourier_terms = 10; 

%Going back to psi, E, and rho based on the approximation coefficients
%calcluated above. This time we can use an arbirary x
[psi0, E0, rho0]= get_approx(x,Cm,K,N_Fourier_terms, Cp);

% Plotting the original data (dashed) and the approximation (line).
% In general, psi looks good, while E and rho are pretty bad. For Gkyl the
% most important variables are E and psi (if rho looks a bit bad, it
% dosen't matter that much). 
% Uncomment the data set you want to view. 
hold on
%plot(x,psi0, X,psi,'--')
plot(x,E0, X,E,'--')
%plot(x,rho0, X,rho,'--')

% Creating a save name with the same info as the data name, but without the
% "-E-rho" part, and also beginning with "COEFS-" for clarity.
save_name = regexp(data_name, '-E-rho','split');
save_name = strcat(DATA_path,'COEFS-',save_name{:});
%Method for saving the data in the right format.
save_coefs(save_name,Cm(1:N_Fourier_terms),K,Cp);



