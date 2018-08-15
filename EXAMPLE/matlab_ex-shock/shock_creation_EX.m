%% Shock creation example
% Examlpe scrip on how to use the shock package Shock_pkg_new to create a
% shock and then using the shock object to calcualte psi of x.
% NOTE: the Shock_pkg_new is with the new normalization, hence the use of
% psi as the variable name for the electrostatic potential.
%
% (c) Andréas Sundström, 2018
clc; clear
% 
SH_path=strcat(getenv('SHOCKLIB'), '/matlabLib');
addpath(SH_path)
import Shock_pkg.*

% The shock input parameters
Z=[1];      % Ion charge number
m=[1];      % Ion mass
n=[1];      % Unperturbed ion density
tau=200;    % Electron to ion temp ratio: T_e/T_i (no Z_i here)
Mach=1.3;   % Mach number
tol=1e-9;   % Wanted numerical tolerance
% For multi-ion species Z, m, and n can be vectors containing the values of
% each species.

psimaxmin_in=[0.6775 0.4929]; % The initial guess of psimax and psimin,
% In this case, the initial guess is spot on, for speed of evaluation.

% Creating the shock ocbject. It is in this process that the actual values
% of psimax and psi min
Sh0=Shock_MB(Z,m,n, tau, Mach, psimaxmin_in, tol);

disp(Sh0)%Displaying the result




%% solving for psi(x)
clc
% The range in which the psi(x) should be calcuated (not too big!)
xmin=-10;
xmax=10;

% This method numerically solves Poisson's equation to give psi(X) at some
% discrete points X, we also get E(X) an rho(X).
[X,psi,E,rho]=Sh0.find_psi(xmin,xmax);

%Then we plot the result
% Plotting psi(x) in figure 1
figure(1), clf 
plot(X,psi), xlabel('x'), ylabel('\psi(x)')

% Plotting E(x) in figure 2
figure(2), clf
plot(X,E), xlabel('x'), ylabel('E(x)')

% Plotting rho(x) in figure 3
figure(3), clf
plot(X,rho), xlabel('x'), ylabel('\rho(x)')




%% Saving the data
% You do not nessecarily need to save this data. You can just as well use 
% X, psi, E directly into the find_approx routine.
clc

%Path to where the data should be saved
DATA_path=strcat(getenv('SHOCKS'), '/scripts/Gkyl_matlab-lib/EXAMPLE/in-data/');

%Gather all the data into one matrix, where each column contains the data
%for X, psi(X), E(X), and rho(X), respectively.
data=real([X,psi,E,rho]); 

% Name of the file the data should be saved to, contains the value of 
% M and tau
save_name=sprintf('psi-E-rho_M%1.3f_tau%d_MB.tsv',Mach,tau);
save_name=regexprep(save_name,'(\d)\.(\d)','$1-$2'); %Changes decimal point "." to "-"
fprintf('Data will be saved in\n%s%s\n',DATA_path,save_name) %Prints out the full path

%Daves the data in the data path and the file name above
save([DATA_path,save_name],'data','-ascii','-double')



