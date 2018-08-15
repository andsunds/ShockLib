function [Cm, K, Cp] = find_approx(X,psi,E,rho,N_sample, exp_threshold,poly_deg)
% Function that calculates the approximation coefficients, Cm, K, and Cp,
% according to the approximation schemes described further down in this
% code. 
% This function takes in the true values of psi, E, and rho, this is the
% recommended way of running this function. However, this function can be
% run with only psi as input.
%
% (c) Andréas Sundström, 2018


% Finding all points beloning to eith negative or positive X values.
Im=find(X<=0); % All points with neg. X
Ip=find(X>=0); % All points with pos. X
% Some times there are two identical datapoints with X=0, if that is the
% case we remove one of them.
if X(Im(end-1))==X(Im(end)); Im=Im(1:end-1); end
if X(Ip(2))==X(Ip(1)); Ip=Ip(2:end); end

% Init of some humman-readable variable names
Xm=X(Im); psim=psi(Im);
Xp=X(Ip); psip=psi(Ip);

% Checks the number of input arguments, to see wheter all of psi, E, and
% rho, were given.
if nargin == 7
    Em=E([Im(1), Im(end)]);
    rho0=rho(Ip(1));
    drho0=( rho(Ip(2))-rho(Ip(1)) )/( Xp(2)-Xp(1) );
    use_exact_low_coef=true;
elseif nargin == 5
    % If only psi was given as the input
    % Approximaive E-filed, the important part is that Em(end)=0;
    Em=-[(psim(2)-psim(1))/(Xm(2)-Xm(1)); 0]*0.9;
    % Approximate rho
    rho0=0;%( psim(end-1)-2*psip(1)+psip(2) )/( (Xm(end)-Xm(end-1))*(Xp(2)-Xp(1)) )
    drho0=0;
    use_exact_low_coef=false;
    % Shifting down the other variables, if E and rho are not given.
    poly_deg=N_sample;
    exp_threshold=rho;
    N_sample=E;
end

%%%%%%%%%%%%%%% Downstream (X<=0) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Since the DS oscillation is periodic, we should be able to fit a Fourier
% series to it. We also know that psi'(0)=0, which leaves us with only
% cosine terms (real part of complex F-coefficients). However, as the data
% points are unevenly spaced, and not aligened with the exact period of
% oscillation, we first create a cubic spline interpolation. Then we have a
% dedicated function, spline_extrema(SPL), that finds the analytical
% min/max points of a cubic spline. We use this function to determine the
% wavelength of the first (closest to x=0) DS oscillation (the first since
% we don't want to include numerical errors introduced by the ODE solver
% that gave us the psi data). Then we create a densly sampled version of
% psi(x) in this first wavelength, using the spline interpolation, and
% find the F-coefficients with FFT. Finally we convert the complex
% F-transforn into a cosine transformation by only choosing the real part
% of the complex coefficients, and also doubling the non-zer frequency
% terms.


psimSP=spline(Xm,[-Em(1); psim; -Em(end)]); %Spline with defined end slopes
x0=spline_extrema(psimSP); %Finds the min/max of the spline
x=linspace(x0(end-2),x0(end),N_sample)'; %densly sampled gripoints.
psift=fft( ppval(psimSP,x) )/N_sample; % F-transf of the densly sampled data. 

fn=1:floor(N_sample/2); % The cosine transf only uses half of the coeeficients,
Cm=real(psift(fn)); % and also only the real part of them
% The non-zero frequency components have to be doubled due to the ones we
% trowed away.
Cm(2:end)=Cm(2:end)*2;
%CT(2:end)=CT(2:end)*(psim(end)-CT(1))/sum(CT(2:end));

% This is the base frequency used in the CT.
K=2*pi/(x(end)-x(1));


%%%%%%%%%%%%%%% Upstream (X>=0) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The idea behind this analytical approximatino is to fit a rational
% function, r(x)=p(x)/q(x), to log(psi(x)). We know that 
% psi(x) ~ exp(C1+C2*x) for x>>1 (C2<0). This means that if we choose 
% q(x) = x^(n-1)+1
% and
% p(x) = C2*x^n + C1*x^(n-1) + log(phi_max) + \sum_{k=2}^{m} A_k\mcode{E}, and
% then we would get the desired asymptotic behavior of r(x), as well as
% giving r(0) and r'(0) the proper values (r'(0)=0). We then determine A_k
% by least square fitting r(x) to log(psi(x)).

%Huamn-readable variable names
pmax=psip(1);
lpmax=log(pmax);

n=poly_deg(1); % Highest order term in p(x)
if length(poly_deg)==2
    m=poly_deg(2); % The highest order fitting term
else
    % If it's not given, we take the fitting terms to go all the way up to
    % the fixed terms of p(x).
    m=n-2; 
end

% Linear fitting to the asymptotic exponential decay
I=find(Xp>exp_threshold);
B=[ones(size(I)), Xp(I)]\log(psip(I));

alpha=2:m;
% p=@(x,A) C(2)*x.^n+C(1)*x.^(n-1)+lpmax + (x.^alpha)*A;
% q=@(x) (x.^(n-1)+1);
% r=@(x,A) p(x,A)./q(x);
r=@(x,A) ( B(2)*x.^n+B(1)*x.^(n-1)+lpmax + (x.^alpha)*A )./( x.^(n-1)+1 );


% Approximating the starting points with the first two non-zero Taylor
% coefficients (x^0 not included) of log(phi(x))*q(x).
A0=zeros(size(alpha))'; %init
a2=-rho0/pmax/2;
a3=-drho0/pmax/6;
A0(1:2)=[a2;a3];

% We only fit the lower order coeficients to the part of psi, not clearly
% just exponential. 
II=find(Xp<exp_threshold); 
if use_exact_low_coef
    % Least-square fitting the coefficients:
    AA=fminsearch(@(a) sum(( log(psip(II))-r(Xp(II),[A0(1:2);a]) ).^2), A0(3:end),...
        optimset('tolX',1e-4, 'tolFun',1e-8, 'maxFunEvals',50000));
    % Here we force r''(0) and r'''(0) to have the right values.
    AA=[A0(1:2);AA];
else
    % Least-square fitting the coefficients:
    AA=fminsearch(@(a) sum(( log(psip(II))-r(Xp(II),a) ).^2), A0,...
        optimset('tolX',1e-4, 'tolFun',1e-8, 'maxFunEvals',50000));
end

% Here, we're defining the coefficients for p(x) for all powers of x, even
% thoug some of them are known to be 0. This helps poring this to other
% functions.
% p(x) = \sum_{i=0}^{n} Cp_{i+1}*x^i.
Cp=zeros(n+1,1);%init
Cp(1)=lpmax;
Cp(end-1:end)=[B(1);B(2)];
Cp( (2:m)+1 )=AA;


end

