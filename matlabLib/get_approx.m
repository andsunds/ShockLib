function [psi,E,rho] = get_approx(X,Cm,K,nf, Cp)
% A function that shows the vaule of the approximation in the points X,
% given the approximation coefficients Cm and Cp. K is the frequency of the
% downstream oscillations, and nf is the number or Fourier modes to use in
% the approximation (usually Cm is 1024 elements long).
%
% (c) Andréas Sundström, 2018

X=reshape(X,[],1); %Reshape to ensure column vector

psi=zeros(size(X)); %inint

Im=find(X<=0); Ip=find(X>0);

%%%%%%%%%%%%%%% Downstream (X<=0) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The DS approximation is a cosine series with base frequency K, and
% coefficients Cm.
if ~isempty(Im) % Check that there are negative values
    k= K*(0:(nf-1)); % vector of frequencies
    % Summing the cosine series
    psi(Im) = cos(X(Im).*k)*Cm(1:nf); 
end

%%%%%%%%%%%%%%% Upstream (X>=0) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The US approximation consista of a rational function, r=p/q, fitted to
% the log of psi. The approximation coeffiitens are however only for p(x),
% since q(x) is simply x^{n-1}+1, where n is the degree of p.
if ~isempty(Ip)% Check that there are positive values
    Cp=reshape(Cp,[],1); %Reshape to ensure column vector
    n=length(Cp)-1; % Highest power term in the numerator polynomial
    a=0:n; %All the terms of the numerator polynomial
    
    % log(psi(x)) ~ r(x) = p(x)/q(x)
    p=X(Ip).^(a)*Cp;    % p(X)
    q=X(Ip).^(n-1)+1;   % q(X)
    r=p./q;             % r(X)
    
    psi(Ip) = exp(r);   % psi(X)
end


%%%%%%%%%%%%%%% E and rho %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If wanted, this functyion can also calcuate the analytical first and
% second derivatives of the approximation functions. 
if nargout>=2 % If asked for, also give E(X)
    E=zeros(size(X));
    if ~isempty(Im) % Check that there are negative values
        % plus sing on the sine since E=-dpsi/dx
        E(Im)=sin(X(Im).*k)*(Cm(1:nf).*k');
    end
    if ~isempty(Ip) % Check that there are positive values
        b=0:(n-1); %Short cut for the powers of dp
        dp=X(Ip).^(b)*(Cp(2:end).*a(2:end)');   % dp = p'(X)
        dq=(n-1)*X(Ip).^(n-2);                  % dq = q'(X)
        dr=r.*(dp./p-dq./q);                      % dr = r'(X)
        % - here since  E=-dpsi/dx
        E(Ip)=-psi(Ip).*dr;                     % E(X)
    end
end
if nargout>=3 % If asked for, also give rho(X)
    rho=zeros(size(X));
    if ~isempty(Im) % Check that there are negative values
        % plus sign on the cosine since rho=-d2psi/dx2
        rho(Im)=cos(X(Im).*k)*(Cm(1:nf).*k.^2');
    end
    if ~isempty(Ip) % Check that there are positive values
        c=0:(n-2);%Short cut for the powers of ddp
        ddp=X(Ip).^(c)*(Cp(3:end).*b(2:end)'.*a(3:end)');     % ddp = p''(X)
        ddq=(n-1)*(n-2)*X(Ip).^(n-3);                         % ddq = q''(X)
        ddr=(ddp-dr.*dq-r.*ddq)./q + (r.*dq.^2-dp.*dq)./q.^2; % ddr = r''(X)
        
        rho(Ip)=E(Ip).*dr -psi(Ip).*ddr;                      % rho(X)
    end
end

end


