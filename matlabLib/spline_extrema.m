function [x0,y0] = spline_extrema(pp_spline)
% A function that finds all the exema of a spline object by analytically
% calculating the roots of the derivative of each spline, and then testing
% whether that root lies within the spline range x_n to x_{n+1}.
%
% Each spline has the form of a cubic polynomial, which means that the
% roots of its derivatives are given by the quadratic formula. This gives
% all the pseudo-roots. For a pseudo-root to be a proper root it must also
% lie within the range of its polynomial.
%
% (c) Andréas Sundström, 2018

% shorthand variables
Coefs=pp_spline.coefs;
xn=reshape(pp_spline.breaks,[],1); % This must be a column vector to be 
                                   % added together with the columnvecotrs
                                   % of pseudo-roots, later on.
% The pseudo-roots of the spline derivative
pseudo_x0 = ...
    (-Coefs(:,2)+sqrt(Coefs(:,2).^2-3*Coefs(:,1).*Coefs(:,3)).*[1,-1])./...
    (3*Coefs(:,1));

% Checking for proper roots
L=size(pseudo_x0,1);
for i=1:L
    % Removing all pseudo-roots to the left of their range
    pseudo_x0(i, pseudo_x0(i,:)<0 )=NaN; 
    % Removing all pseudo-roots to the right of their range
    pseudo_x0(i, pseudo_x0(i,:)>xn(i+1)-xn(i) )=NaN;
    % Removing any complex pseudo-roots
    pseudo_x0(i, imag(pseudo_x0(i,:))~=0 )=NaN;
end

% Finds all the proper roots. Note that pseudo_x0 is an L-by-2 matrix,
% since there are two roots of a quadratic equation.
I_proper=find(~isnan(pseudo_x0));
% To get the corresponding indices in the xn vector (length L), we have to
% take I_proper mod L:
I_proper_xn=mod(I_proper,L);
I_proper_xn(I_proper_xn==0)=L;% However we still want L to be L, and not 0.

% Finally each proper root is set in the range of its polynomial, and then
% all roots are sorted.
x0=sort( pseudo_x0(I_proper) + xn(I_proper_xn) );

% If asked for, this function also provides the spline values at the
% extema.
if nargout==2
    y0=ppval(pp_spline, x0);
end

end % end function

