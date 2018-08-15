function [save_data] = save_coefs(save_path,Cm,K,Cp)
% A function for saving the approximation coefficients  in the right format
% for importing into Gkyl. This format is of the form:
%   -------------
%   |K      Cp_0|
%   |Cm_0   Cp_1|
%   |Cm_1   Cp_2|
%   |...    ... |
%   |0      Cp_n|
%   -------------
% which is what the Lua method requires. Notice that the first column
% starts with the frequency, and then the cosine coefficients.
%
% You should specify the full path (including file name) to where the coefs
% should be saved, in the input argument save_path. The DS coefficient
% vector, Cm, should only contain the coefficients that you want to save.
%
%
% (c) Andréas Sundström, 2018

% the sizes of the coefficient vectors
Nm=length(Cm);
Np=length(Cp);
% Initializing the data matrix to be saved. The length of the matrix is
% the longest of the two different coefficient vectors.
save_data=zeros(max([Nm+1, Np]), 2);
% the first column is for the downstream
save_data(1,1)=K; % first element is the frequency
save_data(2:(Nm+1),1)=Cm; % the rest of the elements are the Fourier coefs
% second column for the upstream
save_data(1:Np,2)=Cp;

%Prints out the full path to the file
fprintf('Approximation coefficients will be saved in\n%s\n\n',...
    save_path) 
%Asks if you want to save.
reply = input('Do you want to save? y/n [n]:','s');

%%if isempty(reply) || reply=='n'
if reply == 'y'
    % Saves the data if y
    save(strcat(save_path),'save_data','-ascii','-double');
    fprintf('Data saved!\n')
else
    % Otherwise no
    fprintf('Data not saved!\n')
end

end

