function [M_lim, phi_maxmin_out] = find_Mlim(Sh_handle,args, Tau,dM,N_tries, MaxMin, stepBack, output_file)


if nargin<=6
    stepBack=0;
end
if nargin<=7
    output_file='data/Mlim-tmp.tsv';
end



%Clean up of old tmp file
if exist(output_file, 'file')==2
    delete(output_file)
end
log_file=fopen([output_file(1:end-3) 'log'],'w');


fprintf(log_file,'\n   Scan started at %s\n\n', datestr(now,'HH:MM:SS.FFF'));
fprintf(log_file,'Scanning in Tau=%d:%1.2f:%d, for %s-limit, with stepBack=%1.2f\n\n',...
    Tau(1),Tau(2)-Tau(1),Tau(end),MaxMin,stepBack);


Mach_index=find(strcmp(args,'Mach'))+1;
M0=args{Mach_index};
phim_maxmin_in_index=Mach_index+1;
phi_maxmin_in=args{phim_maxmin_in_index};

N=length(Tau);


if strcmp(MaxMin,'min')
    i0=N_tries;
    first_try=1;
    UpDown=+1;
elseif strcmp(MaxMin,'max')
    i0=1;
    first_try=N_tries;
    UpDown=-1;
else
    warning('MaxMin must be either ''max'' or ''min''\n')
    M_lim=[]; phi_maxmin_out=[];
    return
end

%inint
M_lim=zeros(N,1);
phi_maxmin_out=zeros(N,2);



for i=1:N
    fprintf(log_file,'=============================================================\n');
    fprintf(log_file,'*************************************************************\n');
    fprintf(log_file,'=============================================================\n');
    tau=Tau(i);
    
    args{4}=tau;
    args{Mach_index}=M0+stepBack*dM*UpDown;
    args{phim_maxmin_in_index}=phi_maxmin_in;
    
    fprintf(log_file,'tau = %d,   \t M = %1.5f \n',tau,M0);
    Sh0=Sh_handle(args{:});
    while isnan(Sh0.phimax)
        args{Mach_index}=args{Mach_index}+dM*UpDown;
        %args{phim_maxmin_in_index}(2)=args{phim_maxmin_in_index}(2)*0.99;
        fprintf(log_file,'Couldn''t find shock at this Mach #, now trying M = %1.4f\n',args{Mach_index});
        Sh0=Sh_handle(args{:});
    end
    
    M0=args{Mach_index};
    Sh_cell=Shock_pkg.Mach_sweep(Sh_handle, Sh0, i0, N_tries, M0, dM, '/dev/null');
    
    if ~isempty(Sh_cell{first_try}) && ~isnan(Sh_cell{first_try}.phimax)
        M0=Sh_cell{first_try}.Mach;
        phi_maxmin_in=[Sh_cell{first_try}.phimax, Sh_cell{first_try}.phimin];
    else
        try
            j=first_try;
            while isempty(Sh_cell{j}) || isnan(Sh_cell{j}.phimax)
                j=j+1*UpDown
            end
            M0=Sh_cell{j}.Mach;
            phi_maxmin_in=[Sh_cell{j}.phimax, Sh_cell{j}.phimin];
        catch
            warning('Cells seems empty.');
            fprintf(log_file,'Cells seems empty.');
            M0=NaN;phi_maxmin_in=[];
        end
    end
    M_lim(i)=M0;
    phi_maxmin_out(i,:)=phi_maxmin_in;
    data=[tau,M0,phi_maxmin_in];
    save(output_file, 'data', '-ascii','-append')
    
    fprintf(log_file,'Scan at tau=%d ended at at %s\n',tau, datestr(now,'HH:MM:SS.FFF'));
end

fprintf(log_file,'=============================================================\n');
fprintf(log_file,'*************************************************************\n');
fprintf(log_file,'=============================================================\n');
fprintf(log_file,'\n   Scan of all tau ended at %s\n', datestr(now,'HH:MM:SS.FFF'));
end

