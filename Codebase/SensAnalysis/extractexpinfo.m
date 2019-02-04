function [ tyout, V_0, sgn, NT, experiment, diagout, absout, R, Q, E, I,STS,STSnorm, sens ] =...
    extractexpinfo( name, path, normalize_sens_factor,Np,mu)
%given a file name, this function extracts the info from the sensitivities
%calculated by Dong Zhang
%Author:Dylan kato
%10/26/18
    p.volt_max = 4.2;
    p.volt_min = 3;
    
    load(horzcat(path,name));
    index=str2double(name(1:end-4));
    NT=length(I);
    sgn=0;
    if mean(I>0)==1
        sgn=1;
    end
    if mean(I<0)==1
        sgn=-1;
    end
        S3 = S_CASADI;
    for k = 1:length(S_CASADI)
        S3(k,:) = normalize_sens_factor .* S3(k,:);
    end
    
        [Q,R,E] = qr(S3,0);
    % Extract Diagonal
    Di = diag(R);
    for j = 1:Np
         D(j) = Di(find(E == j));
    end
    STS=S3'*S3;
    ttt=inv(diag(sqrt(diag(STS))));
    STSnorm=abs(ttt*STS*ttt);
    
    experiment = index*(abs(D)>(mu*sqrt(NT)));
    tyout=T_amb(1);
    diagout=diag(S3'*S3);
    absout=abs(D);
    sens=S3;
end


% load(t);
%     NT=length(I);
%     sgn=0;
%     if mean(I>0)==1
%         sgn=1;
%     end
%     if mean(I<0)==1
%         sgn=1;
%     end
%     %
%     T(i-2) = T_amb(1);
%     V_0_mat(i-2)=V_0;
%     Sign(i-2)=sgn;
%     NT_mat(i-2)=NT;
%     %
%     S3 = S_CASADI;
%     for k = 1:length(S_CASADI)
%         S3(k,:) = normalize_param .* S3(k,:);
%     end
%     
%     [Q,R,E] = qr(S3,0);
%     % Extract Diagonal
%     Di = diag(R);
%     for j = 1:Np
%          D(j) = Di(find(E == j));
%     end
%     experiment = index*(abs(D)>(mu*sqrt(NT)));
%     A(:,i-2) = experiment;
%     Sens_Mag(:,i-2) = diag(S3'*S3);
%     Sens_orth(:,i-2) = abs(D); 
%     AllR{i-2} = R;
%     AllQ{i-2} = Q;
%     AllE{i-2} = E;
%     AllI{i-2} = I;
