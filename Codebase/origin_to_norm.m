function [varargout] = origin_to_norm(mode,selected_param, bounds, Sel_k)

    % By Saehong Park (2017.10.29)
    % Edited by Zach Gima (2018-11-9)
    % Convert original parameters to normalized parameters for simulation
    % or normalize sensitivity
    
    % originS_to_normS functionality integrated into origin_to_norm 
    
    % Sister function: norm_to_origin.m

    %%

    if(length(selected_param) ~= sum(nonzeros(Sel_k)))
        error('Check the number of parameters')
    end

    origin_param = NaN(length(Sel_k),1);
    origin_param(find(Sel_k)) = selected_param;

    emptyblk = NaN;

    if strcmp(mode,'param') %normalizing parameters
        origin2norm = [ (log10(origin_param(1)) - log10(bounds.min(1))) / (log10(bounds.max(1))-log10(bounds.min(1))); ... %Dsn (log)
                         (log10(origin_param(2)) - log10(bounds.min(2))) / (log10(bounds.max(2))-log10(bounds.min(2))); ... %Dsp (log)
                         (origin_param(3) - bounds.min(3)) / (bounds.max(3) - bounds.min(3)); ... %Rsn (Min/Max)
                         (origin_param(4) - bounds.min(4)) / (bounds.max(4) - bounds.min(4)); ... %Rsp (Min/Max)
                         emptyblk;
                         emptyblk;
                         (log10(origin_param(7)) - log10(bounds.min(7))) / (log10(bounds.max(7))-log10(bounds.min(7))); ... %sig_n (log)
                         (log10(origin_param(8)) - log10(bounds.min(8))) / (log10(bounds.max(8))-log10(bounds.min(8))); ... %sig_p (log)
                         (origin_param(9) - bounds.min(9)) / (bounds.max(9) - bounds.min(9));   ... %D_e (Min/Max)
                         (origin_param(10) - bounds.min(10)) / (bounds.max(10) - bounds.min(10));   ... %eps_e_n (Min/Max)
                         (origin_param(11) - bounds.min(11)) / (bounds.max(11) - bounds.min(11));   ... %eps_e_s (Min/Max)
                         (origin_param(12) - bounds.min(12)) / (bounds.max(12) - bounds.min(12));   ... %eps_e_p (Min/Max)
                         (origin_param(13) - bounds.min(13)) / (bounds.max(13) - bounds.min(13));   ... %kappa (Min/Max)
                         (origin_param(14) - bounds.min(14)) / (bounds.max(14) - bounds.min(14));   ... %t_plus (Min/Max)
                         (origin_param(15) - bounds.min(15)) / (bounds.max(15) - bounds.min(15));   ... %dactivity (Min/Max)
                         (log10(origin_param(16)) - log10(bounds.min(16))) / (log10(bounds.max(16))-log10(bounds.min(16))); ... %k_n (log)
                         (log10(origin_param(17)) - log10(bounds.min(17))) / (log10(bounds.max(17))-log10(bounds.min(17))); ... %k_p (log)
                         (origin_param(18) - bounds.min(18)) / (bounds.max(18) - bounds.min(18)); ... %Rfn (Min/Max)
                         (origin_param(19) - bounds.min(19)) / (bounds.max(19) - bounds.min(19)); ... %Rfp (Min/Max)
                         emptyblk;
                         (origin_param(21) - bounds.min(21)) / (bounds.max(21) - bounds.min(21)); %ce_0 (Min/Max)                     
                         (origin_param(22) - bounds.min(22)) / (bounds.max(22) - bounds.min(22)); %E.Dsn (Min/Max)                     
                         (origin_param(23) - bounds.min(23)) / (bounds.max(23) - bounds.min(23)); %E.Dsp (Min/Max)       
                         (origin_param(24) - bounds.min(24)) / (bounds.max(24) - bounds.min(24)); %E.kn (Min/Max)       
                         (origin_param(25) - bounds.min(25)) / (bounds.max(25) - bounds.min(25)); %E.kp (Min/Max)       
                         ];
                         % check (1) min locations (2) max locations (3) numbers (4)
%                          % normalization scheme correct
         norm_param = origin2norm(find(Sel_k));
        varargout{1} = norm_param;
       % varargout{1} = origin2norm;
    elseif strcmp(mode,'sens') %if normalizing sensitivity
        origin2norm = [ (1/log10(exp(1)))*log10(bounds.max(1) / bounds.min(1))*origin_param(1); ...
                (1/log10(exp(1)))*log10(bounds.max(2) / bounds.min(2))*origin_param(2); ...
                bounds.max(3) - bounds.min(3);%(3)
                bounds.max(4) - bounds.min(4);%(4)
                emptyblk;%(5)
                emptyblk;%(6)
                (1/log10(exp(1)))*log10(bounds.max(7)/bounds.min(7)) * origin_param(7); ... %(7)
                (1/log10(exp(1)))*log10(bounds.max(8)/bounds.min(8)) * origin_param(8); ... %(8)
                (bounds.max(9) - bounds.min(9)); ... %(9)
                bounds.max(10) - bounds.min(10); ... %(10)
                bounds.max(11) - bounds.min(11); ... %(11)
                bounds.max(12) - bounds.min(12); ... %(12)
                (bounds.max(13) - bounds.min(13)) ; ... %(13)
                (bounds.max(14) - bounds.min(14)) ; ... %(14)
                (bounds.max(15) - bounds.min(15)); ... %(15)
                (1/log10(exp(1)))*log10(bounds.max(16) / bounds.min(16)) * origin_param(16); ... %(16)
                (1/log10(exp(1)))*log10(bounds.max(17) / bounds.min(17)) * origin_param(17); ... %(17)
                bounds.max(18) - bounds.min(18); ... %(18)
                bounds.max(19) - bounds.min(19); ... %(19)
                emptyblk; %(20)
                bounds.max(21) - bounds.min(21); %(21)   
                bounds.max(22) - bounds.min(22); %(22)                  
                bounds.max(23) - bounds.min(23); %(23)                  
                bounds.max(24) - bounds.min(24); %(24)                  
                bounds.max(25) - bounds.min(25); %(25)                  
                ];
            
%         norm_paramS = origin2norm(find(Sel_k));
%         varargout{1} = norm_paramS; % 22x1
                varargout{1} = origin2norm; % 25x1
    else
        error('Specify either ''params'' or ''sens'' as the first argument to normalize parameters or sensitivity respectively')
    end
end