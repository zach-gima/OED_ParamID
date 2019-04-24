function [origin_param] = norm_to_origin(selected_param, bounds, Sel_k)

    % By Saehong Park (2017.10.29)
    % Convert normalized parameters to original parameters for simulation
    % Sister function: origin_to_norm.m

    %% 

    if(length(selected_param) ~= sum(nonzeros(Sel_k)))
        error('Check the number of parameters')
    end

    norm_param = NaN(length(Sel_k),1);

    norm_param(find(Sel_k)) = selected_param;

    norm2origin = [10^(norm_param(1)*log10(bounds.max(1)/bounds.min(1))+log10(bounds.min(1)));  ... %Dsn (log)
                   10^(norm_param(2)*log10(bounds.max(2)/bounds.min(2))+log10(bounds.min(2)));  ... %Dsp (log)
                   norm_param(3)*(bounds.max(3) - bounds.min(3)) + bounds.min(3);  ... %Rsn (Min/Max)
                   norm_param(4)*(bounds.max(4) - bounds.min(4)) + bounds.min(4);  ... %Rsp (Nominal)
                   NaN; ... %eps_s_n
                   NaN; ... %eps_s_p
                   10^(norm_param(7)*log10(bounds.max(7)/bounds.min(7))+log10(bounds.min(7))); %sig_n (log)
                   10^(norm_param(8)*log10(bounds.max(8)/bounds.min(8))+log10(bounds.min(8))); %sig_p (log)
                   norm_param(9)*(bounds.max(9) - bounds.min(9)) + bounds.min(9); %D_e (Min/Max)
                   norm_param(10)*(bounds.max(10) - bounds.min(10)) + bounds.min(10); %eps_e_n (Min/Max)
                   norm_param(11)*(bounds.max(11) - bounds.min(11)) + bounds.min(11); %eps_e_s (Min/Max)
                   norm_param(12)*(bounds.max(12) - bounds.min(12)) + bounds.min(12);%eps_e_p (Min/Max)
                   norm_param(13)*(bounds.max(13) - bounds.min(13)) + bounds.min(13); %kappa (Min/Max)\
                   norm_param(14)*(bounds.max(14) - bounds.min(14)) + bounds.min(14); %t_plus (Min/Max)
                   norm_param(15)*(bounds.max(15) - bounds.min(15)) + bounds.min(15); %t_plus (Min/Max)
                   10^(norm_param(16)*log10(bounds.max(16)/bounds.min(16))+log10(bounds.min(16)));%k_n (log)
                   10^(norm_param(17)*log10(bounds.max(17)/bounds.min(17))+log10(bounds.min(17))); %k_p (log)
                   norm_param(18)*(bounds.max(18) - bounds.min(18)) + bounds.min(18); %Rfn (Min/Max)
                   norm_param(19)*(bounds.max(19) - bounds.min(19)) + bounds.min(19); %Rfp (Min/Max)
                   NaN;
                   norm_param(21)*(bounds.max(21) - bounds.min(21)) + bounds.min(21) %c_e (Min/Max)
                   norm_param(22)*(bounds.max(22) - bounds.min(22)) + bounds.min(22) %E.Dsn (Min/Max) 
                   norm_param(23)*(bounds.max(23) - bounds.min(23)) + bounds.min(23) %E.Dsp (Min/Max) 
                   norm_param(24)*(bounds.max(24) - bounds.min(24)) + bounds.min(24) %E.kn (Min/Max) 
                   norm_param(25)*(bounds.max(25) - bounds.min(25)) + bounds.min(25) %E.kp (Min/Max) 
                   ];


%     origin_param = norm2origin(find(Sel_k));
    origin_param = norm2origin;

end
