function [norm_paramS] = originS_to_normS(selected_param, bounds, Sel_k)

% By Saehong Park (2017.10.29)
% Normalize parameter sensitivity values for simulation
% Note: This function calculates the normalization factor that will be
% multipled against the calculated sensitivity values

%%

if(length(selected_param) ~= sum(nonzeros(Sel_k)))
    error('Check the number of parameters')
end

origin_param = NaN(length(Sel_k),1);
origin_param(find(Sel_k)) = selected_param;

emptyblk = NaN;

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
                bounds.max(21) - bounds.min(21)]; %(21)                  
                     
norm_paramS = origin2norm(find(Sel_k));


