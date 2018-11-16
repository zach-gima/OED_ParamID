function [ removed ] = findclusterparams(STSnorm,Sens_mag, threshold,Np,num_exp)
%This function takes the values in STSnorm and a given threshold and it
%will find what parameters must be removed
minSTSnorm=elementwisecellmin(STSnorm,Np,num_exp)-diag(ones(Np,1));
avgSTS = mean(Sens_mag');
removed=[];
while true

% find the maximum minimum element excluding those removed
maxS = maxexc(minSTSnorm,removed);
% check if this element is above the threshold
% if not, end
if maxS<threshold
    break
end
% if so, procede

% find the indexes of this element
[row,col] = find(minSTSnorm==maxS,1);
% compare their sensitivities
senscomp = [avgSTS(row),avgSTS(col)];
% save which is higher and which is lower
if senscomp(1) >= senscomp(2)
    rem=col;
end
if senscomp(1) < senscomp(2)
    rem=row;
end
% append the lower to a list of removed elements
removed = [removed,rem];
% repeat for appended matrix
end
removed = sort(removed);
end

