function [out] = twenty2to25(in,flag)
% Created By Dylan Kato 4/26/19
% Transforms between matrix indicies with and without equilibrium
% parameters
noneqidx = [1:4,7:19,21:25];
if flag == '225'
    out = noneqidx(in);
else
    out = find(noneqidx == in);
end
end

