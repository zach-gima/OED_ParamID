function [ card ] = Cardinality( A )
%Finds the cardinality of the admissable experiment sets using the A matrix
card = sum((A~=0)');

end

