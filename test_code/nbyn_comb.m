function [li_aa] = nbyn_comb(a,n)
% A is an n by n square matrix 
% - represents an adjacency matrix
% - (all possible ways of combining 1-n elements)

% inputs:
% a contains a subset of the elements of a
% n is the number of elements

% output: 
% li_aa, linear indices of all possible combinations of subset elements

numele=length(a);
% make sure a is a row vector
if(size(a,1)==1)
    a=a';
end

% all possible combinations of subset of elements
arows=a*ones(1,numele);
arows=arows(:);
acols=a*ones(1,numele);
acols=acols';
acols=acols(:);

% linear indices of combinations
li_aa=n*(acols-1)+arows;