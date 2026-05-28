function estx = postest(x, wts, flag)

%{
This function will return an estimate of a distribution, based on what flag is passed.
So, for e.g. if flag passed is 1 and x is a matrix, the function will compute the mean along
the longer dimension
if flag is 2 then the function will compute the max along the longer dimension... and so on.

ok. so there's a check now that the longer dimension has to be columns!!
%}

switch flag
    case 1
        estx = mean(x,2);
    case 2
        [val idx] = max(wts);
        estx = x(:,idx);
    case 3
        estx = mode(x,2);
    otherwise
        error('brrrr!!');
end
    