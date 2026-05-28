function [t br] =id2tbr(id, factor)
% function [t b]=id2loc(id, factor)
%
% function to give location if id is given
% factor is a power of 10 


br=round((id/factor-floor(id/factor))*factor);

t=floor(id/factor);

