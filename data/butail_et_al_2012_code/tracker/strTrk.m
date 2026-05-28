function Trk_arr=strTrk(siz)
% function Trk_arr=strTrk(siz)
% if siz is a number then Trk_arr is a [1 x siz ]vector of arrays 

Xi=strXi;

Trk=struct('state', zeros(Xi.nX,1), ...
           'k', 0, ...
           'prob', 0);

if numel(siz) ==1
    Trk_arr(1:siz)=Trk;
else
    Trk_arr(1:siz(1),1:siz(2))=Trk;
end
