function Xh = terminate_track(Xh, id, k, Xi, dir)

if dir % leading edge
    Xh(snipx(Xi.cX, id), k+1:end)=0;
else % trailing edge
    Xh(snipx(Xi.cX, id), 1:k-1)=0;
end
