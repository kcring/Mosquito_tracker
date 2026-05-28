function C_arr=strClstr(siz)

C=struct('z_id', [], ...
         'scan', [], ...
         'kgate', []);
if numel(siz)==1     
	C_arr(1:siz)=C;
else
    C_arr(1:siz(1),1:siz(2))=C;
end