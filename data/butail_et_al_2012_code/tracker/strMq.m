function Mq_arr=strMq(siz)
%function fish_arr=strMq(siz)
% if siz is a number then Mq_arr is a [1 x siz ]vector of arrays 

Mq=struct('position', [100 100 900]', ...
          'velocity', 10*ones(3,1), ...
	  'acceleration', zeros(3,1), ...
	  'flags', zeros(2,1), ...
	  'id', 0);
        
if numel(siz) == 1
    Mq_arr(1:siz)=Mq;
else
    Mq_arr(1:siz(1), 1:siz(2))=Mq;
end
