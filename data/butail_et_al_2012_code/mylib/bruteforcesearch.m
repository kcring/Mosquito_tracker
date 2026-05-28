function [X fval] = bruteforcesearch(obj_fun, X0, res0, resd, nelm, options)
%function [X fval] = bruteforcesearch(obj_fun, X0, res0, resd, nelm, options)
%
% brute force search
% obj_fun is the function handle which needs to be evaluated
% X0 is the initial start
% res0 is initial resolution vector of size [nx x 1]
% resd is the desired resolution vector this is also [nx x 1] 
% Think of resd=res0./[some vector];
%
% nelm is the number of elements per dimension [nx x 1]
%
% options.MINMAX==1 implies minimization, 2 is for maximization 
% options.display==1 will show the status as the algorithm proceeds
%
%
%%%%%%%%%%%%%%%% Example (nx=11) %%%%%%%
% rv(:,1)=[1.25 1.25 .5 .25 .25 .25 .02 10 .0625 .0625 .0039]';
% % desired resolution vector
% drv=rv(:,1)./[4 4 4 8 8 8 64 1 2 2 2]';
% % number of elements per dimension in the begining
% % nelm=[2 2 2 5 5 5 6 1 1 1 1]' ;
% nelm=[2 2 2 3 3 3 8 1 1 1 1]' ;

iter=0;
resd_true=0;
Xh_max_wt=X0';
val=0;

while(~resd_true)
    [p_gen N resd_true]=get_particles(res0(:,1)/2^(iter-1), Xh_max_wt, resd, nelm, options);
    if(~resd_true)
        if options.Display
            disp(sprintf('iter %d, N %d, wt(%d) %.1f', iter, N, iter-1, val));
        end
        iter=iter+1;
        if options.MINMAX ==1
            wts=ones(1,N)*10000;
        elseif options.MINMAX==2
            wts=zeros(1,N);
        end
        for ii=1:N
            wts(ii)=obj_fun(p_gen(:,ii)');
        end
        if options.MINMAX==2
            [val idx]=max(wts);
        elseif options.MINMAX==1
            [val idx]=min(wts);
        end
        
        Xh_max_wt=p_gen(:,idx);
    end
    fval=val;
    X=Xh_max_wt';
end


function [p_gen N resd_true]=get_particles(res0, Xh_max_wt, resd, nelm, options)


vl=size(Xh_max_wt,1);
show_rvals=zeros(max(nelm),vl);

if(res0<=resd)
    resd_true=1;
    p_gen=Xh_max_wt;
    N=1;
else

    for r_ii=1:size(res0,1)
        if(res0(r_ii)>resd(r_ii))
            rvals(r_ii).v=linspace(Xh_max_wt(r_ii)-res0(r_ii)/2, Xh_max_wt(r_ii)+res0(r_ii)/2, nelm(r_ii));
        else
            rvals(r_ii).v=Xh_max_wt(r_ii);
        end
        
        show_rvals(1:size(rvals(r_ii).v,2), r_ii)=rvals(r_ii).v';
    end

    switch vl
        case 2
            [Rvals(1).v Rvals(2).v]=ndgrid(rvals(1).v, rvals(2).v);
        case 3
            [Rvals(1).v Rvals(2).v Rvals(3).v]=ndgrid(rvals(1).v, ...
                                                rvals(2).v, rvals(3).v);
        case 4
            [Rvals(1).v Rvals(2).v Rvals(3).v Rvals(4).v] = ndgrid( ...
                                                rvals(1).v, rvals(2).v, rvals(3).v, ...
                                                rvals(4).v);            
        case 5
            [Rvals(1).v Rvals(2).v Rvals(3).v Rvals(4).v Rvals(5).v] = ...
                                            ndgrid( rvals(1).v, rvals(2).v,...
                                            rvals(3).v, rvals(4).v, rvals(5).v);            
        case 6
            [Rvals(1).v Rvals(2).v Rvals(3).v Rvals(4).v Rvals(5).v ...
                Rvals(6).v ] =  ndgrid( rvals(1).v, rvals(2).v,  rvals(3).v, ...
                                        rvals(4).v, rvals(5).v,  rvals(6).v);      
        case 7
            [Rvals(1).v Rvals(2).v Rvals(3).v Rvals(4).v Rvals(5).v ...
                Rvals(6).v Rvals(7).v] = ndgrid( rvals(1).v, rvals(2).v, ...
                                            rvals(3).v, rvals(4).v, ...
                                            rvals(5).v,  rvals(6).v, ...
                                            rvals(7).v);
        case 8
        case 9
        case 10
        case 11
            [Rvals(1).v Rvals(2).v Rvals(3).v Rvals(4).v Rvals(5).v ...
                Rvals(6).v Rvals(7).v Rvals(8).v Rvals(9).v Rvals(10).v ...
                Rvals(11).v]    = ...
                                   ndgrid( rvals(1).v, rvals(2).v,  rvals(3).v, ...
                                        rvals(4).v, rvals(5).v,  rvals(6).v, ...
                                        rvals(7).v, rvals(9).v,  rvals(10).v, ...
                                        rvals(11).v); 
            
        case 12
        otherwise
            disp('[!] Vector length should be between 2-12')
    end
    if(options.Display)
        disp(show_rvals);
    end

    N=numel(Rvals(1).v);
    p_gen=zeros(N,vl);
    for ii=1:vl
        Rvals(ii).v=reshape(Rvals(ii).v,N,1);
        p_gen(:,ii)=Rvals(ii).v;
    end
	p_gen=p_gen';
    
    resd_true=0;
end

