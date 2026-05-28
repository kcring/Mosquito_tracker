function [X fval] = sim_anneal(cost, X0, options)
%function [X fval] = sim_anneal(cost, X0, options)
%
% simulated annealing search
% t_freeze=options.t_freeze;
% peturb=options.perturb;
% cooling_schedule=options.cooling_schedule;
% E_min=options.E_min;
% T=options.hot;
% max_eval_at_T=options.max_eval_at_T;
% max_rej_at_T=options.max_rej_at_T;
%

% From Wikipedia
% s ← s0; e ← E(s)                                // Initial state, energy.
% sbest ← s; ebest ← e                            // Initial "best" solution
% k ← 0                                           // Energy evaluation count.
% while k < kmax and e > emax                     // While time left & not good enough:
%   snew ← neighbour(s)                           // Pick some neighbour.
%   enew ← E(snew)                                // Compute its energy.
%   if enew < ebest then                          // Is this a new best?
%     sbest ← snew; ebest ← enew                  // Save 'new neighbour' to 'best found'.
%   if P(e, enew, temp(k/kmax)) > random() then   // Should we move to it?
%     s ← snew; e ← enew                          // Yes, change state.
%   k ← k + 1                                     // One more evaluation done
% return sbest  
%
% References
% http://www1bpt.bridgeport.edu/sed/projects/449/Fall_2000/fangmin/chapter2.htm
% http://www.mathworks.com/matlabcentral/fileexchange/10548
% General simulated annealing algorithm
% by Joachim Vandekerckhove

t_freeze=options.t_freeze;
perturb1=options.perturb;
cooling_schedule=options.cooling_schedule;
E_min=options.E_min;
sampling_distr=options.sampling_distr;

T=options.hot;
max_rej_at_T=options.max_rej_at_T;

E0=cost(X0);

E=E0;
terminate=0;
X1=X0;

feval=0;
success=0;
failure=0;        

while(~terminate)
    
    % test 
    Xt=perturb1(X1);
    Et=cost(Xt);
    feval=feval+1;
    
    dE=Et-E;

    % A high value of T0 or options.hot means that you are ready to jump
    % around easily initially but will not do so later
	eval=min(1, sampling_distr(T, dE));
%     eval=min(1, sampling_distr(T, success));
	if rand <= eval
		X1=Xt; E=Et;
		success=success+1;
        failure=0;
        T=cooling_schedule(T);
    else
		failure=failure+1;
    end

    if( Et <= E_min || T <= t_freeze || ...
            failure >= max_rej_at_T)
        terminate=1;
    end

end

X=X1;
fval=Et;

fprintf('Total function evaluations = %d\n', feval);
fprintf('Total successful changes in state = %d\n', success);
fprintf('Rejections at current temperature = %d\n', failure);
fprintf('Initial cost = %.2f\n', E0);
fprintf('Final cost = %.2f\n', Et);
fprintf('Final temperature = %.5f\n', T);
fprintf('------------------------------------\n');


