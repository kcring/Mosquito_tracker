function [mqid mqids instr] =t1_get_mqid(Xh, k, Xi, type)

% find the mosquitoes that are being tracked manually at all times
m_mqid_all=show_tracked_mq(Xh, Xi, 0,  0);

% find the mosquitoes that are being tracked manually at this time-step
m_mqid=show_tracked_mq(Xh, Xi, k, 0);
mqid_max=max(m_mqid_all);

nt_init=min(find(~ismember(1:mqid_max, m_mqid_all)));
if isempty(nt_init)
    nt_init=mqid_max+1;
end

% first time start up
if isempty(mqid_max)
    nt_init=1;
end

switch type
    
    case 0 % new track
        mqid_new=input(sprintf('Assign a new id []=%d, or other=user_defined: ', nt_init));
        if isempty(mqid_new)
            mqid=nt_init;
        else
            mqid=mqid_new;
        end
        mqids=mqid;
        fprintf('Creating new id=%d...\n', mqid);  
        
    case 1 % existing track
        
        if(~isempty(m_mqid))
            mqids=m_mqid;
            fprintf('Listing mosquito ids at current time-step...\n');
            mqid_str=sprintf('%d ...', mqids);
            mqid=input(sprintf('%s\nSelect mosquito id. []=new, #=existing: ', mqid_str));
            if(isempty(mqid))
                mqid_new=input(sprintf('Assign a new id []=%d, or other=user_defined: ', nt_init));
                if isempty(mqid_new)
                    mqid=nt_init;
                else
                    if ~ismember(m_mqid_all, mqid_new)
                        mqid=mqid_new;
                    else
                        fprintf('[!] This id is already in use...Try again\n');
                        [mqid mqids]=t1_get_mqid(Xh,k, Xi, 1);
                    end
                end
                mqids=[mqids; mqid];
                % set instructions
            else
                if ~ismember(mqid, mqids)
                    fprintf('[!] Choose an existing id at current time step ..\n');
                    [mqid mqids]=t1_get_mqid(Xh,k, Xi, 1);
                end
            end
        else
            [mqid mqids]=t1_get_mqid(Xh,k, Xi, 0);
        end


%         switch optTrack.trackone.alg
%             case 'pf'
%                 r=snipx(Xi.cX, mqid);
%                 c=snipy(Xi.cX, k);
%                 if Xh(r(1), c)
%                     p=[ randn(3, optTrack.trackone.N)*optTrack.auto.jitter ; 
%                                 randn(3, optTrack.trackone.N).*(optTrack.auto.sigma_dv*ones(1,optTrack.trackone.N))]...
%                                 + Xh(r(1:6),c)*ones(1,optTrack.trackone.N);
%                 else
%                     p=[ randn(3, optTrack.trackone.N)*optTrack.auto.jitter*10 ; 
%                                 randn(3, optTrack.trackone.N).*(optTrack.auto.sigma_dv*ones(1,optTrack.trackone.N))*5] +...
%                                 [100,100,900, 10, 10, 10]'*ones(1,optTrack.trackone.N);
%                 end
%         end
end
instr='[A] Click on a mosquito in the left frame';
