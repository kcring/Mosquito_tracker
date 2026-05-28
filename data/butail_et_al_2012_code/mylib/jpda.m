function Bjt = jpda(m, zh_, S, gamma, nt, Pd, V)
%function Bjt = jpda(m, zh_, S, gamma, nt, Pd, V) 
%     This function computes association probabilities
%     for a Joint Probabilistic Data Association Filter (Bar-Shalom 1987)
%     Inputs
%         m       : measurement vector
%         zh_     : estimated measurements h(x_)
%         S       : measurement covariance
%         gamma   : gate value
%         nt      : number of targets
%         Pd      : Probability of detection
%         V       : Volume of region
% 
%     Outputs
%         Bjt     : A matrix with its elements as probabilities of jth measurement having come from 
%                   t-th target
%
% Jan 2010, Sachit Butail

[mvec nm]=size(m);

vm=zeros(nm, nt+1);

% measurement can always be from clutter 
vm(:,1)=1;


for i=1:nm % measurements
    for j=1:nt % estimated measurements corresponding to targets
        inov(:,j)=m(:,i)-zh_(:,j);
        vmtest=inov(:,j)'*inv(S(:,:,j))*inov(:,j);
        if(vmtest<=gamma)
            vm(i,j+1)=1;
        end
    end
end

% extract feasible events from VM
Omega = vm;
[rows cols]=size(Omega);

% pick all possible rows for a row where there is only one 1 in each row
for j=1:rows
    onesum=sum(Omega(j,:));
    onevec=find(Omega(j,:)==1);
    eval(['rm' sprintf('%.2d', j) ' = zeros(onesum,cols);']);
    for i=1:onesum
        eval(['rm' sprintf('%.2d',j) '(i,onevec(i))=1;'])
    end
end

% make a vector with matrix num as the first digit and row num as the
% second digit
jj=1;
for j=1:rows
    onesum=sum(Omega(j,:));
    for i=1:onesum
        vv(jj)=str2double(sprintf('%.2d%.2d', j,i));
        jj=jj+1;
    end
end

% how many combinations are possible. note that there will be many but we
% don't want repititions of rows from the same matrix which we throw away
% next
try
fms=nchoosek(vv,rows); % feasible matrices
catch
    keyboard;
end
[fmsnum jj]=size(fms);
jj=1;
for ii=1:fmsnum
    vvcheck=fms(ii,:);
    for kk=1:rows
        fc1(kk)=floor(vvcheck(kk)/(10^2));
    end
    if(prod(fc1)==factorial(rows))
        fmsgood(jj,:)=fms(ii,:);
        jj=jj+1;
    end
end
      
% assign the values to feasible events theta
[fmnum jj]=size(fmsgood);

for ii=1:fmnum
    for jj=1:rows
        id=sprintf('%.4d',fmsgood(ii,jj));
        rmid=id(1:2);
        rowid=fmsgood(ii,jj)-floor(fmsgood(ii,jj)/(10^2))*(10^2);
        eval(['thetat(jj,:,ii)=rm' rmid '(rowid,:);']);
    end
end

% pick valid thetas, i.e sum of cols is 1 for all except first col. NOTE that sum of rows is already 1 due to
% our choice of matrices
ii=1;
for jj=1:fmnum
    colsum=sum(thetat(:,:,jj),1);
    if(find(colsum(2:end)>1)) % column condition for validity -- One target cannot produce more than one measurement
        %
    else
        theta(:,:,ii)=thetat(:,:,jj);
        ii=ii+1;
    end
end

% for each theta found we can now compute the binary indicators

[jj kk thetanum]=size(theta);

Bjt=zeros(rows,cols);

% find the probability of each feasible event
for jj=1:thetanum
    for ii=1:rows % measurement association indicator for each theta
        tau(ii)=sum(theta(ii,2:end,jj));
    end
    
    % number of unassociated measurements
    phi(jj)=rows-sum(tau);
    
    % target detection indicator (ignore t=1 which is clutter only)
    tj=zeros(cols-1,1); % tj associates the measurement to target
    
    Pdprod=1;
    for kk=2:cols
        delta(kk-1)=sum(theta(:,kk,jj)); % kk-1 because the 1st col is clutter
        
        % which measurement is associated with this target
        if(delta(kk-1)), tj(kk-1)=find(theta(:,kk,jj)==1); end
        
        Pdprod=Pdprod*Pd^delta(kk-1)*(1-Pd)^(1-delta(kk-1));
    end
    
    % equation 9-45 in Bar-shalom
    Ntj=ones(rows,1);
    for mm=1:rows
        if(tau(mm)) % if there is a target associated with this measurement use "find"
            ati=find(tj==mm); % associated target index
            Ntj(mm)=(1/(2*pi)^(mvec/2)/abs(det(S(:,:,ati)))*...
                exp(-1/2*(m(:,mm)-zh_(:,ati))'*inv(S(:,:,ati))*(m(:,mm)-zh_(:,ati))))^tau(mm);
        end
    end
    
    PthetaZ(jj)=factorial(phi(jj))/V^phi(jj)*prod(Ntj)*Pdprod;
end

% normalize. Probability of each feasible event
PthetaZ=PthetaZ/sum(PthetaZ);

% equation 9-46 in Bar-Shalom
for jj=1:rows % for each measurement
    for tt=1:cols % for every target including clutter
        for ii=1:thetanum
            % the probability that a measurement came from a target is
            % summed over all the feasible association events 
            Bjt(jj,tt)=Bjt(jj,tt)+PthetaZ(ii)*theta(jj,tt,ii);
        end
    end
end