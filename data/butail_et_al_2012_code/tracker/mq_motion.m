function xp1=mq_motion(xp, mm, dt)
Np=size(xp,2);
nx=size(xp,1);

xp1=xp;


switch mm
    case 'wienera'
        dW=randn(3,Np)*20000;
        % wiener process noise acceleration (integral of white noise)
        xp1(1:3,:)=xp(1:3,:) + xp(4:6,:)*dt + xp(7:9,:)*dt^2/2 +dW*dt^2/2;
        xp1(4:6,:)=xp(4:6,:) + xp(7:9,:)*dt + dW*dt;
        xp1(7:9,:)=xp(7:9,:) + dW*dt;

    case 'wna'
        dW=randn(3,Np)*10000;
        xp1(1:3,:)=xp(1:3,:) + xp(4:6,:)*dt + xp(7:9,:)*dt^2/2 +dW*dt^2/2;
        xp1(4:6,:)=xp(4:6,:) + xp(7:9,:)*dt+ dW*dt;
        xp1(7:9,:)= dW;       
        
    case 'cv'
        
        dW=randn(3,Np)*7000;
        xp1(1:3,:)=xp(1:3,:) + xp(4:6,:)*dt +dW*dt^2/2;
        xp1(4:6,:)=xp(4:6,:) + dW*dt;
%         xp1(7:8,:)=0;

    case 'sp1'
        
        % self-propelled particle in world frame
        xp1(1:3,:)=xp(1:3,:)+(ones(3,1)*xp(9,:)).*xp(4:6,:)*dt;
        xp1(4:6,:)=normv(cross(xp(4:6,:), [xp(9,:)*0; xp(7:8,:)]));
        xp1(7:8,:)=xp(7:8,:)+randn(2,Np)*.1;
        xp1(9,:)=xp(9,:) + randn(1,Np)*100;
        xp1(10:nx,:)=xp(10:nx,:);
        
    case 'sp2'
        
        % self-propelled particle in body frame
        se3m1(1:3,1,:)=xp(4:6,:); se3m1(1:3,2,:)=xp(7:9,:); se3m1(1:3,3,:)=xp(10:12,:);
        se3m1(1:3,4,:)=xp(1:3,:); 
        se3m1(4,4,:) = ones(Np,1);
        s=xp(13,:)+randn*100;

        q=xp(14,:)+randn(1,Np)*1;
        h=xp(15,:)+randn(1,Np)*1;
        w=randn(1,Np)*.1;
        
        for c = 1:Np
            xi=[[0 -q(c) -h(c); q(c) 0 -w(c); h(c) w(c) 0], [s(c) 0 0]'; 0 0 0 0];
            expxidt = expm(xi*sqrt(dt));
            se3(:,:,c) = se3m1(:,:,c)*expxidt;

            xp1(1:3,c) = se3(1:3,4,c);
            xp1(4:6,c) = se3(1:3,1,c);
            xp1(7:9,c) = se3(1:3,2,c);
            xp1(10:12,c) = se3(1:3,3,c);
            xp1(13,c)=s(c);
            xp1(14:15,c)=[q(c); h(c)];
        end
        xp1(16:nx,:)=xp(16:nx,:);
        
end
