function X=optimal_triangulation(z1,z2,F, P1, P2)
%function X=optimal_triangulation(z1,z2,F)
%
%
% ref: Algorithm 12.1 Hartley and Zisserman - Multiple view geometry
%
% jump to the last step if z1'*F*z2=0.00

[z1h z2h]=find_optimal_correspondence(z1,z2,F);

% fprintf('z1\n')
% z1
% z1h
% 
% fprintf('z2\n')
% z2
% z2h

A=[ z1(1)*P1(:,3)'-P1(:,1)'; 
    z1(2)*P1(:,3)'-P1(:,2)';
    z2(1)*P2(:,3)'-P2(:,1)'; 
    z2(2)*P2(:,3)'-P2(:,2)'];

