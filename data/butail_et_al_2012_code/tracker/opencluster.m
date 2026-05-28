function [Zk, scan]=opencluster(CL, Z)
Zk=Z.u(:,CL.z_id)';
scan=CL.scan;
