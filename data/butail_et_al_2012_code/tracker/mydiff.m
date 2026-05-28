function vdot=mydiff(v,shift)

vdot=circshift(v,[0,-shift])-v;
vdot(:,end)=vdot(:,end-1);