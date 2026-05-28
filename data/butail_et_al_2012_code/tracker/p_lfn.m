function wts=p_lfn(Zk, p_, cam, Xi, optTrack)

% the standard deviation in end points is based on location of an end-point
% pixel of a fixed grain of rice on a pendulum in a plane parallel to the camera plane
% that is filmed at 25 frames per second. Since we know the length of the
% grain of rice to be fixed, we assume that any variation in length is due
% to noise in end-point computation. Therefore assuming Gaussian noise,
% noise(length) = noise(e1) + noise(e2). Assuming both are the same,
% std(length) = sqrt(std(ep))
%
% The noise is center position is a function of the streak length. Since
% the bounding ellipse also represents a normal distribution around the
% center of the streak, we can equate the major and minor axis length as
% diag([1/a^2, 1/b^2]) = Sigma^-1. Which implies sigma_x=a, sigma_y=b; BUT
% since the ellipse is not necessarily in the image plane frame we perform
% a transformation [cos(t) -si

wts=p_pos(Zk.u, p_(Xi.ri,:), cam, diag(Zk.sigma)).* ...
    p_mq_velocity(Zk, p_, cam, Xi, optTrack).*...
    pdf('unif', sqrt(sum(p_(Xi.rdi,:).^2)), 100, 4000);

%     
%               diag(abs(imstream(cc).Zk(assoc_t(t_ii,cc)).v)/2)).*...
%     pdf('unif', p_(r(Xi.ri(3)),:), optTrack.auto.swarm_boundaries(1), optTrack.auto.swarm_boundaries(2));