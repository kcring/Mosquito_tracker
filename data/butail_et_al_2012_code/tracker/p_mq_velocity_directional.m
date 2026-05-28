function wts = p_mq_velocity_directional(luv, tp, sigma_v, cam, te, Xi)

vp=w2cam_rdot(tp, cam, te, Xi);

wts= normpdf(vp(1,:), luv(1), sigma_v(1)).*normpdf(vp(2,:), luv(2), sigma_v(2));