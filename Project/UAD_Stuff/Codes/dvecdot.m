function ddvec = dvecdot(t,dvec,h,alp,eps1,eps2,u,b)
    ddvec = zeros(size(dvec));
    ddvec(1) = h(t) - (eps1*u/b)*dvec(1);
    ddvec(2) = h(t) - (eps2*u/b)*dvec(2);
    ddvec(3) = alp(t) - (eps1*u/b)*dvec(3);
    ddvec(4) = alp(t) - (eps2*u/b)*dvec(4);
end