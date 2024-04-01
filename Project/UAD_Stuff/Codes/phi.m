function out = phi(t,u,b,p1,p2,e1,e2)
%     psi1 = 0.165;
%     psi2 = 0.335;
%     eps1 = 0.0455;
%     eps2 = 0.3;
    out = 1 - p1*exp(-e1*u*t/b) - p2*exp(-e2*u*t/b);
end