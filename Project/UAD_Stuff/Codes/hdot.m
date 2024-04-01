function out = hdot(t,a,w,phi)
%     a = 1;
%     w = 2*pi*0.25;
    out = -a*w*sin(w*t + phi);
%     out = 0;
end