function out = alpdot(t,a,w,phi)
%     a = 1;
%     w = 2*pi*0.25;
    out = -a*w*sin(w*t + phi);
%     if t<10
%         out = deg2rad(10)/10;
%     else
%         out = 0;
%     end
end