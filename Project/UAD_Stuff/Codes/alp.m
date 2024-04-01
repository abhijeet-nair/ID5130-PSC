function out = alp(t,a,w,phi)
%     a = 1;
%     w = 2*pi*0.25;
    out = a*cos(w*t + phi);
% %     if t<2 && t>4
% %         out = 0;
% %     else
% %         out = deg2rad(10);
% %     end
%     if t<10
%         out = deg2rad(10)*(t/10);
%     else
%         out = deg2rad(10);
%     end
end