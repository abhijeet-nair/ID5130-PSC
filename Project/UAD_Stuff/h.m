function out = h(t,f1,f2,a)
    T1 = 1/f1;
    T2 = 1/f2;
    T = 0.5*(T1 + T2);
    nt = length(t);
    out = zeros(size(t));
    for i = 1:nt
        t1 = rem(t(i),T);
%         fprintf('t1 = %.4f, T1 = %.4f, T2 = %.4f\n',t1,T1,T2)
    
        if t1 < 0.5*T1
            out(i) = a*cos(2*pi*f1*t1);
%             fprintf('\tT1\n')
        else
            out(i) = a*cos(2*pi*f2*(t1 - 0.5*T1) + pi);
%             fprintf('\tT2\n')
        end
    end
end