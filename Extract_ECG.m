function extract_paper_code(x,val,L,k,T)
Fs = 360;
N = length(val);

x = fliplr(x);

ex_wm =[];
a_start = k + ceil((L-2*k)/2);
a_end = N - k - floor((L-2*k)/2);

b1_start = 3;
b1_end = a_start - 1;

b2_start = a_end + 1;
b2_end = N - 2;

X = zeros(L-2*k,2*k);

for i = b1_start:b1_end
    x_tilda = mean([x(i-2) x(i-1) x(i+1) x(i+1)]);
    x_cap = round(x_tilda);
    ei = x(i)-x_cap;
    
    if (x(i)>T-1)&&(x(i)<4095-T)
        if (ei<=-2*T+1) || (ei>=2*T)
            if ei>=2*T
                x(i)=x(i)-T;
            end
            if ei<-2*T+1
                x(i)=x(i)+T-1;
            end
        elseif (ei >-2*T+1) && (ei<2*T)
            b = ei - 2*floor(ei/2);
            ex_wm = [ex_wm b];
            x(i) = (x(i)+x_cap -b)/2;
        end
    elseif (x(i)>=0 && x(i)<=T-1)||(x(i)>=4095-T && x(i)<=4095)
        fb = ex_wm(end);
        ex_wm = ex_wm(1:end-1);
        if fb == 0
            x(i) = x(i);
        else
            if (ei<=-2*T+1) || (ei>=2*T)
                if ei>=2*T
                    x(i)=x(i)-T;
                end
                if ei<-2*T+1
                    x(i)=x(i)+T-1;
                end
            elseif (ei >-2*T+1) && (ei<2*T)
                b = ei - 2*floor(ei/2);
                ex_wm = [ex_wm b];
                x(i) = (x(i)+x_cap -b)/2;
            end
        end
        
        
    end
end


X = zeros(L-2*k,2*k);
for i = a_start:a_end
    xi = [x((i-k):(i-1)) x((i+1):(i+k))];
    x_tilda = mean([x(i-2) x(i-1) x(i+1) x(i+1)]);
    
    xi_orig = x(i);
    
    x(i) = x_tilda;
    
    n = floor((L-2*k)/2);
    m = (i-n):(i+n);
    y = x(m)';
    for j = 1:length(y)
        X(j,:) = [x((m(j)-k):(m(j)-1)) x((m(j)+1):(m(j)+k))];
    end
    v = ((X'*X)^-1)*X'*y;
    x_cap =round(sum(v'.*xi));

    x(i)=xi_orig;
    
    ei = x(i)-x_cap;
    
    if (x(i)>T-1)&&(x(i)<4095-T)
        if (ei<=-2*T+1) || (ei>=2*T)
            if ei>=2*T
                x(i)=x(i)-T;
            end
            if ei<-2*T+1
                x(i)=x(i)+T-1;
            end
        elseif (ei >-2*T+1) && (ei<2*T)
            b = ei - 2*floor(ei/2);
            ex_wm = [ex_wm b];
            x(i) = (x(i)+x_cap -b)/2;
        end
    elseif (x(i)>=0 && x(i)<=T-1)||(x(i)>=4095-T && x(i)<=4095)
        fb = ex_wm(end);
        ex_wm = ex_wm(1:end-1);
        if fb == 0
            x(i) = x(i);
        else
            if (ei<=-2*T+1) || (ei>=2*T)
                if ei>=2*T
                    x(i)=x(i)-T;
                end
                if ei<-2*T+1
                    x(i)=x(i)+T-1;
                end
            elseif (ei >-2*T+1) && (ei<2*T)
                b = ei - 2*floor(ei/2);
                ex_wm = [ex_wm b];
                x(i) = (x(i)+x_cap -b)/2;
            end
        end
        
        
    end
end

for i = b2_start:b2_end
    x_tilda = mean([x(i-2) x(i-1) x(i+1) x(i+1)]);
    x_cap = round(x_tilda);
    ei = x(i)-x_cap;
    
    ei = x(i)-x_cap;
    
    if (x(i)>T-1)&&(x(i)<4095-T)
        if (ei<=-2*T+1) || (ei>=2*T)
            if ei>=2*T
                x(i)=x(i)-T;
            end
            if ei<-2*T+1
                x(i)=x(i)+T-1;
            end
        elseif (ei >-2*T+1) && (ei<2*T)
            b = ei - 2*floor(ei/2);
            ex_wm = [ex_wm b];
            x(i) = (x(i)+x_cap -b)/2;
        end
    elseif (x(i)>=0 && x(i)<=T-1)||(x(i)>=4095-T && x(i)<=4095)
        fb = ex_wm(end);
        ex_wm = ex_wm(1:end-1);
        if fb == 0
            x(i) = x(i);
        else
            if (ei<=-2*T+1) || (ei>=2*T)
                if ei>=2*T
                    x(i)=x(i)-T;
                end
                if ei<-2*T+1
                    x(i)=x(i)+T-1;
                end
            elseif (ei >-2*T+1) && (ei<2*T)
                b = ei - 2*floor(ei/2);
                ex_wm = [ex_wm b];
                x(i) = (x(i)+x_cap -b)/2;
            end
        end 
    end
end

x = fliplr(x);

res = [val; x; val - x];

length(ex_m)