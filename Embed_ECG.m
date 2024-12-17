clc
clear
close all
L = 15;
k = 3;
load('100m.mat')
Fs = 360;
N = length(val);
T = 4;
wm = double(rand(1,3600)>0.5);
S =[];
ind = 1;
val = val + 2048;
x = val;

a_start = k + ceil((L-2*k)/2);
a_end = N - k - floor((L-2*k)/2);

b1_start = 3;
b1_end = a_start - 1;

b2_start = a_end + 1;
b2_end = N - 2;

for i = b1_start:b1_end
    x_tilda = mean([x(i-2) x(i-1) x(i+1) x(i+1)]);
    x_cap = round(x_tilda);
    
    ei = x(i)-x_cap;
    
    if abs(ei)<T
        ei_mod = 2*ei+wm(1);
        S = [S wm(1)];
        wm =wm(2:end);
        ind = ind+1;
    elseif ei>=T
        ei_mod = ei+T;
    elseif ei<=-T
        ei_mod = ei-T+1;
    end
    
    xi_mod = x_cap + ei_mod;
    
    if xi_mod>=0 && xi_mod<=4095
        x(i) = xi_mod;
        ind = ind -1;
    end
    
    if (xi_mod>=0 && xi_mod<=T-1)||(xi_mod>=4095-T && xi_mod<=4095)
        wm =[1 wm];
        ind = ind -1;
    end
    
    if xi_mod<0||xi_mod>4095
        x(i)=x(i);
        wm = [0 wm];
    end
    
end


X = zeros(L-2*k,2*k);
for i = a_start:a_end
    xi = [x((i-k):(i-1)) x((i+1):(i+k))] ;
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
    
    ei = xi_orig-x_cap;
    
    if abs(ei)<T
        ei_mod = 2*ei+wm(1);
        S = [S wm(1)];
        wm =wm(2:end);
        ind = ind + 1;
    elseif ei>=T
        ei_mod = ei+T;
    elseif ei<=-T
        ei_mod = ei-T+1;
    end
    
    xi_mod = x_cap + ei_mod;
    
    if xi_mod>=0 && xi_mod<=4095
        x(i) = xi_mod;
    end
    
    if (xi_mod>=0 && xi_mod<=T-1)||(xi_mod>=4095-T && xi_mod<=4095)
        wm =[1 wm];
        ind = ind -1;
    end
    
    if xi_mod<0||xi_mod>4095
        x(i)=xi_orig;
        wm = [0 wm];
        ind = ind -1;
    end
end

for i = b2_start:b2_end
    x_tilda = mean([x(i-2) x(i-1) x(i+1) x(i+1)]);
    x_cap = round(x_tilda);
    ei = x(i)-x_cap;
    
    if abs(ei)<T
        ei_mod = 2*ei+wm(1);
        S = [S wm(1)];
        wm =wm(2:end);
        ind = ind +1;
    elseif ei>=T
        ei_mod = ei+T;
    elseif ei<=-T
        ei_mod = ei-T+1;
    end
    
    xi_mod = x_cap + ei_mod;
    
    if xi_mod>=0 && xi_mod<=4095
        x(i) = xi_mod;
    end
    
    if (xi_mod>=0 && xi_mod<=T-1)||(xi_mod>=4095-T && xi_mod<=4095)
        wm =[1 wm];
        ind = ind -1;
    end
    
    if xi_mod<0||xi_mod>4095
        x(i)=x(i);
        wm = [0 wm];
        ind = ind -1;
    end
end

% val = val - 2048;
% x = x - 2048;

new_ecg_data = x;
PRD = 100 *((sum((val - new_ecg_data).^2)/sum(val.^2)).^0.5);

res = [val;x];
extract_paper_code(x,val,L,k,T);