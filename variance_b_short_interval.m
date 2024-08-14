function compute_empirical_variance_short_interval()
    clear;
    t = 10^4;
    %generate an array of b(i) for i=t^2+1, ..., 3*t^2.
    s = ceil(t*3^0.5);
    a(1,1:(s+1)) = (0:s).^2;
    temp = repmat(a,(s+1),1);
    clear a;
    sumOfSquares = temp + temp';
    clear temp;
    u = sort(unique(sumOfSquares));
    clear sumOfSquares;
    b1 = find(u>3*t^2,1);
    b2 = find(u>t^2,1);
    u = u(b2:b1-1);
    b_array = zeros(2*t^2,1);
    for i=1:length(u)
        b_array(u(i)-t^2) = 1;
    end
    clear u;
    
    %goal: compute \bar{B}(i)-\bar{B}(i-1) for i=t^2+1, ..., 3*t^2, and subtract it from b_array
    
    %constants related to contribution of singularity at s=1:
    f_val1 = ffunc_at_one(50,50);
    fun_top = @(x) 1./(x.*sqrt(1-x));
    q_top = integral(fun_top,0.75,1,'RelTol',0,'AbsTol',1e-12);
    
    m = 1000; %should be even
    %compute f(s) at (m+1) points s between 3/4 and 1. 
    %(excluding last point, 1.)
    pts = 0:1:m;
    s_pts_top = (1/(4*m))*pts  + (3/4);
    s_pts_top(end)=[];
    f_vals_top = zeros(1,m);
    for i=1:m
        f_vals_top(i) = ffunc(s_pts_top(i),50,50);
    end
    clear pts;
    
    nums_m = 1:m;
    weights_top = (horzcat([1], 2*(mod(nums_m(2:end)+1, 2)+1)))'; % 1,4,2,4,2,... - Simpson weights
    K = 2000; %K - controls memory usage. Larger K = less memory.
    L = ceil(2*t^2 / K);
    f_val1_vec = f_val1 +zeros(1,m);
    s_denom = sqrt(1-s_pts_top).*s_pts_top;
    s_pts_upper_rep = repelem(s_denom,[L],[1]);
    for i=1:K
        min_idx = L*(i-1)+1;
        max_idx = min(L*i,2*t^2);
        xs = ((min_idx + t^2-1):(max_idx+t^2-1))';
        if i==K
            s_pts_upper_rep = repelem(s_denom,[max_idx-min_idx+1],[1]);
        end
        values_top = ((f_vals_top.*(xs+1).^s_pts_top - f_val1_vec.*(xs+1)) -(f_vals_top.*(xs).^s_pts_top - f_val1_vec.*xs))./s_pts_upper_rep;
        
        b_array(min_idx:max_idx) = b_array(min_idx:max_idx) - (mtimes(values_top, weights_top)/(12*pi*m));
    end
    clear weights_top; clear f_vals_top; clear s_pts_top; clear values_top; clear f_val1_vec; clear s_denom;
    
    %subtract contribution of singular component at [3/4,1]
    b_array = b_array - f_val1*q_top/pi;
    
    %compute f(s) at m+1 points s between 1/2 and 3/4 (excluding first one, 1/2).
    pts = 0:1:m;
    s_pts_bottom = (1/(4*m))*pts  + (1/2);
    s_pts_bottom = s_pts_bottom(2:end);
    f_vals_bottom = zeros(1,m);
    for i=1:m
        f_vals_bottom(i) = ffunc(s_pts_bottom(i),50,50);
    end
    clear pts;
    
    %constants related to contribution of singularity at s=1/2:
    f_val2 = ffunc_at_half(50,50);
    q_bottom = 2^(7/4)/3; %integral of sqrt(8)(2s-1)^(-1/4) from 1/2 to 3/4.
    
    nums_m = 2:(m+1);
    weights_bottom = (horzcat( 2*(mod(nums_m(1:end-1)+1, 2)+1), [1]))'; % 4,2,4,2,...,1 - Simpson weights
    clear nums_m; 
    
    f_val_half_vec = f_val2 + zeros(1,m);
    s_denom = sqrt(1-s_pts_bottom).*s_pts_bottom;
    
    s_denom_rep = repelem(s_denom,[L],[1]);
    for i=1:K
        min_idx = L*(i-1)+1;
        max_idx = min(L*i,2*t^2);
        xs = ((min_idx + t^2-1):(max_idx+t^2-1))';
        if i==K
            s_denom_rep = repelem(s_denom,[max_idx-min_idx+1],[1]);
        end
        values_bottom = (f_vals_bottom.*(xs).^s_pts_bottom)./s_denom_rep - f_val_half_vec.*(2*s_pts_bottom-1).^(-1/4).*sqrt(8*xs);
        values_bottom = (f_vals_bottom.*(xs+1).^s_pts_bottom)./s_denom_rep -  f_val_half_vec.*(2*s_pts_bottom-1).^(-1/4).*sqrt(8*(xs+1)) -  values_bottom;       
        b_array(min_idx:max_idx) = b_array(min_idx:max_idx) - (mtimes(values_bottom, weights_bottom)/(12*pi*m));
    end
    clear s_denom_rep; clear weights_bottom; clear f_vals_bottom; 
    clear s_pts_bottom; clear values_bottom; clear f_val_half_vec; 
    clear xs; clear s_denom;
 
    %subtract contribution of singular component at [1/2,3/4]
    nums = (1:2*t^2)';
    b_array = b_array - f_val2*q_bottom./(pi*(sqrt(nums+t^2)+sqrt(nums+t^2-1))); %contribution of singular component at [1/2,3/4]
    clear nums; 
    
    %generate values H, from 1 to X, with (roughly) constant spacing on a logarithmic scale
    hs_number = floor(log(t^2)/log(1.1));
    hs = zeros(hs_number,1);
    for i=1:hs_number
        hs(i) = floor(1.1^i);
    end
    hs = sort(unique(vertcat((1:10)', hs))); %make sure to include smallest H: H=1,2,...,10
    
    %compute variance in short intervals
    hcount = zeros(length(hs),1);
    for k = 1:length(hs)
        window_sums = zeros(t^2,1);
        window_sums(1) = sum(b_array(1:hs(k)));
        for i = 2:(t^2)
            window_sums(i) = window_sums(i-1) + b_array(i+hs(k)-1)-b_array(i-1);
        end
        normh = hs(k)/sqrt(log(t^2));
        hcount(k) = mean(window_sums.^2,1)/normh;
    end
    x_axis = (log(hs)/log(t^2));
    save('results_short_interval.mat', 'x_axis', 'hcount');
end

function chi_val = chi(x)
    %evaluates the non-trivial character modulo 4 at x
    if mod(x,2) == 0
        chi_val = 0;
    else
        chi_val = (-1)^((x-1)/2);
    end
end

function l_val = lfunc(s,N)
    %approximates L(s,chi_4) for s>0, by taking an Euler transform of the
    % series sum(k=1,...) chi(k)/k^s.
    % parameter N controls precision.
    l_val = 0;
    for i = 0:(N-1)
        for k = 0:i
            l_val = l_val + nchoosek(i,k)* (-1)^k * chi(k+1) * (k+1)^(-s) / 2^(i+1);
        end
    end
end

function z_val = zfunc(s,N)
    %approximates zeta(s) for s>0, by taking an Euler transform of 
    % the eta function (1-2^(1-s))*zeta(s) = sum(k=1,...) (-1)^(k-1)/k^s,
    % and then dividing by (1-2^(1-s)).
    %parameter N controls precision.
    z_val = 0;
    for i = 0:(N-1)
        for k = 0:i
            z_val = z_val + nchoosek(i,k)* (-1)^k * (k+1)^(-s) / 2^(i+1);
        end
    end
    z_val = z_val / (1-2^(1-s));
end

function f_val = ffunc_at_one(N,n)
    %approximates f(1), where f appears in the integrand of \bar{B}(x)
    %N, n - control precision
    v = 1;
    for k=1:n
        v = v * (zfunc(2^k,N)*(1-2^(-2^k))/lfunc(2^k,N))^(1/2^(k+1));
    end
    f_val = (pi/2)^(1/2) * v;
end

function f_val = ffunc_at_half(N,n)
    %approximates lim_{s->1/2}f(s)(2s-1)^(1/4), where f appears in the integrand of \bar{B}(x)
    %N, n - control precision
    v1 = (lfunc(1/2,N)/(1-2^(-1/2)))^(1/2);
    v2 = (-zfunc(1/2,N)/2)^(1/2);
    v3 = (2/pi)^(1/4);
    for k=2:n
        v3 = v3 * (zfunc(2^(k-1),N)*(1-2^(-2^(k-1)))/lfunc(2^(k-1),N))^(1/2^(k+1));
    end
    f_val = v1*v2*v3;
end

function f_val = ffunc(s,N,n)
    %approximates the function f(s), appearing in the integrand of \bar{B}(x)
    %N, n - control precision
    v1 = (lfunc(s,N)/(1-2^(-s)))^(1/2);
    v2 = (zfunc(s,N)*(s-1))^(1/2);
    v3 = 1;
    for k=1:n
        v3 = v3 * (zfunc(2^k*s,N)*(1-2^(-2^k*s))/lfunc(2^k*s,N))^(1/2^(k+1));
    end
    f_val = v1*v2*v3;
end 