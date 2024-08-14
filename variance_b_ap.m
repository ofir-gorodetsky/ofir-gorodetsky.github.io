function compute_empirical_variance()
    clear;
    
    t = 3*10^4;
    %generate all sums of two squares up to X := t^2.
    a(1,1:(t+1)) = (0:t).^2;
    temp = repmat(a,(t+1),1);
    clear a;
    sumOfSquares = temp + temp';
    clear temp;
    u2 = sort(unique(sumOfSquares));
    clear sumOfSquares;
    b1 = find(u2>t^2,1);
    u1 = u2(1:b1-1);
    clear u2;

    %generate prime moduli, from 1 to X, with (roughly) constant spacing on a logarithmic scale
    q1 = primes(t^2)';
    q2 = zeros(floor(log(t^2)/log(1.1)),1);
    for j=1:floor(log(t^2)/log(1.1))
        q2(j,1) = q1( find( q1 > round((1.1)^j), 1 ) );
    end
    q = unique(q2);
    clear q2;
    
    %generate prime moduli, from X/2 to X, with (roughly) constant spacing on a logarithmic scale
    q2 = zeros(100,1);
    for j=1:100
        q2(j,1) = q1( find( q1 > round(t^2/2^(j/100)), 1 ) );
    end
    clear q1;
    q_large = unique(q2);
    
    qcount = zeros(size(q,1),1);
    qcount_large = zeros(size(q_large,1),1);
    
    %compute variance in arithmetic progressions
    for k = 1:size(q,1)
        buckets = zeros(q(k,1),1);
        for a = 1:size(u1,1)
            buckets(mod(u1(a,1), q(k,1))+1,1) = buckets(mod(u1(a,1), q(k,1))+1,1) + 1;
        end
        b1 = buckets(2:end,1);
        clear buckets;
        normq = t^2/((sqrt(log(t^2)))*(q(k,1)));
        qcount(k,1) = var(b1)/normq;
        clear b1;
    end
    x_axis = 1 - (log(q)/(log(t^2)));
    save('results.mat', 'x_axis', 'qcount')
    
    %compute variance in arithmetic progressions for large prime moduli
    for k = 1:size(q_large,1)
        buckets = zeros(q_large(k,1),1);
        for a = 1:size(u1,1)
            buckets(mod(u1(a,1), q_large(k,1))+1,1) = buckets(mod(u1(a,1), q_large(k,1))+1,1) + 1;
        end
        b1 = buckets(2:end,1);
        clear buckets;
        normq = t^2/((sqrt(log(t^2)))*(q_large(k,1)));
        qcount_large(k,1) = var(b1)/normq;
        clear b1;
    end
    x_axis = 1 - (log(q_large)/(log(t^2)));
    save('results_large_p.mat', 'x_axis', 'qcount_large')
end
