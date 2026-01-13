function [the_matrix,T1,T2] = generate_linear_system(inlier, outlier, tol_lower,tol_upper,F)

W = 100; H = 100;
N = inlier + outlier;
the_matrix = zeros(N, 9);
T1 = zeros(N,3);
T2 = zeros(N,3);

% Fix scale so tol is meaningful
F = F / norm(F,'fro');

%% -------- inliers --------
cnt = 1;
while cnt <= inlier
    a  = W*rand();
    b  = H*rand();
    x1 = [a; b; 1];
    l  = F*x1;   % l(1)*u + l(2)*v + l(3) = 0

    success = false;

    for tries = 1:10000
        if abs(l(1)) > abs(l(2)) && abs(l(1)) > 1e-12
            v = H*rand();
            u = -(l(2)*v + l(3)) / l(1);
        elseif abs(l(2)) > 1e-12
            u = W*rand();
            v = -(l(1)*u + l(3)) / l(2);
        else
            break
        end

        if (u>=0 && u<=W && v>=0 && v<=H)
            success = true;
            break
        end
    end

    if ~success
        continue   % 重抽一个新的 (a,b)
    end
    x2 = [u; v; 1];   
    T1(cnt,:) = x1';
    T2(cnt,:) = x2';

    the_matrix(cnt,:) = [u*a, u*b, u,  v*a, v*b, v,  a, b, 1];
    cnt = cnt + 1;
end

%% -------- outliers --------
for j = 1:outlier
    r = 0;
    success = false;

    for tries = 1:100000
        a  = W*rand();
        b  = H*rand();
        u  = W*rand();
        v  = H*rand();
        x1 = [a; b; 1];
        x2 = [u; v; 1];

        r = abs(x2' * F * x1);
        if r >= tol_lower && r<=tol_upper
            success = true;
            break
        end
    end

    if ~success
        error('Failed to generate outlier with residual >= tol. Try smaller tol.');
    end
    T1(inlier + j,:) = x1';
    T2(inlier + j,:) = x2';
    the_matrix(inlier + j, :) = [u*a, u*b, u,  v*a, v*b, v,  a, b, 1];
end
perm = randperm(N);
the_matrix = the_matrix(perm,:);
T1 = T1(perm,:);
T2 = T2(perm,:);

end

