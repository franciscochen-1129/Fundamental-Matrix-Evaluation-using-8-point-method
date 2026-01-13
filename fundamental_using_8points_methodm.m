clear;
clc;

seed = 1;
rng(seed);


inlier = 10;
outlier = 100;
total = inlier + outlier;
tol_lower = 0.1;
tol_upper = 10;

[F,vector_F] = generate_random_f();
[A,T1,T2] = generate_linear_system(inlier, outlier, tol_lower,tol_upper,F);


%%8 points method%%
mu1 = mean(T1(:,1:2), 1);
T1_new = T1;
T1_new(:,1:2) = T1(:,1:2) - mu1;
T1_new(:,3) = 1;
mu2 = mean(T2(:,1:2), 1);
T2_new = T2;
T2_new(:,1:2) = T2(:,1:2) - mu2;
T2_new(:,3) = 1;

d_T1 = 0;
d_T2 = 0;
for i = 1:total
    di1 = sqrt(T1_new(i,1)^2+T1_new(i,2)^2);
    di2 = sqrt(T2_new(i,1)^2+T2_new(i,2)^2);
    d_T1 = d_T1 + di1;
    d_T2 = d_T2 + di2;
end
d_T1 = d_T1/total;
d_T2 = d_T2/total;
s_T1 = sqrt(2)/d_T1;
s_T2 = sqrt(2)/d_T2;
T1_new(:,1:2) = s_T1 * T1_new(:,1:2);
T2_new(:,1:2) = s_T2 * T2_new(:,1:2);
final_T1 = [s_T1, 0 , -s_T1*mu1(1);...
    0, s_T1, -s_T1*mu1(2);...
    0, 0, 1];
final_T2 = [s_T2, 0 , -s_T2*mu2(1);...
    0, s_T2, -s_T2*mu2(2);...
    0, 0, 1];

A_norm = zeros(total,9);
for j = 1:total
    u1 = T1_new(j,1);
    v1 = T1_new(j,2);
    u2 = T2_new(j,1);
    v2 = T2_new(j,2);
    x = [u1*u2,v1*u2,u2,u1*v2,v1*v2,v2,u1,v1,1];
    A_norm(j,:) = x;
end
[~,~,V] = svd(A_norm);
estimated_f = V(:,end);
estimated_f = reshape(estimated_f,3,3)';
[U,S,V2] = svd(estimated_f);
S(3,3) = 0;
estimated_f = U * S * V2';
new_F = final_T2'*estimated_f*final_T1;
new_F = new_F./norm(new_F,'fro');


fgt  = F(:);     fgt  = fgt / norm(fgt);
fhat = new_F(:); fhat = fhat / norm(fhat);

theta = acos( min(1, abs(fgt' * fhat)) ); 
theta_deg = theta * 180/pi;


alpha = sum(sum(F .* new_F)) / sum(sum(new_F .* new_F));
relerr = norm(F - alpha*new_F, 'fro') / norm(F, 'fro');

fprintf('Angle error: %.6f deg\n', theta_deg);
fprintf('Relative Frobenius error: %.6e\n', relerr);
