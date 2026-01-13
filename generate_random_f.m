function [F,vector_F] = generate_random_f()
    q = randn(4,1);
    q = q / norm(q);
    w = q(1); x = q(2); y = q(3); z = q(4);

    R = [1-2*(y^2+z^2), 2*(x*y - z*w), 2*(x*z + y*w);
         2*(x*y + z*w), 1-2*(x^2+z^2), 2*(y*z - x*w);
         2*(x*z - y*w), 2*(y*z + x*w), 1-2*(x^2+y^2)];
    t = randn(3,1);
    t = t / norm(t);

    tx = [  0   -t(3)  t(2);
           t(3)   0   -t(1);
          -t(2)  t(1)   0 ];

    F = tx * R;
    F = F / norm(F, 'fro');
    vector_F = F(:);
end