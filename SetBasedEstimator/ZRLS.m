function [C_next, G_next, P_next, theta_next, K] = ZRLS(C, G, P, phi, y, Q_v, lambda,sigma,m)

n = size(P, 1);
I_n = eye(n);
K = P * phi' * inv(phi * P * phi' + lambda*Q_v*m);
C_next = C + K * (y - phi * C);
P_next =( (I_n - K * phi) * P )/lambda;

G_propagated_cell = cell(size(G));
for i = 1:length(G)
    G_propagated_cell{i} = (I_n - K * phi) * G{i} / sqrt(lambda);
end

q_v = zeros(1,m);
G_v = zeros(n, m);
P_v = zeros(n,n);
Gcell = {};
for j = 1:m
    q_v(1,j) = sigma;
    G_v = -K  * q_v;
    Gcell{j} = G_v;

    P_v = P_v + G_v*G_v';
    q_v = zeros(1,m);
    G_v = zeros(n, m);

end


G_next = [G_propagated_cell, Gcell];
for i = 1:length(G_next)

    G_NEXT{i} = G_next{i}';
end

theta_next = matZonotope(C_next', G_NEXT);
end