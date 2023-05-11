function [A,B,C,D,x0, sv]=po_moesp(u,y,s,n)
    u_col = u(1:2*s);
    y_col = y(1:2*s);
    
    u_row = u(2*s:end);
    y_row = y(2*s:end);
    
    U_hankel = hankel(u_col, u_row);
    Y_hankel = hankel(y_col, y_row);

    U2 = [U_hankel(s+1:end, :); U_hankel(1:s, :)];
    input = [U2; Y_hankel];
    [q, r] = qr(input');
    L = transpose(r);
    L_32 = L(3*s+1:4*s, s+1:3*s);

    l = 1;
    [U, S, V] = svd(L_32, 'econ');
    C = U(1, 1:n);
    A = U(1:s*l-l, 1:n)\U(1+l:end, 1:n);

    N = length(y);
    B_curr = zeros(N,n);
    for k= 1:n
        temp_B = zeros(n,1);
        temp_B(k, 1)= 1;
        sys = ss(A,temp_B,C,0,1);
        B_curr(:,k) = lsim(sys, u);
    end
    x0_temp = zeros(N,n);
    for k= 1:N
        x0_temp(k,:) = [C*A^(k-1)];
    end
    phi_temp = [x0_temp B_curr u'];
    size(phi_temp)
    size(y)
    phi = phi_temp\y;
    x0= phi(1:n);
    B = phi(n+1:2*n);
    D = phi(2*n+1);
    sv = svd(L_32, 'econ');
  end