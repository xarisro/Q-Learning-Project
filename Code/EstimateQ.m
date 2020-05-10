function [H, K] = EstimateQ(ro,Zinv,x,u,n,m,l,M)
    %Initial estimation
    H = zeros(n+m,n+m);
    K = zeros(m,n);
    d = zeros(l,1);
    e = 42;
    
    j = 0;
    %Main loop
    while e > 0.0001 || j < 10
        %First calculate d
        for i = 1:l
            d(i) = x(:,i)'*M*x(:,i) + u(i)'*ro*u(i) + x(:,i+1)'*[eye(n) K']*H*[eye(n); K]*x(:,i+1);
        end
        
        W = Zinv * d;
        
        %New estimation of H and K
        H = [W(1) W(5) W(6) W(7)
             W(5) W(2) W(8) W(9)
             W(6) W(8) W(3) W(10)
             W(7) W(9) W(10) W(4)];
        
         Knew = -inv(H(4,4))*H(4,1:3);
         
         e = (Knew-K)*(Knew-K)';
         K = Knew;
         
         j = j+1;
    end
end

