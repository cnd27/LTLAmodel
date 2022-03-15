function bb = getBetaInit(p,agemix,SS)

    F = @(beta,SS1,SS2,SS3)[0 beta*agemix(1,1)*SS1 beta*p.tau*agemix(1,1)*SS1 0 beta*agemix(1,2)*SS1 beta*p.tau*agemix(1,2)*SS1 0 beta*agemix(1,3)*SS1 beta*p.tau*agemix(1,3)*SS1;
        zeros(2,9);
        0 beta*agemix(2,1)*SS2 beta*p.tau*agemix(2,1)*SS2 0 beta*agemix(2,2)*SS2 beta*p.tau*agemix(2,2)*SS2 0 beta*agemix(2,3)*SS2 beta*p.tau*agemix(2,3)*SS2;
        zeros(2,9);
        0 beta*agemix(3,1)*SS3 beta*p.tau*agemix(3,1)*SS3 0 beta*agemix(3,2)*SS3 beta*p.tau*agemix(3,2)*SS3 0 beta*agemix(3,3)*SS3 beta*p.tau*agemix(3,3)*SS3;
        zeros(2,9)];

    V = [p.epsilon 0 0 zeros(1,6);
         -p.d(1,1)*p.epsilon p.gamma 0 zeros(1,6);
         -(1-p.d(1,1))*p.epsilon 0 p.gamma zeros(1,6)
        zeros(1,3) p.epsilon 0 0  zeros(1,3);
        zeros(1,3)  -p.d(2,1)*p.epsilon p.gamma 0 zeros(1,3);
        zeros(1,3)  -(1-p.d(2,1))*p.epsilon 0 p.gamma zeros(1,3)
        zeros(1,6) p.epsilon 0 0 ;
        zeros(1,6)  -p.d(3,1)*p.epsilon p.gamma 0;
        zeros(1,6)  -(1-p.d(3,1))*p.epsilon 0 p.gamma];

    M = @(beta,SS1,SS2,SS3)F(beta,SS1,SS2,SS3)*inv(V);

    lval = @(beta,SS1,SS2,SS3)eig(M(beta,SS1,SS2,SS3));
    
    bb = zeros(p.Nltlas,1);
    for i = 1:p.Nltlas
        lval2 = @(beta)(abs(lval(beta,SS(i,1),SS(i,2),SS(i,3))-p.R0));
        getFirst = @(M,i) M(i);
        lval3 = @(beta)getFirst(lval2(beta),1);
        bb(i) = fminsearch(@(beta)(lval3(beta)),1);
    end
    bb = abs(bb);

end