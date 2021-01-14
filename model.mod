    // Model 11. Home production wiht gender implications for Mexico
    // Dynare code
    // File: model.mod
    // Marco Ramos based on José L. Torres. University of Málaga (Spain)
    // Endogenous variables
    var Y, Cm, Ch, I, K, Lam, Lah, Lbm, Lbh, Wa, Wb, R, A, B,C, Oa, Ob;
    // Exogenous variables
    varexo e, u,x;
    // Parameters
    parameters alpha, beta, delta, gamma, omega, eta, theta, rho1, rho2, i, z, j, p, H,n;
    // Calibration
    alpha = 0.35;
    beta  = 0.98;
    delta = 0.05;
    gamma = 0.40;
    omega = 0.45;
    eta   = 0.80;
    theta = 0.80;
    rho1  = 0.95;
    rho2  = 0.95;
    i  = 0.55;
    z  = 0.70;
    j  = 0.45;
    p  = 0.70;
    H  = 1.00;
    n  = .50;
    


    // Equations of the model economy
    model;
    gamma*omega*(Cm^(eta-1))/(omega*Cm^eta+(1-omega)*C*Ch^eta)
       =n*(1-gamma)/(Wa*(H-Lam-Lah));
 gamma*omega*(Cm^(eta-1))/(omega*Cm^eta+(1-omega)*C*Ch^eta)
       =(1-n)*(1-gamma)/(Wb*(H-Lbm-Lbh));

gamma*(1-omega)*(C*Ch^(eta-1))/(omega*Cm^eta+(1-omega)*C*Ch^eta)
       =n*(1-gamma)/(theta*B*(((i*(Lah^z)+(1-i)*(Lbh^z))^((theta-z)/z))
)*(i*(Lah^(z-1)))*(H-Lam-Lah));

gamma*(1-omega)*(C*Ch^(eta-1))/(omega*Cm^eta+(1-omega)*C*Ch^eta)
       =(1-n)*(1-gamma)/(theta*B*(((i*(Lah^z)+(1-i)*(Lbh^z))^((theta-z)/z))
)*((1-i)*(Lbh^(z-1)))*(H-Lbm-Lbh));

    ((Cm^(eta-1))/(omega*Cm^eta+(1-omega)*C*Ch^eta))/((Cm(+1)^(eta-1))
       /(omega*Cm(+1)^eta+(1-omega)*C*Ch(+1)^eta))=beta*(R(+1)+1-delta);

    Y = A*(K(-1)^alpha)*((j*(Lam^p)+(1-j)*(Lbm^p))^((1-alpha)/p));
    Ch = B*((i*(Lah^z)+(1-i)*(Lbh^z))^(theta/z));
    K = (Y-Cm)+(1-delta)*K(-1);
    I = Y-Cm;
    Wa = (1-alpha)*A*(K(-1)^alpha)*((j*(Lam^p)+(1-j)*(Lbm^p))^((1-alpha-p)/p))*j*(Lam^(p-1));
    Wb = (1-alpha)*A*(K(-1)^alpha)*((j*(Lam^p)+(1-j)*(Lbm^p))^((1-alpha-p)/p))*(1-j)*(Lbm^(p-1));
    R = alpha*A*(K(-1)^(alpha-1))*((j*(Lam^p)+(1-j)*(Lbm^p))^((1-alpha)/p));
    log(A) = rho1*log(A(-1))+e;
    log(B) = rho2*log(B(-1))+u;
    log(C) = rho2*log(C(-1))+x;

    Oa=H-Lam-Lah;
    Ob=H-Lbm-Lbh;
   

    end;

    // Initial values
    initval;
     Y = 1; 
    Cm = 0.75;
    Ch = 0.2; 
    Lam = 0.15;
    Lah = 0.05;
    Lbm = 0.15;
    Lbh = 0.05; 
    Oa=H-Lam-Lah;
    Ob=H-Lbm-Lbh;
    K = 3.5;
    I = 0.25;
    Wa = j*(Lam^(p-1)*(1-alpha)*Y/((j*(Lam^p)+(1-j)*(Lbm^p))^(1/p)));
    Wb = (1-j)*(Lbm^(p-1)*(1-alpha)*Y/((j*(Lam^p)+(1-j)*(Lbm^p))^(1/p)));
    R = alpha*Y/K;
    A = 1;
    B = 1;
    e = 0;
    u = 0;
    C = 1;
    x = 0;

     end;
     // Steady state
    steady;
    // Blanchard-Kahn conditions
    check;

shocks;

var x; stderr 0.01;

end;

stoch_simul;