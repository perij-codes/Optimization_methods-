% Lagrange method equations derivation

function [t,fx_G,fy_G,fx_end,fy_end,ftheta,fdtheta,fv_G] = p48_solve(tspan,p)
    % Initialize symbolic variables
    syms t 
    theta = sym([]);
    theta_p = sym([]);
    dtheta_p = sym([]);
    er = sym([]);
    et = sym([]);
    rG_O = sym([]);
    rG_A = sym([]);
    rA_B = sym([]);
    v_G = sym([]);
    Ek = sym();
    Ep = sym();
    eqns = sym([]);

    % Distance and velocity vectors, kinetic and potential energy
    for i = 1:p.N
        theta(i) = symfun(sym(sprintf('theta_%d(t)', i)), t);
        theta_p(i) = sym(sprintf('theta_p_%d', i));
        dtheta_p(i) = sym(sprintf('dtheta_p_%d', i));
        er(i,:) = [cos(theta_p(i)) sin(theta_p(i)) 0];
        et(i,:) = cross(k_hat,er(i,:));
        rA_B(i,:) = p.l(i)*er(i,:);
        rG_A(i,:) = p.d(i)*er(i,:);
        rG_O(i,:) = rG_A(i,:);
        if (i > 1)
            v_G(i,:) = v_G(i-1,:) + dtheta_p(i-1)*(p.l(i-1)-p.d(i-1))*et(i-1,:)...
                + dtheta_p(i)*p.d(i)*et(i,:);
            for j = 1:(i-1)
                rG_O(i,:) = rG_O(i,:) + rA_B(j,:);
            end
        else
            v_G(1,:) = dtheta_p(1)*p.d(1)*et(1,:);
        end
        Ek = Ek + 1/2*(p.m(i)*sum(v_G(i,:).^2) + p.I_G(i)*dtheta_p(i)^2);
        Ep = Ep - p.m(i)*p.g*rG_O(i,:)*i_hat.';
    end

    % write Lagrange equations
    L = Ek - Ep;
    for i = 1:p.N
        eqns(i) = subs(diff(L,theta_p(i)),[theta_p dtheta_p],[theta diff(theta)])...
        - diff(subs((diff(L,dtheta_p(i))),[theta_p dtheta_p],[theta diff(theta)]),t);
    end

    % Integrate the equations in ODE45
    vars = theta;
    [new_eqns,new_vars] = reduceDifferentialOrder(eqns,vars);
    [M,b] = massMatrixForm(new_eqns,new_vars);

    q.all_vars = new_vars;
    q.M = M;
    q.b = b;
    q.offset = 0;

    q.num_vars = numel(q.all_vars);

    % Functions odeprog.m and odeabort.m created by Tim Franklin
    % Downloaded from MATLAB Central File Exchange
    opts = odeset('OutputFcn',@odeprog,'Events',@odeabort);

    q.type = 1;
    [t,z] = ode45(@p48_RHS,tspan,p.icv,opts,q);

    % Store the variables
    ftheta = zeros(numel(t),p.N);
    fdtheta = zeros(numel(t),p.N);
    fx_G = zeros(numel(t),p.N);
    fy_G = zeros(numel(t),p.N);
    fx_end = zeros(numel(t),p.N);
    fy_end = zeros(numel(t),p.N);
    fer = zeros(numel(t),3,p.N);
    fet = zeros(numel(t),3,p.N);
    fv_G = zeros(numel(t),3,p.N);

    % Angles of links
    for i = 1:p.N
        ftheta(:,i) = z(:,i);
        fdtheta(:,i) = z(:,i+p.N);
    end

    % Positions and velocities of link ends and centers of mass
    for i = 1:p.N
        fx = 0;
        fy = 0;
        for j = 1:i-1
            fx = fx + p.l(j)*cos(ftheta(:,j));
            fy = fy + p.l(j)*sin(ftheta(:,j));
        end
        fx_end(:,i) = fx + p.l(i)*cos(ftheta(:,i));
        fy_end(:,i) = fy + p.l(i)*sin(ftheta(:,i));
        fer(:,:,i) = [cos(ftheta(:,i)) sin(ftheta(:,i)) zeros(numel(t),1)];

        fx_G(:,i) = fx + p.d(i)*cos(ftheta(:,i));
        fy_G(:,i) = fy + p.d(i)*sin(ftheta(:,i));
        for j = 1:numel(t)
            fet(j,:,i) = cross(k_hat,fer(j,:,i));
            if (i == 1)
                fv_G(j,:,1) = fdtheta(j,1)*p.d(1)*fet(j,:,1);
            else
                fv_G(j,:,i) = fv_G(j,:,i-1)...
                    + fdtheta(j,i-1)*(p.l(i-1)-p.d(i-1))*fet(j,:,i-1)...
                    + fdtheta(j,i)*p.d(i)*fet(j,:,i);
            end
        end
    end
end
