% clear; close all; clc;
% sruti vutukury last editted: 04/14/2022

% Uniform pendulums of unit mass and unit length
% no rigid body motion
% neat multi DOF behavior for small thetas
% chaotic behavior for larger thetas
% references: 
% (S.Widnall) MIT Lecture Notes: Lecture L20 - Energy Methods: Lagrange Equations,
% (M. Buche) Cornell University MAE 4730 Intermediate Dynamics Manual,
% (Yesilyurt, Boran ) "Equations of Motion Formulation of a Pendulum Containing N-point Masses"
% the EOMS are verified, but animation looks suspicious, debug....

p.g = -9.81;
tspan = [0 5];

case_num = 1;

for nlinks = [3,5,7,10]
  
    p.n = nlinks;
   
    % ICs = [positions velocities];
    % all_sets = [pi/6 -pi/6 pi/6 1 0.001 -0.001];  %this one I liked
    all_sets = [pi/4*ones(1,p.n), zeros(1,p.n);
                pi/6*ones(1,p.n), 0.001*ones(1,p.n);
                pi/6*ones(1,p.n), -0.001*ones(1,p.n);
                -0.5 + (0.5+0.5)*rand(1,p.n), -0.01 + (0.01+0.01)*rand(1,p.n);
                pi/16*ones(1,p.n), 0.01*ones(1,p.n)];
    num_ics = size(all_sets);

    for i = 1:num_ics(1)
        p.case_num = case_num;
        try
            p.n = nlinks;
            ICs = all_sets(i,:);
    
            p.h = ones(1,p.n);
            p.m = ones(1,p.n);
    
            %derive equations of motion
            [eqns,p] = derive_lagrange(p);
%             [eqns_check, p] = check_eqns2(p); %check double pendulum equations
%             [eqns_check, p] = check_eqns3(p); %check triple pendulum equations
            disp(sprintf('derived equations for case: %d',p.case_num))
            p.eqns = eqns;
    
            %simulate
            opts.reltol = 1e-8; opts.abstol = 1e-8;
            [t_out,z_out]=ode45(@RHSfunction,tspan,ICs,opts,p);
%             p.eqns = eqns_check; %checking double and triple pendulum EOMs
%             [t_out_check,z_out_check]=ode45(@RHSfunction,tspan,ICs,opts,p);
            disp(sprintf('solved for case: %d',p.case_num))
    
            %get positions of ends of pendulums
            th_out = z_out(:,1:p.n);
            x_end = zeros(numel(t_out),p.n);
            y_end = zeros(numel(t_out),p.n);
            x_end(:,1) = p.h(1)*sin(th_out(:,1));
            y_end(:,1) = -p.h(1)*cos(th_out(:,1));
            
            for j = 2:p.n %iterate over links
                x_end(:,j) = x_end(:,j-1) + p.h(j)*sin(th_out(:,j));
                y_end(:,j) = y_end(:,j-1) - p.h(j)*cos(th_out(:,j));
            end
            
            save(sprintf("data_case%d",p.case_num))
    
            %plot positions
            plotting(t_out,x_end,y_end,th_out,p);
            disp(sprintf('plotted for case: %d',p.case_num))
    
            %animations
            animate_pendulum(t_out,x_end,y_end,p)
            disp(sprintf('animated and saved for case: %d',p.case_num))
    
            disp(sprintf('completed case: %d',p.case_num))
            case_num = case_num + 1;
        catch
            disp(sprintf('couldnt complete case: %d',p.case_num))
            case_num = case_num + 1;
            continue
        end
    end
end

function [eqns,p] = derive_lagrange(p)
    %derive the equations of motion with respect to minimal coordinates
    %lagrange method
    syms t;
    th = sym([]);
    thdot = sym([]);
    eqns = sym([]);
    n = p.n;

    for i = 1:n
        th(i) = symfun(sprintf('th_%d(t)', i), t);
        thdot(i) = sym(sprintf('thdot_%d', i));
    end
    x = sym([]);
    y = sym([]);
    xdot = sym([]);
    ydot = sym([]);
    for i = 1:n
        if i == 1
            x(1) = p.h(1)*sin(th(1));
            y(1) = -p.h(1)*cos(th(1));
            xdot(1) = p.h(1)*cos(th(1))*thdot(1);
            ydot(1) = p.h(1)*sin(th(1))*thdot(1);
        else
            x(i) = x(i-1) + p.h(i)*sin(th(i));
            y(i) = y(i-1) -p.h(i)*cos(th(i));
            xdot(i) = xdot(i-1) + p.h(i)*cos(th(i))*thdot(i);
            ydot(i) = ydot(i-1) + p.h(i)*sin(th(i))*thdot(i);
        end
    end
    PE_total = 0;
    KE_total = 0;
    for i = 1:n
        PE = p.m(i)*p.g*y(i);
        PE_total = PE_total+PE;
        KE = 0.5*p.m(i)*((xdot(i))^2 + (ydot(i))^2);
        KE_total = KE_total + KE;
    end
    L = KE_total - PE_total; %substitute all parts with th and thdot to reduce to minimal coords
    for i = 1:n
        %we want to reduce to minimal coordinates: th, thdot
        eqns(i) = subs(diff(L,th(i)),[th thdot],[th diff(th)])...
        - diff(subs((diff(L,thdot(i))),[th thdot],[th diff(th)]),t);
    end

    %reduce system of higher-order DEs to equivalent system of first-order DEs
    vars = th;
    [new_eqns,new_vars] = reduceDifferentialOrder(eqns,vars);
    p.all_vars = new_vars;
    p.num_vars = numel(p.all_vars);
    eqns = new_eqns;
end

function dz = RHSfunction(~,z,p)
    N = length(z)/2; %minimal coordinates for each pendulum: th,thdot
    dz = zeros(2*N,1);
    dz(1:N) = z(N+1:end);
    %subsititude thddot to calculate acceleration
    for i = 1:N
        %propogate state
        eqn = p.eqns(i);
        dz(N+i)=subs(eqn,p.all_vars,z);
    end
end

function plotting(t_out,x_end,y_end,th_out,p)
    figure()
    for i = 1:p.n
        txt = ['link ',num2str(i)];
        plot(x_end(:,i),y_end(:,i),'DisplayName',txt)
        hold on;
    end
    xlabel('x_{end}')
    ylabel('y_{end}')
    hold off
    legend show
    saveas(gcf,sprintf('y_vs_x_case%d.jpg',p.case_num))
    close all
       
    figure()
    for i = 1:p.n
        txt = ['link ',num2str(i)];
        plot(t_out,th_out(:,i),'DisplayName',txt)
        hold on;
    end
    xlabel('t')
    ylabel('theta')
    hold off
    legend show
    saveas(gcf,sprintf('theta_vs_t_case%d.jpg',p.case_num))
    close all
end

function animate_pendulum(t_out,x_end,y_end,p)
    % setup graphic
    p.animation_timescale = 1; %increase to speed up
    p.filename = sprintf('nlinkpend_case%d.gif',p.case_num);
    figure();
    grid on
    grid minor
    axis equal
    axis([-1 1 -1 1]*1.1*sum(p.h))
    hold on
    set(gca,'fontsize',10);
    timer_text = annotation(gcf,'textbox',[0.825 0.75 0.1 0.1],...
            'string',{'$t = 0.000$'},'fontsize',10,'Interpreter','LaTeX');

    % links' and masses' initial position
    link = zeros(1,p.n);
    mass = zeros(1,p.n);
    link_width = 2;
    mass_size = 20;

    link(1) = plot([0 x_end(1,1)],[0 y_end(1,1)],'k','linewidth',link_width);
    plot(0,0,'+','markersize',20); %origin
    mass(1) = plot(x_end(1,1),y_end(1,1),'.','markersize',mass_size);
    for i = 2:p.n
        link(i) = plot([x_end(1,i-1) x_end(1,i)],[y_end(1,i-1) y_end(1,i)],...
            'k','linewidth',link_width);
        mass(i) = plot(x_end(1,i),y_end(1,i),'.','markersize',mass_size);
    end

    % update positions of the links and masses in real time
    pause(1)
    xint = zeros(1,p.n);
    yint = zeros(1,p.n);
    time = 0; tic
    counter = 0;
    while time < t_out(end)/p.animation_timescale
        for i = 1:p.n
            %vq = interp1(x,v,xq)
            xint(i) = interp1(t_out,x_end(:,i),time*p.animation_timescale);
            yint(i) = interp1(t_out,y_end(:,i),time*p.animation_timescale);
            if i == 1
                set(mass(i),'xdata',xint(1),'ydata',yint(1))
                set(link(i),'xdata',[0 xint(1)],'ydata',[0 yint(1)]);
            else
                set(mass(i),'xdata',xint(i),'ydata',yint(i))
                set(link(i),'xdata',[xint(i-1) xint(i)],'ydata',[yint(i-1) yint(i)])
            end
        end
        drawnow;
        time = toc;

        % update timer
        set(timer_text,'string',sprintf('$t = %.3f$',time*p.animation_timescale));
        xlabel('$x(t)$','Interpreter','LaTeX');
        ylabel('$y(t)$','Interpreter','LaTeX');
        title(sprintf('animation with sped up: %dx',p.animation_timescale))

        frame = getframe();
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if counter == 0
            imwrite(imind,cm,p.filename,'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,p.filename,'gif','WriteMode','append');
        end
        counter = counter + 1;
        
    end
    close all
end

function [eqns,p] = check_eqns2(p)
    %checking against double pendulum equations from MIT lecture notes
    % (S.Widnall) MIT Lecture Notes: Lecture L20 - Energy Methods: Lagrange’s Equations
    syms t;
    th = sym([]);
    thdot = sym([]);
    eqns = sym([]);
    n = p.n;

    for i = 1:n
        th(i) = symfun(sprintf('th_%d(t)', i), t);
        thdot(i) = sym(sprintf('thdot_%d', i));
    end

    KE_total = 0.5*p.m(1)*p.h(1)^2*thdot(1)^2 + 0.5*p.m(2)*(p.h(1)^2*thdot(1)^2 + p.h(2)^2*thdot(2)^2+ ...
                    2*p.h(1)*p.h(2)*thdot(1)*thdot(2)*cos(th(1)-th(2)));
    PE_total = -(p.m(1)+p.m(2))*p.g*p.h(1)*cos(th(1)) - p.m(2)*p.g*p.h(2)*cos(th(2));

    L = KE_total - PE_total; %substitute all parts with th and thdot to reduce to minimal coords
    for i = 1:p.n
        %we want to reduce to minimal coordinates: th, thdot
        eqns(i) = subs(diff(L,th(i)),[th thdot],[th diff(th)])...
        - diff(subs((diff(L,thdot(i))),[th thdot],[th diff(th)]),t);
    end

    %Reduce system of higher-order DEs to equivalent system of first-order DEs
    vars = th;
    [new_eqns,new_vars] = reduceDifferentialOrder(eqns,vars);
    p.all_vars = new_vars;
    p.num_vars = numel(p.all_vars);
    eqns = new_eqns;
end

function [eqns,p] = check_eqns3(p)
    %checking against triple pendulum equations from research gate paper
    % (Yesilyurt, Boran ) "Equations of Motion Formulation of a Pendulum Containing N-point Masses"
    syms t;
    th = sym([]);
    thdot = sym([]);
    eqns = sym([]);
    n = p.n;

    for i = 1:n
        th(i) = symfun(sprintf('th_%d(t)', i), t);
        thdot(i) = sym(sprintf('thdot_%d', i));
    end
    term1 = 0.5*p.m(1)*p.h(1)^2*thdot(1)^2; 
    term2 = p.m(1)*p.g*p.h(1)^2*cos(th(1)); 
    term3 = 0.5*p.m(2)*(p.h(1)^2*thdot(1)+p.h(2)^2*thdot(2)^2+ ...
        2*p.h(1)*p.h(2)*cos(th(1)-th(2))*thdot(1)*thdot(2));
    term4 = p.m(2)*p.g*(p.h(1)*cos(th(1)) +p.h(2)*cos(th(2)));
    L12 = term1 + term2 + term3 + term4;
    
    K3 = 0.5*p.m(3)*(p.h(1)^2*thdot(1)^2 + p.h(2)^2*thdot(2)^2+ ...
        p.h(3)^2*thdot(3)^2 + 2*p.h(1)*p.h(2)*cos(th(1)-th(2))*thdot(1)*thdot(2) + ...
        2*p.h(1)*p.h(3)*cos(th(1)-th(3))*thdot(1)*thdot(3) + ...
        2*p.h(2)*p.h(3)*cos(th(2)-th(3))*thdot(2)*thdot(3));

    U3 = -p.m(3)*p.g*(p.h(1)*cos(th(1)) + p.h(2)*cos(th(2)) + p.h(3)*cos(th(3)));
    L3 = K3-U3;

    L = L12 + L3; %substitute all parts with th and thdot to reduce to minimal coords
    for i = 1:p.n
        %we want to reduce to minimal coordinates: th, thdot
        eqns(i) = subs(diff(L,th(i)),[th thdot],[th diff(th)])...
        - diff(subs((diff(L,thdot(i))),[th thdot],[th diff(th)]),t);
    end

    %Reduce system of higher-order DEs to equivalent system of first-order DEs
    vars = th;
    [new_eqns,new_vars] = reduceDifferentialOrder(eqns,vars);
    p.all_vars = new_vars;
    p.num_vars = numel(p.all_vars);
    eqns = new_eqns;
end