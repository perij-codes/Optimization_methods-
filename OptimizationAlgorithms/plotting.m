% written by Sruti Vutukury
% ouput generation and plotting script
addpath(pwd,"Project_Problems_MATLAB")
addpath(pwd,"Methods")
addpath(pwd,"Outputs")

%% individual problem plots, unconstrained
clc; close all; clear all;
figure()
for p = 1:12
    for m = 1:11
        runname = strcat("problem",string(p),"_method",string(m),".mat");
        fullFileName = fullfile(pwd,"Outputs",runname);
        D = load(fullFileName);
    
        fhold = D.outputs.fhold;
        time = D.outputs.time;

        k = length(fhold);
        curvename = strcat("Method ",string(m));
        
        plot(fhold,'*','LineWidth',1,DisplayName=curvename);
        hold on
        %set(gca, 'YScale', 'log') 
        set(gca, 'XScale', 'log') 
    end
end
xlabel("Iterate",'FontSize',14)
ylabel("fvalue",'FontSize',14)
legend show
figname = strcat("Problem ",string(p));
title(figname,'FontSize',14)
gcf
savefig(figname)

%% individual problem plots, constrained
clc; close all; clear all;
figure()
for p = 13:14
    for m = 12:13
        runname = strcat("problem",string(p),"_method",string(m),".mat");
        fullFileName = fullfile(pwd,"Outputs",runname);
        D = load(fullFileName);
    
        fhold = D.outputs.fhold;
        time = D.outputs.time;

        k = length(fhold);
        curvename = strcat("Method ",string(m));
        
        plot(fhold,'*','LineWidth',1,DisplayName=curvename);
        hold on
        %set(gca, 'YScale', 'log') 
        set(gca, 'XScale', 'log') 
    end
end
xlabel("Iterate",'FontSize',14)
ylabel("fvalue",'FontSize',14)
legend show
figname = strcat("Problem ",string(p));
title(figname,'FontSize',14)
gcf
savefig(figname)

%% output all tables, unconstrained

T1 = table(); %f* table
T2 = table(); %iterations table
T3 = table(); %CPU time table
T4 = table(); %fevals table
T5 = table(); %gevals table
T6 = table(); %hevals table
for p = 1:12
    for m = 1:11
        runname = strcat("problem",string(p),"_method",string(m),".mat");
        fullFileName = fullfile(pwd,"/Outputs", runname);
        load(fullFileName)
        T1(p,m) = {f};
        T2(p,m) = {outputs.k};
        T3(p,m) = {outputs.fc};
        T4(p,m) = {outputs.gc};
        T5(p,m) = {outputs.hc};
    end
end
disp(T1)
disp(T2)
disp(T3)
disp(T4)
disp(T5)
disp(T6)

%% output all tables, constrained

T1 = table(); %f* table
T2 = table(); %iterations table
T3 = table(); %CPU time table
T4 = table(); %fevals table
T5 = table(); %gevals table
T6 = table(); %hevals table
for p = 13:14
    for m = 12:13
        runname = strcat("problem",string(p),"_method",string(m),".mat");
        fullFileName = fullfile(pwd,"/Outputs", runname);
        load(fullFileName)
        T1(p,m) = {f};
        T2(p,m) = {outputs.k};
        T3(p,m) = {outputs.fc};
        T4(p,m) = {outputs.gc};
        T5(p,m) = {outputs.hc};
    end
end
disp(T1)
disp(T2)
disp(T3)
disp(T4)
disp(T5)
disp(T6)





