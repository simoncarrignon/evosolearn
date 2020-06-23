% Stabilizng selection in asexual population
%
% Feb.29. 2020

close all
clear all

if 1
    rng('shuffle');
    seed=rng;
    save seed;
else
    load seed;
    rng(seed)
end
ones(10)*ones(10);

% Parameters
% theta=0;          % optimum trait
b=2;                % average number of offspring
T=10000;  
skip=2000;          % number of transient generations to skip
Runs=40;

Delta=[0 .1 .2 .4];

save_figures=1;

Sigma=[1 2 4 10000];
Mu=[0.001 0.01];
M=[0.05 0.1 0.2];
C=[500 1000 2000];

% sigma=1;        % strength of stabilizing selection
% K=1000;         % carrying capacity
% mu=0.01;        % probability of mutation
% m=0.2;          % st. dev. of mutation


for i0=1:1
    delta=Delta(i0);

    % Useful
    N=nan(T+1,1);
    Mean_x=nan(T+1,1);
    Var_x=nan(T+1,1);
    Mean_v=nan(T+1,1);

    NN=nan(Runs,1);
    MX=nan(Runs,1);
    VX=nan(Runs,1);
    MV=nan(Runs,1);

    Z =nan(length(Mu),length(Sigma),length(M),length(C),2);    % variance
    Z1=nan(length(Mu),length(Sigma),length(M),length(C));      % error bars

    figure('Position',[50 50 1500 900],'Visible','on');

    tic
    for i1=1:length(Mu)
        mu=Mu(i1);
        for i2=1:length(Sigma)
            sigma=Sigma(i2);
            alpha=1/2/sigma^2;
            for i3=1:length(M)
                m=M(i3);
                for i4=1:length(C)
                    K=C(i4);               
                    for run=1:Runs
                        run
                        n=K;                    % population size
                        x=2*rand(n,1)-1;        % initial state
                        theta=delta*randn(T,1);  % precompute optimum values

                        for t=1:T

                            % statistics for this generation
                            Mean_x(t)=mean(x);
                            Var_x(t)=var(x);
                            N(t)=n;

                            % selection
                            w=exp(-alpha*(x-theta(t)).^2);      % condition
                            v=1./(1+(b-1)/K*n./w);              % viability
                            Mean_v(t)=mean(v);
                            parents=find(v>rand(n,1));          % indices of survivors
                            no_of_parents=length(parents);      % their number

                            % reproduction
                            R=poissrnd(b,no_of_parents,1);      % number of offspring for each surviving parent
                            n=sum(R);                           % total number of offspring
                            offspring=rude(R,x(parents))';      % their traits before mutation 

                            % mutation
                            mutants=find(rand(n,1)<mu);         % mutants
                            offspring(mutants)=offspring(mutants)+m*randn(length(mutants),1);   % their traits

                            x=offspring;                        % new generation
                        end

                    % data for each run
                    NN(run)=mean(N(skip:T));
                    MX(run)=mean(Mean_x(skip:T));
                    VX(run)=mean(Var_x(skip:T));
                    MV(run)=mean(Mean_v(skip:T));

                    end

                    gam=m^2/(2*sigma^2);
                    hermisson=m^2 *(gam*mean(NN)+1)/(4*gam*mean(NN))*...
                        (sqrt(1+2*gam*mean(NN)*4*mu*mean(NN)/(gam*mean(NN)+1)^2)-1);

                    Z(i1,i2,i3,i4,1)=mean(VX);    % observed
                    Z(i1,i2,i3,i4,2)=hermisson;   % predicted

                    Z1(i1,i2,i3,i4) =std(VX);      % error bars

                    [mean(NN)  mean(MV) mean(MX) mean(VX) hermisson]

                end
            end
        end
    end
    toc

    I=1;
    A=[1 2 3; 4 5 6; 7 8 9];
    for i1=1:length(Mu)
        for i2=1:length(Sigma)
            subplot(length(Mu),length(Sigma),I);
    %         plot(A,squeeze(Z(i1,i2,:,:,1)),'o')
            errorbar(A,squeeze(Z(i1,i2,:,:,1)),squeeze(Z1(i1,i2,:,:)),'o','MarkerSize',8)
            hold on
            plot(A,squeeze(Z(i1,i2,:,:,2)),'xk','MarkerSize',8)
            hold off
            xlim([0 10])
            ylim([0 inf])
            if I==1
                title({['$\delta$=',num2str(delta)];['$\mu$=',num2str(Mu(i1)),', $\sigma$=',num2str(Sigma(i2))]},'Interpreter','Latex','fontweight','bold','fontsize',16);
            else
                title(['$\mu$=',num2str(Mu(i1)),', $\sigma$=',num2str(Sigma(i2))],'Interpreter','Latex','fontweight','bold','fontsize',16);
            end
            set(gca,'XTick',[2 5 8]);
            set(gca,'XTickLabel',M);
            xlabel('m');
            if I==1
                legend('K=500','K=1000','K=2000','Location','NorthWest')
            end
            I=I+1;
        end
    end

    if save_figures
        eps_file=sprintf('delta.%2.1f.eps',delta);
        print('-depsc2',eps_file);
        png_file=sprintf('delta.%2.1f.png',delta);
        print('-dpng',png_file);    
    end

end

% figure
% subplot(2,2,1)
% plot(N)
% ylabel('N');
% subplot(2,2,3)
% plot(Mean_x)
% ylabel('X')
% subplot(2,2,2)
% plot(Mean_v)
% ylabel('$\overline{v}$','Interpreter','Latex')
% subplot(2,2,4)
% plot(Var_x)
% ylabel('var')

