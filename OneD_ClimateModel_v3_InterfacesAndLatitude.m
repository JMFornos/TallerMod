clear all
close all

latitude=discretize('varas');
N_zones=length(latitude);

%% PARAMETERS %% 
DEBUG_STAT=0; %% Enable for debug plots (stationary)
DEBUG_DYN=0; %% Enable for debug plots (dynamic)

global S epsilon alpha backRadCoef h_AtmEarthConvec F G delta eps
delta=0.25;  %% Absolute tolerance
eps=0.1;  %% Relative Tolerance for convergence control

S=zeros(N_zones,1);
for i=1:1:N_zones
    S(i)=calc_inso_by_lat(latitude(i),180);
end

%S=calc_inso_by_lat_ez(latitude,300); %% Incoming radiation profile

epsilon=ones(N_zones,1)*0.8; %% Emmissivity of the atmosphere (per zone)
alpha=ones(N_zones,1)*0.25; %% Albedo of each zone
backRadCoef=0.75; %% Fracción de la radiación atmosférica emitida a la tierra
h_AtmEarthConvec=0.1; %% Término calor convección - Coeficiente de película entre atmósfera y tierra
% Vector F and G represent the exchange in each INTERFACE. Thus the value
% of F(i) describes the flux between zones i and i+1. F(i-1) then describes
% the flux between zones i and i-1. 
F=ones(N_zones-1,1)*0.5; 
G=zeros(N_zones-1,1);

%%%% SOLVE STATIONARY PROBLEM
maxit_stat=100; % Number of iterations for newtons method
run_stat=1; 
if run_stat==1
    Temps_init_stat=ones(2*N_zones,1)*273;
    [temps_cell_stat, max_q_reached]=calc_stationary(Temps_init_stat,maxit_stat);
    meantemps(temps_cell_stat,max_q_reached) % We calculate the 
    plot_stationary(latitude, temps_cell_stat,max_q_reached,DEBUG_STAT) 
end

%% Parameters specific to the dynamic problem 
ts=200; %% Timestep
numof_ts=1000;  %% Number of timesteps to calculate
t0=0; %% Initial time
time=zeros(numof_ts,1);
time(1)=t0;
Cp=[ones(N_zones,1)*10000;ones(N_zones,1)*100000]; % Prepared to consider the different masses and
%inertias of earth and atmosphere: Top half is atmosphere, bottom half is earth
maxit_dyn= 20; % Maximum iteration number for the newtons method
th=0.5; %% Theta method theta parameter value


run_dynamic=1; %% Enable the calculation of the dynamic problem
if run_dynamic==1
    
    % Data structure: Each row is an iteration, each column is a timestep
    temps_dyn=cell(maxit_dyn,numof_ts);  
    %Distribución inicial de temperaturas
    temps_dyn{1,1}=ones(2*N_zones,1)*273;
    iters_used=zeros(numof_ts,1);
    %% NEWTON
    %Iterate until either max iterations are calculated or the convergence criterion is met
    for t=1:1:numof_ts
        if t>1
            time(t)=time(t-1)+ts;
        end
        %Define the seed of the iteration process using forward euler
        T_iter = temps_dyn{1,t} +(ts)*balance(N_zones,temps_dyn{1,t})./Cp;

        %Iterate with newton
        for q=1:1:maxit_dyn
            funcF_j=balance(N_zones,temps_dyn{q,t})./Cp;
            funcF_j_1=balance(N_zones,T_iter)./Cp;
            %Define F vector
            vectF=T_iter-temps_dyn{q,t}-(ts)*(th*funcF_j+(1-th)*funcF_j_1);
            %Define DF. Evaluate the value of the Jacobian matrix for that
            DF=eye(2*N_zones)-(ts)*(1-th)*jacobian(N_zones,T_iter);
            %Solve the linear system DF(y(q))*w=-F
            deltaTemp_iter=linsolve(DF,-vectF);
            %Add w to yq to obtain yq+1
            T_iter=T_iter+deltaTemp_iter;
            %Evaluate convergence criteria, break for loop if both met
            %Recalculate F with the updated value to use in the convergence
            %criteria
            funcF_j_1_updt=balance(N_zones,T_iter)./Cp;
            vectF=T_iter-temps_dyn{q,t}-ts*(th*funcF_j+(1-th)*funcF_j_1_updt);
            if (norm(deltaTemp_iter)/1+norm(T_iter))<eps && delta>norm(vectF)
                disp("Timestep: "+str(time(t)))
                disp("Convergence achieved in iteration " + str(q))
                break
            end
            if q<maxit_dyn
                temps_dyn{q+1,t}=T_iter;
            end
        end
        iters_used(t)=q;
        if (DEBUG_DYN==1 && mod(t,1000)==0)
            plot_iterations_dynamics(latitude, temps_dyn,t,iters_used(t))
        end
        if t<numof_ts
            temps_dyn{1,t+1}=temps_dyn{q,t};
        end
    end
    plot_dynamic(latitude,temps_dyn,iters_used,numof_ts)
end

function sf_boltz=boltz()
%% This function takes no inputs and returns the value of stefan-boltzmanns constant
    sf_boltz=5.67*10^-8;
end

function latitude=discretize(varargin)
%% This function calculates the mean latitude of each zone taking the total number of zones as an input
    if (strcmpi(varargin{1},'varas')==1)
        disp('Discretizing as described in class: polar caps, equatorial zone, zone centered around each tropic and 2 further zones')
        lat_north=[78.37;58.57;42.79;23.26];
        lat_south=-1*flip(lat_north);
        latitude=[lat_north;0;lat_south];        
    elseif isinteger(varargin{1})
        N=varargin{1};
        disp('Discretizing in '+str(N)+' regular intervals')
        latitude=ones(N,1); %% Latitud media de cada zona
        latitude(1)=90-180/(2*N);
        for n=2:1:N
            latitude(n)=latitude(n-1)-(180/N);
        end
    else
        disp("Bad argument. Discretizing with 11 regular zones")
        N=11;
        latitude=ones(N,1); %% Latitud media de cada zona
        latitude(1)=90-180/(2*N);
        for n=2:1:N
            latitude(n)=latitude(n-1)-(180/N);
        end
    end
end

%% GRUPO SONSOLES + DIEGO %%
function Sdia=calc_inso_by_lat(latitude,dia)
    % Sdia = calc_inso_by_lat(latitude,dia)
    % 
    % Sonsoles Mateo, Jose María Reinoso, Diego Zapardiel
    %
    disp("latitude")
    disp(latitude)
    
    delta = 23.5*cosd(360*(dia-172)/365);         % Declinación solar
    D = dia-81;                                   % Di
    n = 1.496e13;                                 % Distancia Sol-Tierra media
    d = n*(1-0.017*sind(0.9865.*D));              % Distancia Sol-Tierra dia
    h = 90 - latitude + delta;                    % Angulo solar
    
    if (latitude==-90)
        omega = -90;
        if(h > 180)
            N = 0;
        else 
            N = 24;
        end
    elseif (latitude==90)
        omega = 90;
        if(h < 0)
            N = 0;
        else 
            N = 24;
        end
    else 
        omega = real(90 - acosd(-tand(latitude)*tand(delta)));
        amanecer  = 12 - 1/15*omega;       
        atardecer = 12 + 1/15*omega;
        N         = 12 + amanecer - atardecer; % Horas en la latitud
    end
    
    % Cálculo de la radiación
    S0    = 1367*(n/d)^2;
    alpha = asind(sind(latitude)*sind(delta) + cosd(latitude)*cosd(delta)*cosd(omega));
    S     = S0*sind(alpha);
    Sdia  = S*3600*(sind(h)*2*N/pi);
end

function S=calc_inso_by_lat_ez(latitude,dia)
% This function calculates the value of Si from the mean latitude of the zone
    N=size(latitude,1);
    S=zeros(N,1);
    for i=1:1:N
        %S(i)=1362*(90^2-abs(latitude(i))^2)/90^2;
        S(i)=1362*(-0.5/(90*90)*latitude(i)^2+1);
    end
end

function jacMatrix=jacobian(N,temps)
    %% This function defines the Jacobian matrix for Newton's Method
    global epsilon backRadCoef h_AtmEarthConvec F G
    
    % Dividimos el vector de temperaturas en dos mitades, atmosfera y tierra
    temps_atm=temps(1:N);
    temps_surf=temps(N+1:end);
    
    % Inicializamos la matriz jacobiana con el tamaño adecuado, en cuartos
    Q1=zeros(N,N);
    Q2=Q1;
    Q3=Q1;
    Q4=Q1;
    %% Calculamos las submatrices
    for i=1:1:N % Filas
        for j=1:1:N % Columnas
            if (i==1 && j==1)
                Q1(i,j)=h_AtmEarthConvec-F(i)-4*boltz*temps_atm(i)^3;
                Q2(i,j)=4*boltz*epsilon(i)*temps_surf(i)^3-h_AtmEarthConvec;
                Q3(i,j)=4*boltz*backRadCoef*temps_atm(i)^3-h_AtmEarthConvec;
                Q4(i,j)=h_AtmEarthConvec-G(i)-4*boltz*temps_surf(i)^3;
            elseif (i==N && j==N)
                Q1(i,j)=h_AtmEarthConvec-F(i-1)-4*boltz*temps_atm(i)^3;
                Q2(i,j)=4*boltz*epsilon(i)*temps_surf(i)^3-h_AtmEarthConvec;
                Q3(i,j)=4*boltz*backRadCoef*temps_atm(i)^3-h_AtmEarthConvec;
                Q4(i,j)=h_AtmEarthConvec-G(i-1)-4*boltz*temps_surf(i)^3;
            else
                if i==j % Diagonal
                    Q1(i,j)=h_AtmEarthConvec-F(i)-F(i-1)-4*boltz*temps_atm(i)^3;
                    Q2(i,j)=4*boltz*epsilon(i)*temps_surf(i)^3-h_AtmEarthConvec;
                    Q3(i,j)=4*boltz*backRadCoef*temps_atm(i)^3-h_AtmEarthConvec;
                    Q4(i,j)=h_AtmEarthConvec-G(i)-G(i-1)-4*boltz*temps_surf(i)^3;
                elseif i==(j-1) % Superdiagonal
                    Q1(i,j)=F(i);
                    Q4(i,j)=G(i);
                elseif i==(j+1) % Subdiagonal 
                    Q1(i,j)=F(i-1);
                    Q4(i,j)=G(i-1);
                end
            end
        end
    end
    % Ensamblamos las submatrices 
    
    % | Q1 | Q2 |
    % | Q3 | Q4 |
    
    jacMatrix=[Q1,Q2;Q3,Q4];
end

function funcF=balance(N,temps)
    global S epsilon alpha backRadCoef h_AtmEarthConvec F G

    %This function calculates the vectorial balance equations.
    topHalf=zeros(N,1); %% Ecuaciones balance atmosfera
    botHalf=zeros(N,1); %% Ecuaciones balance tierra

    % Dividimos el vector de temperaturas en dos mitades, atmosfera y tierra
    temps_atm=temps(1:N);
    temps_surf=temps(N+1:end);

    for i=1:1:N
        % i=1, i=N son zonas polares con unicamente un intercambio, de ahi
        % el tratamiento diferenciado. 
        if i==1
            %disp("i equals one")
            topHalf(i)= h_AtmEarthConvec*(temps_atm(i)-temps_surf(i)) ...
                        +epsilon(i)*boltz*temps_surf(i)^4 ...
                        -boltz*temps_atm(i)^4 ...
                        +F(i)*(temps_atm(i+1)-temps_atm(i));

            botHalf(i)= (1-alpha(i))*S(i)/4 ...
                        +backRadCoef*boltz*temps_atm(i)^4 ...
                        -boltz*temps_surf(i)^4 ...
                        +G(i)*(temps_surf(i+1)-temps_surf(i)) ...
                        -h_AtmEarthConvec*(temps_atm(i)+temps_surf(i));
        elseif i==N
            %disp("i equals N")
            topHalf(i)= h_AtmEarthConvec*(temps_atm(i)-temps_surf(i)) ...
                        +epsilon(i)*boltz*temps_surf(i)^4 ...
                        -boltz*temps_atm(i)^4 ... 
                        -F(i-1)*(temps_atm(i)-temps_atm(i-1));

            botHalf(i)= (1-alpha(i))*S(i)/4 ...
                        +backRadCoef*boltz*temps_atm(i)^4 ...
                        -boltz*temps_surf(i)^4 ...
                        -G(i-1)*(temps_surf(i)-temps_surf(i-1)) ...
                        -h_AtmEarthConvec*(temps_atm(i)+temps_surf(i));
        else
            %disp("i equals neither N nor one")
            topHalf(i)= h_AtmEarthConvec*(temps_atm(i)-temps_surf(i)) ...
                        +epsilon(i)*boltz*temps_surf(i)^4 ...
                        -boltz*temps_atm(i)^4 ...
                        +F(i)*(temps_atm(i+1)-temps_atm(i)) ...
                        -F(i-1)*(temps_atm(i)-temps_atm(i-1));

            botHalf(i)=(1-alpha(i))*S(i)/4 ...
                        +backRadCoef*boltz*temps_atm(i)^4 ...
                        -boltz*temps_surf(i)^4 ...
                        +G(i)*(temps_surf(i+1)-temps_surf(i)) ...
                        -G(i-1)*(temps_surf(i)-temps_surf(i-1)) ...
                        -h_AtmEarthConvec*(temps_atm(i)+temps_surf(i));
        end
    end
    % Ensamblamos los vectores
    funcF=[topHalf;botHalf];
end

function [temps_cell_stat, last_iter_stat]=calc_stationary(Temps_init,max_iter)
    global eps delta

    N_zones=length(Temps_init)/2; % Calculo numero de zonas
    temps_cell_stat=cell(max_iter,1); % Inicializamos la cell
    temps_cell_stat{1}=Temps_init; % Distribución inicial 0ºC Constante
    
    %% Newtons method
    for q=2:1:max_iter
        jac=jacobian(N_zones,temps_cell_stat{q-1,1});
        funcF=balance(N_zones,temps_cell_stat{q-1,1});
        deltaTemps=linsolve(jac,-funcF);
        temps_cell_stat{q,1}=temps_cell_stat{q-1,1}+deltaTemps;
        %Convergence control
        funcF_updt=balance(N_zones,temps_cell_stat{q,1});
        if (norm(deltaTemps)/(1+norm(temps_cell_stat{q,1})))<eps && delta>norm(funcF_updt)
            %disp("Convergence achieved in iteration " + str(q))
            break
        end
    end
    last_iter_stat=q; % Save the last iteration calculated (for plotting when convergence is achieved before)
end

function meantemps(temps_cell_stat,lastiter)
    N_zones=length(temps_cell_stat{1,1})/2;
    %Calculate mean temperatures from the cell structure
    mean_temp_atmo=mean(temps_cell_stat{lastiter}(1:N_zones));
    mean_temp_earth=mean(temps_cell_stat{lastiter}(N_zones+1:end));
    absolute_mean=mean(temps_cell_stat{lastiter});
    disp("Stationary results:")
    disp("Mean atmospheric temperature: ")
    disp(mean_temp_atmo)
    disp("Mean earth temperature: ")
    disp(mean_temp_earth)
    disp("Mean temperature: ")
    disp(absolute_mean)
end

function plot_stationary(latitude, temps_cell_stat,varargin)
    N_zones=length(temps_cell_stat{1,1})/2;
    %% Plot final results for stationary case: 
    figure("Name","Temperaturas por zona, Estacionario")
    hold on
    scatter(latitude,temps_cell_stat{varargin{1},1}(1:N_zones,1),"DisplayName","Temperaturas atmósfera","MarkerFaceColor","flat");
    scatter(latitude,temps_cell_stat{varargin{1},1}(N_zones+1:end,1),"DisplayName","Temperaturas terrestres","MarkerFaceColor","flat");
    xlabel("Latitud [º]")
    xlim([-90,90])
    ylabel("Temperatura; [K]")
    legend show
    set(gcf, 'Windowstyle', 'docked')
    
    if nargin==3
    %% First variable length argument enables the plot of the iterations, the second one that must be passed in this case has to be the last iteration calculated
        if varargin{2}==1
            %% Plot iterative calculus for debug
            figure("Name","Temperaturas por zona - Iteraciones")
            subplot(2,1,1)
            hold on
            for k=1:50:varargin{1}
                scatter(latitude,temps_cell_stat{k,1}(1:N_zones,1),"DisplayName","Temp atm. Iteration"+num2str(k),"MarkerFaceColor","flat");
            end
            xlabel("Latitud [º]")
            xlim([-90,90])
            ylabel("Temperatura; [K]")
            legend show
            set(gcf, 'Windowstyle', 'docked')
            subplot(2,1,2)
            hold on
            for k=1:50:varargin{1}
                scatter(latitude,temps_cell_stat{k,1}(N_zones+1:end,1),"DisplayName","Temperaturas terrestres, iteration"+num2str(k),"MarkerFaceColor","flat");
            end
            xlabel("Latitud [º]")
            xlim([-90,90])
            ylabel("Temperatura; [K]")
            legend show
            set(gcf, 'Windowstyle', 'docked')
        end
    end
end

function plot_dynamic(latitude,temps_dyn,iters_used, timesteps)
    N_zones=length(temps_dyn{1,1})/2;
    figure("Name","Temporal Evolution")
    subplot(2,1,1)
    title("Atmosphere")
    xlabel("Latitud [º]")
    xlim([-90,90])
    ylabel("Temperatura [K]")
    hold on
    for l=1:timesteps/5:timesteps
        plot(latitude,temps_dyn{iters_used(l),l}(1:N_zones,1),'DisplayName','timestep='+string(l),'LineWidth',1.5)
    end
    subplot(2,1,2)
    title("Earth")
    xlabel("Latitud [º]")
    xlim([-90,90])
    ylabel("Temperatura [K]")
    hold on
    for l=1:timesteps/5:timesteps
        plot(latitude,temps_dyn{iters_used(l),l}(N_zones+1:end,1),'DisplayName','timestep='+string(l),'LineWidth',1.5)
    end
    legend show
    set(gcf, 'Windowstyle', 'docked')
end

function plot_iterations_dynamics(latitude, temps_dyn, timestep, lastiter)
    N=length(temps_dyn{1,1})/2;
    % For debug, Iteration plots
    figure("Name","Iteration Evolution")
    subplot(2,1,1)
    title("Atmosphere")
    xlabel("Latitud [º]")
    xlim([-90,90])
    ylabel("Temperatura [K]")
    hold on
    for l=1:10:lastiter
        plot(latitude,temps_dyn{l,timestep}(1:N,1),'DisplayName','t='+string(l),'LineWidth',1.5)
    end
    subplot(2,1,2)
    title("Earth")
    xlabel("Latitud [º]")
    xlim([-90,90])
    ylabel("Temperatura [K]")
    hold on
    for l=1:10:lastiter
        plot(latitude,temps_dyn{l,timestep}(N+1:end,1),'DisplayName','t='+string(l),'LineWidth',1.5)
    end
    legend show
    set(gcf, 'Windowstyle', 'docked')
end

%%% BACKUP CODE - I WAS TRYING STUFF %%%
% function latitude=discretize_regular(N)
% %% This function calculates the mean latitude of each zone taking the total number of zones as an input
%     latitude=ones(N,1); %% Latitud media de cada zona
%     latitude(1)=90-180/(2*N);
%     for n=2:1:N
%         latitude(n)=latitude(n-1)-(180/N);
%     end
% end

