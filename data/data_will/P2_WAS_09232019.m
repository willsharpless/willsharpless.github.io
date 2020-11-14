
%% Simulation of autonomous UAV's - E150 Project2 - Will Sharpless
% September 24, 2019

%% Define Constants
%Dynamics and Integration Params
A=1; C=0.25; m=10; v_a=zeros(1,3); rho_a=1.225; dt=0.2; t_f=60; F_pi = 200;

%Object & Interaction Params
agent_sight=5; crash_range=2;
N_m=15; %15;
N_T=100; %100
N_o=25; %25

%Genetic Algorithm Params
parents=6; %6
S=20; %20
G=100; %100
w1=70; w2=10; w3=20;
M_star = zeros(G,S);
T_star = zeros(G,S);
L_star = zeros(G,S);
Pi = zeros(G,S);
Orig = zeros(G,S);

%Initialize Target and Obstacle locations (uniform random)
Txy = 200*rand(2,N_T) - 100 ;
Tz = 20*rand(1,N_T) - 10;
T = [Txy; Tz];
Oxy = 200*rand(2,N_o) - 100; Oz = 20*rand(1,N_o) - 10;
O = [Oxy; Oz];

%Initialize Swarm Members and Positions in 3 planes x=-110,-130,-150 in X formation (min d = 14)
rx = [-110*ones(1,5) -130*ones(1,5) -150*ones(1,5)];
ry = [10 10 0 -10 -10 10 10 0 -10 -10 10 10 0 -10 -10];
rz = [10 -10 0 10 -10 10 -10 0 10 -10 10 -10 0 10 -10];

r_0=[rx; ry; rz]; v_0=zeros(3,N_m);
N_t = t_f/dt;
r_0(:,:,2:N_t)=inf; %3D mat to plot paths over time

%Generate random lambda_i variables (S x 15 array)
lambda = 2*rand(S,15); 
score = nan(S,1); %for economic computation
lambda = [lambda score];

%Activity Matricies
T0_active = ones(N_T,1); m0_active = ones(N_m,1);

%plotting or debugging
flag=0; %flag variable for breaking loops

%% Iterations

for g = 1:G %generation loop
  if g==1
      fprintf('Generation 1 \n')
  else
      fprintf('Generation %d, best score: %.2f \n',g,Pi(g-1,1))
  end
  
  
  for s = 1:S %String simulation loop
    
    if ~(isnan(lambda(s,16))) %string already simulated
        continue
    end
    
    fprintf('    Simulating string %d \n',s)
    
    crash_loc = [200; 200; 200]; map_loc = [200; 200; 200]; %for plotting crashes and target mappings
    
    [Wmt, Wmo, Wmm, wt1, wt2, wo1, wo2, wm1, wm2, a1, a2, b1, b2, c1, c2] = ...
    deal(lambda(s,1),lambda(s,2),lambda(s,3),lambda(s,4),lambda(s,5),lambda(s,6), ...
    lambda(s,7),lambda(s,8),lambda(s,9),lambda(s,10),lambda(s,11),lambda(s,12),lambda(s,13),lambda(s,14),lambda(s,15)); %Parse lambda
    

    %Reset activity states of targets and drones
    T_active = T0_active; m_active = m0_active;
    r = r_0; v = v_0;
  
    t=1;
    while t<=N_t %Time loop                    
        
        for i=1:N_m  %Agent loop
            if m_active(i) == 0 %self is dead agent
                continue
            end
            
            N_mt = 0; N_mm = 0; N_mo = 0;
                    
            for j=1:N_T
                if T_active(j) == 0 %mapped target
                    continue
                end
                
                d_mt_ij = norm(T(:,j)-r(:,i,t));
                
                if d_mt_ij < agent_sight %target mapped
%                    fprintf('       Target mapped! \n')
                   map_loc = [map_loc T(:,j)];
                   T_active(j) = 0;
                else
                    n_mt_ij = (T(:,j)-r(:,i,t)) ./ d_mt_ij;
                    n_mt_ij_hat = (wt1*exp(-a1*d_mt_ij)-wt2*exp(-a2*d_mt_ij))*n_mt_ij;
                    N_mt = N_mt + n_mt_ij_hat;
                end
            end %target loop
            
            for j=1:N_m %agent-agent interactions
                if m_active(j) == 0 %other is dead agent
                    continue
                end
                if i==j %dont compare w self
                    continue
                end
                
                d_mm_ij = norm(r(:,j,t)-r(:,i,t));
                
                if d_mm_ij < crash_range %drones crash
%                    fprintf('       two agents collided \n')
                   crash_loc = [crash_loc r(:,i,t)];
                   m_active(j) = 0;
                   m_active(i) = 0;
                   flag=1;
                   break
                else
                    n_mm_ij = (r(:,j,t)-r(:,i,t)) ./ d_mm_ij;
                    n_mm_ij_hat = (wm1*exp(-b1*d_mm_ij)-wm2*exp(-b2*d_mm_ij))*n_mm_ij;
                    N_mm = N_mm + n_mm_ij_hat;
                end
            end %agent loop
            
            if flag==1 %agent crashed
                flag=0;
                continue
            end
            
            for j=1:N_o %loop over obstacles
                d_mo_ij = norm(O(:,j)-r(:,i,t));
                
                if d_mo_ij < crash_range %agent crashed into obstacle
%                    fprintf('       agent crashed \n')
                   crash_loc = [crash_loc r(:,i,t)];
                   m_active(i) = 0;
                   flag=1;
                   break
                else
                    n_mo_ij = (O(:,j)-r(:,i,t)) ./ d_mo_ij;
                    n_mo_ij_hat = (wo1*exp(-c1*d_mo_ij)-wo2*exp(-c2*d_mo_ij))*n_mo_ij;
                    N_mo = N_mo + n_mo_ij_hat;
                end
            end %obstacle loops
            
            if flag==1 %agent crashed
                flag=0;
                continue
            end
            
            %calculating Psi
            N_tot = Wmt*N_mt + Wmm*N_mm + Wmo*N_mo;
            n_prop = N_tot ./ norm(N_tot);
            F_p = F_pi*n_prop;
            F_d = 0.5*rho_a*C*A*norm(-v(:,i))*(-v(:,i));
            Psi = F_p + F_d;
            
            %updating v and r
            v(:,i) = v(:,i) + (dt/m)*Psi;
            
            if t+1<N_t
                
                r(:,i,t+1) = r(:,i,t) + dt*v(:,i);
                
                if abs(r(1,i,t+1))>150 || abs(r(2,i,t+1))>150 || abs(r(3,i,t+1))>60 %out of bounds = death
%                     fprintf('       agent flew away \n')
                    crash_loc = [crash_loc r(:,i,t+1)];
                    m_active(i) = 0;
                continue
                end
            end
        
        end %agent loop (iteration for each agent)
        
        t = t + 1;
        t_stop = (t-1)*dt; %for T_star calculation
        
        if sum(T_active) == 0 || sum(m_active) == 0 %end sim if targets or agents done
            break
        end
    end %time loop (single simulation of all agents)
    
    %store individual arrays
    M_star(g,s) = sum(T_active)/N_T;
    T_star(g,s) = t_stop/t_f;
    L_star(g,s) = (N_m-sum(m_active))/N_m;
   
    %compute Pi(lambda, L*, M*, T*, store individual arrays
    lambda(s,16) = w1*M_star(g,s) + w2*T_star(g,s) + w3*L_star(g,s);
       
  end %string simulation scoring
 
  [Pi(g,:), Orig(g,:)] = sort(lambda(:,16)'); %sorting the scores
  lambda = lambda(Orig(g,:),:); %resorting lambda
  
  for k = 1:2:parents
    phi1 = rand(1,15);
    phi2 = rand(1,15);
    child_1 = (lambda(k,1:15).*phi1) + (lambda(k+1,1:15).*(1-phi1));
    child_2 = (lambda(k,1:15).*(1-phi2)) + (lambda(k+1,1:15).*phi2);
    lambda(parents+k,:) = [child_1 NaN];
    lambda(parents+k+1,:) = [child_2 NaN];
  end
 
  lambda((2*parents+1):end,:) = [2*rand(S-2*parents, 15) nan(S-2*parents,1)];
  %now lambda = [parents; children; new guesses]
  
  %storing parent values of M*, L*, S* that will not be simualted
  if g~=G
      for l = 1:6
        M_star(g+1,1:6) = M_star(g,Orig(g,l));
        T_star(g+1,1:6) = T_star(g,Orig(g,l));
        L_star(g+1,1:6) = L_star(g,Orig(g,l));
      end
  end
 
end %generation

%% Plotting Convergence of Cost

% Pi Convergence
  figure
    semilogy(1:size(Pi,1),Pi(:,1));
    hold on
    semilogy(1:size(Pi,1),mean(Pi(:,1:parents),2));
    semilogy(1:size(Pi,1),mean(Pi(:,:),2));
    title('Cost per Generation','FontSize',17)
    xlabel('Generation')
    ylabel('Pi')
    legend('Best','Mean Parents','Mean Population','FontSize',14)
    grid on
    hold off
savefig('Cost per Generation');
%% Plotting Convergence of M_star, T_star, L_star 

 % M_star Convergence
  figure
    plot(1:size(M_star,1),M_star(:,1),'*-');
    hold on
    plot(1:size(M_star,1),mean(M_star(:,1:parents),2),'o-');
    plot(1:size(M_star,1),mean(M_star(:,:),2));
    title('Fraction of Targets Missed per Generation','FontSize',17)
    xlabel('Generation')
    ylabel('M star')
    legend('Best','Mean Parents','Mean Population','FontSize',14)
    grid on
    hold off

  % T_star Convergence
  figure
    plot(1:size(T_star,1),T_star(:,1),'*-');
    hold on
    plot(1:size(T_star,1),mean(T_star(:,1:parents),2),'o-');
    plot(1:size(T_star,1),mean(T_star(:,:),2));
    title('Fraction of Time Spent per Generation','FontSize',17)
    xlabel('Generation')
    ylabel('T star')
    legend('Best','Mean Parents','Mean Population','FontSize',14)
    grid on
    hold off

    
  % L_star Convergence
  figure
    plot(1:size(L_star,1),L_star(:,1),'*-');
    hold on
    plot(1:size(L_star,1),mean(L_star(:,1:parents),2),'o-');
    plot(1:size(L_star,1),mean(L_star(:,:),2));
    title('Fraction of Crashed Agents per Generation','FontSize',17)
    xlabel('Generation')
    ylabel('L star')
    legend('Best','Mean Parents','Mean Population','FontSize',14)
    grid on
    hold off

save('info','Orig', 'lambda', 'Pi')

