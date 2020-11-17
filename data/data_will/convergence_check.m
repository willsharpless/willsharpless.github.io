% Checking the desolve convergence of several communities for Kyle and Skylar

traj = readtable('/Users/willsharpless/Downloads/PSC_1305_desolve_timecourse_smaller.csv');
params = readtable('/Users/willsharpless/Downloads/25mar20_KillerComm_1305_stode_final_parameters.csv');

%%

figure('color','w') %,'Units','normalized'
hold on
for i = 1:100
    
if i == 29 || i ==81  % somethin fcked abt these sims (contain 'N/A')
    continue
end

traji = cellfun(@str2num, traj{traj.comm == i,3:end});
mui = params{i,1:11}'; Ai = reshape(params{1,12:end},11,11)';
fpi = Ai\(-mui);
dist = sqrt(sum((traji - fpi').^2,2));

if max(max(dist)) > 2
    strcat({'Simulation '},num2str(i),{' has norm > 2'})
end

plot(traj(traj.comm == 1,:).time,dist);

end

title('Convergence of Community Trajectories')
ylim([0 0.5])
xlabel('Time (hours)','FontSize',10)
ylabel('||state - equilibrium||')