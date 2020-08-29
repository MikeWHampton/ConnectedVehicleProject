%For loop method is very slow
%alpha01=0.3; beta01=0.4; gamma01=0.3; beta02=0.3; gamma02 = 0.2; % Designed Parameters]
alpha01Vec = 0;
beta01Vec = 0:0.01:0.05;
beta02Vec = 0:0.01:0.1;
beta03Vec = 0:0.01:0.3;
gamma01Vec = 0:0.01:0.1;
gamma02Vec = 0:.05:0.1;
gamma03Vec = 0:.05:0.2;

energyData = []
counter = 1
numIterations = length(alpha01Vec)*length(beta01Vec)*length(beta02Vec)*length(beta03Vec)*length(gamma01Vec)*length(gamma02Vec)*length(gamma03Vec);
for alpha01=alpha01Vec
    for beta01=beta01Vec
        for beta02 = beta02Vec
            for beta03 = beta03Vec
                for gamma01 = gamma01Vec
                    for gamma02 = gamma02Vec
                        for gamma03 = gamma03Vec
                            energy = EConsumpFunc4([alpha01, beta01, beta02, beta03, gamma01, gamma02, gamma03]);
                            energyData(counter, :) = [alpha01, beta01, beta02, beta03, gamma01, gamma02, gamma03, energy];
                            counter = counter+1;
                            disp(['Percent Complete: ', num2str(counter/(numIterations+1)*100)]);
                        end
                    end
                end
            end
        end
    end
end
[minEnergy, minEnergyIdx] = min(energyData(:,8));
disp(energyData(minEnergyIdx,:))