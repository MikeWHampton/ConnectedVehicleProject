clc; clearvars;
%For loop method is very slow
alpha01Vec=0.034;
beta01Vec=0.068;
beta02Vec=0.092;
beta03Vec=0.085;
beta04Vec=0;
beta05Vec = 0;
beta06Vec = 0;
gamma01Vec = 0;
gamma02Vec = 0;
gamma03Vec = 0.2; 
gamma04Vec = 0:0.001:0.2;
gamma05Vec = 0.05;
gamma06Vec = 0.24;

counter = 1;
numIterations = length(alpha01Vec)*length(beta01Vec)*length(beta02Vec)*length(beta03Vec)...
    *length(beta04Vec)*length(beta05Vec)*length(beta06Vec)...
    *length(gamma01Vec)*length(gamma02Vec)*length(gamma03Vec)...
    *length(gamma04Vec)*length(gamma05Vec)*length(gamma06Vec)
energyData = zeros(numIterations,14);
for alpha01=alpha01Vec
    for beta01=beta01Vec
        for beta02=beta02Vec
            for beta03=beta03Vec
                for beta04 = beta04Vec
                    for beta05 = beta05Vec
                        for beta06 = beta06Vec
                            for gamma01 = gamma01Vec
                                for gamma02 = gamma02Vec
                                    for gamma03 = gamma03Vec
                                        for gamma04 = gamma04Vec
                                            for gamma05 = gamma05Vec
                                                for gamma06 = gamma06Vec
                                                    energy = EConsumpFunc8Try2([alpha01, beta01, beta02, beta03, beta04, beta05, beta06,...
                                                        gamma01, gamma02, gamma03, gamma04, gamma05, gamma06]);
                                                    energyData(counter, :) = [alpha01, beta01, beta02, beta03, beta04, beta05, beta06,...
                                                        gamma01, gamma02, gamma03, gamma04, gamma05, gamma06, energy];
                                                    counter = counter+1;
                                                    disp(['Percent Complete: ', num2str(counter/(numIterations+1)*100)])
                                                    
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
[minEnergy, minEnergyIdx] = min(energyData(:,14));
bestCase = (energyData(minEnergyIdx,:))