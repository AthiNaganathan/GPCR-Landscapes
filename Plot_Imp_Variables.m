% Plotting differential scanning calorimetry profiles, free energy profiles,
% residue folding probabilities as a function of the reaction coordinate,
% and residue-residue coupling matrices.

%% Choosing a GPCR
gpcr = 'gpcr1'; % Enter the GPCR for which plots are to be generated
state = 'i';    % Enter the required state (i: inactive, a: active)

%% Loading the data
eval(['load ',gpcr,state,'.mat;']);
eval(['Cpd=Cpd_',gpcr,state,';']);
eval(['fes=fes_310_',gpcr,state,';']);
eval(['fes2D=fes2D_310_',gpcr,state,';']);
eval(['Fpath=Fpath_310_',gpcr,state,';']);
eval(['CouplingMat=CouplingMat_310_',gpcr,state,';']);

%% Plotting the differential scanning calorimetry profile
figure;
plot(273:373, Cpd, 'k-', 'LineWidth', 2);
title('Differential Scanning Calorimetry Profile');
xlabel('Temperature (K)')
ylabel('C_p (kJ mol^{-1} K^{-1})')
xlim([273 373])
set(gca,'LineWidth',2);

%% Plotting the 1D free energy profile
figure;
plot(linspace(0,1,length(fes)), fes, 'k-', 'LineWidth', 2);
title('1D Free Energy Profile')
xlabel('Fraction of Sructured Blocks')
ylabel('FE (kJ mol^{-1})')
set(gca,'LineWidth',2);

%% Plotting the 2D free energy profile
figure;
surf(fes2D);
view(150,50)
colormap jet;
shading interp;
title('2D Free Energy Profile')
xlabel('n_C')
ylabel('n_N')
zlabel('FE (kJ mol^{-1})')

%% Plotting the residue folding probabilities as a function of the reaction coordinate
figure;
pcolor(linspace(0,1,length(fes)), 1:size(Fpath,1), Fpath);
colormap jet;
shading interp;
hcb1 = colorbar;
colorTitleHandle1 = get(hcb1,'Title');
set(colorTitleHandle1 ,'String','Probability');
xlabel('Fraction of Sructured Blocks')
ylabel('Residue Number')
title('Residue Folding Probabilities')

%% Plotting the residue-residue coupling matrix
figure;
pcolor(CouplingMat);
colormap jet;
shading interp;
hcb2 = colorbar;
colorTitleHandle2 = get(hcb2,'Title');
set(colorTitleHandle2 ,'String','FE (kJ mol^{-1})');
xlabel('Residue Number')
ylabel('Residue Number')
title('Coupling Matrix')
