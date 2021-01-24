% revise_distance_metric.m
% Nathanael Kazmierczak, 03/22/2020

% Nondimensionalized fitting quantities
trig_prefactors = ones(1,24);
displacement_vector = [0,1.4];
R12_true = extendedTrigFittingFunctions(displacement_vector,trig_prefactors,0);
% R12_true = trig_prefactors;  % implicit definition of AA; work on this later.

% Set the neighborhood distance we are interested in probing
epsilon = 0.1;

% randomly sample N points in R12. This construction will always include
% the entirety of each neighborhood because the uniform distribution gives
% hypercubes, which is equivalent to L_infinity norm.
N = 1e3;
R12_all = zeros(N,12);
useRand = false;
useLatin = false;
if useRand
    for i = 1:12
        R12_all(:,i) = unifrnd(R12_true(i)-epsilon,R12_true(i)+epsilon,N,1);
    end
elseif useLatin
    R12_all = lhsdesign(N,12);
end


%% Try the mapping
% Build a scalar objective function for each norm
% weight_vector_local = ones(1,12);
weight_vector_unweighted = [1*ones(1,6),1*ones(1,6),1*ones(1,6),0*ones(1,6)];
weight_vector_local = [1*ones(1,6),0*ones(1,6),0*ones(1,6),0*ones(1,6)];
weight_vector_local2 = [0*ones(1,6),1*ones(1,6),0*ones(1,6),0*ones(1,6)];
weight_vector_local3 = [0*ones(1,6),0*ones(1,6),1*ones(1,6),0*ones(1,6)];
weight_vector_local4 = [0*ones(1,6),0*ones(1,6),0*ones(1,6),1*ones(1,6)];
objfcn_generic = @(DSCvector,k,weight_vector) (sum(((abs(R12_true - extendedTrigFittingFunctions(DSCvector,trig_prefactors,0))).^k).*weight_vector)).^(1/k);
% objfcn_generic_unweighted = @(DSCvector,k) (sum(((abs(R12_true - extendedTrigFittingFunctions(DSCvector,trig_prefactors,0))).^k))).^(1/k);
% objfcn_L0p01 = @(DSCvector) objfcn_generic(DSCvector,0.01);
% objfcn_L1 = @(DSCvector) objfcn_generic(DSCvector,1);
objfcn_L2_unweighted = @(DSCvector) objfcn_generic(DSCvector,2,weight_vector_unweighted);
objfcn_L2_firstring = @(DSCvector) objfcn_generic(DSCvector,2,weight_vector_local);
objfcn_L2_secondring = @(DSCvector) objfcn_generic(DSCvector,2,weight_vector_local2);
objfcn_L2_thirdring = @(DSCvector) objfcn_generic(DSCvector,2,weight_vector_local3);
objfcn_L2_fourthring = @(DSCvector) objfcn_generic(DSCvector,2,weight_vector_local4);

% objfcn_L3 = @(DSCvector) objfcn_generic(DSCvector,3);
% objfcn_L10 = @(DSCvector) objfcn_generic(DSCvector,10);
% objfcn_L2_wpenalty = @(DSCvector) objfcn_generic(DSCvector,2) + max(abs(R12_true - trigFittingFunctions(DSCvector,trig_prefactors,0)));
% objfcn_L2_wpenalty3 = @(DSCvector) objfcn_generic(DSCvector,2) + 3*max(abs(R12_true - trigFittingFunctions(DSCvector,trig_prefactors,0)));
% objfcn_L2_w2penalty = @(DSCvector) objfcn_generic(DSCvector,2) + getTopN(abs(R12_true - trigFittingFunctions(DSCvector,trig_prefactors,0)),4)*[1;1;1;4];

xbase = -1.5:0.05:1.5;
ybase = 0:0.05:1.5;
[xspace,yspace] = meshgrid(xbase,ybase);
% RMSR_L0p01 = zeros(size(xspace));
% RMSR_L1 = zeros(size(xspace));
RMSR_L2_ring1 = zeros(size(xspace));
RMSR_L2_ring2 = zeros(size(xspace));
RMSR_L2_ring3 = zeros(size(xspace));
RMSR_L2_ring4 = zeros(size(xspace));
RMSR_L2_unweighted = zeros(size(xspace));
% RMSR_L3 = zeros(size(xspace));
% RMSR_L10 = zeros(size(xspace));
% RMSR_L2_wpenalty = zeros(size(xspace));
% RMSR_L2_wpenalty3 = zeros(size(xspace));
% RMSR_L2_w2penalty = zeros(size(xspace));
for i = 1:numel(ybase)
    i
    for j = 1:numel(xbase)
        this_DSC = [xspace(i,j),yspace(i,j)];
%         result_L0p01 = objfcn_L0p01(this_DSC);
%         result_L1 = objfcn_L1(this_DSC);
        RMSR_L2_ring1(i,j) = objfcn_L2_firstring(this_DSC);
        RMSR_L2_ring2(i,j) = objfcn_L2_secondring(this_DSC);
        RMSR_L2_ring3(i,j) = objfcn_L2_thirdring(this_DSC);
        RMSR_L2_ring4(i,j) = objfcn_L2_fourthring(this_DSC);
        RMSR_L2_unweighted(i,j) = objfcn_L2_unweighted(this_DSC);
%         result_L3 = objfcn_L3(this_DSC);
%         result_L10 = objfcn_L10(this_DSC);
%         result_L2_wpenalty = objfcn_L2_wpenalty(this_DSC);
%         result_L2_wpenalty3 = objfcn_L2_wpenalty3(this_DSC);
%         result_L2_w2penalty = objfcn_L2_w2penalty(this_DSC);
% %         % Norms in R2, displacement space
%         RMSR_L0p01(i,j) = result_L0p01;
%         RMSR_L1(i,j) = result_L1;
%          = result_L2;
%          = result_L2_unweighted;
%         RMSR_L3(i,j) = result_L3;
%         RMSR_L10(i,j) = result_L10;
%         RMSR_L2_wpenalty(i,j) = result_L2_wpenalty;
%         RMSR_L2_wpenalty3(i,j) = result_L2_wpenalty3;
%         RMSR_L2_w2penalty(i,j) = result_L2_w2penalty;
    end
end

%% Plot the distance results.
% % figure
% % contourf(xbase,ybase,RMSR_L0p01/(max(max(RMSR_L0p01))),30);
% % title('L0.01 Norm distance from AA','FontSize',12)
% % axis equal
% % colormap jet
% % colorbar
% % % caxis([0.5,3]);
% % 
% % figure
% % contourf(xbase,ybase,RMSR_L1/(max(max(RMSR_L1))),30);
% % title('L1 Norm distance from AA','FontSize',12)
% % colormap jet
% % axis equal 
% % colorbar
% % % caxis([0.5,3]);
figure
subplot(2,2,1);
contourf(xbase,ybase,RMSR_L2_unweighted/(max(max(RMSR_L2_unweighted))),20);
title('L2 Norm unweighted distance from AA','FontSize',12)
colormap jet
axis equal 
colorbar

subplot(2,2,2);
contourf(xbase,ybase,RMSR_L2_ring1/(max(max(RMSR_L2_ring1))),20);
title('L2 Norm first ring distance from AA','FontSize',12)
colormap jet
axis equal 
colorbar

subplot(2,2,3);
contourf(xbase,ybase,RMSR_L2_ring2/(max(max(RMSR_L2_ring2))),20);
title('L2 Norm second ring distance from AA','FontSize',12)
colormap jet
axis equal 
colorbar

subplot(2,2,4);
contourf(xbase,ybase,RMSR_L2_ring3/(max(max(RMSR_L2_ring3))),20);
title('L2 Norm third ring distance from AA','FontSize',12)
colormap jet
axis equal 
colorbar

set(gcf,'Position',[0,0,1200,600]);

figure 
contourf(xbase,ybase,RMSR_L2_ring4/(max(max(RMSR_L2_ring4))),15);
title('L2 Norm fourth ring distance from AA','FontSize',12)
colormap jet
axis equal 
colorbar


% % figure
% % contourf(xbase,ybase,RMSR_L3/(max(max(RMSR_L3))),30);
% % title('L3 Norm distance from AA','FontSize',12)
% % colormap jet
% % axis equal 
% % colorbar
% % 
% % figure
% % contourf(xbase,ybase,RMSR_L10/(max(max(RMSR_L10))),30);
% % title('L10 Norm distance from AA','FontSize',12)
% % colormap jet
% % axis equal 
% % colorbar
% % 
% % figure
% % contourf(xbase,ybase,RMSR_L2_wpenalty/(max(max(RMSR_L2_wpenalty))),30);
% % title('L2 Norm with penalty distance from AA','FontSize',12)
% % colormap jet
% % axis equal 
% % colorbar
% % 
% % figure
% % contourf(xbase,ybase,RMSR_L2_wpenalty3/(max(max(RMSR_L2_wpenalty3))),30);
% % title('L2 Norm with triple weight penalty distance from AA','FontSize',12)
% % colormap jet
% % axis equal 
% % colorbar
% % 
% % figure
% % contourf(xbase,ybase,RMSR_L2_w2penalty/(max(max(RMSR_L2_w2penalty))),30);
% % title('L2 Norm with doubled penalty distance from AA','FontSize',12)
% % colormap jet
% % axis equal 
% % colorbar

% % Define the neighborhoods
% radii_function = @(k) sum((abs(R12_all - repmat(R12_true,N,1))).^k,2).^(1/k);
% L2_radii = radii_function(2);
% R12_L2 = R12_all(L2_radii < epsilon,:);
% size(R12_L2)
% 
% L1_radii = radii_function(1);
% R12_L1 = R12_all(L1_radii < epsilon,:);
% size(R12_L1)
% 
% L3_radii = radii_function(1);
% R12_L3 = R12_all(L3_radii < epsilon,:);
% size(R12_L3)

