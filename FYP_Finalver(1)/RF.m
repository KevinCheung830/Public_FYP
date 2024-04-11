%%this file is hybrid precoding
close all
clear all
%Simulation parameter
implementation = 500;


%Channel Parameter
N_s=2;                                 % N_s= 1 or 2
N_T = 64;                              % Number of tx antennas
N_R = 16;                               % Number of rx antennas

Nrf=6;
Nr_RF=6;


Nc  =8; %Cluster in paper fig 3
Np=10; %No of path in cluster

% Antenna separation distance in wavelengths
Fc = 28e9; %freq
lamda= 3e8/ Fc; %length of carrier wave
d=lamda/2;
%d=0.5;

P = 1;                                  % Power of the channel gain
%                                  % number of multi-paths
% theta_t = rand(Np,1)*pi;                 % Angle of departure in radians
% theta_r = rand(Np,1)*pi;                 % Angle of arrival in radians
rho=1;
%%%%%% Constructing the channel matrix%%%%%
                  
nR = (0:(N_R-1))';
nT = (0:(N_T-1))';



BS1 = pi*1/3;
Max_Angle_BaseStation =BS1;
Max_Usr=pi; %range of angle user

snr = 1.0002:0.0001:1.0008;
R_RF = 4:1:10;
%ASSS =7.5/2;%7.5;
ASSS =7.5/2;%7.5;
angleSpread = ASSS * pi / 180;%standard/regularize
coeff = 2 * pi * d / lamda;
total_Np = Nc * Np; %Cluster * Path

K=1; %Channel Generation time
H = zeros(N_R, N_T,K); 


% at_BS=zeros(N_T, Np);
% ar_BS=zeros(N_R, Np);

%Lapacian
mu=0;
b=1;
%TempVar
N_I=min(N_T,N_R);
%%%%%% End of Constructing the channel matrix%%%%%

%%%%%% Trial in 3Feb on file in 31 Jan %%%%%
spectralEffSVD= zeros(implementation,numel(R_RF));
SVD_SqErr = zeros(implementation,numel(R_RF));
spectralEff_BS = zeros(implementation,numel(R_RF));



for index = 1:numel(R_RF)
    
    snr_i = snr(index);
    snr_linear = 10^(snr_i/10);
    sigma2 = 1/snr_linear;
    R_RF_i=R_RF(index);
    Nrf=R_RF_i;
    for impl = 1:implementation
        %%%%%%%%%%%%%%%%%%   The channel matrix %%%%%%%%%%%%%%%%%
         for i = 1 : K
            impl
%             LapcianRan1=mu-b*sign(rand(Nc,1)-0.5).*log(1-2*abs(rand(Nc,1)-0.5));
%             LapcianRan2=mu-b*sign(rand(Nc,Np)-0.5).*log(1-2*abs(rand(Nc,Np)-0.5));
            %%Transmit
            nT = (0 : (N_T - 1))';   
            Lacpcain_phi = (Max_Angle_BaseStation - angleSpread) * rand(Nc, 1) * ones(1, Np) - angleSpread * sign(rand(Nc,Np)-0.5).*log(1-2*abs(rand(Nc,Np)-0.5));
            Lapcian_distri_Angle=Lacpcain_phi;
            Lacpcain_phi = coeff * sin(reshape(Lacpcain_phi',[1, total_Np]));
            at_Temp = nT * Lacpcain_phi;
            H_At = exp(1i * at_Temp) / sqrt(N_T);
            
            % Receive
            nR = (0 : (N_R - 1))';   
            Lacpcain_theta = (Max_Usr - angleSpread) * rand(Nc, 1) * ones(1, Np) - angleSpread * sign(rand(Nc,Np)-0.5).*log(1-2*abs(rand(Nc,Np)-0.5));
            Lacpcain_theta = coeff * (reshape(Lacpcain_theta',[1, total_Np]));
            ar_Temp = nR * Lacpcain_theta;
            H_Ar = exp(1i * ar_Temp) / sqrt(N_R);
   
            alpha = sqrt(N_T * N_R / total_Np/2) * diag(randn(1,total_Np)+randn(1,total_Np)*1j) ;
            H  = H_Ar * alpha * H_At';
        end
        

        %H = (randn(N_R,N_T)+1j*randn(N_R,N_T))/sqrt(2);
        %%%%%%%%%%%%%%%%%%   End of The channel matrix %%%%%%%%%%%%%%%%%

        
        %%%%%%%%%%%%%%%%%%   Optimal Unconstraint SVD %%%%%%%%%%%%%%%%%
        [U,S,V] = svd(H); %H=U*S*V'
        F_opt = V(:,1:N_s);  % V=[V1 V2] which V1 size is  Nt X Ns
        %W_opt = (1/sqrt(snr_linear)) * (F_opt'* H'*H*F_opt+N_s/snr_linear*eye(N_s,N_s))\(F_opt'*H');
        %DecentMat = S;
        W_opt = U';
%       W_opt = W_opt(1:N_s,:);
        DecentMat = S;
     
        I=eye(N_I);
        Temp_SVD = W_opt*H*F_opt;
        Temp_SVD = Temp_SVD*(Temp_SVD');
        R_SVD = W_opt'*W_opt*sigma2;


        P_waterfilling=water_filling(S,P*N_s,N_s);
        P_alloc = zeros(N_R,N_R);
        P_alloc(1:N_s,1:N_s) = P_waterfilling;
        SVDEff = log2(real(det( I + (rho/N_s)* P_alloc*Temp_SVD*inv(R_SVD)))); 
        %SVDEff = log2(
        % real(det( I + (rho/N_s)*Temp_SVD*inv(R_SVD))));     

        spectralEffSVD(impl,index)=SVDEff;

        SVD_SqErr(impl,index) = log2(real(det( I +H*H'/sigma2)));
        %%%%%%%%%%%%%%%%%% End of  Optimal unconstraint Var SVD %%%%%%%%%%%%%%%%%
  


        %%%%%%%%%%%%%%%%%%  Hybrid Precoding %%%%%%%%%%%%%%%%%
        %%Here i  have reconstruct the array response since above H_At is
        %%the product  of array response with lapacian ditributed,i reconstrc a
        %%uniform distri At

        At = zeros(N_T,N_T);
        Ar = zeros(N_R,N_R);    
        index_BS = 0;
        index_BSR = 0;
        for theta_i = -1:2/N_T:1-2/N_T
            index_BS = index_BS+1;
            At(index_BS,:) = exp(1j*pi.*nT*(theta_i));
        end

        for theta_i = -1:2/N_R:1-2/N_R
            index_BSR = index_BSR+1;
            Ar(index_BSR,:) = exp(1j*pi.*nR*(theta_i));
        end
        At=At';
        Ar=Ar';

        % At=H_At;
        % Ar=H_Ar;

        %         [NT,~]=size(At); %initialisation with size of array response of transmission
        %         F_RF = zeros(N_T,Nrf); %initialisation with NT X NRF
        waterfill=P_alloc(N_s,N_R,:);
        F_RF=[]; %Paper Step1 Empty Matrix for F_RF
        F_res=F_opt*waterfill;
        %F_res=F_opt;  %Paper Step2 declare Residual Precoding matrix F_res
        for i=1:Nrf  %Paper Step 3 i<=Nrf_t
            psi=At'*F_res; %Paper Step 4 Psi=At
            [~,k]=max(diag(psi*psi')); %Paper Step 5,arg maxof(psipsi*)Np,Np diagonally
            F_RF=[F_RF At(:,k)];      % Paper Step 6(Not sure how to append here,Use NP.hstack? or can just do [A B] in for lp,Append the selected vector in step 5 aboveline  to F_RF
            F_BB=((inv(F_RF'*F_RF))*F_RF')*F_opt; %Paper step 7
            F_res=(F_opt-F_RF*F_BB)/norm(F_opt-F_RF*F_BB,'fro'); %Paper Step8,Cal F_res,with Frobenius norm,aka Euclidean norm, aim to measure the magnitude of a matrix
            %             norm(F_RF*F_BB-F_opt,'fro')
        end  % Paper Step 9 ,end of For loop
        F_BB = sqrt(N_s)*F_BB/norm(F_RF*F_BB,'fro');   %Paper Step 10 Normalised

        Eyy=(rho/N_s)*H*F_RF*F_BB*F_BB'*F_RF'*H'+(sigma2)*eye(N_R,N_R);
        TemoWMMSE=F_BB'*F_RF'*H'*H*F_RF*F_BB+ ((sigma2)*N_s/rho);
        W_MMSEH=(1/sqrt(rho))*inv(TemoWMMSE)*F_BB'*F_RF'*H';
        W_MMSE = W_MMSEH';
        %         [N_R,~] = size(Ar); %initialisation with size of array response of receive
        %         W_RF = zeros(N_R,Nr_RF);  %initialisation with NR X Nr_RF
        W_RF=[];
        W_res = W_MMSE; % Paper Step 2
        for i = 1 : Nr_RF %Paper Step 3 i <= Nr_RF
            psi_r = Ar'*Eyy*W_res; %Paper Step 4 find the contribution???
            [~,k] =  max(diag(psi_r*psi_r'));  %Paper Step 5,arg max of(psi_r psi_r*)Np,Np diagonally
            W_RF = [W_RF Ar(:,k)]; %Paper Steps 6 append the selected array response to W_RF
            W_BB = (inv(W_RF'*Eyy*W_RF))*W_RF'*Eyy*W_MMSE; %Paper step 7
            W_res = (W_MMSE-W_RF*W_BB)/norm(W_MMSE-W_RF*W_BB,'fro');      %Paper Step 8 Elimated the selected atom's contribution to the  res
            %             norm(W_RF*W_BB-W_MMSE,'fro')
        end
        W_BB = sqrt(N_s)*W_BB/norm(W_RF*W_BB,'fro');

        %         SNR = rho / (sigma2); % in paper page 1509

        R_n = W_BB'*W_RF'*W_RF*W_BB*sigma2;  %R_noise is the noise covariance matrix suggested in the paper expression (3) in page 1501

        temp_OMP = (rho/N_s)* (inv(R_n)) * W_BB' *W_RF' * H * F_RF * F_BB* F_BB' * F_RF' * H' * W_RF * W_BB ; % for clear representation ,Cross product in the expression (3)

        Temp_SpectralEfficiency_OMP = log2(abs(det( eye(N_s,N_s) + temp_OMP)));  %Paper 1501 expression(3)
        spectralEff_OMP(impl,index)=Temp_SpectralEfficiency_OMP;
        %%%%%%%%%%%%%%%%%%  End Of Hybrid Precoding %%%%%%%%%%%%%%%%%



        %%%%%%%%%%%%%%%%%% Beam Steering  %%%%%%%%%%%%%%%%%
        %%find the  best path for data stream and add to dicitionary
        W_opt = W_opt(1:N_s,:);
        F_BS=zeros(N_T,Nrf);
        W_BS=zeros(Nrf,N_R);
        Peaksave=[];
        grid_ = Ar'*H*At;
        grid_Save = Ar'*H*At;
        for iiii = 1:Nrf
            
            [value_,index_] = max(abs(grid_),[],'all');
            [row_,col_] = ind2sub(size(grid_),index_);
           
            W_BS(iiii,:) = Ar(:,row_)';
            F_BS(:,iiii) = At(:,col_);
            %Peaksave(iiii,:)=[row_ col_ value_];
            %Peaksave
            %abs(Ar(:,row_)'*H*At(:,col_))-value_
            %Prevent reselection of the Peak which already add to W_BS and F_BS
            
            grid_(row_,:) = 0;
            grid_(:,col_) = 0;
            
            %abs(Ar(:,row_)'*H*At(:,col_))-value_

        end
% 
%         [indexArrResponseTx,indexArrResponseRx] =findSteeringVector(H,at_BS,ar_BS,N_s);
%         for n=1:N_s
%             F_BS=[F_BS At(:,indexArrResponseTx)];   
%             W_BS=[W_BS Ar(:,indexArrResponseRx)];
%         end

        F_BB_bs = inv(F_BS'*F_BS)*F_BS'*F_opt;
        F_BB_bs = sqrt(N_s)*F_BB_bs/norm(F_BS*F_BB_bs,'fro');

        W_BB_bs = W_opt*W_BS'*inv(W_BS*W_BS');
        W_BB_bs = sqrt(N_s)*W_BB_bs/norm(W_BS'*W_BB_bs','fro');
        
        RnBS = (W_BB_bs*W_BS)*(W_BS'*W_BB_bs')*sigma2;
        
        spect_BS =log2(abs(det(  eye(N_s,N_s) +  (rho/N_s) * (RnBS\(W_BB_bs*W_BS*H*F_BS*F_BB_bs)*(F_BB_bs'*F_BS'*H'*W_BS'*W_BB_bs') ))));
        spectralEff_BS(impl,index)=spect_BS;

        %%%%%%%%%%%%%%%%%% Beam Steering  %%%%%%%%%%%%%%%%%

    end

end

%SE_ave = SE/implementation;

spectralEff_OMP=mean(spectralEff_OMP,1);
spectralEffSVD=mean(spectralEffSVD,1);
spectralEff_BS=mean(spectralEff_BS,1);

%SVD_SqErr = mean(SVD_SqErr,1);

figure
hold on
l1 = plot(R_RF,spectralEffSVD,'-o');; %SVD with blue dotted x
l2 = plot(R_RF,spectralEff_OMP,'-*');    %Hybrid Precoding with  orange line
l3 = plot(R_RF,spectralEff_BS,'c'); 
title('Spectral Efficiency achieved by various precoding solutions')
xlabel('Number of radio Frequency chain used') 
ylabel('Spectrual Efficiency (bits/s)/Hz ') 
legend('Optimal Unconstrained Precoding','Orthogonal Match Pursuit','Beam Steering')

figure
surf(abs(grid_))
title('Spectral Efficiency achieved by various precoding solutions')

figure
surf(abs(DecentMat))
title('Opt Precoding SVD-Σ Matrix ')
xlabel('No of transmit Antenna') 
ylabel('No of Receive Antenna') 

figure
surf(abs(V*V'))
title('Opt Precoding SVD-V Matrix ')
xlabel('No of transmit Antenna') 
ylabel('No of transmit Antenna') 

figure
surf(abs(U*U'))
title('Opt Precoding SVD-U Matrix ')
xlabel('No of Receive Antenna') 
ylabel('No of Receive Antenna') 

figure
surf(abs(H))
title('H Channel Matrix ')
xlabel('Total Number of paths(No  of Cluster X MultiPath)') 
ylabel('No of Receive Antenna') 
zlabel('Peak_Gain') 

figure;
histogram(Lapcian_distri_Angle)
title('CDF of AOE')
% xlabel('No of transmit Antenna') 
% ylabel('No of Receive Antenna') 

