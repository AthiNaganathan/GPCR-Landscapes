% Single Sequence Approximation

% clear; clc;
% tic;

%% Looping PDB files

for run=1
%     clearvars -except run pdb ene;
    pdb=char('pdbfilename'); %modify to the name of the relevant PDB file
    
    %% Contact Maps obtained through PDB
    
    aa=pdb(run,:);
    eval(['load contactmapmatElecB',aa,'.dat;']);
    eval(['cmapmask=contactmapmatElecB',aa,';']);
    eval(['load contdistElecB',aa,'.dat;']);
    eval(['contdist=contdistElecB',aa,';']);
    eval(['load BlockDet',aa,'.dat;']);
    eval(['BlockDet=BlockDet',aa,';']);
    nres=length(cmapmask);
    nres1 = length(unique(BlockDet(:,1)));
    
    Mw=nres1*110;
    
    %% Variables
    
    ene=-48.2/1000; %modify to the required vdW interaction energy value
    DS=-10/1000;
    DCp=-0.3579/1000;
    DDS=6.0606;
    
    int = 1.7;
    slope = 6.7/1000;
    
    %WT
    disr = [];
    ppos=[];
    hh9=fopen('struct.txt','rt'); %secondary structure prediction file generated using STRIDE
    lkk=fgetl(hh9);
    while lkk > 0
        if length(lkk)> 25 && strcmp(lkk(1:3),'ASG')
            if strcmp(lkk(6:8),'PRO')
                ppos(end+1)=str2double(lkk(13:15));
            elseif strcmp(lkk(25),'H')==0 && strcmp(lkk(25),'E')==0 && strcmp(lkk(25),'G')==0
                disr(end+1)=str2double(lkk(13:15));
            end
        end
        lkk=fgetl(hh9);
    end
    fclose(hh9);
    ppos2= [];
    disr2 = []';
    
    T=[310]';
    %T=263;
    IS=0.1;
    C=0;
    
    %% Constants
    R=0.008314;
    Tref = 385;
    
    zval=exp(DS./R);
    zvalc=exp((DS-(DDS/1000))./R);
    zjj=zval.*ones(nres1,1);
    zjj(disr)=zvalc;
    zjj(ppos)=1;
    zjj(ppos2)=exp((DS-(18/1000))./R);
    %zjj(disr2)=exp((DS-(6.1/1000))./R);
    
    zvec = ones(BlockDet(end,2),1);
    for i=1:length(BlockDet)
        zvec(BlockDet(i,2)) = zvec(BlockDet(i,2)) * zjj(BlockDet(i,1));
    end
    
    pert=(nres+1)';
    npert=length(pert);
    permag=0.5;
    
    %% calculate Zfin
    ResProb = zeros(nres,length(pert));
    fesmat = zeros(nres,length(pert));
    Zfin=zeros(length(pert));
    pResProb = zeros(nres,nres,length(pert));
    pResProbSD = zeros(nres,nres,length(pert));
    pResProbSDL = zeros(nres,nres,length(pert));
    
    for kk=1:npert
        
        if pert(kk)<nres
            cmapmaskx=cmapmask;
            cmapmaskx(kk,:)=cmapmask(kk,:).*permag;
        else
            cmapmaskx=cmapmask;
        end
        
%         kk
        for ll=1:length(C)
            %% Elec. energy involving IS contribution
            ISfac=5.66*sqrt(IS/T)*sqrt(80/4);
            emapmask=zeros(nres);
            for i=1:nres
                for iin=i:nres
                    x1=find(contdist(:,1)==i & contdist(:,2)==iin);
                    emapmask(i,iin)=sum(contdist(x1,5).*exp(-ISfac.*contdist(x1,3)));
                end
            end
            pepval=zeros(nchoosek(nres+1,2)+2*nchoosek(nres+1,4),7);
            k=1;
            %% Generating Combinations for SSA
            for i=1:nres
                for iin=1:nres-i+1
                    stabE=0; eneE=0; v=0;
                    ncont=sum(sum(cmapmaskx(i:i+iin-1,i:i+iin-1)));
                    stabE=(ncont*(ene+DCp*(T-Tref)-T*DCp*log(T/Tref)));
                    eneE=sum(sum(emapmask(i:i+iin-1,i:i+iin-1)));
                    sw=exp(-(stabE+eneE)/(R*T))*prod(zvec(i:i+iin-1));
                    pepval(k,:)=[(iin) sw i iin 0 0 1];
                    pept1 = zeros(nres); pept1(i:i+iin-1,i:i+iin-1)=1;
                    pResProb(:,:,kk) =  pResProb(:,:,kk) + sw*pept1;
                    pResProbSD(:,:,kk) =  pResProbSD(:,:,kk) + sw*pept1;
                    pResProbSDL(:,:,kk) =  pResProbSDL(:,:,kk) + sw*pept1;
                    k=k+1;
                end
                %disp(i);
            end
            
            %% Generating Combinations for DSA
            for i=1:nres
                for iin=1:nres-i+1
                    for j=i+iin+1:nres
                        for jin=1:nres-j+1
                            stabE1=0; eneE1=0; stabE2=0; eneE2=0;
                            ncont1=sum(sum(cmapmaskx(i:i+iin-1,i:i+iin-1)));
                            stabE1=(ncont1*(ene+DCp*(T-Tref)-T*DCp*log(T/Tref)));
                            eneE1=sum(sum(emapmask(i:i+iin-1,i:i+iin-1)));
                            ncont2=sum(sum(cmapmaskx(j:j+jin-1,j:j+jin-1)));
                            stabE2=(ncont2*(ene+DCp*(T-Tref)-T*DCp*log(T/Tref)));
                            eneE2=sum(sum(emapmask(j:j+jin-1,j:j+jin-1)));
                            stabE=stabE1+stabE2;
                            eneE=eneE1+eneE2;
                            ncont=ncont1+ncont2;
                            sw=exp(-(stabE+eneE)/(R*T))*prod(zvec(i:i+iin-1))*prod(zvec(j:j+jin-1));
                            pepval(k,:)=[(iin+jin) sw i iin j jin 2];
                            pept1 = zeros(nres); vv = [(i:i+iin-1) (j:j+jin-1)]; pept1(vv,vv)=1;
                            pResProb(:,:,kk) =  pResProb(:,:,kk) + sw*pept1;
                            pResProbSD(:,:,kk) =  pResProbSD(:,:,kk) + sw*pept1;
                            k=k+1;
                        end
                    end
                end
                %disp(i);
            end
            
            %% Generating Combinations for DSAw/L
            for i=1:nres
                for iin=1:nres-i+1
                    for j=i+iin+1:nres
                        for jin=1:nres-j+1
                            ncont1=0; ncont2=0; stabE1=0; eneE1=0; stabE2=0; eneE2=0;  stabE3=0; eneE3=0;
                            vv = [(i:i+iin-1) (j:j+jin-1)];
                            ncont=sum(sum(cmapmaskx(vv,vv)));
                            stabE=(ncont*(ene+DCp*(T-Tref)-T*DCp*log(T/Tref)));
                            eneE=sum(sum(emapmask(vv,vv)));
                            if sum(sum(cmapmaskx(i:i+iin-1,j:j+jin-1)))~=0 || sum(sum(emapmask(i:i+iin-1,j:j+jin-1)))~=0
                                sw=exp(-(stabE+eneE)/(R*T))*prod(zvec(i:i+iin-1))*prod(zvec(j:j+jin-1))*zvalc^(j-(i+iin));
                            else
                                sw = 0;
                            end
                            pepval(k,:)=[(iin+jin) sw i iin j jin 3];
                            pept1 = zeros(nres); vv = [(i:i+iin-1) (j:j+jin-1)]; pept1(vv,vv)=1;
                            pResProb(:,:,kk) =  pResProb(:,:,kk) + sw*pept1;
                            pResProbSDL(:,:,kk) =  pResProbSDL(:,:,kk) + sw*pept1;
                            k=k+1;
                        end
                    end
                end
                %disp(i);
            end
%             toc;
            Zfin(kk) = sum(pepval(:,2));
            pResProb(:,:,kk) = pResProb(:,:,kk)./Zfin(kk);
            pResProbSD(:,:,kk) = pResProbSD(:,:,kk)./Zfin(kk);
            pResProbSDL(:,:,kk) = pResProbSDL(:,:,kk)./Zfin(kk);
            pepval = pepval(pepval(:,2)~=0,:);
            
            %% 1D free energy profile
            % The variable fes contains the 1D free energy profile
            fes = zeros(nres,1);
            for i=1:nres
                fes(i) = sum(pepval(pepval(:,1)==i,2));
            end
            fes = fes./sum(fes);
            fes = -R*T*log(fes);
            
            fesmat(:,kk)=fes;
            
            pepvalx(:,:,kk)=pepval;
            
            %% Residue folding probability as a function of the RC and 2D free-energy profile
            % The variable fes2D contains the 2D free energy profile
            % The variable Fpath contains the residue folding probability as a function of the RC
            hv = round(nres/2);
            Fpath = zeros(nres1,nres);
            FpathZ = zeros(nres,1);
            conv = zeros(nres,2);
            fes2D = zeros(hv+1,nres-hv+1);
            fes2DResProb = zeros(hv+1,nres-hv+1,nres);
            for i=1:nres
                f = find(BlockDet(:,2)==i);
                conv(i,:) = [BlockDet(f(1),1) BlockDet(f(end),1)];
            end
            
            for i=1:length(pepval)
                pept = zeros(nres1,1);
                pept(conv(pepval(i,3),1):conv(pepval(i,3)+pepval(i,4)-1,2))=1;
                if pepval(i,end)>1; pept(conv(pepval(i,5),1):conv(pepval(i,5)+pepval(i,6)-1,2))=1; end
                Fpath(:,pepval(i,1)) = Fpath(:,pepval(i,1)) + pepval(i,2)*pept;
                FpathZ(pepval(i,1)) = FpathZ(pepval(i,1)) + pepval(i,2);
                pept = zeros(nres,1);
                pept(pepval(i,3):pepval(i,3)+pepval(i,4)-1)=1;
                if pepval(i,end)>1; pept(pepval(i,5):pepval(i,5)+pepval(i,6)-1)=1; end
                ResProb(:,kk) = ResProb(:,kk) + pepval(i,2)*pept; % To calculate residue folding probability
                fes2D(sum(pept(1:hv))+1,sum(pept(hv+1:end))+1) = fes2D(sum(pept(1:hv))+1,sum(pept(hv+1:end))+1)+pepval(i,2);
                fes2DResProb(sum(pept(1:hv))+1,sum(pept(hv+1:end))+1,:) = fes2DResProb(sum(pept(1:hv))+1,sum(pept(hv+1:end))+1,:)+reshape(pepval(i,2)*pept,1,1,length(pept));
            end
            ResProb(:,kk) = ResProb(:,kk)./Zfin(kk);
            for i=1:nres
                Fpath(:,i) = Fpath(:,i)./FpathZ(i);
            end
            for i=1:length(fes2D(:,1))
                for j=1:length(fes2D(1,:))
                    fes2DResProb(i,j,:) = fes2DResProb(i,j,:)./fes2D(i,j);
                end
            end
            fes2D = fes2D./sum(sum(fes2D));
            fes2D = -R*T(kk)*log(fes2D);
            
        end
    end
end

%% Coupling calculation

x=squeeze(pepval(:,:,1));
nmic=length(squeeze(pepvalx(:,1,1)));
pv=zeros(nmic,npert);
respa=zeros(nres,npert);
respn=zeros(nres,npert);

pjfkf=zeros(nres,nres,npert);
pjfku=zeros(nres,nres,npert);
pjukf=zeros(nres,nres,npert);
pjuku=zeros(nres,nres,npert);

for kk=1:npert
    
%     kk
    
    Z(kk,1)=sum(squeeze(pepvalx(:,2,kk)));
    pv(:,kk)=squeeze(pepvalx(:,2,kk))./Z(kk);
    
    for j=1:nres
%         disp(j)
        a=find( (x(:,4)+x(:,3)-1)>=j & x(:,3)<=j | ((x(:,6)+x(:,5)-1)>=j & x(:,5)<=j) );
        respa(j,kk)=sum(squeeze(pepvalx(a,2,kk)))./Z(kk); % all states in which res j is structured
        
        n=setdiff(1:length(x),a,'stable'); % all states in which j in unstructured
        respn(j,kk)=sum(squeeze(pepvalx(n,2,kk)))./Z(kk); %sum(x(n,2))./Z;
        
        for k=1:nres
            
            a1=find( (x(:,4)+x(:,3)-1)>=k & x(:,3)<=k | ((x(:,6)+x(:,5)-1)>=k & x(:,5)<=k));
            pjfkf(j,k,kk)=sum(pv(intersect(a,a1,'stable'),kk));
            
            n1=setdiff(1:length(x),a1,'stable');
            a2=intersect(a,n1,'stable');
            pjfku(j,k,kk)=sum(pv(a2,kk));
            
            a3=intersect(n,a1,'stable');
            pjukf(j,k,kk)=sum(pv(a3,kk));
            
            a4=intersect(n,n1,'stable');
            pjuku(j,k,kk)=sum(pv(a4,kk));
            
            clear a1 a2 a3 a4 n1
            
        end
        
    end
end

%% Coupling Analysis

sthx=zeros(nres,nres);
sthx2=zeros(nres,nres);
sthx3=zeros(nres,nres,nres);
sthx4=zeros(nres,nres,nres);
sthx5=zeros(nres,nres,nres);

dGwtcp=R*T*log( (pjuku(:,:,end).*pjfkf(:,:,end))./ (pjfku(:,:,end).*pjukf(:,:,end))); % effective coupling free energy
chipluswt=R*T*log(squeeze(pjfkf(:,:,end))./(pjukf(:,:,end))); % positive coupling free energy matrix
chiminuswt=R*T*log(squeeze(pjfku(:,:,end))./(pjuku(:,:,end))); %negative coupling free energy matrix

for i=1:npert-1
        
    chiplus=R*T*log(squeeze(pjfkf(:,:,i))./(pjukf(:,:,i)));
    sthx=sthx+ chiplus; 
    
    chiminus=R*T*log(squeeze(pjfku(:,:,i))./(pjuku(:,:,i)));
    sthx2=sthx2+ chiminus;
    
    sthx3(:,:,i)=R*T*log((squeeze(pjuku(:,:,i)).*squeeze(pjfkf(:,:,i)))./(squeeze(pjfku(:,:,i)).*squeeze(pjukf(:,:,i))));
    
    sthx4(:,:,i)= (chiplus-chiminus); % is equivalent to sthx3
    
    sthx5(:,:,i)= (chiplus-chiminus)-(chipluswt-chiminuswt); 
    
    totpert(:,i)= nansum(squeeze(sthx5(:,:,i)));
    mutper(:,i)= nansum(squeeze(sthx3(:,:,i)));
    
end

a=(chipluswt+chipluswt')./2;
a2=(chiminuswt+chiminuswt')./2;

a(isinf(a))=NaN;
a2(isinf(a2))=NaN;

as=nansum(a);
a2s=nansum(a2);

am=nanmean(a);
a2m=nanmean(a2);

ast=nanstd(a);
a2st=nanstd(a2);

save('FEData.mat')