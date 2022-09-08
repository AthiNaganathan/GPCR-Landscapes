% Single Sequence Approximation

clear; clc;
tic;

%% Looping PDB files

for run=1
    clearvars -except run;
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
    
    %% Variables
    ene=-48.2/1000; %modify to the required vdW interaction energy value
    DS=-10/1000;
    DCp=-0.3579/1000;
    
    int = 1.6;
    slope = 6.7/1000;
    
    Mw=nres1*110;
    disr = [];
    ppos=[];
    
    T=(273:1:373);
    IS=0.1; % in molar units
    C=0;
    
    hh9=fopen('struct.txt','rt');
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
    
    %% Constants
    R=0.008314;
    Tref=385;
    zval=exp(DS./R);
    zvalc=exp((DS-(6.0606/1000))./R);
    zjj=zval.*ones(nres1,1);
    zjj(disr)=zvalc;
    zjj(ppos)=1;
    zvec = ones(BlockDet(end,2),1);
    for i=1:length(BlockDet);
        zvec(BlockDet(i,2)) = zvec(BlockDet(i,2)) * zjj(BlockDet(i,1));
    end
    
    %% calculate Zfin
    Zfin=zeros(length(T),1);
    parfor kk=1:length(T)
        for ll=1:length(C)
            %% Elec. energy involving IS contribution
            ISfac=5.66*sqrt(IS/T(kk))*sqrt(80/4);
            emapmask=zeros(nres);
            for i=1:nres
                for iin=i:nres
                    x1=find(contdist(:,1)==i & contdist(:,2)==iin);
                    emapmask(i,iin)=sum(contdist(x1,5).*exp(-ISfac.*contdist(x1,3)));
                end
            end
            k=1;
            %% Generating Combinations for SSA
            for i=1:nres
                for iin=1:nres-i+1
                    stabE=0; eneE=0; v=0;
                    ncont=sum(sum(cmapmask(i:i+iin-1,i:i+iin-1)));
                    stabE=(ncont*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                    eneE=sum(sum(emapmask(i:i+iin-1,i:i+iin-1)));
                    sw=exp(-(stabE+eneE)/(R*T(kk)))*prod(zvec(i:i+iin-1));
                    Zfin(kk) = Zfin(kk)+sw;
                    k=k+1;
                end
            end
            
            %% Generating Combinations for DSA
            for i=1:nres
                for iin=1:nres-i+1
                    for j=i+iin+1:nres
                        for jin=1:nres-j+1
                            stabE1=0; eneE1=0; stabE2=0; eneE2=0;
                            
                            ncont1=sum(sum(cmapmask(i:i+iin-1,i:i+iin-1)));
                            stabE1=(ncont1*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                            eneE1=sum(sum(emapmask(i:i+iin-1,i:i+iin-1)));
                            ncont2=sum(sum(cmapmask(j:j+jin-1,j:j+jin-1)));
                            stabE2=(ncont2*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                            eneE2=sum(sum(emapmask(j:j+jin-1,j:j+jin-1)));
                            stabE=stabE1+stabE2;
                            eneE=eneE1+eneE2;
                            ncont=ncont1+ncont2;
                            sw=exp(-(stabE+eneE)/(R*T(kk)))*prod(zvec(i:i+iin-1))*prod(zvec(j:j+jin-1));
                            Zfin(kk) = Zfin(kk)+sw;
                            k=k+1;
                        end
                    end
                end
            end
            
            %% Generating Combinations for DSAw/L
            for i=1:nres
                for iin=1:nres-i+1
                    for j=i+iin+1:nres
                        for jin=1:nres-j+1
                            vv = [(i:i+iin-1) (j:j+jin-1)];
                            ncont=sum(sum(cmapmask(vv,vv)));
                            stabE=(ncont*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                            eneE=sum(sum(emapmask(vv,vv)));
                            if sum(sum(cmapmask(i:i+iin-1,j:j+jin-1)))~=0 || sum(sum(emapmask(i:i+iin-1,j:j+jin-1)))~=0
                                sw=exp(-(stabE+eneE)/(R*T(kk)))*prod(zvec(i:i+iin-1))*prod(zvec(j:j+jin-1))*zvalc^(j-(i+iin));
                            else
                                sw = 0;
                            end
                            Zfin(kk) = Zfin(kk)+sw;
                            k=k+1;
                        end
                    end
                end
            end
        end
        disp(kk);
    end
    logZ=log(Zfin);
    Tint=(T(1)-10:0.1:T(end)+10)';
    logZint=interp1(T,logZ,Tint,'spline','extrap');
    
    der1d=diff(logZint)./diff(Tint);
    der1df=interp1(Tint(1:end-1),der1d,T,'spline');
    der1d2=interp1(T,der1df,Tint,'spline','extrap');
    der2d=diff(der1d2)./diff(Tint);
    der2df=interp1(Tint(1:end-1),der2d,T,'spline');
    
    Cpd=2*R*T.*der1df+R*T.^2.*der2df;
    Cpfin=Cpd + (int+slope*(T-273.15))*Mw./1000;
%     plot(T,Cpfin,'b');
%     hold on
%     plot(T,Cpd,'r');
end
save('DSCdata.mat')