% Computes VdW Contact Map and Electrostatic energetics

% clear; clc;
tic;

%% Looping files
for run=1
%     clearvars -except run;
    pdblist=char('pdbfilename'); %modify to the name of the relevant PDB file
    BlockSize = 4;
    cdata=cat(2,deblank(pdblist(run,:)),'.pdb');
    pdbr=fopen(cdata,'rt'); % Raw PDB data
    pdbc=fopen('buffer1.txt','wt'); % Curated PDB data
    
    %% Remove Hydrogens and hetero atoms including nucleotides
    % gets the raw PDB file and generates the buffer1.txt required for
    % calculations
    line=fgetl(pdbr);
    while feof(pdbr)~=1
        if(length(line)>4)
            if(strcmp(line(1:4),'ATOM')==1 && strcmp(line(13),'H')==0 && strcmp(line(13),'Q')==0 && strcmp(line(14),'H')==0 && strcmp(line(18),' ')==0)
                % if(strcmp(line(1:4),'ATOM')==1 && strcmp(line(13),'H')==0 && strcmp(line(13),'Q')==0 && strcmp(line(14),'H')==0 && strcmp(line(18),' ')==0 && strcmp(line(13:14),'CA')==1)
                fprintf(pdbc,'%s\n',line);
            end
            if(strcmp(line(1:3),'TER')==1)
                break;
            end
        end
        line=fgetl(pdbr);
    end
    fclose(pdbr);
    fclose(pdbc);
    
    %% Get atmoic details with their corresponding residue number
    % atomn: Atom name, resno: Residue number, x,y,z: coordinates
    
    pdbc=fopen('buffer1.txt','rt');
    i=1; j=1; prevres=0; nnegres=0; nposres=0; nhphilres=0; iin=1; iim=1;
    aares={'GLY','ALA','VAL','LEU','ILE','MET','PHE','TYR','TRP','SER','ASP','ASN','THR','GLU','GLN','HIS','LYS','ARG','PRO','CYS'};
    aacode={'G','A','V','L','I','M','F','Y','W','S','D','N','T','E','Q','H','K','R','P','C'};
    charres={'HIS','LYS','ARG','GLU','ASP'};
    posres=[16 17 18]; % [HIS LYS ARG]
    negres=[11 14]; % [ASP GLU]
    hphilres=[8 10 12 13 15 20]; % [TYR SER ASN THR GLN CYS]
    atomc={'NE ','NH1','NH2','NZ ','OD1','OD2','OE1','OE2','ND1','NE2'};
    atombb={'C','N','CA'};
    charmag7=[0.33 0.33 0.33 1 -0.5 -0.5 -0.5 -0.5 0 0]'; % pH 7
    charmag5=[0.33 0.33 0.33 1 -0.5 -0.5 -0.5 -0.5 0.5 0.5]'; % pH 5
    charmag3p5=[0.33 0.33 0.33 1 -0.25 -0.25 -0.25 -0.25 0.5 0.5]'; % pH 3.5
    charmag2=[0.33 0.33 0.33 1 0 0 0 0 0.5 0.5]'; % pH 2
    
    while feof(pdbc)~=1
        line=fgetl(pdbc);
        if(length(line)>4)
            atomn(i,:)=line(14:16);
            resno(i,1)=str2double(line(23:26));
            x(i,1)=str2double(line(32:38));
            y(i,1)=str2double(line(40:46));
            z(i,1)=str2double(line(48:54));
            if strcmp(line(14:16),'N  ')
                nx(iin,1)=str2double(line(32:38));
                ny(iin,1)=str2double(line(40:46));
                nz(iin,1)=str2double(line(48:54));
                iin=iin+1;
            end
            if strcmp(line(14:16),'C  ')
                cx(iim,1)=str2double(line(32:38));
                cy(iim,1)=str2double(line(40:46));
                cz(iim,1)=str2double(line(48:54));
                iim=iim+1;
            end
            % Assign charge magnitude at specified pH
            charmag(i,1)=0;
            if sum(strcmp(charres,line(18:20)))
                [~,x1] = find(strcmp(atomc,atomn(i,:)));
                if isempty(x1)==0
                    charmag(i,1)=charmag7(x1); %-----------------> To be changed <------------------
                end
            end
            
            % Calculate Seq Composition
            if (resno(i)-prevres~=0)
                
                % Get sequence from PDB
                [~,x2]=find(strcmp(line(18:20),aares));
                protseq(1,j)=aacode(x2);
                prevres=resno(i);
                
                % Amino Acid Composition
                if sum(hphilres==x2)
                    nhphilres=nhphilres+1;
                elseif sum(posres==x2)
                    nposres=nposres+1;
                elseif sum(negres==x2)
                    nnegres=nnegres+1;
                end
                j=j+1;
            end
        end
        i=i+1;
    end
    
    %% Compute Short and Long-Range interaction map
    % Short_Range: no. of atomic contacts and h-bonds b/w residues
    % Long_Range: details for electrostatic interactions
    
    srcutoff=5.0; % Cutoff dist for VdW interaction
    ecutoff=1000; % Cutoff dist for electrostatic interaction
    seqsepcutoff=0;
    %seqsepcutoff=3; % Sequence seperation cutoff to classify local and non local contacts
    
    %-----Initializarion-----------------------------------------------------------------------
    atomt=length(x(:,1));
    nres=resno(atomt,1)-resno(1,1)+1;
    loccont=zeros(nres); locconte=zeros(nres);
    nloccont=zeros(nres); nlocconte=zeros(nres);
    srcont=zeros(nres);
    contbetrese=zeros(nres); contbetres=zeros(nres);
    srtcont=0; etcont=0;
    
    %------------------------------------------------------------------------------------------
    
    startres=resno(1,1);
    resno(:,1)=resno(:,1)-startres+1; %Changing res. no. according to the PDB file
    k=1;
    
    for i=1:(atomt-1)
        for j=(i+1):atomt
            dist=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2+(z(i)-z(j))^2);
            
            % Calculate Short-Range Interactions
            if((resno(j)-resno(i))>=1 && (charmag(i)==0 || charmag(j)==0))
                dist=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2+(z(i)-z(j))^2);
                
                if(dist<=srcutoff)
                    srcont(resno(i),resno(j))=srcont(resno(i),resno(j))+1;
                    
                    % Classify local/non-local SR interaction
                    if(abs(resno(i)-resno(j))<=seqsepcutoff);
                        loccont(resno(i),resno(j))=loccont(resno(i),resno(j))+1;
                    else
                        nloccont(resno(i),resno(j))=nloccont(resno(i),resno(j))+1;
                    end
                    
                    srtcont=srtcont+1;
                    contbetres(resno(i),resno(j))=contbetres(resno(i),resno(j))+1;
                    
                end
            end
            % Calculate Long-Range Interactions and associated Energy
            if ((resno(j)-resno(i))>=1 && (charmag(i)~=0 && charmag(j)~=0))
                diste=sqrt(((x(i)-x(j))^2)+((y(i)-y(j))^2)+((z(i)-z(j))^2));
                if(diste<=ecutoff)
                    intene=332*4.184*charmag(i)*charmag(j)/(4*diste); %Elec. Contibution
                    
                    % Classify local/non-local LR interaction
                    if(abs(resno(i)-resno(j))<=seqsepcutoff);
                        locconte(resno(i),resno(j))=locconte(resno(i),resno(j))+1;
                    else
                        nlocconte(resno(i),resno(j))=nlocconte(resno(i),resno(j))+1;
                    end
                    etcont=etcont+1;
                    contbetrese(resno(i),resno(j))=contbetrese(resno(i),resno(j))+1;
                    
                    % Elec. Info in matrix form
                    ElecMat(k,:)=[resno(i) resno(j) diste abs(resno(i)-resno(j)) intene];
                    k=k+1;
                end
            end
        end
    end
    
    %% Get SecStr Info
%     cmd = cat(2,'stride ',cdata,' -fstruct.txt');
%     system(cmd);
    
    hh9=fopen('struct.txt','rt'); i=1; j=1; %secondary structure prediction file generated using STRIDE
    lkk=fgetl(hh9);
    while lkk > 0
        if length(lkk)> 25 && strcmp(lkk(1:3),'ASG')
            h(i,1) = 0;
            if (strcmp(lkk(25),'H') || strcmp(lkk(25),'E') || strcmp(lkk(25),'G'))
                h(i,1) = 1;
            end
            i=i+1;
        end
        lkk=fgetl(hh9);
    end
    fclose(hh9);
    
    h1=[[0;h] [h;0]];
    hbeg = find(h1(:,1)==0 & h1(:,2)==1);
    hend = find(h1(:,1)==1 & h1(:,2)==0);
    StructBlock = startres; hbeg = hbeg+startres-1; hend=hend+startres-1;
    for i=1:length(hbeg)
        StructBlock = [StructBlock;hbeg(i);hend(i)];
    end
    StructBlock = [StructBlock;resno(end)+startres];
    residual = rem(StructBlock(2:end)-StructBlock(1:end-1),BlockSize);
    nBlocks = floor((StructBlock(2:end)-StructBlock(1:end-1))/BlockSize);
    
    k=1; BlockID = 0; uresno = unique(resno);
    BlockMat = zeros(length(uresno),2); flag=0;
    for i=2:length(StructBlock)
        for j=1:nBlocks(i-1)
            BlockID = BlockID+1;
            BlockMat(k:k+BlockSize-1,:) = [uresno(k:k+BlockSize-1) BlockID*ones(BlockSize,1)];
            k=k+BlockSize; flag=1;
        end
        if residual(i-1) == 1
            if flag==0; BlockID = BlockID+1; end
            BlockMat(k,:) = [uresno(k) BlockID];
            k=k+1;
        elseif residual(i-1) > 1
            BlockID = BlockID+1;
            BlockMat(k:k+residual(i-1)-1,:) = [uresno(k:k+residual(i-1)-1) BlockID*ones(residual(i-1),1)];
            k=k+residual(i-1);
        end
    end
    
    contactMapB = zeros(BlockID);
    for i=1:nres
        for j=i+1:nres
            contactMapB(BlockMat(i,2),BlockMat(j,2)) = contactMapB(BlockMat(i,2),BlockMat(j,2)) + srcont(i,j);
        end
    end
    
    contdistMapB = zeros(length(ElecMat(:,1)),5);
    for i=1:length(ElecMat(:,1))
        contdistMapB(i,:) = [BlockMat(ElecMat(i,1),2) BlockMat(ElecMat(i,2),2) ElecMat(i,3:end)];
    end
    
    %% Printing Regular Comtact Map files
    eval(['save contactmapmatElec',deblank(pdblist(run,:)),'.dat srcont -ascii']);
    eval(['save contdistElec',deblank(pdblist(run,:)),'.dat ElecMat -ascii']);
    eval(['save contdistElecB',deblank(pdblist(run,:)),'.dat contdistMapB -ascii']);
    eval(['save contactmapmatElecB',deblank(pdblist(run,:)),'.dat contactMapB -ascii']);
    eval(['save BlockDet',deblank(pdblist(run,:)),'.dat BlockMat -ascii']);
end
