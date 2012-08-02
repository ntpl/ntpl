%function phonopy_ld
%--------------------------------------------------------------------------
%-------------------------INPUT--------------------------------------------
%--------------------------------------------------------------------------
constant.kb = 1.3806E-23;                    %aJ/k (1.3806E-23 J/K)
constant.hbar = 1.054E-34;                %J/s
constant.i = sqrt(-1);
constant.c = 29979245800.00019;      %cm/s
constant.s2ps = 10E-13;

constant.mass.C = 12.0/(6.022 * 10^23)/(10^3);   %kg
constant.mass.Si = 28.0855/(6.022 * 10^23)/(10^3);   %kg

constant.eV2J = 1.60217646E-19;
constant.Ang2m = 1E10;
constant.avo = 6.023E23;

%--------------------------------------------------------------------------
            [status,str.main] = system('pwd');
%--------------------------------------------------------------------------

%str.main = strcat(str.main,'/2x');

%--------------------------------------------------------------------------
%This is required to get the kpt + dk to properly input to GULP
            format long
%--------------------------------------------------------------------------

%KPTLIST
 	LD.kptlist(:,1:3) = load(strcat(str.main,'/kptlist.dat')); 
	[LD.NUM_KPTS, blank] = size(LD.kptlist(:,1:3));

	LD.dk = 0.00001;

%LD PARAMETERS
	%tmp = load(strcat(str.main,'/ld.param'));
	LD.NUM_POSCAR = 1; 
    LD.SPOSCAR.NUM_ATOMS = 16; 
    LD.NUM_ATOM_UCELL = 2;
   

%--------------------------------------------------------------------------
%This is required to get the kpt + dk to properly input to GULP
            format long
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%READ SPOSCAR LATVEC
%--------------------------------------------------------------------------
       str.str1 = 'grep -A '; 
	str.str2 = int2str(4);
	str.str3 = ' "Si" ';
	str.read = strcat(str.main,'/SPOSCAR');
	str.str4 = ' > sposcar.latvec';
	str.str_cmd = [str.str1 str.str2 str.str3 str.read str.str4];
	[status,tmp] = system(str.str_cmd);
	str.str_cmd = ('sed ''s/Si//g'' sposcar.latvec > sposcar.latvec2');
	[status,tmp] = system(str.str_cmd);
	LD.SPOSCAR.latvec = dlmread(strcat(str.main,'/sposcar.latvec2'));

%--------------------------------------------------------------------------	
%READ SPOSCAR DIRECT
%--------------------------------------------------------------------------
    str.str1 = 'grep -A '; 
	str.str2 = int2str(LD.SPOSCAR.NUM_ATOMS);
	str.str3 = ' "Direct" ';
	str.str4 = ' > sposcar.direct';
	str.str_cmd = [str.str1 str.str2 str.str3 str.read str.str4];
	[status,LD.POSCAR.direct] = system(str.str_cmd);
	str.str_cmd = ('sed ''s/Direct//g'' sposcar.direct > sposcar.direct2');
	[status,tmp] = system(str.str_cmd);
	LD.SPOSCAR.direct = load(strcat(str.main,'/sposcar.direct2'));

%--------------------------------------------------------------------------
%USE INPUT TO PROCESS
%--------------------------------------------------------------------------

	LD.SPOSCAR.cart(:,1) = LD.SPOSCAR.direct(:,1)*LD.SPOSCAR.latvec(2,1)...
        + LD.SPOSCAR.direct(:,2)*LD.SPOSCAR.latvec(3,1)...
        + LD.SPOSCAR.direct(:,3)*LD.SPOSCAR.latvec(4,1) ;
	LD.SPOSCAR.cart(:,2) = LD.SPOSCAR.direct(:,1)*LD.SPOSCAR.latvec(2,2)...
        + LD.SPOSCAR.direct(:,2)*LD.SPOSCAR.latvec(3,2)...
        + LD.SPOSCAR.direct(:,3)*LD.SPOSCAR.latvec(4,2) ;
	LD.SPOSCAR.cart(:,3) = LD.SPOSCAR.direct(:,1)*LD.SPOSCAR.latvec(2,3)...
        + LD.SPOSCAR.direct(:,2)*LD.SPOSCAR.latvec(3,3)...
        + LD.SPOSCAR.direct(:,3)*LD.SPOSCAR.latvec(4,3) ;
    
    LD.m(1:size(LD.SPOSCAR.cart,1)) = constant.mass.Si;

%Create a "cloud" to find nearest-image in dynam
    LD.CLOUD = [];
    for N1 = -1:1:1
        for N2 = -1:1:1
            for N3 = -1:1:1
%                 if N1==0 & N2==0 & N3==0
%                 else
                    tmp_pos(:,1) = N1*LD.SPOSCAR.latvec(2,1) +...
                        N2*LD.SPOSCAR.latvec(3,1) +...
                        N3*LD.SPOSCAR.latvec(4,1) +...
                        LD.SPOSCAR.cart(:,1);
                    
                    tmp_pos(:,2) = N1*LD.SPOSCAR.latvec(2,2) +...
                        N2*LD.SPOSCAR.latvec(3,2) +...
                        N3*LD.SPOSCAR.latvec(4,2) +...
                        LD.SPOSCAR.cart(:,2);
                    
                    tmp_pos(:,3) = N1*LD.SPOSCAR.latvec(2,3) +...
                        N2*LD.SPOSCAR.latvec(3,3) +...
                        N3*LD.SPOSCAR.latvec(4,3) +...
                        LD.SPOSCAR.cart(:,3);
                    
                LD.CLOUD = [LD.CLOUD;tmp_pos];
                    
%                 end
            end
        end
    end
    
plot3(LD.SPOSCAR.cart(:,1),LD.SPOSCAR.cart(:,2),LD.SPOSCAR.cart(:,3),'.',...
LD.CLOUD(1:16:size(LD.CLOUD,1),1),LD.CLOUD(1:16:size(LD.CLOUD,1),2),...
LD.CLOUD(1:16:size(LD.CLOUD,1),3),'.')

%--------------------------------------------------------------------------
%READ FORCE CONSTANTS
%--------------------------------------------------------------------------

%Linux
str.read=strcat(str.main,'/FORCE_CONSTANTS');
dummy = dlmread(str.read); 
%resize dummy to chop off first line  
dummy = dummy(2:size(dummy,1),:); 
    FC.id_cnt=1;
    for i1=1:4:length(dummy)
        FC.id(FC.id_cnt,1) = dummy( i1,1);         %i
        FC.id(FC.id_cnt,2) = dummy( i1,2);         %j
        FC.FC(FC.id_cnt,1:3,1:3) = dummy( (i1+1) : (i1+3),1:3)*...
            (constant.eV2J)*(constant.Ang2m^2)  ;      %eV/Ang^2
    FC.id_cnt=FC.id_cnt+1;
    end

%--------------------------------------------------------------------------
%pause
%-------------------------------------------------------------------------- 


    LD.alat(1:3) = 5.43E-10/2;
    kappa = [0.25 0.0 0.0];
%convert LD.SPOSCAR.cart to meters
    LD.SPOSCAR.cart = LD.SPOSCAR.cart/constant.Ang2m; 
    LD.CLOUD = LD.CLOUD/constant.Ang2m;
    rtol = 1E-5;
    
    cnt_Ir = 1; cnt_single=1; r_cutoff = 7E-10;

%[freq,eigVsorted] = phonopy_dynam_matrix(LD,FC,constant,kappa);

%FUNCTION: find dynamical matrix
DYNAM(1:3*LD.NUM_ATOM_UCELL,1:3*LD.NUM_ATOM_UCELL) = 0.0;
 
for iucell = 1:LD.NUM_ATOM_UCELL  
%identify atom in the ucell, phononpy FC.FC list has redundancy
    Ifc = find(FC.id(:,1)==iucell);
    iposcar = (iucell-1)*(LD.SPOSCAR.NUM_ATOMS/LD.NUM_ATOM_UCELL) +1;
    for iFC=1:size(Ifc,1)  
%find nearest image by "cloud" convention
%people.virginia.edu/~lz2n/mse627/notes/Boundary.pdf
        iimages = FC.id(Ifc(iFC),2):LD.SPOSCAR.NUM_ATOMS:size(LD.CLOUD,1);    
        r_images = LD.CLOUD(iimages,1:3);
        rij = bsxfun(@minus,...
             r_images, LD.SPOSCAR.cart( iposcar ,1:3) );
        r = sqrt(rij(:,1).^2 + rij(:,2).^2 + rij(:,3).^2);
        [Y,Irij] = min(r);
%find all the minimum images, there are sometimes more than 1
        Ir = find( abs(r - min(r))/max(r) < rtol );
%plot(r*1E10)
        plot3(LD.SPOSCAR.cart( iposcar ,1),...
            LD.SPOSCAR.cart( iposcar ,2),...
            LD.SPOSCAR.cart( iposcar ,3),...
            '.',...
            LD.CLOUD(iimages(Ir),1),...
            LD.CLOUD(iimages(Ir),2),...
            LD.CLOUD(iimages(Ir),3),...
            '.',...
            LD.SPOSCAR.cart(FC.id(Ifc(iFC),2),1),...
            LD.SPOSCAR.cart(FC.id(Ifc(iFC),2),2),...
            LD.SPOSCAR.cart(FC.id(Ifc(iFC),2),3),...
            '.')
%         LD.CLOUD(iimages,1),...
%         LD.CLOUD(iimages,2),...
%         LD.CLOUD(iimages,3),...
%         '.')

%if r(Ir(1)) < r_cutoff



%sort(r)

%size(Ir,1)

%--------------------------------------------------------------------------
%pause
%-------------------------------------------------------------------------- 

%sum over all the images
invmexpikR = 0;
    for ir = 1:size(Ir,1)
    dot_product =...
        (2*pi)*(...
        kappa(1)*rij(Ir(ir),1)/LD.alat(1) + ...
        kappa(2)*rij(Ir(ir),2)/LD.alat(2) + ...
        kappa(3)*rij(Ir(ir),3)/LD.alat(3) ); 
    invmexpikR = invmexpikR +...
        exp(sqrt(-1)*dot_product) ; 
    end
    
%the copies of iucell=1 are in first half of SPOSCAR   
    ii = (FC.id(Ifc(iFC),1)-1)*3 + 1 : (FC.id(Ifc(iFC),1))*3;
    if FC.id(Ifc(iFC),2)<=LD.SPOSCAR.NUM_ATOMS/LD.NUM_ATOM_UCELL
        jj = (1-1)*3 + 1 : (1)*3;
    else
        jj = (2-1)*3 + 1 : (2)*3;
    end  
%possible sign change?
    if FC.id(Ifc(iFC),1)==FC.id(Ifc(iFC),2)
        Phi(1:3,1:3) = invmexpikR*FC.FC(Ifc(iFC),1:3,1:3);
    else
        Phi(1:3,1:3) = -invmexpikR*FC.FC(Ifc(iFC),1:3,1:3);
    end
    
    DYNAM( ii , jj  ) = DYNAM(  ii , jj ) +...
        Phi/size(Ir,1)*...
        (1/sqrt(LD.m(FC.id(Ifc(iFC),1))*LD.m(FC.id(Ifc(iFC),2)))) ;

%else
%end
    end
end
 
DYNAM = ( DYNAM + conj(transpose(DYNAM))) / 2;

%Calculate freqs, eigV, and v_g
    [eigV,W2]=eig(DYNAM); W=sqrt(W2);
%CREATE FREQ VECTOR          
    for i=1:length(W) w(i) = W(i,i); end   
%SORT FREQS AND EIGV
	[freq,I]=sort(w);
	eigVsorted(1:3*LD.NUM_ATOM_UCELL,1:3*LD.NUM_ATOM_UCELL) = eigV(:,I);  


freq/(2*pi)


plot(1:length(r_single),r_single,...
1:length(r_multiple),r_multiple)


%     dis_image = sqrt(...
%     (LD.SPOSCAR.cart( iposcar ,1) - LD.CLOUD(iimages(Ir),1) )^2 +...
%     (LD.SPOSCAR.cart( iposcar ,2) - LD.CLOUD(iimages(Ir),2) )^2 +...
%     (LD.SPOSCAR.cart( iposcar ,3) - LD.CLOUD(iimages(Ir),3) )^2 )
% 
%     dis_SPOSCAR = sqrt(...
%     (LD.SPOSCAR.cart( iposcar ,1) - LD.SPOSCAR.cart(FC.id(Ifc(iFC),2),1))^2 +...
%     (LD.SPOSCAR.cart( iposcar ,2) - LD.SPOSCAR.cart(FC.id(Ifc(iFC),2),2))^2 +...
%     (LD.SPOSCAR.cart( iposcar ,3) - LD.SPOSCAR.cart(FC.id(Ifc(iFC),2),3))^2 )
%--------------------------------------------------------------------------
%pause
%-------------------------------------------------------------------------- 









%--------------------------------------------------------------------------
%OUPUT PHONON PROPS--------------------------------------------------------
%--------------------------------------------------------------------------
%         str_write_vec=strcat(str.main,'/eigvec.dat');
%         str_write_freq=strcat(str.main,'/freq.dat');
%         str_write_vel=strcat(str.main,'/vel.dat');        
% 
%             for ikpt = 1:LD.NUM_KPTS
%                 LD.kptlist(ikpt,1:3)
%                 [freq,eigVsorted,velocity] = phonopy_LD(LD,FC,constant,LD.kptlist(ikpt,1:3));
%                     dlmwrite(str_write_vec,eigVsorted,'-append','delimiter',' ','precision',8);
%                     dlmwrite(str_write_freq,freq,'-append','delimiter',' ','precision',8);
%                     dlmwrite(str_write_vel,velocity,'-append','delimiter',' ','precision',8);
% 
%             end

%--------------------------------------------------------------------------
%END MAIN------------------------------------------------------------------
%-------------------------------------------------------------------------- 



% function [freq,eigVsorted,velocity] = phonopy_LD(LD,FC,constant,kappa)
% 
% %CALC DYNAM MATRIX
%             [freq,eigVsorted] = phonopy_dynam_matrix(LD,FC,constant,kappa);
%             
% %CALCULATE VG numerically              
%          
%             for idim = 1:3
%                 if kappa(idim)==0.5
%                     [HLDpks1] = phonopy_dynam_matrix(LD,FC,constant,kappa);
%                     kappa(idim) = kappa(idim) - LD.dk;
%                     [HLDpks2] = phonopy_dynam_matrix(LD,FC,constant,kappa);
%                     velocity(:,idim) = (( HLDpks1 - HLDpks2 )/ LD.dk );
%                     kappa(idim) = kappa(idim) + LD.dk;
%                 elseif kappa(idim)==-0.5
%                     [HLDpks1] = phonopy_dynam_matrix(LD,FC,constant,kappa);
%                     kappa(idim) = kappa(idim) + LD.dk;
%                     [HLDpks2] = phonopy_dynam_matrix(LD,FC,constant,kappa);
%                     velocity(:,idim) = (( HLDpks1 - HLDpks2 )/LD.dk);
%                     kappa(idim) = kappa(idim) - LD.dk;
%                 elseif kappa(idim)==0.0
%                     kappa(idim) = kappa(idim) + LD.dk;
%                     [HLDpks1] = phonopy_dynam_matrix(LD,FC,constant,kappa);
%                     kappa(idim) = kappa(idim) - LD.dk;
%                     [HLDpks2] = phonopy_dynam_matrix(LD,FC,constant,kappa);
%                     velocity(:,idim) = (( HLDpks1 - HLDpks2 )/LD.dk);
%                 else
%                     kappa(idim) = kappa(idim) + LD.dk;
%                     [HLDpks1] = phonopy_dynam_matrix(LD,FC,constant,kappa);
%                     kappa(idim) = kappa(idim) - 2*LD.dk;
%                     [HLDpks2] = phonopy_dynam_matrix(LD,FC,constant,kappa);
%                     velocity(:,idim) = (( HLDpks1 - HLDpks2 )/(2*LD.dk));
%                     kappa(idim) = kappa(idim) + LD.dk;
%                 end
%             end
% 
% end  
% %--------------------------------------------------------------------------
% %END FUNCTION--------------------------------------------------------------
% %--------------------------------------------------------------------------            
%             
%             
% function [freq,eigVsorted] = phonopy_dynam_matrix(LD,FC,constant,kappa)
% %FUNCTION: find dynamical matrix
%  DYNAM(1:3*LD.NUM_ATOMS_UCELL,1:3*LD.NUM_ATOMS_UCELL) = 0.0;
%     for i1=1:length(FC.id)  
% 	id.i = FC.id(i1,1)-floor( (FC.id(i1,1)-1)/LD.NUM_ATOMS_UCELL)*LD.NUM_ATOMS_UCELL;
% 	id.j = FC.id(i1,2)-floor( (FC.id(i1,2)-1)/LD.NUM_ATOMS_UCELL)*LD.NUM_ATOMS_UCELL;
% 
% %	FC.id(i1,1)
% %	FC.id(i1,2)
% 
% %	id.i
% %	id.j
% 
% 	rij = LD.SPOSCAR.cart(FC.id(i1,1),1:3) - LD.SPOSCAR.cart(FC.id(i1,2),1:3);
%         dot_product = (2*pi)*( kappa(1)*rij(1)/LD.alat(1) + kappa(2)*rij(2)/LD.alat(2) + kappa(3)*rij(3)/LD.alat(3) ); 
%         invmexpikR = exp(sqrt(-1)*dot_product)*(1/sqrt(LD.m(FC.id(i1,1))*LD.m(FC.id(i1,1)))); 
%         if FC.id(i1,1)==FC.id(i1,2)
%             Phi(1:3,1:3) = invmexpikR*FC.FC(1:3,1:3,i1);
%         else
%             Phi(1:3,1:3) = invmexpikR*FC.FC(1:3,1:3,i1);
%         end
%         ii = (id.i-1)*3 + 1 : (id.i)*3;
%         jj = (id.j-1)*3 + 1 : (id.j)*3;
% %--------------------------------------------------------------------------
% %pause
% %-------------------------------------------------------------------------- 
%     DYNAM( ii , jj  ) = DYNAM(  ii , jj ) + Phi ; %real(Phi)+sqrt(-1)*imag(Phi);
%     end
% %Calculate freqs, eigV, and v_g
%             [eigV,W2]=eig(DYNAM); W=sqrt(abs(W2));
% %CREATE FREQ VECTOR          
%             for i=1:length(W) w(i) = W(i,i); end   
% %SORT FREQS AND EIGV
% 	[freq,I]=sort(abs(w));
% 	eigVsorted(1:3*LD.NUM_ATOMS_UCELL,1:3*LD.NUM_ATOMS_UCELL) = eigV(:,I);  
% end
% %--------------------------------------------------------------------------
% %END FUNCTION--------------------------------------------------------------
% %-------------------------------------------------------------------------- 
% 
% 
% end
%--------------------------------------------------------------------------
%END MAIN------------------------------------------------------------------
%-------------------------------------------------------------------------- 

