function x0 = lj_superlattice_create( x0 )

        x0 = lj_superlattice_unit_cell( x0 );
        x0.pos(:,3:5) =...
            m_build_supercell(...
        x0.superlattice.direct, x0.superlattice.period, x0.Nx, x0.Ny, x0.Nz );
	x0.pos(:,3:5) = x0.alat(2).*x0.pos(:,3:5);
	%size instead length
	x0.pos(:,2)=repmat(x0.m',length(x0.pos(:,2))/length(x0.m),1);
	x0.pos(:,1)=(1:1:length(x0.pos(:,1)))';
	if strcmp(x0.type(2).str,'GAMMA')
		x0.kptlist=[0, 0, 0];
		x0.ucell.cart= x0.alat(2).*x0.pos(:,3:5);		
	else
		%x0.kptlist=create_kptlist(  x0.latvec/x0.superlattice.period, x0.Nx , x0.Ny, x0.Nz );
		x0.kptlist=create_kptlist(  x0.latvec, x0.Nx , x0.Ny, x0.Nz );
		%x0.kptlist(:,1)=zeros(length(x0.kptlist(:,1)),1)
		x0.kptlist(:,1)=x0.Nx.*x0.kptlist(:,1);%*x0.superlattice.period(1,1);
		x0.kptlist(:,2)=x0.Ny.*x0.kptlist(:,2);
		x0.kptlist(:,3)=x0.Nz.*x0.kptlist(:,3);
		x0.ucell.cart= x0.alat(2).*x0.superlattice.direct;
		dlmwrite('kptlist.dat',x0.kptlist,'delimiter',' ');
	end
end
