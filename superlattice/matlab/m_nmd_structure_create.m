function NMD = nmd_structure_create( NMD )
%--------------------------------------------------------------------------
%Unit cell and lattice vectors
%--------------------------------------------------------------------------

if strcmp(NMD.system(1).str,'LJ')==1

    if strcmp(NMD.system(2).str,'alloy')==1
        
        NMD.x0 = lj_alloy_create( NMD.x0 );
               
    elseif strcmp(NMD.system(2).str,'superlattice')==1

        NMD.x0 = lj_superlattice_create(NMD.x0);
    
    end

%elseif NMD.system(1).str == 'Si'
end
    

%NMD.kptlist = create_kptlist( NMD.x0.latvec, NMD.x0.Nx, NMD.x0.Ny, NMD.x0.Nz );

%output_universal( NMD );

end


