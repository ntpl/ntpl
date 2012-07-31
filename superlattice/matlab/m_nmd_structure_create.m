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

end

end


