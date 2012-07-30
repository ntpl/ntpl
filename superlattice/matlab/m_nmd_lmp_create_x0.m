function nmd_lmp_create_x0( NMD )

%output lammps    
        str.orig = 'NUM_ATOMS';
        str.change = [int2str(length(NMD.x0.pos(:,1)))];
        str.cmd1 = ['-e ''s/\<' str.orig '\>/' str.change '/g'' '];
   
        str.orig = 'NUM_ATOMS_TYPE';
        str.change = [int2str(NMD.NUM_ATOMS_TYPE)];
        str.cmd2 = ['-e ''s/\<' str.orig '\>/' str.change '/g'' '];
        
        str.orig = 'LX';
        str.change = [num2str((NMD.x0.Nx)*NMD.x0.alat(1))]% *NMD.x0.superlattice.period(1,1)) ];
        str.cmd3 = ['-e ''s/\<' str.orig '\>/' str.change '/g'' '];
        str.orig = 'LY';
        str.change = [num2str((NMD.x0.Ny)*NMD.x0.alat(2))]% *NMD.x0.superlattice.period(2,2))];
        str.cmd4 = ['-e ''s/\<' str.orig '\>/' str.change '/g'' '];
        str.orig = 'LZ';
        str.change = [num2str((NMD.x0.Nz)*NMD.x0.alat(3))]% *NMD.x0.superlattice.period(3,3))];
        str.cmd5 = ['-e ''s/\<' str.orig '\>/' str.change '/g'' '];
        str.orig = 'ATOM_MASS_1';
        str.change = ['1 ' num2str(NMD.x0.mass(1))];
        str.cmd6 = ['-e ''s/\<' str.orig '\>/' str.change '/g'' '];
        
    	str.orig = 'ATOM_MASS_2';
    	str.change = ['2 ' num2str(NMD.x0.mass(2))];
    	str.mass2 = ['-e ''s/\<' str.orig '\>/' str.change '/g'' '];
        
        str.cmd8 =...
        ['lmp.in.x0.temp > ./' int2str(NMD.seed.superlattice ) '/lmp.in.x0.' NMD.system(2).str];
        
        str.cmd = ['sed ' str.cmd1 str.cmd2 str.cmd3 str.cmd4 str.cmd5...
            str.cmd6 str.mass2 str.cmd8];
        
        system(str.cmd);
        
        output = [NMD.x0.pos(:,1:5)];
        str.write=['./' int2str(NMD.seed.superlattice ) '/lmp.in.x0.' NMD.system(2).str];
        dlmwrite(str.write,output,'-append','delimiter','\t');



end
