clear
    NMD.x0.LJ.eps = 1.67E-21;              %aJ (1.67E-21 Joules) aJ=1E-18 J
    NMD.x0.LJ.sigma = 3.4E-10;                 %Angstroms 3.4E-10 meters
    NMD.x0.LJ.a_0 = 5.2686E-10/NMD.x0.LJ.sigma; %the lattice constant of Ar: http://www.infoplease.com/periodictable.php?id=18
    NMD.x0.LJ.mass = 6.6326E-26;               %1E-28 kg (6.6326E-26 kg)
    NMD.x0.LJ.tau = sqrt((NMD.x0.LJ.mass*(NMD.x0.LJ.sigma^2))/NMD.x0.LJ.eps);
    NMD.x0.constant.c = 29979245800.00019;      %cm/s
    NMD.x0.constant.s2ps = 1E-12;
    NMD.x0.constant.ang2m = 1E-10;
    NMD.x0.constant.eV2J = 1.60217646E-19;
%--------------------------------------------------------------------------
%STRUCTURE-----------------------------------------------------------------    
%--------------------------------------------------------------------------  
    [tmp,NMD.str.main]=system('pwd');
    NMD.title = '1';
%--------------------------------------------------------------------------
    NMD.system(1).str = 'LJ'; NMD.system(2).str = 'superlattice';
%--------------------------------------------------------------------------
    NMD.x0.mass(1) = 1.0; NMD.x0.mass(2) = 3.0; NMD.NUM_ATOMS_TYPE = 2;
%    NMD.x0.alat= 6.31/4;
    if strcmp(NMD.system(2).str,'alloy')
    	NMD.x0.alloy_conc = 0.5;
    	NMD.x0.vm =...
        ( (1-NMD.x0.alloy_conc)*NMD.x0.mass(1) +...
        NMD.x0.alloy_conc*NMD.x0.mass(2) );
    	NMD.seed.alloy = 1;
    
    elseif strcmp(NMD.system(2).str,'superlattice')
    	NMD.x0.superlattice.period = [12.0  0.0  0.0 
                                      0.0  1.0  0.0 
                                      0.0  0.0  1.0];
	NMD.x0.atype(1).str = 'Ar1' ;  NMD.x0.atype(2).str = 'Ar2' ;
	NMD.x0.amass(1) = 39.948;      NMD.x0.amass(2) = 119.844;
	NMD.seed.superlattice = 1;
        NMD.x0.alat(1) = 6.2529/4*NMD.x0.superlattice.period(1,1);
        NMD.x0.alat(2) = 6.2529/4*NMD.x0.superlattice.period(2,2);
        NMD.x0.alat(3) = 6.2529/4*NMD.x0.superlattice.period(3,3);
    else    
    end
%--------------------------------------------------------------------------  
    NMD.env = 'gilgamesh'
	 
    if strcmp(NMD.env,'gilgamesh') 
    	NMD.lmp.path = '/home/jason/lammps/lammps-2Nov10/src';
    	NMD.lmp.exe = 'lmp_generic';
    	NMD.lmp.walltime = 24; 
    	NMD.lmp.cpu = 8; 
    	NMD.lmp.mem = 2;
    	NMD.matlab.path = '/opt/mcgaugheygroup/matlab_R2011a/bin';
    	NMD.matlab.lib = '/opt/mcgaugheygroup/matlab_R2011a/work';
    	NMD.gulp.matlab.lib = '/home/jason/Documents/MATLAB/nmd';
    	NMD.matlab.exe = 'matlab'
    	NMD.matlab.walltime = 24; 
    	NMD.matlab.cpu = 1; 
    	NMD.matlab.mem = 2;
	NMD.gulp.path = '/home/jason/GULP/gulp.4.0/Src';
    	NMD.gulp.exe = 'gulp_ifort';
    elseif strcmp(NMD.env,'ntpl')
    	NMD.lmp.path = '/home/jason/lammps/lammps-10Oct10/src';
    	NMD.lmp.exe = 'lmp_openmpi';
    	NMD.lmp.walltime = 24; 
    	NMD.lmp.cpu = 1; 
    	NMD.lmp.mem = 2;
    	NMD.matlab.path = '/home/jason/matlab/matlab_R2010a/bin';
    	NMD.matlab.lib = '/home/jason/matlab/matlab_R2010a/work';
    	%NMD.gulp.matlab.lib = '/home/jason/Documents/MATLAB/nmd';
	NMD.gulp.matlab.lib = NMD.str.main;
    	NMD.matlab.exe = 'matlab'
    	NMD.matlab.walltime = 24; 
    	NMD.matlab.cpu = 1; 
    	NMD.matlab.mem = 8;
	NMD.gulp.path = '/home/jason/GULP/Src';
    	NMD.gulp.exe = 'gulp';
	NMD.gulp.str.main = NMD.str.main;
    end 
%--------------------------------------------------------------------------
    NMD.x0.Nx = 4; NMD.x0.Ny = 6; NMD.x0.Nz = 6;
%--------------------------------------------------------------------------
    NMD.seed.initial = 1:1:5;
    NMD.seed.NUM_SEEDS = size(NMD.seed.initial,2);
%--------------------------------------------------------------------------
%SED PARAMETERS------------------------------------------------------------    

    NMD.x0.type(1).str ='SED'; NMD.x0.type(1).str ='XCORR';   
    NMD.x0.type(2).str ='DISP'; NMD.x0.type(2).str ='nGAMMA';
%--------------------------------------------------------------------------    
    NMD.t_total = 2^20; NMD.t_fft = 2^16; NMD.t_step = 2^5; NMD.dt = 0.002;
%--------------------------------------------------------------------------  
    NMD.NUM_MODESLICES = NMD.x0.Nx*NMD.x0.Ny*NMD.x0.Nz; 
%--------------------------------------------------------------------------
%Run GULP->LAMMPS->NMD
    NMD = m_nmd_structure_create( NMD );
%GULP
    NMD.gulp.matlab.lib = NMD.str.main;
    NMD.gulp.str.main = NMD.str.main;
    NMD.gulp = m_gulp_disp (NMD.gulp, NMD.x0);
    %NMD = nmd_gulp_submit(NMD);
%LAMMPS    
    m_nmd_lmp_create_sh(NMD)
    m_nmd_lmp_create_x0(NMD)
    m_nmd_lmp_create_in(NMD)
%NMD
    m_nmd_nmd_create(NMD)
%fit
    m_fit_create(NMD)


