function gulp = m_gulp_disp( gulp, x0 ) 


%--------------------------------------------------------------------------
gulp.dk = 10E-5;
%--------------------------------------------------------------------------

%if exist([gulp.str.main '/eigvec.dat'], 'file') ~= 0
    %str.cmd=['rm -f ' gulp.str.main '/eigvec.dat'];
    str.cmd=['rm -f eigvec.dat'];
    system(str.cmd);
    %str.cmd=['rm -f ' gulp.str.main '/freq.dat'];
    str.cmd=['rm -f freq.dat'];
    system(str.cmd);
    %str.cmd=['rm -f ' gulp.str.main '/vel.dat'];
    str.cmd=['rm -f vel.dat'];
    system(str.cmd);
    str.cmd=['rm -f gulpkpt.dat'];
    system(str.cmd);
%end    

for ikpt=1:size(x0.kptlist,1)       
%only works for orthogonal lattice vector systems!
    gulp.kpt = x0.kptlist(ikpt,:)./[x0.Nx x0.Ny x0.Nz];
    %gulp.kpt = gulp.kpt./[x0.superlattice.period(1,1) x0.superlattice.period(2,2) x0.superlattice.period(3,3)]    
    dlmwrite('gulpkpt.dat',gulp.kpt,'-append','delimiter',' ');
    gulp = m_gulp_prepare( gulp , x0 );

eigvec =...
    m_gulp_grep_eig ( gulp, x0 );
freq =...
    m_gulp_grep_freq ( gulp, x0 );
vel =...
    m_gulp_grep_vel ( gulp, x0 );
%Output formatted eigvec		
    str.write=strcat(gulp.str.main,'/eigvec.dat');
    dlmwrite(str.write,eigvec,'-append','delimiter',' ');
%Output formatted freqs
    str.write=strcat(gulp.str.main,'/freq.dat');
    dlmwrite(str.write,freq,'-append','delimiter',' ');
%Output velocities
    str.write=strcat(gulp.str.main,'/vel.dat');
    dlmwrite(str.write,vel,'-append','delimiter',' ');
end

%end
