function gulp = m_gulp_superlattice_stringchanges( gulp, x0 )

gulp.filename = strcat(gulp.matlab.lib,'/gulp_lj_superlattice.gin');

gulp.orig(1).str = 'EPSILON';
gulp.change(1).str = num2str(x0.LJ.eps);
gulp.orig(2).str = 'SIGMA';
gulp.change(2).str = num2str(x0.LJ.sigma);
gulp.orig(3).str = 'ATOM_MASS_1';
gulp.change(3).str = num2str(x0.amass(1));
gulp.orig(4).str = 'ATOM_MASS_2';
gulp.change(4).str = num2str(x0.amass(2));
if strcmp(x0.type(2).str,'GAMMA')
  gulp.orig(5).str = 'ALATX';
  gulp.change(5).str = num2str(x0.LJ.sigma.*(1e10).*x0.alat(1)*x0.superlattice.period(1,1)*x0.Nx);
  gulp.orig(6).str = 'ALATY';
  gulp.change(6).str = num2str(x0.LJ.sigma.*(1e10).*x0.alat(1)*x0.superlattice.period(2,2)*x0.Ny);
  gulp.orig(7).str = 'ALATZ';
  gulp.change(7).str = num2str(x0.LJ.sigma.*(1e10).*x0.alat(1)*x0.superlattice.period(3,3)*x0.Nz);
else
  gulp.orig(5).str = 'ALATX';
  gulp.change(5).str = num2str(x0.LJ.sigma.*(1e10).*x0.alat(1));%*x0.superlattice.period(1,1));
  gulp.orig(6).str = 'ALATY';
  gulp.change(6).str = num2str(x0.LJ.sigma.*(1e10).*x0.alat(2));%*x0.superlattice.period(2,2));
  gulp.orig(7).str = 'ALATZ';
  gulp.change(7).str = num2str(x0.LJ.sigma.*(1e10).*x0.alat(3));%*x0.superlattice.period(3,3));
end

gulp.orig(8).str = 'KPT';
gulp.change(8).str = num2str(gulp.kpt);


gulp = m_coords2gulp( gulp, x0 )

gulp.orig(9).str = 'COORDS';
gulp.change(9).str = gulp.coords;
end
