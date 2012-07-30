function gulp = m_gulp_prepare (gulp, x0)

if isfield(x0,'alloy')
    if x0.alloy.conc == 0.0
    gulp = gulp_disp_crystal_alloy_single_stringchanges( gulp, x0 );
    change_file_strings(...
        gulp.filename , gulp.orig, strcat(gulp.str.main,'/disp.gin'), gulp.change)
    else
        
    end
elseif isfield(x0,'superlattice')
    %gulp = gulp_disp_crystal_alloy_single_stringchanges( gulp, x0 );
    gulp = m_gulp_superlattice_stringchanges( gulp, x0 );
	gulp.filename 
	gulp.orig
	strcat(gulp.str.main,'/disp.gin') 
	gulp.change
    m_change_file_strings(...
        gulp.filename , gulp.orig, strcat(gulp.str.main,'/disp.gin') , gulp.change)


elseif isfield(x0,'crystal')

end
    
end
