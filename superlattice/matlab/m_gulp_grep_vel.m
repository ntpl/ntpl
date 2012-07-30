function vel = m_gulp_grep_vel( gulp , x0)

gulp.dk = 10E-5;

%2) Input dk kpts to measure group velocities
	vel = zeros(3*size(x0.ucell.cart,1),3);		

for idim = 1:3
    if gulp.kpt(idim)==0.5
        gulp = m_gulp_prepare( gulp , x0 )
        freq = m_gulp_grep_freq( gulp , x0 );
        if idim==1
		gulp.kpt(idim) = gulp.kpt(idim) - x0.superlattice.period(1,1)*gulp.dk;
        else
		gulp.kpt(idim) = gulp.kpt(idim) - gulp.dk;
        end
        gulp = m_gulp_prepare( gulp , x0 )
        freq_mdk = m_gulp_grep_freq( gulp , x0 );
        vel(:,idim) = (( freq - freq_mdk )/(gulp.dk))/4;
        %Put kpt back to orig
        if idim==1
                gulp.kpt(idim) = gulp.kpt(idim) +  x0.superlattice.period(1,1)*gulp.dk;
        else
                gulp.kpt(idim) = gulp.kpt(idim) + gulp.dk;
        end        

    elseif gulp.kpt(idim)==-0.5
        gulp = m_gulp_prepare( gulp , x0 )
        freq = m_gulp_grep_freq( gulp , x0 );
        if idim==1
                gulp.kpt(idim) = gulp.kpt(idim) +  x0.superlattice.period(1,1)*gulp.dk;
        else
                gulp.kpt(idim) = gulp.kpt(idim) + gulp.dk;
        end
        gulp = m_gulp_prepare( gulp , x0 )
        freq_pdk = m_gulp_grep_freq( gulp , x0 );
        vel(:,idim) = (( freq_pdk - freq )/(gulp.dk))/4;
        %Put kpt back to orig
        if idim==1
                gulp.kpt(idim) = gulp.kpt(idim) -  x0.superlattice.period(1,1)*gulp.dk;
        else
                gulp.kpt(idim) = gulp.kpt(idim) - gulp.dk;
        end

    elseif gulp.kpt(idim)==0.0
        gulp = m_gulp_prepare( gulp , x0 )
        freq = m_gulp_grep_freq( gulp , x0 );
        if idim==1
                gulp.kpt(idim) = gulp.kpt(idim) +  x0.superlattice.period(1,1)*gulp.dk;
        else
                gulp.kpt(idim) = gulp.kpt(idim) + gulp.dk;
        end
        gulp = m_gulp_prepare( gulp , x0 )
        freq_pdk = m_gulp_grep_freq( gulp , x0 );
        vel(:,idim) = (( freq_pdk - freq )/(gulp.dk))/4;
        %Put kpt back to orig
        if idim==1
                gulp.kpt(idim) = gulp.kpt(idim) - x0.superlattice.period(1,1)*gulp.dk;
        else
                gulp.kpt(idim) = gulp.kpt(idim) - gulp.dk;
        end
        
    else
        %freq = m_gulp_grep_freq( gulp , x0 );
        if idim==1
                gulp.kpt(idim) = gulp.kpt(idim) + x0.superlattice.period(1,1)*gulp.dk;
        else
                gulp.kpt(idim) = gulp.kpt(idim) + gulp.dk;
        end
        gulp = m_gulp_prepare( gulp , x0 )
        freq_pdk = m_gulp_grep_freq( gulp , x0 );
        if idim==1
                gulp.kpt(idim) = gulp.kpt(idim) -x0.superlattice.period(1,1)*2*gulp.dk;
        else
                gulp.kpt(idim) = gulp.kpt(idim) -2*gulp.dk;
        end
        gulp = m_gulp_prepare( gulp , x0 )
        freq_mdk = m_gulp_grep_freq( gulp , x0 );
      	vel(:,idim) = (( freq_pdk - freq_mdk )/(gulp.dk))/8;
    	%Put kpt back to orig
    	if idim==1
                gulp.kpt(idim) = gulp.kpt(idim) +  x0.superlattice.period(1,1)*gulp.dk;
        else
                gulp.kpt(idim) = gulp.kpt(idim) + gulp.dk;
        end
    end
end

end

