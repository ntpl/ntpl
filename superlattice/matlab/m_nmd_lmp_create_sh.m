function m_nmd_lmp_create_sh(NMD)


%for NMD.seed.superlattice =1:size(NMD.seed.superlattice ,2)

    str.cmd = ['rm -r ./' int2str(NMD.seed.superlattice )];
    system(str.cmd);

    str.cmd = ['mkdir -p ./' int2str(NMD.seed.superlattice ) '/NMD'];
    system(str.cmd);

    str.cmd = ['mkdir -p ./' int2str(NMD.seed.superlattice )];
    system(str.cmd);
    
    %str.cmd = ['cp sh.nmd_lmp.sh ./' int2str(NMD.seed.superlattice ) '/lmp_submit.sh'];
    %system(str.cmd);    


%loops over initial seeds
    for iseed=1:size(NMD.seed.initial,2)
        
%lmp_ISEED.sh------------------------------------------------------        
        str.orig = 'lmp.sh.temp';
        str.change = ['lmp' int2str(iseed) '.sh'];
        str.cmd1 = ['-e ''s/\<' str.orig '\>/' str.change '/g'' '];
        str.orig = 'runpath';
        str.change = strcat(NMD.str.main,'/',int2str(NMD.seed.superlattice ));
        str.temp = strcat('-e ''s|',str.orig,'|',str.change);
        str.cmd2 = [str.temp '|g'' '];
        str.orig = 'LMP.TEMP';
        str.change = ['lmp.in.sed.' NMD.system(2).str '.' int2str(iseed)];
        str.cmd3 = ['-e ''s/\<' str.orig '\>/' str.change '/g'' '];
        str.orig = 'lmp_temp';
        str.change = ['lmp' int2str(iseed)];
        str.cmd4 = ['-e ''s/\<' str.orig '\>/' str.change '/g'' '];
    
       str.cmd5 = ['sh.nmd_lmp.sh > ./' int2str(NMD.seed.superlattice ) '/lmp' int2str(iseed) '.sh'];
    
    str.cmd = ['sed ' str.cmd1 str.cmd2 str.cmd3 str.cmd4 str.cmd5];
    
        system(str.cmd);
               

        
 	%lmp_submit.sh-------------------------------------------------------------
     if strcmp(NMD.env,'gilgamesh') 
    	output =...
        ['qsub -l walltime=' int2str(NMD.lmp.walltime)...
        ':00:00 -l nodes=1:ppn=' int2str(NMD.lmp.cpu)...
        ' lmp' int2str(iseed) '.sh'];
    	str.write = strcat(NMD.str.main,'/',int2str(NMD.seed.superlattice ),'/lmp_submit.sh');
    	dlmwrite(str.write,output,'-append','delimiter','');
     elseif strcmp(NMD.env,'ntpl')
	str.tmp=strcat(NMD.str.main,'/')
 	output =...
        ['qsub ' str.tmp int2str(NMD.seed.superlattice ) '/lmp' int2str(iseed) '.sh'];
    	str.write = strcat(NMD.str.main,'/',int2str(NMD.seed.superlattice ),'/lmp_submit.sh');
    	dlmwrite(str.write,output,'-append','delimiter','');
     end
    end 
%end
end
