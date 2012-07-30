function nmd_lmp_create_in( NMD )

	%seeds1{:}=num2str(ceil(abs(100000*randn(1,NMD.seed.NUM_SEEDS))));
	%seeds = sprintf('%s \t',seeds1{:})
     for l=1:1:NMD.seed.NUM_SEEDS
	 seeds1{:}=num2str(ceil(abs(100000*randn(1,1))));
	 seeds = sprintf('%s \t',seeds1{:})
%lmp.in.sed.iseed
        str.orig = 'IN.X0';
        str.change = ['lmp.in.x0.' NMD.system(2).str];
        str.cmd1 = ['-e ''s/\<' str.orig '\>/' str.change '/g'' '];
        str.orig = 'LMP_TEMP';
        str.change = ['lmp.in.sed.' NMD.system(2).str '.' int2str(l)];
        str.cmd2 = ['-e ''s/\<' str.orig '\>/' str.change '/g'' '];
        str.orig = 'IX0';
        str.change = [NMD.system(2).str];
        str.cmd3 = ['-e ''s/\' str.orig '>/' str.change '/g'' '];
        str.orig = 'ISEED_TMP';
        str.change = [int2str(l)];
        str.cmd4 = ['-e ''s/\<' str.orig '\>/' str.change '/g'' '];
        str.orig = 'SEED_TMP';
        str.change = [seeds];
        str.cmd5 = ['-e ''s/\<' str.orig '\>/' str.change '/g'' '];
        str.orig = 'T_STEP';
        str.change = [num2str(NMD.t_step)];
        str.cmd6 = ['-e ''s/\<' str.orig '\>/' str.change '/g'' '];
        str.orig = 'T_FFT';
        str.change = [num2str(NMD.t_fft)];
        str.cmd7 = ['-e ''s/\<' str.orig '\>/' str.change '/g'' '];
        str.orig = 'T_TOTAL';
        str.change = [num2str(NMD.t_total)];
        str.cmd8 = ['-e ''s/\<' str.orig '\>/' str.change '/g'' '];
        
        str.cmd9 =...
            ['lmp.in.sed.temp > ./' int2str(NMD.seed.superlattice ) '/lmp.in.sed.' NMD.system(2).str '.' int2str(l)];
        
        str.cmd =...
            ['sed ' str.cmd1 str.cmd2 str.cmd3 str.cmd4 str.cmd5...
            str.cmd6 str.cmd7 str.cmd8 str.cmd9 ];       
        system(str.cmd);
     end

    str.cmd = ['chmod +x ' int2str(NMD.seed.superlattice) '/lmp_submit.sh'];
    system(str.cmd);
    str.cmd = [int2str(NMD.seed.superlattice) '/lmp_submit.sh'];
    system(str.cmd);
end
