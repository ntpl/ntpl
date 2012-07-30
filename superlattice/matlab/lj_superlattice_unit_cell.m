function x0 = lj_superlattice_unit_cell( x0 )

x0.latvec =        [1.0  0.0  0.0 
                    0.0  1.0  0.0 
                    0.0  0.0  1.0];

x0.ucell.direct =      [0.0  0.0  0.0
                        0.5  0.5  0.0
                        0.5  0.0  0.5
                        0.0  0.5  0.5];

x0.superlattice.direct = m_build_supercell(x0.ucell.direct,x0.latvec,...
			x0.superlattice.period(1,1), x0.superlattice.period(2,2),...
			x0.superlattice.period(3,3))

%x0.ucell.cart= x0.alat.*x0.superlattice.direct

%Set m vector to create mass difference interface
Imass = find(x0.superlattice.direct(:,1)<(x0.superlattice.period(1,1)/2));
x0.m(Imass) = 1;
Imass = find(x0.superlattice.direct(:,1)>=(x0.superlattice.period(1,1)/2));
x0.m(Imass) = 2;

x0.latvec_rec =     [1.0 0.0  0.0 
                    0.0  1.0  0.0 
                    0.0  0.0  1.0];
end
