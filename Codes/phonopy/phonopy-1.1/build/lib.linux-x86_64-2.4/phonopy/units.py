# Copyright (C) 2011 Atsushi Togo
#
# This file is part of phonopy.
#
# Phonopy is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Phonopy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with phonopy.  If not, see <http://www.gnu.org/licenses/>.

from math import pi, sqrt

kb_J = 1.3806504e-23 # [J/K]
PlanckConstant = 4.13566733e-15 # [eV s]
Hbar = PlanckConstant/(2*pi) # [eV s]
Avogadro = 6.02214179e23
SpeedOfLight = 299792458 # [m/s]
AMU = 1.6605402e-27 # [kg]
Newton = 1.0        # [kg m / s^2]
Joule = 1.0         # [kg m^2 / s^2]
EV = 1.60217733e-19 # [J]
Angstrom = 1.0e-10  # [m]
THz = 1.0e12        # [/s]
Mu0 = 4.0e-7 * pi
Epsilon0 = 1.0 / Mu0 / SpeedOfLight**2
Me = 9.10938215e-31

Bohr = 4e10 * pi * Epsilon0 * Hbar**2 / Me  # Bohr radius [A] 0.5291772
Hartree = Me * EV / 16 / pi**2 / Epsilon0**2 / Hbar**2 # Hartree [eV] 27.211398
Rydberg = Hartree / 2 # Rydberg [eV] 13.6056991

THzToEv = PlanckConstant * 1e12 # [eV]
Kb = kb_J / EV  # [eV/K] 8.6173383e-05
THzToCm = 1.0e12 / (SpeedOfLight * 100) # [cm^-1] 33.356410
CmToEv = THzToEv / THzToCm # [eV] 1.2398419e-4
VaspToEv = sqrt(EV/AMU)/Angstrom/(2*pi)*PlanckConstant # [eV] 6.46541380e-2
VaspToTHz = sqrt(EV/AMU)/Angstrom/(2*pi)/1e12 # [THz] 15.633302
VaspToCm =  VaspToTHz * THzToCm # [cm^-1] 521.47083
EvTokJmol = EV / 1000 * Avogadro # [kJ/mol] 96.4853910
Wien2kToTHz = sqrt(Rydberg/1000*EV/AMU)/(Bohr*1e-10)/(2*pi)/1e12 # [THz] 3.44595837
