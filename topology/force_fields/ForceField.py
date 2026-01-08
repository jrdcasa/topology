import periodictable as pt


class ForceField(object):

    # ===========================================================================
    def __init__(self, ffname, log=None):

        self._logger = log
        self._allowed_forcefields = ["OPLSAA", "TRAPPEAA"]
        self._atomtypes_dict = None
        self._ff_param_filename = None
        self._ffname = ffname

        if ffname.upper() not in self._allowed_forcefields:
            m = "\n\t\t The force field {0:s} does not exist\n".format(ffname)
            m += "\t\t Available force fields {}\n".format(", ".join(self._allowed_forcefields))
            print(m) if self._logger is None else self._logger.error(m)
            exit()
        elif ffname.upper() == "OPLSAA":
            self._oplsaa()

    # ===========================================================================
    def __opls_parameters(self):

        # Number of exclusions
        self._nrexcl = 3

        # 1-4 interactions
        self._nbfunctiontype = 1   # 1 (Lennard-Jones) or 2 (Buckingham)
        self._combrule = 3         # 1 and 3 (Geometric rule) and 2 (Arithmetic rule)
        self._genpairs = "yes"     # Generate pairs (yes/no). If yes generates pairs not present in [pairtypes]
        self._fudge_lj = 0.5       # Factor to multiply LJ 1-4 interctions, default 1. Not used if genpairs=no
        self._fudge_coulomb = 0.5  # Factor to multiply electrostatic 1-4 interctions, default 1. Always used.

        # Name of the general parameters for the force field
        self._general_itp_name = "ff_oplsaa.itp"

        # Atomtype: [element, bond_type_name, mass[g/mol], charge, sigma(A), epsilon(kJ/mol)
        self._atomtypes_dict = {
            "opls_122": ["C" , "CT",  6,  12.01100,  0.248, 3.80, 0.209200],
            "opls_123": ["Cl", "Cl",  17, 35.45300, -0.062, 3.47, 1.112950],
            "opls_135": ["C" , "CT",  6,  12.01100, -0.180, 3.50, 0.276144],
            "opls_136": ["C" , "CT",  6,  12.01100, -0.120, 3.50, 0.276144],
            "opls_137": ["C" , "CT",  6,  12.01100, -0.060, 3.50, 0.276144],
            "opls_139": ["C" , "CT",  6,  12.01100,  0.000, 3.50, 0.276144],
            "opls_140": ["H" , "HC",  1,   1.00800,  0.060, 2.50, 0.125520],
            "opls_145": ["C" , "CA",  6,  12.01100, -0.115, 3.55, 0.292880],
            "opls_146": ["H" , "HA",  1,   1.00800,  0.115, 2.42, 0.125520],
            "opls_151": ["Cl", "Cl",  17, 35.45300, -0.200, 3.40, 1.255200],
            "opls_152": ["C" , "CT",  6,  12.01100, -0.006, 3.50, 0.276144],
            "opls_153": ["H" , "HC",  1,   1.00800,  0.103, 2.50, 0.125520],
            "opls_164": ["F" ,  "F",  9,  18.99840, -0.206, 2.94, 0.255224],
            "opls_263": ["C" , "CA",  6,  12.01100,  0.180, 3.55, 0.292880],
            "opls_264": ["Cl" ,"Cl", 17,  35.45300, -0.180, 3.40, 1.255200],
            "opls_280": ["C" , "C_2", 6,  12.01100, -0.470, 3.75, 0.439320],
            "opls_281": ["O" , "O_2", 8,  15.99940, -0.470, 2.96, 0.878640],
            "opls_722": ["Br", "Br", 35,  79.90400, -0.220, 3.47, 0.196648],
        }
        # Bonds types:
        #  ["1 Harmonic bond",  "b0 (nm)" , "kb (kJmol-1nm-2"]
        self._bondtypes_dict = {
            "CT-HC"  : [1,  0.10900, 284512.0 ],    # CT HC 1 0.10900 284512.0 ; CHARMM 22 parameter file
            "CT-Cl"  : [1,  0.17810, 205016.0 ],    # CT Cl 1 0.17810 205016.0 ; wlj- from MM2 (Tet 31, 1971 (75))
            "CT-CT"  : [1,  0.15290, 224262.4 ],    # CT CT 1 0.15290 224262.4 ; CHARMM 22 parameter file
            "CT-F"   : [1,  0.13320, 307105.6 ],    # CT  F 1 0.13320 307105.6 ; PAK CT-F for CHF3 (emd 5-09-94)
            "CT-Br"  : [1,  0.19450, 205016.0 ],    # CT Br 1 0.19450 205016.0 ; wlj
            "CT-C_2" : [1,  0.15220, 265265.6 ],    # C_2 CT 1 0.15220 265265.6;
            "C_2-O_2": [1,  0.12290, 476976.0 ],    # C_2 O_2 1 0.12290 476976.0;
            "CA-CA"  : [1,  0.14000, 392459.2 ],    # CA CA 1 0.14000 392459.2; TRP,TYR,PHE
            "CA-HA"  : [1,  0.10800, 307105.6 ],    # CA HA 1 0.10800 307105.6; PHE, etc.
            "CA-Cl"  : [1,  0.17250, 251040.0 ],    # CA Cl 1 0.17250 251040.0; wlj
        }

        # Angle types:
        #  ["1 Harmonic bond",  "theta0 (º)" , "kb (kJmol-1rad-2"]
        self._angletypes_dict = {
            "HC-CT-HC"  : [1, 107.800, 276.144 ],    # HC CT HC 1 107.800 276.144 ; CHARMM 22 parameter file
            "HC-CT-Cl"  : [1, 107.600, 426.768 ],    # HC CT Cl 1 107.600 426.768 ;
            "HC-CT-F"   : [1, 107.000, 334.720 ],    # HC CT F  1 107.000 334.720 ; wlj
            "HC-CT-Br"  : [1, 107.600, 426.768 ],    # HC CT Br 1 107.600 426.768 ; wlj
            "F-CT-Cl"   : [1, 110.580, 432.207 ],    # From LigParGen
            "F-CT-Br"   : [1, 110.580, 432.207 ],    # From LigParGen
            "Cl-CT-Br"  : [1, 110.580, 432.207 ],    # From LigParGen
            "Cl-CT-Cl"  : [1, 111.700, 652.704 ],    # Cl CT Cl 1 111.700 652.704 ;
            "CT-CT-HC"  : [1, 110.700, 313.800 ],    # CT CT HC 1 110.700 313.800 ; CHARMM 22 parameter file
            "CT-CT-CT"  : [1, 112.700, 488.273 ],    # CT CT CT 1 112.700 488.273 ; CHARMM 22 parameter file
            "CT-CT-Cl"  : [1, 109.800, 577.392 ],    # CT CT Cl 1 109.800 577.392 ; wlj - from MM2
            "CT-C_2-O_2": [1, 120.400, 669.440 ],    # CT C_2 O_2 1 120.400 669.440;
            "CT-C_2-CT" : [1, 116.000, 585.760 ],    # CT C_2 CT  1 116.000 585.760; wlj 7/96
            "HC-CT-C_2" : [1, 109.800, 577.392 ],    # C_2 CT HC  1 109.500 292.880;
            "CA-CA-CA"  : [1, 120.000, 527.184 ],    # CA CA CA 1 120.000 527.184; PHE(OL)
            "CA-CA-HA"  : [1, 120.000, 292.880 ],    # CA CA HA 1 120.000 292.880;
            "CA-CA-Cl"  : [1, 120.000, 627.600 ],    # CA CA Cl 1 120.000 627.600;
        }

        # Dihedral types:
        #  ["3 Ryckaert-Bellemans dihedral",  "C0, C1, C2, C3, C4, C5 (kJmol-1)"]
        self._dihedraltypes_dict = {
            # CT CT CT Cl 3 0.83680 2.51040 0.00000 -3.34720 0.00000 0.00000 ; alkyl chloride
            "CT-CT-CT-Cl": [3, [0.83680, 2.51040, 0.00000, -3.34720, 0.00000, 0.00000] ],
            # HC CT CT HC 3 0.62760 1.88280 0.00000 -2.51040 0.00000 0.00000 ; hydrocarbon *new* 11/99
            "HC-CT-CT-HC": [3, [0.62760, 1.88280, 0.00000, -2.51040, 0.00000, 0.00000] ],
            # CT CT CT HC 3 0.62760 1.88280 0.00000 -2.51040 0.00000 0.00000 ; hydrocarbon all-atom
            "HC-CT-CT-CT": [3, [0.62760, 1.88280, 0.00000, -2.51040, 0.00000, 0.00000] ],
            # Cl CT CT HC 3 0.83680 2.51040 0.00000 -3.34720 0.00000 0.00000 ; alkyl chloride
            "HC-CT-CT-Cl": [3, [0.83680, 2.51040, 0.00000, -3.34720, 0.00000, 0.00000] ],
            # CT CT CT CT 3 2.92880 -1.46440 0.20920 -1.67360 0.00000 0.00000 ; hydrocarbon all-atom
            "CT-CT-CT-CT": [3, [2.92880, -1.46440, 0.20920, -1.67360, 0.00000, 0.00000] ],
            #CT C_2 CT HC 3 0.57530 1.72590 0.00000 -2.30120 0.00000 0.00000 ; ketone
            "CT-C_2-CT-HC": [3, [0.57530, 1.72590, 0.00000, -2.30120, 0.00000, 0.00000] ],
            # HC CT C_2 O_2 3 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 ; aldehyde
            "HC-CT-C_2-O_2": [3, [0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000] ],
            # X CA CA X 3 30.33400 0.00000 -30.33400 0.00000 0.00000 0.00000; aromatic ring
            "CA-CA-CA-CA": [3, [30.33400, 0.00000, -30.33400, 0.00000, 0.00000, 0.00000] ],
            "CA-CA-CA-HA": [3, [30.33400, 0.00000, -30.33400, 0.00000, 0.00000, 0.00000] ],
            "CA-CA-CA-Cl": [3, [30.33400, 0.00000, -30.33400, 0.00000, 0.00000, 0.00000] ],
            "HA-CA-CA-Cl": [3, [30.33400, 0.00000, -30.33400, 0.00000, 0.00000, 0.00000] ],
            "HA-CA-CA-HA": [3, [30.33400, 0.00000, -30.33400, 0.00000, 0.00000, 0.00000] ],
        }

    # ===========================================================================
    def _oplsaa(self):

        # Load default parameters
        self.__opls_parameters()

        # Wrtie the default parameters
        self._ff_param_filename = "ff_oplsaa.itp"

        m = "\t\t Writting {0:s}".format(self._ff_param_filename+'\n')
        print(m) if self._logger is None else self._logger.info(m)

        lines_itp = "[atomtypes]\n"
        lines_itp += ";Parameters are taken from gromacs/top/oplsaa.ff and https://virtualchemistry.org/ff.php\n"
        lines_itp += "; name  bond_type atomic_number   mass    charge   ptype          sigma(A)      epsilon(kJ/mol)\n"

        for key, item in self._atomtypes_dict.items():
            iline = "{0:10s}".format(key)                      # Atom type
            iline += "{0:2s}".format(item[1])                  # Bond type
            iline += "{0:3d}".format(item[2])                  # Atomic number
            iline += "  {0:10.6f}".format(item[3])             # Mass (g/mol)
            iline += "  {0:10.6f}".format(item[4])             # charge
            iline += "  A    {0:10.6f}".format(item[5]/10.0)   # sigma(nm)
            iline += " {0:10.6f}".format(item[6])              # epsilon(kJ/mol)
            lines_itp += iline+"\n"

        with open(self._ff_param_filename, 'w') as ff_file:
            ff_file.writelines(lines_itp)
