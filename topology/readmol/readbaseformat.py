from collections import defaultdict
import datetime
import numpy as np
import os


class ReadBaseFormat(object):

    __slots__ = ['_fnamepath', '_fname', '_natoms', '_nbonds',
                 '_atom3d_element', '_atom3d_imol', '_atom3d_molname',
                 '_nmols', '_mol_residue_list', '_atom3d_map', '_atom3d_idxsd',
                 '_atom3d_parent', '_atom3d_isbackbone', '_atom3d_residue',
                 '_atom3d_resname', '_atom3d_charge', '_atom3d_mass', '_unitcell',
                 '_boxlength', '_boxangle', '_assign_bo', '_topology', '_spacegroup',
                 '_zvalue', '_nres', '_isthere_boxdimension', '_nx', '_ny', '_nz',
                 '_universe', '_atom3d_xyz', '_atom3d_bfactor', '_atom3d_occupancy',
                 '_heads', '_tails', '_atom3d_kindmolecule', '_dimensions']

    # *************************************************************************
    def __init__(self, filenamepath, assign_bondorders=False):

        """
        Read XSD format files creating some useful information for the
        MolecularGraph

        Args:
            filenamepath (str): Path and name of the file for the xsd file

        Members:
            | ``self._filenamepath`` -- Path and name of the file for the xsd file
            | ``self._file`` -- Name of the xsd file
            | ``self._nmols`` -- Number of molecules in the system
            | ``self._nres`` -- Number of residues in the system
            | ``self._natoms`` -- Number of atoms in the system
            | ``self._molreslist`` -- "molreslist  : [ [mol1 name residues] [mol2 name residues], ...]
                                    for molecular graphs
            | ``self._atomgrouplist`` -- "molreslist  :
                     self._atomgrouplist[index_atom (int)]:
                                          [
                                              Atom index,
                                              Local index in the residue,
                                              Name of the atom in the system,
                                              Type of atom in the FF,
                                              Mass atom,
                                              Atom charge,
                                              imol number,
                                              Global number of the residue,
                                              Local number of the residue,
                                              Is or not backbone atom,
                                          ]
            | ``self.map_idx_to_id_dict`` -- dict  : Map the index atomic number to the ID label in the xsd file.
                This map can be used to get the bonds

            | ``self._unitcell`` -- array [3,3] : Contains the box vectors in angstroms
                ([[ a_xx a_xy a_xz ]
                  [ 0.0  b_yy b_yz ]
                  [ 0.0  0.0  c_zz ]]
            | ``self._boxlengths`` -- array [3] : Lengths of the simulation box in angstroms
            | ``self._boxangles`` -- array [3] : Angles of the simulation box in radians
            | ``self._topology`` -- Topology : Topology of the system

        """

        self._fnamepath = filenamepath
        self._fname = filenamepath.split("/")[-1]

        self._natoms = 0
        self._nbonds = 0
        self._topology = None
        self._universe = None
        self._atom3d_xyz = defaultdict(list)
        self._atom3d_element = defaultdict()
        self._atom3d_imol = defaultdict()
        self._atom3d_molname = defaultdict()
        self._atom3d_map = defaultdict()
        self._atom3d_idxsd = defaultdict()
        self._atom3d_parent = defaultdict()
        self._atom3d_isbackbone = defaultdict()
        self._atom3d_residue = defaultdict()
        self._atom3d_resname = defaultdict()
        self._atom3d_charge = defaultdict(float)
        self._atom3d_mass = defaultdict(float)
        self._nmols = 0
        self._nres = 0
        self._mol_residue_list = []
        self._unitcell = np.zeros([3, 3], dtype=np.float32)
        self._boxlength = np.zeros(3, dtype=np.float32)
        self._boxangle = np.zeros(3, dtype=np.float32)
        self._assign_bo = assign_bondorders
        self._spacegroup = ''
        self._zvalue = 1
        self._isthere_boxdimension = False
        self._nx = 1
        self._ny = 1
        self._nz = 1

    # *************************************************************************
    def get_natoms(self):

        return self._natoms

    # *************************************************************************
    def get_nbonds(self):

        return self._nbonds

    # *************************************************************************
    def get_nmols(self):

        return self._nmols

    # *************************************************************************
    def get_nres(self):

        return self._nres

    # *************************************************************************
    def get_residuelist(self):

        return self._mol_residue_list

    # *************************************************************************
    def get_charges(self):

        """
        Return a dictionary ot atomistic charges

        :return: dictionary
        """

        return self._atom3d_charge

    # *************************************************************************
    def write_xyz(self, filename_xyz="test.xyz"):

        """
        Write a xyz file to check the structure
        """

        with open(filename_xyz, 'w') as fxyz:

            fxyz.write("{}\n".format(int(self.get_natoms())))
            fxyz.write("MOL\n")

            for idx in range(self._natoms):
                line = "{0:s} {1:.5f} {2:.5f} {3:.5f}\n".format(self._atom3d_element[idx],
                                                                self._atom3d_xyz[idx][0],
                                                                self._atom3d_xyz[idx][1],
                                                                self._atom3d_xyz[idx][2])
                fxyz.write(line)

    # *************************************************************************
    def write_gro(self, filename_gro="test.gro"):

        """
        Write a gro file to check the structure
        """

        fmt = {
            'n_atoms': "{0:5d}\n",  # number of atoms
            # coordinates output format, see http://chembytes.wikidot.com/g-grofile
            'xyz': "{resid:>5d}{resname:<5.5s}{name:>5.5s}{index:>5d}{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}\n",
            # unitcell
            'box_orthorhombic': "{0:10.5f} {1:9.5f} {2:9.5f}\n",
            'box_triclinic': "{0:10.5f} {1:9.5f} {2:9.5f} "
                             "{3:9.5f} {4:9.5f} {5:9.5f} "
                             "{6:9.5f} {7:9.5f} {8:9.5f}\n"
        }

        with open(filename_gro, 'w') as fgro:

            fgro.write('Written by MStoLAGRO (J.Ramos) \n')
            fgro.write(fmt['n_atoms'].format(int(self.get_natoms())))

            for idx in range(self._natoms):

                try:
                    resname = "{}".format(self._atom3d_resname[idx][0:3])
                except KeyError:
                    resname = "{}".format(self._atom3d_molname[idx][0:3])

                fgro.write(fmt['xyz'].format(
                    resid=self._atom3d_residue[idx],
                    resname="{}".format(resname),
                    index=idx+1,
                    name=self._atom3d_element[idx],
                    pos=[0.1*i for i in self._atom3d_xyz[idx]]
                ))

            # self._unitcell =                  GROMACS box vectors
            #       ([[ a_xx a_xy a_xz ]  Transpose   ([[ v1(x) 0.0   0.0   ]
            #         [ 0.0  b_yy b_yz ]   =======>     [ v2(x) v2(y) 0.0   ]
            #         [ 0.0  0.0  c_zz ]]               [ v3(x) v3(y) v3[z] ]]
            #
            #    v1(x) --> [0,0]; v2(y)--> [1,1]; v3(z)--> [2,2]
            #    v1(y) --> [1,0]; v1(z)--> [2,0]; v2(x)--> [0,1]
            #    v2(z) --> [2,1]; v3(x)--> [0,3]; v3(z)--> [1,2]
            #   Gromacs v1(y) = v1(z) = v2(z) = 0.0
            if self._unitcell[1, 0] != 0.0 or self._unitcell[2, 0] != 0.0 or self._unitcell[2, 1] != 0.0:
                fgro.write(fmt['box_triclinic'].format(self._unitcell[0, 0]*0.1,
                                                       self._unitcell[1, 1]*0.1,
                                                       self._unitcell[2, 2]*0.1,
                                                       self._unitcell[1, 0]*0.1,
                                                       self._unitcell[2, 0]*0.1,
                                                       self._unitcell[0, 1]*0.1,
                                                       self._unitcell[2, 1]*0.1,
                                                       self._unitcell[0, 2]*0.1,
                                                       self._unitcell[1, 2]*0.1))
            else:
                fgro.write(fmt['box_orthorhombic'].format(self._unitcell[0, 0]*0.1,
                                                          self._unitcell[1, 1]*0.1,
                                                          self._unitcell[2, 2]*0.1))

    # *************************************************************************
    def write_pdb(self, filename_pdb="test.pdb", separate_chains=False,
                  renumber=False, head_idx_atom=0, atom_kind_molecule_label=None):

        """
        Write a pdb file to check the structure.
        Adapted from MDAnalysis software (https://www.mdanalysis.org/)
        """

        fmt = {
            'ATOM': (
                "ATOM  {serial:5d} {name:<4s}{altLoc:<1s}{resName:<4s}"
                "{chainID:1s}{resSeq:4d}{iCode:1s}"
                "   {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}{occupancy:6.2f}"
                "{tempFactor:6.2f}      {segID:<4s}{element:>2s}\n"),
            'HETATM': (
                "HETATM{serial:5d} {name:<4s}{altLoc:<1s}{resName:<4s}"
                "{chainID:1s}{resSeq:4d}{iCode:1s}"
                "   {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}{occupancy:6.2f}"
                "{tempFactor:6.2f}      {segID:<4s}{element:>2s}\n"),
            'REMARK': "REMARK     {0}\n",
            'COMPND': "COMPND    {0}\n",
            'HEADER': "HEADER    {0}\n",
            'TITLE': "TITLE     {0}\n",
            'MODEL': "MODEL     {0:>4d}\n",
            'NUMMDL': "NUMMDL    {0:5d}\n",
            'ENDMDL': "ENDMDL\n",
            'END': "END\n",
            'CRYST1': ("CRYST1{0:9.3f}{1:9.3f}{2:9.3f}"
                       "{3:7.2f}{4:7.2f}{5:7.2f} "
                       "{6:<11s}{7:4d}\n"),
            'CONECT': "CONECT{0}\n"
        }

        spacegroup = "P -1"
        zvalue = 1
        radtodeg = 180/np.pi
        with open(filename_pdb, 'w') as fpdb:

            fpdb.write(fmt['REMARK'].format('Renumbered with topology (J.Ramos)'))
            fpdb.write(fmt['REMARK'].format('Backbone atoms --> Beta field == 0.00'))
            fpdb.write(fmt['REMARK'].format('Branch   atoms --> Beta field >= 1.00'))
            fpdb.write(fmt['REMARK'].format('Head     atoms --> Occupanccy field == 1.00'))
            fpdb.write(fmt['REMARK'].format('Tail     atoms --> Occupanccy field == 2.00'))
            fpdb.write(fmt['REMARK'].format('Middle   atoms --> Occupanccy field == 0.00'))
            fpdb.write(fmt['CRYST1'].format(self._boxlength[0], self._boxlength[1], self._boxlength[2],
                                            self._boxangle[0]*radtodeg,
                                            self._boxangle[1]*radtodeg,
                                            self._boxangle[2]*radtodeg,
                                            spacegroup, zvalue))

            for idx in range(self._natoms):

                try:
                    resname = "{}".format(self._atom3d_resname[idx][0:3])
                except KeyError:
                    resname = "{}".format(self._atom3d_molname[idx][0:3])

                if atom_kind_molecule_label is None:
                    fpdb.write(fmt['HETATM'].format(
                        serial=idx+1,
                        name=self._atom3d_element[idx],
                        altLoc=" ",
                        resName=resname,
                        chainID=" ",
                        resSeq=self._atom3d_residue[idx],
                        iCode=" ",
                        pos=[i for i in self._atom3d_xyz[idx]],
                        occupancy=self._atom3d_occupancy[idx],
                        tempFactor=self._atom3d_isbackbone[idx],
                        segID="    ",
                        element=self._atom3d_element[idx]
                    ))
                else:
                    fpdb.write(fmt['HETATM'].format(
                        serial=idx+1,
                        name=self._atom3d_element[idx],
                        altLoc=" ",
                        resName=resname,
                        chainID=atom_kind_molecule_label[idx],
                        resSeq=self._atom3d_residue[idx],
                        iCode=" ",
                        pos=[i for i in self._atom3d_xyz[idx]],
                        occupancy=self._atom3d_occupancy[idx],
                        tempFactor=self._atom3d_isbackbone[idx],
                        segID="    ",
                        element=self._atom3d_element[idx]
                    ))

            fpdb.write('END\n')

            for idx in range(self._natoms):

                fpdb.write('CONECT{0:>5d}'.format(idx+1))
                for ineigh in self._topology._graphdict[idx]:
                    fpdb.writelines('{0:>5d}'.format(ineigh+1))
                fpdb.writelines("\n")
                if self._topology.elements[idx] == 'H' and len(self._topology._graphdict[idx]) >= 2:
                    print("WARNING: Hydrogen atom with more than 1 atom connected. idx: {} ({})"
                          .format(idx, len(self._topology._graphdict[idx])))
                elif self._topology.elements[idx] == 'C' and len(self._topology._graphdict[idx]) >= 5:
                    print("WARNING: Carbon atom with more than 4 atoms connected. idx: {} ({})"
                          .format(idx, len(self._topology._graphdict[idx])))
                elif len(self._topology._graphdict[idx]) > 5:
                    print("WARNING: {} atom with more than 1 atom connected. idx: {} ({})".
                          format(self._topology.elements[idx], idx, len(self._topology._graphdict[idx])-1))

        if separate_chains:
            self._write_separate_chains(filename_pdb)

    # *************************************************************************
    def _write_separate_chains(self, fpdbname):

        pattern = os.path.splitext(fpdbname)[0]

        fmt = {
            'ATOM': (
                "ATOM  {serial:5d} {name:<4s}{altLoc:<1s}{resName:<4s}"
                "{chainID:1s}{resSeq:4d}{iCode:1s}"
                "   {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}{occupancy:6.2f}"
                "{tempFactor:6.2f}      {segID:<4s}{element:>2s}\n"),
            'HETATM': (
                "HETATM{serial:5d} {name:<4s}{altLoc:<1s}{resName:<4s}"
                "{chainID:1s}{resSeq:4d}{iCode:1s}"
                "   {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}{occupancy:6.2f}"
                "{tempFactor:6.2f}      {segID:<4s}{element:>2s}\n"),
            'REMARK': "REMARK     {0}\n",
            'COMPND': "COMPND    {0}\n",
            'HEADER': "HEADER    {0}\n",
            'TITLE': "TITLE     {0}\n",
            'MODEL': "MODEL     {0:>4d}\n",
            'NUMMDL': "NUMMDL    {0:5d}\n",
            'ENDMDL': "ENDMDL\n",
            'END': "END\n",
            'CRYST1': ("CRYST1{0:9.3f}{1:9.3f}{2:9.3f}"
                       "{3:7.2f}{4:7.2f}{5:7.2f} "
                       "{6:<11s}{7:4d}\n"),
            'CONECT': "CONECT{0}\n"
        }

        spacegroup = "P -1"
        zvalue = 1
        radtodeg = 180/np.pi

        ich_current = -1
        for idx in range(self._natoms):
            ich = self._atom3d_imol[idx]
            if ich != ich_current:
                ich_current = ich
                idx_stride = idx
                try:
                    fpdb.write('END\n')
                    fpdb.write(line_connect)
                    line_connect = ""
                    idx_local = 0
                    fpdb.close()
                except UnboundLocalError:
                    idx_local = 0
                    line_connect = ""
                    pass
                fpdb = open("{0:s}_{1:04d}.pdb".format(pattern, ich_current), 'w')
                fpdb.write(fmt['REMARK'].format('Written by MStoLAGRO (J.Ramos). Molecule {0:4d}'.format(ich_current)))
                fpdb.write(fmt['CRYST1'].format(self._boxlength[0], self._boxlength[1], self._boxlength[2],
                                                self._boxangle[0] * radtodeg,
                                                self._boxangle[1] * radtodeg,
                                                self._boxangle[2] * radtodeg,
                                                spacegroup, zvalue))
            fpdb.write(fmt['HETATM'].format(
                serial=idx_local + 1,
                name=self._atom3d_element[idx],
                altLoc=" ",
                resName="{}".format(self._atom3d_molname[idx][0:3]),
                chainID=" ",
                resSeq= 1,
                iCode=" ",
                pos=[i for i in self._atom3d_xyz[idx]],
                occupancy=1.0,
                tempFactor=self._atom3d_isbackbone[idx],
                segID="    ",
                element=self._atom3d_element[idx]
            ))

            line_connect += 'CONECT{0:>5d}'.format(idx_local + 1)
            for ineigh in self._topology._graphdict[idx]:
                line_connect += '{0:>5d}'.format((ineigh -idx_stride) + 1)
            line_connect += "\n"

            idx_local += 1

        # Last chain
        fpdb.write('END\n')
        fpdb.write(line_connect)
        fpdb.close()

    # *************************************************************************
    @staticmethod
    def lengths_and_angles_to_box_vectors(a_length, b_length, c_length, alpha, beta, gamma):

        """
         Convert from the lengths/angles of the unit cell to the box
            vectors (Bravais vectors). The angles should be in degrees.

            Parameters
            ----------
            a_length : scalar or np.ndarray
                length of Bravais unit vector **a** (angstroms)
            b_length : scalar or np.ndarray
                length of Bravais unit vector **b** (angstroms)
            c_length : scalar or np.ndarray
                length of Bravais unit vector **c** (angstroms)
            alpha : scalar or np.ndarray
                angle between vectors **b** and **c**, in radians.
            beta : scalar or np.ndarray
                angle between vectors **c** and **a**, in radians.
            gamma : scalar or np.ndarray
                angle between vectors **a** and **b**, in radians.

        """
        if np.all(alpha > 2*np.pi) and np.all(beta > 2*np.pi) and np.all(gamma > 2*np.pi):
            print('ERROR: All your angles were greater than 2*pi. Did you accidentally give me degrees?')
            exit()
        #
        # alpha = alpha * np.pi / 180.
        # beta = beta * np.pi / 180.
        # gamma = gamma * np.pi / 180.

        a = np.array([a_length, np.zeros_like(a_length), np.zeros_like(a_length)])
        b = np.array([b_length * np.cos(gamma), b_length * np.sin(gamma), np.zeros_like(b_length)])
        cx = c_length * np.cos(beta)
        cy = c_length * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
        cz = np.sqrt(c_length * c_length - cx * cx - cy * cy)
        c = np.array([cx, cy, cz])

        if not a.shape == b.shape == c.shape:
            raise TypeError('Shape is messed up.')

        # Make sure that all vector components that are _almost_ 0 are set exactly
        # to 0
        tol = 1e-6
        a[np.logical_and(a > -tol, a < tol)] = 0.0
        b[np.logical_and(b > -tol, b < tol)] = 0.0
        c[np.logical_and(c > -tol, c < tol)] = 0.0

        return a.T, b.T, c.T

    # *************************************************************************
    def write_psf(self, filename_psf="test.psf", improper_angles=None):

        """
        Write a PSF file using the topology object.
        Improper angles must be given explicitly as a list:
            [[i1, j1, k1, l1], [i2, j2, k2, l2], ...]
        """

        with open(filename_psf, 'w') as f:

            # Title section
            f.write("PSF\n")
            f.write("\n")
            line = str("%8d %s\n" % (2, '!NTITLE'))
            f.write(line)
            # TODO
            line =  "* File created by Topology code *\n"
            line +=  "* by Javier Ramos (IEM-CSIC) *\n"
            now = datetime.datetime.now()
            line += "* DATE: {}\n".format(now)
            f.write(line)
            f.write("\n")

            # Natoms section
            line = str("%8d %s\n" % (self._natoms, '!NATOM'))
            f.write(line)
            for i in range(self._natoms):
                # If self._isbackbone is empty all atoms are labelled as backbones
                # Atom backbone --> 0
                # Atom no backbone --> >0
                try:
                    if self._atom3d_isbackbone[i]:
                        tmp = '           0'
                    else:
                        tmp = '           1'
                except IndexError:
                    tmp = '           0'
                if len(self._topology._types) == 0:
                    line = str("{0:8d} {1:>4d}  {2:<3d} {3:>3s}  {4:<3s}  {5:<3s}  {6:10.6f}      {7:8.4f}{8:12s}\n"
                           .format(i + 1,
                            self._topology._iatch[i]+1,
                            self._topology._iatch[i]+1,
                            self._topology._resname[i][0:3],
                            self._topology.elements[i],
                            self._topology.elements[i],
                            self._topology.charge[i],
                            self._topology.mass[i],
                            tmp))
                else:
                    line = str("{0:8d} {1:>4d}  {2:<3d} {3:>3s}  {4:<3s}  {5:<3s}  {6:10.6f}      {7:8.4f}{8:12s}\n"
                           .format(i + 1,
                            self._topology._iatch[i]+1,
                            self._topology._iatch[i]+1,
                            self._topology._resname[i][0:3],
                            self._topology.elements[i],
                            self._topology._types[i],
                            self._topology.charge[i],
                            self._topology.mass[i],
                            tmp))
                f.write(line)
            f.write("\n")

            # NBonds section ======================================================
            nbonds = len(self._topology.get_edges())
            line = str("{0:8d} {1:s}\n".format(nbonds, '!NBOND: bonds'))
            f.write(line)
            iline = 1
            for i in range(self._natoms):
                ibond = self._topology.find_all_paths_length(i, 1)
                for j in range(len(ibond)):
                    iat0 = ibond[j][0]
                    iat1 = ibond[j][1]
                    if iat0 < iat1:
                        continue
                    line = str("{0:8d}{1:8d}").format(iat0 + 1, iat1 + 1)
                    if iline % 4 == 0:
                        f.write(line+"\n")
                    else:
                        f.write(line)
                    iline += 1
            if nbonds % 4 != 0:
                f.write("\n\n")
            else:
                f.write("\n")

            # NAngles section =====================================================
            iline = 1
            line = ""
            ntheta = 0
            for i in range(self._natoms):
                itheta = self._topology.find_all_paths_length(i, 2)
                for j in range(len(itheta)):
                    iat0 = itheta[j][0]
                    iat1 = itheta[j][1]
                    iat2 = itheta[j][2]
                    # Exchange extremes
                    if iat0 < iat2:
                        continue
                    ntheta += 1
                    line += str("{0:8d}{1:8d}{2:8d}").format(iat0 + 1, iat1 + 1, iat2 + 1)
                    if iline % 3 == 0:
                        line += "\n"
                    iline += 1

            line0 = str("{0:8d} {1:s}\n".format(ntheta, '!NTHETA: angles'))
            f.write(line0)
            f.write(line)
            if ntheta % 3 != 0:
                f.write("\n\n")
            else:
                f.write("\n")

            # NPhi section ========================================================
            nphi = 0
            iphilist = []
            iphich_tmp = []
            ich_prev = -1
            for i in range(self._natoms):
                iphi = self._topology.find_all_paths_length(i, 3)
                for j in range(len(iphi)):

                    iat0 = iphi[j][0]
                    iat1 = iphi[j][1]
                    iat2 = iphi[j][2]
                    iat3 = iphi[j][3]

                    # Remove repeat dihedrals ===================
                    if ich_prev == -1:
                        ich_prev = self._topology._iatch[iat0]
                        ich_curr = ich_prev
                    else:
                        ich_curr = self._topology._iatch[iat0]

                    if ich_curr != ich_prev:
                        iphich_tmp = []
                        ich_prev = ich_curr
                    iphich_tmp.append([iat0, iat1, iat2, iat3])
                    if iphich_tmp.count([iat3, iat2, iat1, iat0]):
                        continue
                    # Remove repeat dihedrals ====================
                    iphilist.append([iat0, iat1, iat2, iat3])
                    nphi += 1

            iline = 1
            line = ""
            nphi = 0
            for idx in iphilist:
                line += str("{0:8d}{1:8d}{2:8d}{3:8d}").format(idx[3] + 1, idx[2] + 1, idx[1] + 1, idx[0] + 1)
                nphi += 1
                if iline % 2 == 0:
                    line += "\n"
                iline += 1

            line0 = str("{0:8d} {1:s}\n".format(nphi, '!NPHI: dihedrals'))
            f.write(line0)
            f.write(line)
            if nphi % 2 != 0:
                f.write("\n\n")
            else:
                f.write("\n")

            # NImp section ========================================================
            if improper_angles is not None:
                nimpr = len(improper_angles)
                line0 = str("{0:8d} {1:s}\n".format(nimpr, '!NIMPHI: impropers'))
                iline = 1
                line = ""
                for iimpr in improper_angles:
                    line += str("{0:8d}{1:8d}{2:8d}{3:8d}").format(iimpr[0], iimpr[1], iimpr[2], iimpr[3])
                    if iline % 2 == 0:
                        line += "\n"
                    iline += 1
                f.write(line0)
                f.write(line)
                if nimpr % 2 != 0:
                    f.write("\n\n")
                else:
                    f.write("\n")

            else:
                line0 = str("{0:8d} {1:s}\n".format(0, '!NIMPHI: impropers'))
                f.write(line0)
                f.write("\n")

            line0 = str("{0:8d} {1:s}\n".format(0, '!NDON: donors'))
            f.write(line0)
            f.write("\n")
            line0 = str("{0:8d} {1:s}\n".format(0, '!NACC: acceptors'))
            f.write(line0)
            f.write("\n")
            line0 = str("{0:8d} {1:s}\n".format(0, '!NNB'))
            f.write(line0)
            iline = 1
            line1 = ""