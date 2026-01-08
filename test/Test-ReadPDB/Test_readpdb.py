import unittest
import utils
import datetime
from copy import copy
from topology.readmol.readpdbformat import ReadPdbFormat
from topology.readmol.readxsdformat import ReadXsdFormat


class TestReadPdb(unittest.TestCase):

    # #########################################################################
    @classmethod
    def setUpClass(cls):

        cls.filelog = "test_readpdb.log"
        cls.log = utils.init_logger("Output", fileoutput=cls.filelog, append=False, inscreen=False)
        m = "\n\t***************** START READPDB TEST *****************"
        m += "\n\t\tTest length is about XX seconds."
        print(m) if cls.log is None else cls.log.info(m)
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        cls.log.info("\t\tStarting: \t {}\n".format(now))

    # #######################################################################################
    def test_01_readfile_polymer_periodic(self):

        m = "\n\t============== START TEST_01 ================================\n"
        m += "\t   Reading a polymer from XSD file created with Materials Studio\n"
        m += "\t   Check basic topology without bond orders \n"
        print(m) if self.log is None else self.log.info(m)

        # Open file
        filename = "../data/PET_3mon.pdb"
        pdb = ReadPdbFormat(filename)

        nmols = pdb.get_nmols()
        nres = pdb.get_nres()
        natoms = pdb.get_natoms()
        mol_reslist = pdb.get_residuelist()
        # #
        self.assertEqual(nmols, 2)
        self.assertEqual(natoms, 136)
        self.assertEqual(nres, 6)
        self.assertEqual(mol_reslist,  ['001', '002', '003',
                                        '004', '005', '006'])

        pdb.write_xyz(filename_xyz="PET_3mon.xyz")
        pdb.write_gro(filename_gro="PET_3mon.gro")
        pdb.write_pdb(filename_pdb="PET_3mon.pdb")

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        self.log.info("\t\tFinishing: \t {}\n".format(now))
        m = "\t============== END   TEST_01 ================================"
        print(m) if self.log is None else self.log.info(m)

    # # #####################################################################################
    def test_02_checktopology_polymer_periodic(self):

        m = "\n\t============== START TEST_02 ================================\n"
        m += "\t   Reading a polymer from PDB file created with Materials Studio \n"
        m += "\t   Check the topology with bond orders \n"
        print(m) if self.log is None else self.log.info(m)

        # Open file
        filename = "../data/PET_3mon.pdb"
        filename2 = "../data/PET_3mon.xsd"
        pdb = ReadPdbFormat(filename, assign_bondorders=True)
        xsd = ReadXsdFormat(filename2, assign_bondorders=True)

        pdb._topology.draw_graph_forest_pygraphviz(title="graphs/test02_pet2ch3mon_")

        # Check natoms
        self.assertEqual(len(pdb._topology.get_vertices()), pdb.get_natoms())

        # Check bonds
        bonds = pdb._topology.get_allbonds()
        edges = pdb._topology.get_edges()
        for item in edges:
            i = sorted(list(item))
            self.assertIn(i, bonds)

        self.assertEqual(pdb.get_nbonds(), len(bonds))

        # Check cycles
        self.assertEqual(pdb._topology.iscyclic(), [True, True])

        # Check molecular charges
        self.assertEqual(pdb._topology._totalcharge_mol, [0, 0])
        self.assertEqual(pdb._topology._totalcharge, 0)

        print(pdb._topology.get_allbonds())
        print(xsd._topology.get_allbonds())
        a = pdb._topology._orderbonds

        # Check bond orders
        bo = {0: 1.0, 1: 1.0,2: 2.0, 7: 1.5, 8: 1.5, 17: 1.0, 4: 1.5, 6: 1.0, 10: 1.5, 9: 1.0,
              3: 1.0, 5: 1.5, 12: 1.5, 13: 1.0, 11: 1.0, 14: 2.0, 15: 1.0, 16: 1.0, 18: 1.0, 19: 1.0,
              20: 1.0, 21: 1.0, 22: 1.0, 23: 1.0, 44: 1.0, 45: 1.0, 46: 1.0, 40: 1.0, 42: 1.0,
              43: 1.0, 39: 1.0, 37: 1.0, 38: 2.0, 29: 1.5, 36: 1.5, 27: 1.0, 28: 1.5, 34: 1.5, 35: 1.0,
              30: 1.0, 31: 1.5, 24: 1.0, 32: 1.5, 25: 1.0, 26: 2.0, 33: 1.0, 41: 1.0, 67: 1.0, 68: 1.0,
              69: 1.0, 63: 1.0, 65: 1.0, 66: 1.0, 62: 1.0, 60: 1.0, 61: 2.0, 52: 1.5, 59: 1.5, 50: 1.0,
              51: 1.5, 57: 1.5, 58: 1.0, 53: 1.0, 54: 1.5, 47: 1.0, 55: 1.5, 48: 1.0, 49: 2.0, 56: 1.0,
              64: 1.0, 70: 1.0, 71: 1.0, 72: 2.0, 77: 1.5, 78: 1.5, 87: 1.0, 74: 1.5, 76: 1.0, 80: 1.5,
              79: 1.0, 73: 1.0, 75: 1.5, 82: 1.5, 83: 1.0, 81: 1.0, 84: 2.0, 85: 1.0, 86: 1.0, 88: 1.0,
              89: 1.0, 90: 1.0, 91: 1.0, 92: 1.0, 93: 1.0, 114: 1.0, 115: 1.0, 116: 1.0, 110: 1.0,
              112: 1.0, 113: 1.0, 109: 1.0, 107: 1.0, 108: 2.0, 99: 1.5, 106: 1.5, 97: 1.0, 98: 1.5,
              104: 1.5, 105: 1.0, 100: 1.0, 101: 1.5, 94: 1.0, 102: 1.5, 95: 1.0, 96: 2.0, 103: 1.0,
              111: 1.0, 137: 1.0, 138: 1.0, 139: 1.0, 133: 1.0, 135: 1.0, 136: 1.0, 132: 1.0, 130: 1.0,
              131: 2.0, 122: 1.5, 129: 1.5, 120: 1.0, 121: 1.5, 127: 1.5, 128: 1.0, 123: 1.0, 124: 1.5,
              117: 1.0, 125: 1.5, 118: 1.0, 119: 2.0, 126: 1.0, 134: 1.0}

        self.assertEqual(pdb._topology._orderbonds, bo)
        #
        # Backbone atoms
        bb = 136*[0]
        self.assertEqual(bb, pdb._topology._isbackbone)

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        self.log.info("\t\tFinishing: \t {}\n".format(now))
        m = "\t============== END   TEST_02 ================================"
        print(m) if self.log is None else self.log.info(m)

    # # #######################################################################################
    def test_03_checktopology_nonpolymer_periodic(self):

        m = "\n\t============== START TEST_03 ================================\n"
        m += "\t   Reading a molecular system from PDB file created with Materials Studio \n"
        m += "\t   Check the topology with bond orders \n"
        print(m) if self.log is None else self.log.info(m)

        # Open file
        filename = "../data/bicarbonate_periodic.pdb"
        pdb = ReadPdbFormat(filename, assign_bondorders=True)

        pdb._topology.draw_graph_forest_pygraphviz(title="graphs/test03_bicarbonate2mol_")

        # Check natoms
        self.assertEqual(len(pdb._topology.get_vertices()), pdb.get_natoms())

        # Check bonds
        bonds = pdb._topology.get_allbonds()
        edges = pdb._topology.get_edges()
        for item in edges:
            i = sorted(list(item))
            self.assertIn(i, bonds)

        # Check cycles
        self.assertEqual(pdb._topology.iscyclic(), [False, False])

        # Check molecular charges
        self.assertNotEqual(pdb._topology._totalcharge_mol, [-1, -1])
        self.assertNotEqual(pdb._topology._totalcharge, -2)

        m = "\t\tCharges are not calculated from PDB files."
        print(m) if self.log is None else self.log.info(m)

        # Check bond orders
        print(pdb._topology.get_allbonds())

        bo = {0: 1.5, 1: 1.5, 2: 1.0, 3: 1.0, 4: 1.5, 5: 1.5, 6: 1.0, 7: 1.0}
        self.assertEqual(pdb._topology._orderbonds, bo)

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        self.log.info("\t\tFinishing: \t {}\n".format(now))
        m = "\t============== END   TEST_03 ================================"
        print(m) if self.log is None else self.log.info(m)

    # #######################################################################################
    def test_04_separatechains_polymer_periodic_real(self):

        m = "\n\t============== START TEST_04 ================================\n"
        m += "\t   Reading a molecular system from PDB file created with Materials Studio \n"
        m += "\t   Separate each chain in a unique pdb file. This can be used to make easier the \n"
        m += "\t   typing of the molecule. \n"
        m += "\t   10 chains PE-SCB C300 4 butyl branches, all atomsitic model. \n"
        print(m) if self.log is None else self.log.info(m)

        # Reading XSD format
        s = datetime.datetime.now()
        m = "\tReading PDB format (9500 atoms) ..."
        print(m) if self.log is None else self.log.info(m)
        filename = "../data/01-Min_Frame_001.pdb"
        pdb = ReadPdbFormat(filename, assign_bondorders=True)
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        # Writing topology
        s = datetime.datetime.now()
        m = "\tWriting PDB and separate chains..."
        print(m) if self.log is None else self.log.info(m)
        pdb.write_pdb(filename_pdb="01-Min_Frame_001.pdb", separate_chains=True)
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        self.log.info("\t\tFinishing: \t {}\n".format(now))
        m = "\t============== END   TEST_08 ================================"
        print(m) if self.log is None else self.log.info(m)

    # #######################################################################################
    def test_05_polymercrystal_supercell(self):

        m = "\n\t============== START TEST_05 ================================\n"
        m += "\t   Reading a molecular system from PDB file created with Materials Studio \n"
        m += "\t   8 Polymer chains without repeat units and periodic \n"
        m += "\t   This system has been created from a cell of unit of the P4HB crystal.\n"
        print(m) if self.log is None else self.log.info(m)

        # Reading XSD format
        s = datetime.datetime.now()
        m = "\tReading PDB format (400 atoms) ..."
        print(m) if self.log is None else self.log.info(m)
        filename = "../data/P4HB_crystal_2x2x2.pdb"
        pdb = ReadPdbFormat(filename, assign_bondorders=True)

        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        # Drawing topology
        s = datetime.datetime.now()
        m = "\tDrawing forest for the XSD graph..."
        print(m) if self.log is None else self.log.info(m)
        pdb._topology.draw_graph_forest_pygraphviz(title="graphs/test10_P4HB_per_")
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        # Writing topology
        s = datetime.datetime.now()
        m = "\tWriting XYZ and GRO formats..."
        print(m) if self.log is None else self.log.info(m)
        pdb.write_xyz(filename_xyz="P4HB_crystal_2x2x2.xyz")
        pdb.write_gro(filename_gro="P4HB_crystal_2x2x2.gro")
        pdb.write_pdb(filename_pdb="P4HB_crystal_2x2x2.pdb")
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        # Check natoms & bonds
        s = datetime.datetime.now()
        m = "\tChecking number of atoms and bonds..."
        print(m) if self.log is None else self.log.info(m)
        self.assertEqual(len(pdb._topology.get_vertices()), pdb.get_natoms())
        # Check bonds
        bonds = pdb._topology.get_allbonds()
        edges = pdb._topology.get_edges()
        for item in edges:
            i = sorted(list(item))
            self.assertIn(i, bonds)
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        # Check cycles and molecular charges
        s = datetime.datetime.now()
        m = "\tChecking cycles and charges..."
        print(m) if self.log is None else self.log.info(m)
        self.assertEqual(pdb._topology.iscyclic(), 8*[False])
        # Check molecular charges
        self.assertEqual(pdb._topology._totalcharge_mol, 8*[0])
        self.assertEqual(pdb._topology._totalcharge, 0)
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        # Check bond orders
        s = datetime.datetime.now()
        m = "\tChecking bond orders..."
        print(m) if self.log is None else self.log.info(m)
        db = [2, 14, 202, 214, 26, 41, 225, 239, 52, 64, 76, 91, 250, 262,
              102, 114, 126, 141, 273, 287, 152, 164, 176, 191, 298, 310, 321,
              335, 346, 358, 369, 383]
        for k, v in pdb._topology._orderbonds.items():
            if k in db:
                self.assertEqual(v, 2.0)
            else:
                self.assertEqual(v, 1.0)
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

    # #######################################################################################
    def test_06_polymercrystal_supercell_noconect(self):

        m = "\n\t============== START TEST_06 ================================\n"
        m += "\t   Reading a molecular system from PDB file created with Materials Studio \n"
        m += "\t   8 Polymer chains without repeat units and periodic \n"
        m += "\t   CONECT is not present in the pdb \n"
        m += "\t   This system has been created from a cell of unit of the P4HB crystal.\n"
        print(m) if self.log is None else self.log.info(m)

        # Reading XSD format
        s = datetime.datetime.now()
        m = "\tReading PDB format (400 atoms) ..."
        print(m) if self.log is None else self.log.info(m)
        filename = "../data/P4HB_crystal_2x2x2_noconect.pdb"
        pdb = ReadPdbFormat(filename, assign_bondorders=True)

        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        # Drawing topology
        s = datetime.datetime.now()
        m = "\tDrawing forest for the XSD graph..."
        print(m) if self.log is None else self.log.info(m)
        pdb._topology.draw_graph_forest_pygraphviz(title="graphs/test10_P4HB_per_")
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        # Writing topology
        s = datetime.datetime.now()
        m = "\tWriting XYZ, PDB and GRO formats..."
        print(m) if self.log is None else self.log.info(m)
        pdb.write_xyz(filename_xyz="P4HB_crystal_2x2x2_noconect.xyz")
        pdb.write_gro(filename_gro="P4HB_crystal_2x2x2_noconect.gro")
        pdb.write_pdb(filename_pdb="P4HB_crystal_2x2x2_noconect.pdb")
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        # Check natoms & bonds
        s = datetime.datetime.now()
        m = "\tChecking number of atoms and bonds..."
        print(m) if self.log is None else self.log.info(m)
        self.assertEqual(len(pdb._topology.get_vertices()), pdb.get_natoms())
        # Check bonds
        bonds = pdb._topology.get_allbonds()
        edges = pdb._topology.get_edges()
        for item in edges:
            i = sorted(list(item))
            self.assertIn(i, bonds)
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        # Check cycles and molecular charges
        s = datetime.datetime.now()
        m = "\tChecking cycles and charges..."
        print(m) if self.log is None else self.log.info(m)
        self.assertEqual(pdb._topology.iscyclic(), 8*[False])
        # Check molecular charges
        self.assertEqual(pdb._topology._totalcharge_mol, 8*[0])
        self.assertEqual(pdb._topology._totalcharge, 0)
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        # Check bond orders
        s = datetime.datetime.now()
        m = "\tChecking bond orders..."
        print(m) if self.log is None else self.log.info(m)
        db = [2, 14, 202, 214, 26, 41, 225, 239, 52, 64, 76, 91, 250, 262,
              102, 114, 126, 141, 273, 287, 152, 164, 176, 191, 298, 310, 321,
              335, 346, 358, 369, 383]
        for k, v in pdb._topology._orderbonds.items():
            if k in db:
                self.assertEqual(v, 2.0)
            else:
                self.assertEqual(v, 1.0)
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

    # #######################################################################################
    def test_07_deepcopy_pdbfile(self):

        from MDAnalysis import Universe, exceptions

        m = "\n\t============== START TEST_07 ================================\n"
        m += "\t   Reading a polymer from PDB file created with Materials Studio\n"
        m += "\t   Check deepcopy \n"
        print(m) if self.log is None else self.log.info(m)

        # Open file
        filename = "../data/PET_3mon.pdb"
        pdb = ReadPdbFormat(filename)
        pdbcopy = copy(pdb)
        self.assertTrue(pdb == pdbcopy)

    # ##################################################################################################################
    @classmethod
    def tearDownClass(cls):

        m = "\n\t***************** END MStoLAGROAM TEST *****************"
        print(m) if cls.log is None else cls.log.info(m)
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        cls.log.info("\t\tFinishing: \t {}\n".format(now))


if __name__ == '__main__':

    unittest.main()
