import unittest
import utils
import datetime
from topology.readmol.readxsdformat import ReadXsdFormat


class TestReadXsd(unittest.TestCase):

    # #########################################################################
    @classmethod
    def setUpClass(cls):

        cls.filelog = "test_readxsd.log"
        cls.log = utils.init_logger("Output", fileoutput=cls.filelog, append=False, inscreen=False)
        m = "\n\t***************** START MStoLAGROAM TEST *****************"
        m += "\n\t\tTest length is about 97 seconds."
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
        filename = "../data/PET_3mon.xsd"
        xsd = ReadXsdFormat(filename)

        nmols = xsd.get_nmols()
        nres = xsd.get_nres()
        natoms = xsd.get_natoms()
        mol_reslist = xsd.get_residuelist()

        self.assertEqual(nmols, 2)
        self.assertEqual(natoms, 136)
        self.assertEqual(nres, 6)

        self.assertEqual(mol_reslist,
                         ['ethylene_terepthalate', 'ethylene_terepthalate', 'ethylene_terepthalate',
                          'ethylene_terepthalate', 'ethylene_terepthalate', 'ethylene_terepthalate'])
        xsd.write_xyz()
        xsd.write_gro()

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        self.log.info("\t\tFinishing: \t {}\n".format(now))
        m = "\t============== END   TEST_01 ================================"
        print(m) if self.log is None else self.log.info(m)

    #####################################################################################
    def test_02_checktopology_polymer_periodic(self):

        m = "\n\t============== START TEST_02 ================================\n"
        m += "\t   Reading a polymer from XSD file created with Materials Studio \n"
        m += "\t   Check the topology with bond orders \n"
        print(m) if self.log is None else self.log.info(m)

        # Open file
        filename = "../data/PET_3mon.xsd"
        xsd = ReadXsdFormat(filename, assign_bondorders=True)

        xsd._topology.draw_graph_forest_pygraphviz(title="graphs/test02_pet2ch3mon_")

        # Check natoms
        self.assertEqual(len(xsd._topology.get_vertices()), xsd.get_natoms())

        # Check bonds
        bonds = xsd._topology.get_allbonds()
        edges = xsd._topology.get_edges()
        for item in edges:
            i = sorted(list(item))
            self.assertIn(i, bonds)

        self.assertEqual(xsd.get_nbonds(), len(bonds))

        # Check cycles
        self.assertEqual(xsd._topology.iscyclic(), [True, True])

        # Check molecular charges
        self.assertEqual(xsd._topology._totalcharge_mol, [0, 0])
        self.assertEqual(xsd._topology._totalcharge, 0)

        # Check bond orders
        bo = {4: 2.0, 5: 1.0, 6: 1.0, 10: 1.5, 12: 1.5, 45: 1.0, 8: 1.5, 11: 1.0, 13: 1.5, 14: 1.0,
              7: 1.5, 9: 1.0, 15: 1.5,
              17: 1.0, 16: 1.0, 18: 2.0, 19: 1.0, 20: 1.0, 1: 1.0, 2: 1.0, 21: 1.0, 0: 1.0,
              3: 1.0, 22: 1.0, 23: 1.0,
              25: 1.0, 26: 1.0, 24: 1.0, 43: 1.0, 44: 1.0, 42: 1.0, 40: 1.0, 41: 2.0, 30: 1.5,
              38: 1.5, 31: 1.5, 32: 1.0,
              36: 1.5, 39: 1.0, 33: 1.5, 34: 1.0, 28: 1.0, 35: 1.5, 27: 2.0, 29: 1.0, 37: 1.0,
              69: 1.0, 46: 1.0, 48: 1.0,
              49: 1.0, 47: 1.0, 66: 1.0, 67: 1.0, 65: 1.0, 63: 1.0, 64: 2.0, 53: 1.5, 61: 1.5,
              54: 1.5, 55: 1.0, 59: 1.5,
              62: 1.0, 56: 1.5, 57: 1.0, 51: 1.0, 58: 1.5, 50: 2.0, 52: 1.0, 60: 1.0, 68: 1.0,
              74: 2.0, 75: 1.0, 76: 1.0,
              80: 1.5, 82: 1.5, 115: 1.0, 78: 1.5, 81: 1.0, 83: 1.5, 84: 1.0, 77: 1.5, 79: 1.0,
              85: 1.5, 87: 1.0, 86: 1.0,
              88: 2.0, 89: 1.0, 90: 1.0, 71: 1.0, 72: 1.0, 91: 1.0, 70: 1.0, 73: 1.0, 92: 1.0,
              93: 1.0, 95: 1.0, 96: 1.0,
              94: 1.0, 113: 1.0, 114: 1.0, 112: 1.0, 110: 1.0, 111: 2.0, 100: 1.5, 108: 1.5,
              101: 1.5, 102: 1.0, 106: 1.5,
              109: 1.0, 103: 1.5, 104: 1.0, 98: 1.0, 105: 1.5, 97: 2.0, 99: 1.0, 107: 1.0, 139: 1.0,
              116: 1.0, 118: 1.0,
              119: 1.0, 117: 1.0, 136: 1.0, 137: 1.0, 135: 1.0, 133: 1.0, 134: 2.0, 123: 1.5,
              131: 1.5, 124: 1.5, 125: 1.0,
              129: 1.5, 132: 1.0, 126: 1.5, 127: 1.0, 121: 1.0, 128: 1.5, 120: 2.0, 122: 1.0, 130: 1.0, 138: 1.0}

        self.assertEqual(xsd._topology._orderbonds, bo)

        # Backbone atoms
        bb = [0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1,
              0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1,
              1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1,
              0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1,
              0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1,
              0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1,
              0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1,
              1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0,
              1, 0, 1, 1, 0, 1, 1, 1]


        xsd.write_xyz(filename_xyz="PET_3mon.xyz")
        xsd.write_gro(filename_gro="PET_3mon.gro")
        xsd.write_pdb(filename_pdb="PET_3mon.pdb")

        self.assertEqual(bb, xsd._topology._isbackbone)

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        self.log.info("\t\tFinishing: \t {}\n".format(now))
        m = "\t============== END   TEST_02 ================================"
        print(m) if self.log is None else self.log.info(m)

    # # #######################################################################################
    def test_03_checktopology_nonpolymer_periodic(self):

        m = "\n\t============== START TEST_03 ================================\n"
        m += "\t   Reading a molecular system from XSD file created with Materials Studio \n"
        m += "\t   Check the topology with bond orders \n"
        print(m) if self.log is None else self.log.info(m)

        # Open file
        filename = "../data/bicarbonate_periodic.xsd"
        xsd = ReadXsdFormat(filename, assign_bondorders=True)

        xsd.write_xyz(filename_xyz="bicarbonate_periodic.xyz")
        xsd.write_gro(filename_gro="bicarbonate_periodic.gro")
        xsd.write_pdb(filename_pdb="bicarbonate_periodic.pdb")

        xsd._topology.draw_graph_forest_pygraphviz(title="graphs/test03_bicarbonate2mol_")

        # Check natoms
        self.assertEqual(len(xsd._topology.get_vertices()), xsd.get_natoms())

        # Check bonds
        bonds = xsd._topology.get_allbonds()
        edges = xsd._topology.get_edges()
        for item in edges:
            i = sorted(list(item))
            self.assertIn(i, bonds)

        # Check cycles
        self.assertEqual(xsd._topology.iscyclic(), [False, False])

        # Check molecular charges
        self.assertEqual(xsd._topology._totalcharge_mol, [-1, -1])
        self.assertEqual(xsd._topology._totalcharge, -2)

        # Check bond orders
        bo = {2: 1.5, 0: 1.5, 1: 1.0, 3: 1.0, 6: 1.5, 4: 1.5, 5: 1.0, 7: 1.}
        self.assertEqual(xsd._topology._orderbonds, bo)

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        self.log.info("\t\tFinishing: \t {}\n".format(now))
        m = "\t============== END   TEST_03 ================================"
        print(m) if self.log is None else self.log.info(m)

    # #######################################################################################
    def test_04_checktopology_nonpolymer_nonperiodic(self):

        m = "\n\t============== START TEST_04 ================================\n"
        m += "\t   Reading a molecular system from XSD file created with Materials Studio \n"
        m += "\t   Non-periodic and Check the topology with bond orders \n"
        print(m) if self.log is None else self.log.info(m)

        filename = "../data/methane_ethane_water.xsd"
        xsd = ReadXsdFormat(filename, assign_bondorders=True)

        xsd.write_xyz(filename_xyz="methane_ethane_water_nonp.xyz")
        xsd.write_gro(filename_gro="methane_ethane_water_nonp.gro")

        xsd._topology.draw_graph_forest_pygraphviz(title="graphs/test04_methaneethanewater_")

        # Check nmols
        nmols = xsd.get_nmols()
        self.assertEqual(nmols, 3)

        # Check nres
        nres = xsd.get_nres()
        self.assertEqual(nres, 3)

        # Check natoms
        self.assertEqual(len(xsd._topology.get_vertices()), xsd.get_natoms())

        # Check bonds
        bonds = xsd._topology.get_allbonds()
        edges = xsd._topology.get_edges()
        for item in edges:
            i = sorted(list(item))
            self.assertIn(i, bonds)

        # Check cycles
        self.assertEqual(xsd._topology.iscyclic(), [False, False, False])

        # Check molecular charges
        self.assertEqual(xsd._topology._totalcharge_mol, [0, 0, 0])
        self.assertEqual(xsd._topology._totalcharge, 0)

        # Check bond orders
        bo = {0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0, 5: 1.0,
              6: 1.0, 7: 1.0, 8: 1.0, 9: 1.0, 10: 1.0, 11: 1.0, 12: 1.0}
        self.assertEqual(xsd._topology._orderbonds, bo)

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        self.log.info("\t\tFinishing: \t {}\n".format(now))
        m = "\t============== END   TEST_04 ================================"
        print(m) if self.log is None else self.log.info(m)

    # #######################################################################################
    def test_05_checktopology_nonpolymer_periodic(self):

        m = "\n\t============== START TEST_05 ================================\n"
        m += "\t   Reading a molecular system from XSD file created with Materials Studio \n"
        m += "\t   Periodic and Check the topology with bond orders \n"
        print(m) if self.log is None else self.log.info(m)

        filename = "../data/met-eth-com_periodic.xsd"
        xsd = ReadXsdFormat(filename, assign_bondorders=True)

        xsd.write_xyz(filename_xyz="met-eth.com_per.xyz")
        xsd.write_gro(filename_gro="met-eth.com_per.gro")

        xsd._topology.draw_graph_forest_pygraphviz(title="graphs/test05_methaneethanecom_")

        # Check natoms
        self.assertEqual(len(xsd._topology.get_vertices()), xsd.get_natoms())

        # Check bonds
        bonds = xsd._topology.get_allbonds()
        edges = xsd._topology.get_edges()
        for item in edges:
            i = sorted(list(item))
            self.assertIn(i, bonds)

        # Check cycles
        self.assertEqual(xsd._topology.iscyclic(), [False, False, True])

        # Check molecular charges
        self.assertEqual(xsd._topology._totalcharge_mol, [0, 0, -1])
        self.assertEqual(xsd._topology._totalcharge, -1)

        # Check bond orders
        bo = {0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 2.0, 5: 1.0, 6: 1.0, 7: 1.0, 8: 1.0,
              13: 1.5, 14: 1.5, 34: 1.0, 9: 1.5, 15: 1.0, 12: 1.5, 18: 1.0, 32: 1.5, 33: 1.5,
              10: 1.5, 16: 1.0, 11: 1.5, 24: 1.0, 17: 1.0, 19: 1.0, 22: 1.0, 23: 1.0, 20: 1.0,
              21: 1.0, 31: 1.0, 25: 1.0, 29: 1.0, 30: 1.0, 26: 1.0, 27: 1.0, 28: 1.0}
        self.assertEqual(xsd._topology._orderbonds, bo)

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        self.log.info("\t\tFinishing: \t {}\n".format(now))
        m = "\t============== END   TEST_05 ================================"
        print(m) if self.log is None else self.log.info(m)
    #
    # #######################################################################################
    def test_06_checktopology_nonpolymer_periodic(self):

        m = "\n\t============== START TEST_06 ================================\n"
        m += "\t   Reading a molecular system from XSD file created with Materials Studio \n"
        m += "\t   Mol and polymer. Periodic and Check the topology with bond orders \n"
        print(m) if self.log is None else self.log.info(m)

        filename = "../data/mol-pol_periodic.xsd"
        xsd = ReadXsdFormat(filename, assign_bondorders=True)

        xsd.write_xyz(filename_xyz="mol-pol_periodic.xyz")
        xsd.write_gro(filename_gro="mol-pol_periodic.gro")

        xsd._topology.draw_graph_forest_pygraphviz(title="graphs/test06_mol-pol_per_")

        # Check natoms
        self.assertEqual(len(xsd._topology.get_vertices()), xsd.get_natoms())

        # Check bonds
        bonds = xsd._topology.get_allbonds()
        edges = xsd._topology.get_edges()
        for item in edges:
            i = sorted(list(item))
            self.assertIn(i, bonds)

        # Check cycles
        self.assertEqual(xsd._topology.iscyclic(), [True, False])

        # Check molecular charges
        self.assertEqual(xsd._topology._totalcharge_mol, [-1, 0])
        self.assertEqual(xsd._topology._totalcharge, -1)

        # Check bond orders
        bo = {0: 1.5, 1: 1.5, 2: 1.5, 3: 1.5, 4: 1.5, 5: 1.5, 6: 1.0, 7: 1.0, 8: 1.0, 9: 1.0,
              10: 1.0, 11: 1.0, 12: 1.0, 13: 1.0, 14: 1.0, 15: 1.0, 16: 1.0, 17: 1.0, 18: 1.0,
              19: 1.0, 20: 1.0, 21: 1.0, 22: 1.0, 23: 1.5, 24: 1.5, 25: 1.0, 26: 1.0, 27: 1.0,
              28: 1.0, 29: 1.0, 30: 1.0, 31: 1.0, 32: 1.0, 33: 1.0, 34: 1.0, 35: 1.0, 36: 1.0,
              37: 1.0, 38: 1.0, 39: 1.0, 40: 1.0, 41: 1.0, 42: 1.0, 43: 1.0, 44: 1.0, 45: 1.0, 46: 1.0,
              47: 1.0, 48: 1.0, 49: 1.0, 50: 1.0, 51: 1.0, 52: 1.0, 53: 1.0, 54: 1.0,  55: 1.0, 56: 1.0,
              57: 1.0, 58: 1.0, 59: 1.0, 60: 1.0, 61: 1.0, 62: 1.0, 63: 1.0, 64: 1.0,  65: 1.0, 66: 1.0}
        self.assertEqual(xsd._topology._orderbonds, bo)

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        self.log.info("\t\tFinishing: \t {}\n".format(now))
        m = "\t============== END   TEST_06 ================================"
        print(m) if self.log is None else self.log.info(m)

    ######################################################################################
    def test_07_checktopology_polymer_periodic_real(self):

        m = "\n\t============== START TEST_07 ================================\n"
        m += "\t   Reading a molecular system from XSD file created with Materials Studio \n"
        m += "\t   5 PE chains and 150 monomers (C300). Periodic and Check the topology with bond orders \n"
        print(m) if self.log is None else self.log.info(m)

        # Reading XSD format
        s = datetime.datetime.now()
        m = "\tReading XSD format (4500 atoms) ..."
        print(m) if self.log is None else self.log.info(m)
        filename = "../data/03-NPT_Frame_001.xsd"
        xsd = ReadXsdFormat(filename, assign_bondorders=True)
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        # Drawing topology
        s = datetime.datetime.now()
        m = "\tDrawing forest for the XSD graph..."
        print(m) if self.log is None else self.log.info(m)
        xsd._topology.draw_graph_forest_pygraphviz(title="graphs/test07_PE05C300_per_")
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        # Writing topology
        s = datetime.datetime.now()
        m = "\tWriting XYZ and GRO formats..."
        print(m) if self.log is None else self.log.info(m)
        xsd.write_xyz(filename_xyz="mol-pol_periodic.xyz")
        xsd.write_gro(filename_gro="mol-pol_periodic.gro")
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        # Check natoms & bonds
        s = datetime.datetime.now()
        m = "\tChecking number of atoms and bonds..."
        print(m) if self.log is None else self.log.info(m)
        self.assertEqual(len(xsd._topology.get_vertices()), xsd.get_natoms())
        # Check bonds
        bonds = xsd._topology.get_allbonds()
        edges = xsd._topology.get_edges()
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
        self.assertEqual(xsd._topology.iscyclic(), [False, False, False, False, False])
        # Check molecular charges
        self.assertEqual(xsd._topology._totalcharge_mol, [0, 0, 0, 0, 0])
        self.assertEqual(xsd._topology._totalcharge, 0)
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        # Check bond orders
        s = datetime.datetime.now()
        m = "\tChecking bondorders..."
        print(m) if self.log is None else self.log.info(m)
        for k, v in xsd._topology._orderbonds.items():
            self.assertEqual(v, 1.0)
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        self.log.info("\t\tFinishing: \t {}\n".format(now))
        m = "\t============== END   TEST_07 ================================"
        print(m) if self.log is None else self.log.info(m)

    #######################################################################################
    def test_08_separatechains_polymer_periodic_real(self):

        m = "\n\t============== START TEST_08 ================================\n"
        m += "\t   Reading a molecular system from XSD file created with Materials Studio \n"
        m += "\t   Separate each chain in a unique pdb file. This can be used to make easier the \n"
        m += "\t   typing of the molecule. \n"
        m += "\t   10 chains PE-SCB C300 4 butyl branches, all atomsitic model. \n"
        print(m) if self.log is None else self.log.info(m)

        # Reading XSD format
        s = datetime.datetime.now()
        m = "\tReading XSD format (9500 atoms) ..."
        print(m) if self.log is None else self.log.info(m)
        filename = "../data/01-Min_Frame_001.xsd"
        xsd = ReadXsdFormat(filename, assign_bondorders=True)
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        # Writing topology
        s = datetime.datetime.now()
        m = "\tWriting PDB and separate chains..."
        print(m) if self.log is None else self.log.info(m)
        xsd.write_pdb(filename_pdb="01-Min_Frame_001.pdb", separate_chains=True)
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        self.log.info("\t\tFinishing: \t {}\n".format(now))
        m = "\t============== END   TEST_08 ================================"
        print(m) if self.log is None else self.log.info(m)

    # #######################################################################################
    def test_09_polymerlinearchain_nonperiodic(self):

        m = "\n\t============== START TEST_09 ================================\n"
        m += "\t   Reading a molecular system from XSD file created with Materials Studio \n"
        m += "\t   Polymer chains with repeat units and not periodic \n"
        m += "\t   1 chains P4HB 100 monomers, all atomsitic model. \n"
        print(m) if self.log is None else self.log.info(m)

        # Reading XSD format
        s = datetime.datetime.now()
        m = "\tReading XSD format (1205 atoms) ..."
        print(m) if self.log is None else self.log.info(m)
        filename = "../data/Poly4HB_100monomers.xsd"
        xsd = ReadXsdFormat(filename, assign_bondorders=False)
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        # Writing topology
        s = datetime.datetime.now()
        m = "\tWriting PDB and separate chains..."
        print(m) if self.log is None else self.log.info(m)
        xsd.write_pdb(filename_pdb="Poly4HB_100monomers.pdb", separate_chains=True)
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        self.log.info("\t\tFinishing: \t {}\n".format(now))
        m = "\t============== END   TEST_08 ================================"
        print(m) if self.log is None else self.log.info(m)

    # #######################################################################################
    def test_10_polymercrystal_supercell(self):

        m = "\n\t============== START TEST_10 ================================\n"
        m += "\t   Reading a molecular system from XSD file created with Materials Studio \n"
        m += "\t   8 Polymer chains without repeat units and periodic \n"
        m += "\t   This system has been created from a cell of unit of the P4HB crystal.\n"
        print(m) if self.log is None else self.log.info(m)

        # Reading XSD format
        s = datetime.datetime.now()
        m = "\tReading XSD format (400 atoms) ..."
        print(m) if self.log is None else self.log.info(m)
        filename = "../data/P4HB_crystal_2x2x2.xsd"
        xsd = ReadXsdFormat(filename, assign_bondorders=True)
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        # Drawing topology
        s = datetime.datetime.now()
        m = "\tDrawing forest for the XSD graph..."
        print(m) if self.log is None else self.log.info(m)
        xsd._topology.draw_graph_forest_pygraphviz(title="graphs/test10_P4HB_per_")
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        # Writing topology
        s = datetime.datetime.now()
        m = "\tWriting XYZ and GRO formats..."
        print(m) if self.log is None else self.log.info(m)
        xsd.write_xyz(filename_xyz="P4HB_crystal_2x2x2.xyz")
        xsd.write_gro(filename_gro="P4HB_crystal_2x2x2.gro")
        xsd.write_pdb(filename_pdb="P4HB_crystal_2x2x2.pdb")
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        # Check natoms & bonds
        s = datetime.datetime.now()
        m = "\tChecking number of atoms and bonds..."
        print(m) if self.log is None else self.log.info(m)
        self.assertEqual(len(xsd._topology.get_vertices()), xsd.get_natoms())
        # Check bonds
        bonds = xsd._topology.get_allbonds()
        edges = xsd._topology.get_edges()
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
        self.assertEqual(xsd._topology.iscyclic(), 8*[False])
        # Check molecular charges
        self.assertEqual(xsd._topology._totalcharge_mol, 8*[0])
        self.assertEqual(xsd._topology._totalcharge, 0)
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

        # Check bond orders
        s = datetime.datetime.now()
        m = "\tChecking bond orders..."
        print(m) if self.log is None else self.log.info(m)
        db = [2, 15, 194, 206, 25, 39, 216, 229, 50, 63, 73, 87, 240, 252,
              98, 111, 121, 135, 146, 262, 275, 159, 169, 183, 286, 298, 308,
              321, 332, 344, 354, 367]
        for k, v in xsd._topology._orderbonds.items():
            if k in db:
                self.assertEqual(v, 2.0)
            else:
                self.assertEqual(v, 1.0)
        e = datetime.datetime.now()
        duration = e - s
        m = "\t\t({0:.3f} seconds)".format(duration.total_seconds())
        print(m) if self.log is None else self.log.info(m)

    # # #######################################################################################
    def test_11_polymercrystal_supercell(self):
        #
        #     m = "\n\t============== START TEST_11 ================================\n"
        #     m += "\t   Reading a molecular system from XSD file created with Materials Studio \n"
        #     m += "\t   8 Polymer chains without repeat units and periodic \n"
        #     m += "\t   This system has been created from a cell of unit of the P4HB crystal.\n"
        #     print(m) if self.log is None else self.log.info(m)

        # Open file
        filename = "../data/P3HB_3mon_01c.xsd"
        xsd = ReadXsdFormat(filename)

        nmols = xsd.get_nmols()
        nres = xsd.get_nres()
        natoms = xsd.get_natoms()
        mol_reslist = xsd.get_residuelist()

        self.assertEqual(nmols, 1)
        self.assertEqual(natoms, 39)
        self.assertEqual(nres, 1)

        self.assertEqual(mol_reslist, ['UNK'])
        xsd.write_xyz(filename_xyz="P3HB_3mon_01c.xyz")
        xsd.write_gro(filename_gro="P3HB_3mon_01c.gro")
        xsd.write_pdb(filename_pdb="P3HB_3mon_01c.pdb")

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        self.log.info("\t\tFinishing: \t {}\n".format(now))
        m = "\t============== END   TEST_01 ================================"
        print(m) if self.log is None else self.log.info(m)

    # ##################################################################################################################
    @classmethod
    def tearDownClass(cls):

        m = "\n\t***************** END MStoLAGROAM TEST *****************"
        print(m) if cls.log is None else cls.log.info(m)
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        cls.log.info("\t\tFinishing: \t {}\n".format(now))


if __name__ == '__main__':

    unittest.main()
