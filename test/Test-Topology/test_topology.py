import unittest
import numpy as np
from copy import copy
import topology as top
import utils
import datetime
from topology.readmol.readpdbformat import ReadPdbFormat
from topology.readmol.readxsdformat import ReadXsdFormat


# noinspection PyUnresolvedReferences
class TopologyTests (unittest.TestCase):

    # #########################################################################
    @classmethod
    def setUpClass(cls):

        cls.filelog = "test_topology.log"
        cls.log = utils.init_logger("Output", fileoutput=cls.filelog, append=False, inscreen=False)
        m = "\n\t***************** START Gecos rdkit TEST *****************"
        print(m) if cls.log is None else cls.log.info(m)
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        cls.log.info("\t\tStarting: \t {}\n".format(now))

    # #########################################################################
    def test_01_detect_connectivity_01(self):

        # Allylbenzene ==================
        m = "\t============== START TEST_01 ================================\n"
        m += "\t   Testing connectivity in allylbenzene molecule (CID:9309)  \n"
        m += "\t   without hydrogens  \n"
        m += "\t============== END   TEST_01 ================================"
        print(m) if self.log is None else self.log.info(m)

        filenamepdb = "../data/bond_orders/Conformer3D_CID_9309_withoutH.pdb"
        s1 = top.Segment(filecoord=filenamepdb)

        self.assertEqual(s1._topology._bonds,
                         [{0, 1}, {0, 2},  {0, 3}, {1, 7}, {2, 4},
                          {3, 5}, {4, 6}, {5, 6}, {7, 8}])

        s2 = top.Segment(filecoord=filenamepdb, filetop=filenamepdb)
        self.assertEqual(s2._topology._bonds,
                         [{0, 1}, {0, 2},  {0, 3}, {1, 7}, {2, 4},
                          {3, 5}, {4, 6}, {5, 6}, {7, 8}])

    # #########################################################################
    def test_02_detect_connectivity_assignbo_noH(self):

        # Allylbenzene ==================
        m = "\n\t============== START TEST_02 ================================\n"
        m += "\t   Testing connectivity in allylbenzene molecule (CID:9309)  \n"
        m += "\t   Assign bond orders without Hs"
        print(m) if self.log is None else self.log.info(m)

        filenamepdb = "../data/bond_orders/Conformer3D_CID_9309_withoutH.pdb"
        s1 = top.Segment(filecoord=filenamepdb)

        self.assertEqual(s1._topology._bonds,
                         [{0, 1}, {0, 2},  {0, 3}, {1, 7}, {2, 4},
                          {3, 5}, {4, 6}, {5, 6}, {7, 8}],)

        s1._topology.assign_bond_orders()
        m = "\tBond orders of Allylbenzene without H (CID: 9309)" \
            " ARE NOT correctly predicted!!. Test02"
        print(m) if self.log is None else self.log.info(m)

        bo = {0: 1.0, 1: 1.5, 2: 1.5, 3: 1.0, 4: 1.5, 5: 1.5, 6: 1.5, 7: 1.5, 8: 2.0}
        self.assertNotEqual(s1._topology._orderbonds, bo)
        m = "\t============== END   TEST_02 ================================"
        print(m) if self.log is None else self.log.info(m)

    # #########################################################################
    def test_03_detect_connectivity_assignbo_noH(self):

        # 1,3-Benzithiazole ==================
        m = "\n\t============== START TEST_03 ================================\n"
        m += "\t   Testing connectivity in 1,3-benzothiazole molecule (CID:7222)  \n"
        m += "\t   Assign bond orders with Hs"
        print(m) if self.log is None else self.log.info(m)

        filenamepdb = "../data/bond_orders/Conformer3D_CID_7222.pdb"
        s1 = top.Segment(filecoord=filenamepdb)

        self.assertEqual(s1._topology._bonds,
                         [{0, 2}, {0, 8},  {1, 3}, {8, 1}, {2, 3},
                          {2, 4}, {3, 5}, {4, 6}, {9, 4}, {5, 7},
                          {10, 5}, {6, 7}, {11, 6}, {12, 7}, {8, 13}])

        s1._topology.assign_bond_orders()
        bo = {0: 1.0, 1: 1.0, 4: 1.5, 5: 1.5, 3: 2.0,
              14: 1.0, 2: 1.0, 6: 1.5, 7: 1.5, 8: 1.0, 9: 1.5,
              10: 1.0, 11: 1.5, 13: 1.0, 12: 1.0}
        self.assertEqual(s1._topology._orderbonds, bo)
        m = "\t============== END   TEST_03 ================================"
        print(m) if self.log is None else self.log.info(m)

    #########################################################################
    def test_04_detect_maximalvalence_02(self):

        # Allylbenzene distorted with more atoms==================
        m = "\n\t============== START TEST_04 ================================\n"
        m += "\t         Testing maximum valence molecule (CID:9309)\n"
        print(m) if self.log is None else self.log.info(m)

        filenamepdb = "../data/bond_orders/Conformer3D_CID_9309_withoutH_badcoordination.pdb"
        s1 = top.Segment(filecoord=filenamepdb)
        m = "\t  An atom in the system exceeds the maximum valence (CID: 9309)" \
            "Check the topology.detect_connectivity method.!!. Test04"
        self.assertEqual(s1._topology._bonds,
                         [{0, 1}, {0, 2}, {0, 3}, {1, 7}, {2, 4},
                          {3, 5}, {4, 6}, {5, 6}, {7, 8}, {7, 10},
                          {7, 11}, {8, 12}, {10, 11}],)
        print(m) if self.log is None else self.log.info(m)

        s2 = top.Segment(filecoord=filenamepdb, filetop=filenamepdb)
        m = "\t  An atom in the system exceeds the maximum valence (CID: 9309)" \
            "Check the topology.detect_connectivity method.!!. with CONNECT Test04"
        self.assertEqual(s2._topology._bonds,
                         [{0, 1}, {0, 2}, {0, 3}, {1, 7}, {2, 4},
                          {3, 5}, {4, 6}, {5, 6}, {7, 8}, {7, 9}, {7, 10}, {7, 11}, {7, 12},
                          {9, 11}, {10, 11}],
                         msg=m)
        print(m) if self.log is None else self.log.info(m)

        m = "\t============== END   TEST_04 ================================"
        print(m) if self.log is None else self.log.info(m)

    ########################################################################
    def test_05_detect_nrings_01(self):

        m = "\n\t============== START TEST_05 ================================\n"
        m += "\t       Testing rings in cyclic molecules (CID:9309)     \n"
        m += "\t       All molecules without rings     "
        print(m) if self.log is None else self.log.info(m)

        # One aromatic ring
        filenamepdb = "../data/bond_orders/Conformer3D_CID_9309_withoutH.pdb"
        s1 = top.Segment(filecoord=filenamepdb, filetop=filenamepdb)
        s1._topology.perception_rings()
        self.assertTrue(s1._topology._cycles == [[0, 2, 4, 6, 5, 3, 0]])

        # No rings
        filenamepdb = "../data/bond_orders/Conformer3D_CID_12204213_withoutH.pdb"
        s2 = top.Segment(filecoord=filenamepdb, filetop=filenamepdb)

        s2._topology.perception_rings()
        self.assertTrue(s2._topology._cycles == [])

        # Two aromatic rings
        filenamepdb = "../data/bond_orders/Conformer3D_CID_2031_withoutH.pdb"
        s3 = top.Segment(filecoord=filenamepdb, filetop=filenamepdb)
        s3._topology.perception_rings()
        self.assertTrue(s3._topology._cycles == [[8, 13, 17, 21, 19, 15, 8],
                                                 [9, 14, 18, 22, 20, 16, 9]])

        # A cyclic molecule
        filenamepdb = "../data/bond_orders/Conformer3D_CID_8135_withoutH.pdb"
        s4 = top.Segment(filecoord=filenamepdb, filetop=filenamepdb)
        s4._topology.perception_rings()
        self.assertTrue(s4._topology._cycles == [[7, 3, 2, 6, 4, 0, 1, 5, 7]])

        # Condensed rings
        filenamepdb = "../data/bond_orders/Conformer3D_CID_3632_withoutH.pdb"
        s5 = top.Segment(filecoord=filenamepdb, filetop=filenamepdb)

        s5._topology.perception_rings()
        self.assertTrue(s5._topology._cycles == [[9, 11, 4, 5, 3, 6, 9],
                                                 [10, 8, 3, 5, 4, 7, 10],
                                                 [10, 7, 4, 11, 9, 6, 3, 8, 10],
                                                 [15, 14, 12, 8, 10, 1, 15],
                                                 [15, 14, 12, 8, 3, 5, 4, 7, 10, 1, 15],
                                                 [15, 14, 12, 8, 3, 6, 9, 11, 4, 7, 10, 1, 15],
                                                 [20, 18, 15, 14, 17, 19, 20],
                                                 [20, 18, 15, 1, 10, 8, 12, 14, 17, 19, 20],
                                                 [20, 18, 15, 1, 10, 7, 4, 5, 3, 8, 12, 14, 17, 19, 20],
                                                 [20, 18, 15, 1, 10, 7, 4, 11, 9, 6, 3, 8, 12, 14, 17, 19, 20]])

        # Condensed rings --> Norbornane
        filenamepdb = "../data/bond_orders/Conformer3D_CID_9233_withoutH.pdb"
        s6 = top.Segment(filecoord=filenamepdb, filetop=filenamepdb)

        s6._topology.perception_rings()
        self.assertTrue(s6._topology._cycles == [[1, 2, 0, 3, 5, 1],
                                                 [1, 2, 0, 4, 6, 1],
                                                 [1, 5, 3, 0, 4, 6, 1]])

        # Condensed rings
        filenamepdb = "../data/bond_orders/Conformer3D_CID_931_withoutH.pdb"
        s7 = top.Segment(filecoord=filenamepdb, filetop=filenamepdb)

        s7._topology.perception_rings()
        self.assertTrue(s7._topology._cycles == [[1, 0, 2, 6, 7, 3, 1],
                                                 [1, 0, 4, 8, 9, 5, 1],
                                                 [1, 3, 7, 6, 2, 0, 4, 8, 9, 5, 1]])

        m = "\t============== END   TEST_05 ================================"
        print(m) if self.log is None else self.log.info(m)

    #########################################################################
    def test_06_copy_topology(self):

        """
        Copy a topology
        Returns
        -------

        """
        m = "\n\t============== START TEST_06 ================================\n"
        m += "\t              Test Copy topology object                     "
        print(m) if self.log is None else self.log.info(m)

        # One aromatic ring --> D-glucose
        filenamepdb = "../data/bond_orders/Conformer3D_CID_5793.pdb"
        s1 = top.Segment(filecoord=filenamepdb, filetop=filenamepdb)
        s1._topology.perception_rings()

        t_original = s1._topology
        t_copy = copy(s1._topology)

        self.assertTrue(t_original == t_copy)

        # _natoms, _nringsCauchy, _undirected
        self.assertEqual(t_original.natoms, t_copy.natoms)
        self.assertEqual(t_original._nringsCauchy, t_copy._nringsCauchy)
        self.assertEqual(t_original._undirected, t_copy._undirected)
        t_copy.natoms += 1
        t_copy._nringsCauchy += 1
        t_copy._undirected = not t_original._undirected
        self.assertNotEqual(t_original.natoms, t_copy.natoms)
        self.assertNotEqual(t_original._nringsCauchy, t_copy._nringsCauchy)
        self.assertNotEqual(t_original._undirected, t_copy._undirected)

        # _bonds ======
        self.assertEqual(t_original._bonds, t_copy._bonds)
        t_copy._bonds.append({11, 12})
        self.assertNotEqual(t_original._bonds, t_copy._bonds)

        # _cycles ======
        self.assertEqual(t_original._cycles, t_copy._cycles)
        t_copy._cycles.append([11, 12, 11])
        self.assertNotEqual(t_original._cycles, t_copy._cycles)

        # _nmols ======
        self.assertEqual(t_original._nmols, t_copy._nmols)
        t_copy._nmols.append([11, 12, 11])
        self.assertNotEqual(t_original._nmols, t_copy._nmols)

        # _graphdict =======
        self.assertEqual(t_original._graphdict, t_copy._graphdict)
        t_copy._graphdict['a'] = [10]
        self.assertNotEqual(t_original._graphdict, t_copy._graphdict)

        # _orderbonds =======
        try:
            np.testing.assert_equal(t_original._orderbonds, t_copy._orderbonds)
            res = True
        except AssertionError as err:
            res = False
            print(err)
        self.assertTrue(res)

        t_copy._orderbonds = np.zeros([1, 1])
        t_copy._orderbonds[0, 0] = 10000
        res = np.array_equal(t_original._orderbonds, t_copy._orderbonds)
        self.assertFalse(res)

        m = "\t============== END   TEST_06 ================================"
        print(m) if self.log is None else self.log.info(m)

    #########################################################################
    def test_07_generate_topology_fromPSF_UA(self):

        m = "\n\t============== START TEST_07 ================================\n"
        m += "\t       Test Generate_topology fromPSF object. United Atom    \n"
        m += "\t       Assign bond order FALSE"
        print(m) if self.log is None else self.log.info(m)

        filepsf = "../data/0003Ch-C020-002br04/namd.psf"
        t = top.Topology(logger=self.log)
        t.get_bonds_topologyPSF(filepsf)

        try:
            import pygraphviz as pgv
            t.draw_graph_forest_pygraphviz(title="../data/0003Ch-C020-002br04/graph01_test05_")

        except ModuleNotFoundError:
            try:
                import networkx as nx
                t.draw_graph_networkx(title="../data/0003Ch-C020-002br04/graph01_test05_")
                t.draw_graph_forest_networkx(title="./data/0003Ch-C020-002br04/graph01_test05_")
            except ModuleNotFoundError:
                m = "\ttest_05 is not run. The pygraphviz and/or networkx module is not installed."
                print(m) if self.log is None else self.log.error(m)
                pass

        ll = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 20, 21],
              [24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 46, 47, 44, 45],
              [48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 70, 71, 68, 69]]

        edges = [{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {8, 7},
                 {8, 9}, {9, 10}, {9, 20}, {10, 11}, {11, 12}, {12, 13}, {13, 14},
                 {14, 15}, {16, 15}, {22, 15}, {16, 17}, {17, 18}, {18, 19}, {20, 21},
                 {22, 23}, {24, 25}, {25, 26}, {26, 27}, {27, 28}, {28, 29}, {29, 30},
                 {30, 31}, {32, 31}, {32, 33}, {33, 34}, {33, 44}, {34, 35}, {35, 36},
                 {36, 37}, {37, 38}, {38, 39}, {40, 39}, {46, 39}, {40, 41}, {41, 42},
                 {42, 43}, {44, 45}, {46, 47}, {48, 49}, {49, 50}, {50, 51}, {51, 52},
                 {52, 53}, {53, 54}, {54, 55}, {56, 55}, {56, 57}, {57, 58}, {57, 68},
                 {58, 59}, {59, 60}, {60, 61}, {61, 62}, {62, 63}, {64, 63}, {70, 63},
                 {64, 65}, {65, 66}, {66, 67}, {68, 69}, {70, 71}]
        self.assertEqual(t._nmols, ll)
        self.assertEqual(t.natoms, 72)
        self.assertEqual(t.get_edges(), edges)

        m = "\t============== END   TEST_07 ================================"
        print(m) if self.log is None else self.log.info(m)

    #########################################################################
    def test_08_generate_topology_fromPSF_UA_bo(self):

        m = "\n\t============== START TEST_08 ================================\n"
        m += "\t       Test Generate_topology fromPSF object. United Atom    \n"
        m += "\t       Assign bond order TRUE"
        print(m) if self.log is None else self.log.info(m)

        filepsf = "../data/0003Ch-C020-002br04/namd.psf"
        t = top.Topology(logger=self.log)
        t.get_bonds_topologyPSF(filepsf, assign_bo=True)

        try:
            import pygraphviz as pgv
            t.draw_graph_forest_pygraphviz(title="./graphs/test08_graph01_")

        except ModuleNotFoundError:
            try:
                import networkx as nx
                t.draw_graph_networkx(title="./graphs/test08_graph01_")
                t.draw_graph_forest_networkx(title="./graphs/test08_graph01_")
            except ModuleNotFoundError:
                m = "\ttest_08 is not run. The pygraphviz and/or networkx module is not installed."
                print(m) if self.log is None else self.log.error(m)
                pass

        ll = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 20, 21],
              [24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 46, 47, 44, 45],
              [48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 70, 71, 68, 69]]

        edges = [{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {8, 7},
                 {8, 9}, {9, 10}, {9, 20}, {10, 11}, {11, 12}, {12, 13}, {13, 14},
                 {14, 15}, {16, 15}, {22, 15}, {16, 17}, {17, 18}, {18, 19}, {20, 21},
                 {22, 23}, {24, 25}, {25, 26}, {26, 27}, {27, 28}, {28, 29}, {29, 30},
                 {30, 31}, {32, 31}, {32, 33}, {33, 34}, {33, 44}, {34, 35}, {35, 36},
                 {36, 37}, {37, 38}, {38, 39}, {40, 39}, {46, 39}, {40, 41}, {41, 42},
                 {42, 43}, {44, 45}, {46, 47}, {48, 49}, {49, 50}, {50, 51}, {51, 52},
                 {52, 53}, {53, 54}, {54, 55}, {56, 55}, {56, 57}, {57, 58}, {57, 68},
                 {58, 59}, {59, 60}, {60, 61}, {61, 62}, {62, 63}, {64, 63}, {70, 63},
                 {64, 65}, {65, 66}, {66, 67}, {68, 69}, {70, 71}]
        self.assertEqual(t._nmols, ll)
        self.assertEqual(t.natoms, 72)
        self.assertEqual(t.get_edges(), edges)

        bo = {0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0, 5: 1.0, 6: 1.0, 7: 1.0, 8: 1.0, 9: 1.0,
              19: 1.0, 10: 1.0, 20: 1.0, 11: 1.0, 12: 1.0, 13: 1.0, 14: 1.0, 15: 1.0, 21: 1.0,
              16: 1.0, 22: 1.0, 17: 1.0, 18: 1.0, 23: 1.0, 24: 1.0, 25: 1.0, 26: 1.0, 27: 1.0,
              28: 1.0, 29: 1.0, 30: 1.0, 31: 1.0, 32: 1.0, 42: 1.0, 33: 1.0, 43: 1.0, 34: 1.0,
              35: 1.0, 36: 1.0, 37: 1.0, 38: 1.0, 44: 1.0, 39: 1.0, 45: 1.0, 40: 1.0, 41: 1.0,
              46: 1.0, 47: 1.0, 48: 1.0, 49: 1.0, 50: 1.0, 51: 1.0, 52: 1.0, 53: 1.0, 54: 1.0,
              55: 1.0, 65: 1.0, 56: 1.0, 66: 1.0, 57: 1.0, 58: 1.0, 59: 1.0, 60: 1.0, 61: 1.0,
              67: 1.0, 62: 1.0, 68: 1.0, 63: 1.0, 64: 1.0}
        self.assertEqual(t._orderbonds, bo)
        self.assertEqual(t._totalcharge_mol, [0, 0, 0])

        m = "\t============== END   TEST_08 ================================"
        print(m) if self.log is None else self.log.info(m)

    #########################################################################
    def test_09_generate_topology_fromPSF_AA_bo(self):

        m = "\n\t============== START TEST_09 ================================\n"
        m += "\t       Test Generate_topology fromPSF object. Atomistic model \n"
        m += "\t       Assign bond order TRUE"
        print(m) if self.log is None else self.log.info(m)

        filepsf = "../data/20-iPP-OPLS-AA/namd_withimproper.psf"
        t = top.Topology(logger=self.log)
        t.get_bonds_topologyPSF(filepsf, assign_bo=True)

        try:
            import pygraphviz as pgv
            t.draw_graph_forest_pygraphviz(title="./graphs/test09_graph01_")

        except ModuleNotFoundError:
            try:
                import networkx as nx
                t.draw_graph_networkx(title="./graphs/test09_graph01_")
                t.draw_graph_forest_networkx(title="./graphs/test09_graph01_")
            except ModuleNotFoundError:
                m = "\ttest_09 is not run. The pygraphviz and/or networkx module is not installed."
                print(m) if self.log is None else self.log.error(m)
                pass

        ll = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 41, 42, 43, 44, 45, 46, 39, 40, 34, 35, 36, 37, 38, 32, 33, 27, 28, 29,
              30, 31, 25, 26, 20, 21, 22, 23, 24, 18, 19, 13, 14, 15, 16, 17, 10, 11, 12],
              [47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 88, 89, 90, 91, 92, 93, 86, 87, 81, 82, 83, 84, 85, 79, 80,
              74, 75, 76, 77, 78, 72, 73, 67, 68, 69, 70, 71, 65, 66, 60, 61, 62, 63, 64, 57, 58, 59]]

        edges = [{0, 1}, {0, 10}, {0, 11}, {0, 12}, {1, 2}, {1, 13}, {1, 17}, {2, 3}, {2, 18}, {2, 19}, {3, 4},
                 {3, 20}, {24, 3}, {4, 5}, {25, 4}, {26, 4}, {5, 6}, {27, 5}, {5, 31}, {6, 7}, {32, 6}, {33, 6}, {8, 7},
                 {34, 7}, {38, 7}, {8, 9}, {8, 39}, {8, 40}, {9, 41}, {9, 45}, {9, 46}, {13, 14}, {13, 15}, {16, 13},
                 {20, 21}, {20, 22}, {20, 23}, {27, 28}, {27, 29}, {27, 30}, {34, 35}, {34, 36}, {34, 37}, {41, 42},
                 {41, 43}, {41, 44}, {48, 47}, {57, 47}, {58, 47}, {59, 47}, {48, 49}, {48, 60}, {48, 64}, {49, 50},
                 {65, 49}, {49, 66}, {50, 51}, {50, 67}, {50, 71}, {51, 52}, {72, 51}, {73, 51}, {52, 53}, {74, 52},
                 {52, 78}, {53, 54}, {53, 79}, {80, 53}, {54, 55}, {81, 54}, {85, 54}, {56, 55}, {86, 55}, {87, 55},
                 {56, 88}, {56, 92}, {56, 93}, {60, 61}, {60, 62}, {60, 63}, {67, 68}, {67, 69}, {67, 70}, {74, 75},
                 {74, 76}, {74, 77}, {81, 82}, {81, 83}, {81, 84}, {88, 89}, {88, 90}, {88, 91}]
        self.assertEqual(t._nmols, ll)
        self.assertEqual(t.natoms, 94)
        self.assertEqual(t.get_edges(), edges)

        bo = {0: 1.0, 9: 1.0, 10: 1.0, 11: 1.0, 1: 1.0, 12: 1.0, 16: 1.0, 2: 1.0, 17: 1.0, 18: 1.0,
              13: 1.0, 14: 1.0, 15: 1.0, 3: 1.0, 19: 1.0, 23: 1.0, 4: 1.0, 24: 1.0, 25: 1.0, 20: 1.0,
              21: 1.0, 22: 1.0, 5: 1.0, 26: 1.0, 30: 1.0, 6: 1.0, 31: 1.0, 32: 1.0, 27: 1.0, 28: 1.0,
              29: 1.0, 7: 1.0, 33: 1.0, 37: 1.0, 8: 1.0, 38: 1.0, 39: 1.0, 34: 1.0, 35: 1.0, 36: 1.0,
              40: 1.0, 44: 1.0, 45: 1.0, 41: 1.0, 42: 1.0, 43: 1.0, 46: 1.0, 55: 1.0, 56: 1.0, 57: 1.0,
              47: 1.0, 58: 1.0, 62: 1.0, 48: 1.0, 63: 1.0, 64: 1.0, 59: 1.0, 60: 1.0, 61: 1.0, 49: 1.0,
              65: 1.0, 69: 1.0, 50: 1.0, 70: 1.0, 71: 1.0, 66: 1.0, 67: 1.0, 68: 1.0, 51: 1.0, 72: 1.0,
              76: 1.0, 52: 1.0, 77: 1.0, 78: 1.0, 73: 1.0, 74: 1.0, 75: 1.0, 53: 1.0, 79: 1.0, 83: 1.0,
              54: 1.0, 84: 1.0, 85: 1.0, 80: 1.0, 81: 1.0, 82: 1.0, 86: 1.0, 90: 1.0, 91: 1.0, 87: 1.0,
              88: 1.0, 89: 1.0}
        self.assertEqual(t._orderbonds, bo)
        self.assertEqual(t._totalcharge_mol, [0, 0])

        m = "\t============== END   TEST_09 ================================"
        print(m) if self.log is None else self.log.info(m)

    #########################################################################
    def test_10_directTopology(self):

        m = "\n\t============== START TEST_10 ================================\n"
        m += "\t             Building a Topology object                      \n"
        print(m) if self.log is None else self.log.info(m)

        t = top.Topology(natoms=4, listbonds=[(0, 1), (2, 3)])

        try:
            import pygraphviz as pgv
            t.draw_graph_pygraphviz(title="graphs/test10_topo01")

        except ModuleNotFoundError:
            try:
                import networkx as nx
                t.draw_graph_networkx(title="graphs/test10_topo01")
            except ModuleNotFoundError:
                print("test_10 is not run. The pygraphviz and/or networkx module is not installed.")
                pass

        t = top.Topology(natoms=6, listbonds=[(0, 1), (1, 5), (1, 6), (2, 3)])
        try:
            import pygraphviz as pgv
            t.draw_graph_pygraphviz(title="graphs/test10_topo02")

        except ModuleNotFoundError:
            try:
                import networkx as nx
                t.draw_graph_networkx(title="graphs/test10_topo02")
            except ModuleNotFoundError:
                print("test_10 is not run. The pygraphviz and/or networkx module is not installed.")
                pass

        m = "\t============== END   TEST_10 ================================"
        print(m) if self.log is None else self.log.info(m)

    #########################################################################
    def test_11_getarraymols(self):

        m = "\n\t============== START TEST_11 ================================\n"
        m += "\t              Checking get_array_mols_neigh                   \n"
        print(m) if self.log is None else self.log.info(m)

        t = top.Topology(natoms=4, listbonds=[(0, 1), (2, 3)])

        nmols_array, l_neigh_array = t.get_array_mols_neigh()
        nmols_array_test = [[0, 1], [2, 3]]
        l_neigh_array_test = [[1], [0], [3], [2]]
        np.testing.assert_array_equal(nmols_array_test, nmols_array)
        np.testing.assert_array_equal(l_neigh_array_test, l_neigh_array)

        t = top.Topology(natoms=6, listbonds=[(0, 1), (1, 5), (1, 6), (2, 3)])
        nmols_array, l_neigh_array = t.get_array_mols_neigh()
        nmols_array_test = [[0, 1, 5], [2, 3, -1], [4, -1, -1]]
        l_neigh_array_test = [[1, -1], [0, 5], [3, -1], [2, -1], [-1, -1], [1, -1]]
        np.testing.assert_array_equal(nmols_array_test, nmols_array)
        np.testing.assert_array_equal(l_neigh_array_test, l_neigh_array)

        m = "\t============== END   TEST_11 ================================"
        print(m) if self.log is None else self.log.info(m)

    #########################################################################
    def test_12_checktopologyPDBconnect_atomistic_bo(self):

        m = "\n\t============== START TEST_12 ================================\n"
        m += "\t              Topology from  PDB CONECT. All atomistic        \n"
        m += "\t              with bond orders                                \n"
        print(m) if self.log is None else self.log.info(m)

        t = top.Topology()
        fnamepdb = "../data/n-hexane.pdb"
        t.get_bonds_topologyCONNECTPDB(filenamePDB=fnamepdb, assign_bo=True)

        t.draw_graph_networkx(title="graphs/test12_topo01")

        ll = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
        edges = [{0, 1}, {0, 3}, {0, 4}, {0, 10}, {1, 2}, {1, 5}, {1, 6}, {2, 7},
                 {8, 2}, {9, 2}, {10, 11}, {10, 13}, {10, 14}, {11, 12}, {11, 15},
                 {16, 11}, {17, 12}, {18, 12}, {19, 12}]
        e = ['C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H']

        lp = sorted(t._nmols[0])
        self.assertEqual(lp, ll)
        self.assertEqual(t.natoms, 20)
        self.assertEqual(t.get_edges(), edges)
        self.assertEqual(t.elements, e)
        self.assertEqual(t._totalcharge_mol, [0])
        for key, val in t._orderbonds.items():
            self.assertEqual(val, 1.0, msg="Order bonds: {}".format(key))
        ib = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
        self.assertEqual(ib, t._isbackbone)

        m = "\t============== END   TEST_12 ================================"
        print(m) if self.log is None else self.log.info(m)

    #########################################################################
    def test_13_checktopologyPDBconnect_UA_bo(self):

        m = "\n\t============== START TEST_13 ================================\n"
        m += "\t              Topology from  PDB CONECT. United atom          \n"
        m += "\t              with bond orders                                \n"
        m += "\t          Bond order cannot be assigned to UA model          \n"
        print(m) if self.log is None else self.log.info(m)

        t = top.Topology()
        fnamepdb = "../data/n-hexane_UA.pdb"
        t.get_bonds_topologyCONNECTPDB(filenamePDB=fnamepdb, assign_bo=True, is_unitedatom=True)

        t.draw_graph_networkx(title="graphs/test13_topo01")

        ll = [0, 1, 2, 3, 4, 5]
        edges = [{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}]
        e = ['C', 'C', 'C', 'C', 'C', 'C']

        lp = sorted(t._nmols[0])
        self.assertEqual(lp, ll)
        self.assertEqual(t.natoms, 6)
        self.assertEqual(t.get_edges(), edges)
        self.assertEqual(t.elements, e)
        self.assertEqual(t._totalcharge_mol, [0])
        # Order bonds cannot be assigned in united atoms models.
        self.assertIsNone(t._orderbonds)
        ib = [1, 1, 1, 1, 1, 1]
        self.assertEqual(ib, t._isbackbone)

        m = "\t============== END   TEST_13 ================================"
        print(m) if self.log is None else self.log.info(m)

    #########################################################################
    def test_14_checktopologyMDAnalysis_UA_tpr_nobo(self):

        m = "============== START TEST_14 ================================\n"
        m += "          Topology from MDAnalysis. UA model                 \n"
        m += "          without bond order                                 \n"
        print(m) if self.log is None else self.log.info(m)

        filename_tpr = "../data/0003Ch-C020-002br04/RUN-001/topol.tpr"

        # Setup topology
        t = top.Topology()
        t.get_bonds_topologyMDAnalysis(filename_tpr, assign_bo=False)
        t.draw_graph_pygraphviz(title="graphs/test14_topo01")

        ll = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 20, 21],
              [24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 46, 47, 44, 45],
              [48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 70, 71, 68, 69]]

        edges = [{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {8, 7},
                 {8, 9}, {9, 10}, {9, 20}, {10, 11}, {11, 12}, {12, 13}, {13, 14},
                 {14, 15}, {16, 15}, {22, 15}, {16, 17}, {17, 18}, {18, 19}, {20, 21},
                 {22, 23}, {24, 25}, {25, 26}, {26, 27}, {27, 28}, {28, 29}, {29, 30},
                 {30, 31}, {32, 31}, {32, 33}, {33, 34}, {33, 44}, {34, 35}, {35, 36},
                 {36, 37}, {37, 38}, {38, 39}, {40, 39}, {46, 39}, {40, 41}, {41, 42},
                 {42, 43}, {44, 45}, {46, 47}, {48, 49}, {49, 50}, {50, 51}, {51, 52},
                 {52, 53}, {53, 54}, {54, 55}, {56, 55}, {56, 57}, {57, 58}, {57, 68},
                 {58, 59}, {59, 60}, {60, 61}, {61, 62}, {62, 63}, {64, 63}, {70, 63},
                 {64, 65}, {65, 66}, {66, 67}, {68, 69}, {70, 71}]

        e = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
             'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
             'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
             'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
             'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']

        self.assertEqual(t._nmols, ll)
        self.assertEqual(t.natoms, 72)
        self.assertEqual(t.get_edges(), edges)
        self.assertEqual(t.elements, e)
        self.assertEqual(t._totalcharge_mol, [0, 0, 0])
        self.assertIsNone(t._orderbonds)
        # self.assertEqual(ib, t._isbackbone)
        # print(t)

        m = "\t============== END   TEST_14 ================================"
        print(m) if self.log is None else self.log.info(m)

    #########################################################################
    def test_15_checktopologyMDAnalysis_UA_tpr_bo(self):

        m = "\t============== START TEST_15 ================================\n"
        m += "\t          Topology from MDAnalysis. UA model                 \n"
        m += "\t          with bond order                                    \n"
        m += "\t          Bond order cannot be assigned to UA model          \n"
        print(m) if self.log is None else self.log.info(m)

        filename_tpr = "../data/0003Ch-C020-002br04/RUN-001/topol.tpr"

        # Setup topology
        t = top.Topology()
        t.get_bonds_topologyMDAnalysis(filename_tpr, assign_bo=True, is_unitedatom=True)
        t.draw_graph_pygraphviz(title="graphs/test15_topo01")

        ll = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 20, 21],
              [24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 46, 47, 44, 45],
              [48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 70, 71, 68, 69]]

        edges = [{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {8, 7},
                 {8, 9}, {9, 10}, {9, 20}, {10, 11}, {11, 12}, {12, 13}, {13, 14},
                 {14, 15}, {16, 15}, {22, 15}, {16, 17}, {17, 18}, {18, 19}, {20, 21},
                 {22, 23}, {24, 25}, {25, 26}, {26, 27}, {27, 28}, {28, 29}, {29, 30},
                 {30, 31}, {32, 31}, {32, 33}, {33, 34}, {33, 44}, {34, 35}, {35, 36},
                 {36, 37}, {37, 38}, {38, 39}, {40, 39}, {46, 39}, {40, 41}, {41, 42},
                 {42, 43}, {44, 45}, {46, 47}, {48, 49}, {49, 50}, {50, 51}, {51, 52},
                 {52, 53}, {53, 54}, {54, 55}, {56, 55}, {56, 57}, {57, 58}, {57, 68},
                 {58, 59}, {59, 60}, {60, 61}, {61, 62}, {62, 63}, {64, 63}, {70, 63},
                 {64, 65}, {65, 66}, {66, 67}, {68, 69}, {70, 71}]

        e = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
             'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
             'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
             'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
             'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']

        self.assertEqual(t._nmols, ll)
        self.assertEqual(t.natoms, 72)
        self.assertEqual(t.get_edges(), edges)
        self.assertEqual(t.elements, e)
        self.assertEqual(t._totalcharge_mol, [0, 0, 0])
        self.assertIsNone(t._orderbonds)

        m = "\t============== END   TEST_15 ================================"
        print(m) if self.log is None else self.log.info(m)

    #########################################################################
    def test_16_assign_backbone_atoms_psf(self):

        m = "\t============== START TEST_16 ================================\n"
        m += "\t          Assign backbone from PSF. UA model                 \n"
        print(m) if self.log is None else self.log.info(m)

        # Setup Trajectory to analyze
        filename_psf = "../data/0003Ch-C020-002br04/namd_out.psf"
        t = top.Topology(logger=self.log)
        t.get_bonds_topologyPSF(filename_psf, assign_bo=True)

        # Backbone atoms from PSF file
        t.assign_backbone_atoms_psf()

        ll = [True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
              True, True, True, False, False, False, False, True, True, True, True, True, True, True, True, True, True,
              True, True, True, True, True, True, True, True, True, True, False, False, False, False, True, True, True,
              True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
              False, False, False, False]

        self.assertEqual(t._isbackbone, ll)

        m = "\t============== END   TEST_16 ================================"
        print(m) if self.log is None else self.log.info(m)

    #########################################################################
    def test_17_canonize_hexane_atomistic(self):

        m = "\n\t============== START TEST_17 ================================\n"
        m += "\t              Topology from  PDB CONECT. All atomistic        \n"
        m += "\t              with bond orders                                \n"
        print(m) if self.log is None else self.log.info(m)

        t1 = top.Topology()
        fnamepdb = "../data/n-hexane.pdb"
        t1.get_bonds_topologyCONNECTPDB(filenamePDB=fnamepdb, assign_bo=True)

        t1.draw_graph_pygraphviz(title="graphs/test17_topo01")

        order1 = t1.canonize_atom_ordering()

        t2 = top.Topology()
        fnamepdb = "../data/n-hexane_other.pdb"
        t2.get_bonds_topologyCONNECTPDB(filenamePDB=fnamepdb, assign_bo=True)

        t2.draw_graph_pygraphviz(title="graphs/test17_topo02")

        order2 = t2.canonize_atom_ordering()

        norder1_t = (7, 8, 9, 17, 18, 19, 5, 6, 15, 16, 3, 4, 13, 14, 2, 12, 1, 11, 0, 10)
        norder2_t = (10, 11, 12, 17, 18, 19, 8, 9, 15, 16, 6, 7, 13, 14, 0, 5, 1, 4, 2, 3)
        #
        self.assertEqual(norder1_t, order1)
        self.assertEqual(norder2_t, order2)

        # Basic test for the canonical order
        for idx in range(len(norder1_t)):
            idx1 = norder1_t[idx]
            idx2 = norder2_t[idx]
            self.assertEqual(t1.elements[idx1], t2.elements[idx2])
            self.assertEqual(len(t1.get_neighbours(idx1)), len(t2.get_neighbours(idx2)))

        m = "\t============== END   TEST_17 ================================"
        print(m) if self.log is None else self.log.info(m)

    #########################################################################
    def test_18_write_psf(self):

        m = "\n\t============== START TEST_18 ================================\n"
        m += "\t          PSF Topology from  PDB. All atomistic        \n"
        m += "\t              with bond orders                                \n"
        print(m) if self.log is None else self.log.info(m)

        t1 = top.Topology()
        fnamepdb = "../data/n-hexane.pdb"
        t1.get_bonds_topologyCONNECTPDB(filenamePDB=fnamepdb, assign_bo=True)

        t1.writepsf(filename_psf="./psf/test_18.psf")

        m = "\t============== END   TEST_18 ================================"
        print(m) if self.log is None else self.log.info(m)

    #########################################################################
    def test_19_write_psf_with_impropers(self):

        m = "\n\t============== START TEST_19 ================================\n"
        m += "\t          PSF Topology from  GROMACS files. CG       \n"
        m += "\t              without bond orders                                \n"
        print(m) if self.log is None else self.log.info(m)

        t1 = top.Topology()
        fnametopo = "../data/0003Ch-C020-002br04/RUN-001/topol.tpr"
        fnamecoord = "../data/0003Ch-C020-002br04/RUN-001/confout.gro"
        t1.get_bonds_topologyMDAnalysis(filenametopo=fnametopo, filecoord=fnamecoord,
                                        assign_bo=False, is_unitedatom=True)

        improper_file = "../data/0003Ch-C020-002br04/improper.ndx"
        improper_list = t1.get_allimpropers_from_file(improper_file)
        t1.writepsf(filename_psf="./psf/test_19.psf", improper_angles=improper_list)

        m = "\t============== END   TEST_19 ================================"
        print(m) if self.log is None else self.log.info(m)

    #########################################################################
    def test_20_write_psf_with_impropers(self):

        m = "\n\t============== START TEST_20 ================================\n"
        m += "\t          PSF Topology from  GROMACS files. CG       \n"
        m += "\t              without bond orders                                \n"
        print(m) if self.log is None else self.log.info(m)

        xtc1 = "../data/0003Ch-C020-002br04/RUN-001/traj_comp.xtc"
        xtc2 = "../data/0003Ch-C020-002br04/RUN-002/traj_comp.part0002.xtc"
        xtc3 = "../data/0003Ch-C020-002br04/RUN-003/traj_comp.part0003.xtc"
        fnametopo = "../data/0003Ch-C020-002br04/namd_out.psf"

        self.trj_small_psf = top.ExtTrajectory([xtc1, xtc2, xtc3], topfile=fnametopo, logger=self.log)

        m = "\t============== END   TEST_20 ================================"
        print(m) if self.log is None else self.log.info(m)

    #########################################################################
    def test_21_renumbering_by_head(self):

        m = "\n\t============== START TEST_21 ================================\n"
        m += "\t          Given an atom number as head renumbering the atoms  \n"
        m += "\t          folowing the graph                                  \n"
        print(m) if self.log is None else self.log.info(m)

        m = "\t\t01 MOLECULE P3HB_3mon_01c.pdb"
        print(m) if self.log is None else self.log.info(m)
        filenamepdb = "../data/P3HB_3mon_01c.pdb"
        pdb = ReadPdbFormat(filenamepdb)
        test = pdb.write_renumber_pdb(head_idx_atom=[4],  tail_idx_atom=[36])
        dd = {0: 4, 1: 37, 2: 0, 3: 6, 4: 3, 5: 9, 6: 10, 7: 11, 8: 1, 9: 7, 10: 8, 11: 2, 12: 5, 13: 28, 14: 24,
              15: 30, 16: 27, 17: 33, 18: 34, 19: 35, 20: 25, 21: 31, 22: 32, 23: 26, 24: 29, 25: 16, 26: 12,
              27: 18, 28: 15, 29: 21, 30: 22, 31: 23, 32: 13, 33: 19, 34: 20, 35: 14, 36: 17, 37: 36, 38: 38}
        self.assertEqual(test, dd)

        m = "\t\t02 MOLECULE iPP_8mon_1Ch.xsd"
        print(m) if self.log is None else self.log.info(m)
        filenamexsd = "../data/iPP_8mon_1Ch.xsd"
        xsd = ReadXsdFormat(filenamexsd)
        xsd.write_pdb(filename_pdb="iPP_8mon_1Ch.pdb")
        pdb = ReadPdbFormat("./iPP_8mon_1Ch.pdb")
        test = pdb.write_renumber_pdb(head_idx_atom=[0], tail_idx_atom=[58])
        dd = {0: 0, 1: 9, 2: 10, 3: 11, 4: 1, 5: 12, 6: 8, 7: 25, 8: 26, 9: 72,
              10: 2, 11: 13, 12: 14, 13: 4, 14: 18, 15: 3, 16: 15, 17: 16, 18: 17,
              19: 5, 20: 19, 21: 20, 22: 7, 23: 24, 24: 6, 25: 21, 26: 22, 27: 23,
              28: 35, 29: 52, 30: 53, 31: 28, 32: 39, 33: 27, 34: 36, 35: 37, 36: 38,
              37: 29, 38: 40, 39: 41, 40: 31, 41: 45, 42: 30, 43: 42, 44: 43, 45: 44,
              46: 32, 47: 46, 48: 47, 49: 34, 50: 51, 51: 33, 52: 48, 53: 49, 54: 50,
              55: 59, 56: 70, 57: 71, 58: 55, 59: 63, 60: 54, 61: 60, 62: 61, 63: 62,
              64: 56, 65: 64, 66: 65, 67: 58, 68: 69, 69: 73, 70: 57, 71: 66, 72: 67, 73: 68}
        self.assertEqual(test, dd)

        m = "\t\t03 MOLECULE PS_3mon_1Ch.xsd"
        print(m) if self.log is None else self.log.info(m)
        filenamexsd = "../data/PS_3mon_1Ch.xsd"
        xsd = ReadXsdFormat(filenamexsd)
        xsd.write_pdb(filename_pdb="PS_3mon_1Ch.pdb")
        pdb = ReadPdbFormat("./PS_3mon_1Ch.pdb")
        test = pdb.write_renumber_pdb(head_idx_atom=[13], tail_idx_atom=[44])
        print(test)
        dd = {0: 13, 1: 14, 2: 15,  3: 16, 4: 11,  5: 12, 6: 0, 7: 9, 8: 7, 9: 5, 10: 3,
              11: 1,  12: 2, 13: 10, 14: 4, 15: 8, 16: 6, 17: 30, 18: 31, 19: 32, 20: 28,
              21: 29,  22: 17, 23: 26, 24: 24, 25: 22, 26: 20, 27: 18, 28: 19, 29: 27, 30: 21,
              31: 25, 32: 23, 33: 47, 34: 48, 35: 49, 36: 44, 37: 45, 38: 46, 39: 33, 40: 42, 41: 40,
              42: 38, 43: 36, 44: 34, 45: 35, 46: 43, 47: 37, 48: 41, 49: 39}
        self.assertEqual(test, dd)

        m = "\t============== END   TEST_21 ================================"
        print(m) if self.log is None else self.log.info(m)

    #########################################################################
    def test_22_renumbering_by_head_severalchains(self):

        m = "\n\t============== START TEST_22 ================================\n"
        m += "\t          Given a list of heads renumbering the atoms         \n"
        m += "\t          following the graph                                 \n"
        print(m) if self.log is None else self.log.info(m)

        filenamepdb = "../data/P4HB_crystal_3Ch.pdb"
        pdb = ReadPdbFormat(filenamepdb)

        head_idx_atomlist = [16, 64, 100]
        tail_idx_atomlist = [84, 132, 24]

        test = pdb.write_renumber_pdb(head_idx_atom=head_idx_atomlist, tail_idx_atom=tail_idx_atomlist)
        print(test)

        dd = {0: 16, 1: 144, 2: 3, 3: 10, 4: 11, 5: 2, 6: 8, 7: 9, 8: 1, 9: 6, 10: 7, 11: 0, 12: 5, 13: 4,
              14: 15, 15: 22, 16: 23, 17: 14, 18: 20, 19: 21, 20: 13, 21: 18, 22: 19, 23: 12, 24: 17, 25: 88,
              26: 75, 27: 82, 28: 83, 29: 74, 30: 80, 31: 81, 32: 73, 33: 78, 34: 79, 35: 72, 36: 77, 37: 76,
              38: 87, 39: 94, 40: 95, 41: 86, 42: 92, 43: 93, 44: 85, 45: 90, 46: 91, 47: 84, 48: 147, 49: 89,
              50: 64, 51: 146, 52: 51, 53: 58, 54: 59, 55: 50, 56: 56, 57: 57, 58: 49, 59: 54, 60: 55, 61: 48,
              62: 53, 63: 52, 64: 63, 65: 70, 66: 71, 67: 62, 68: 68, 69: 69, 70: 61, 71: 66, 72: 67, 73: 60,
              74: 65, 75: 136, 76: 123, 77: 130, 78: 131, 79: 122, 80: 128, 81: 129, 82: 121, 83: 126, 84: 127,
              85: 120, 86: 125, 87: 124, 88: 135, 89: 142, 90: 143, 91: 134, 92: 140, 93: 141, 94: 133, 95: 138,
              96: 139, 97: 132, 98: 149, 99: 137, 100: 100, 101: 148, 102: 111, 103: 118, 104: 119, 105: 110,
              106: 116, 107: 117, 108: 109, 109: 114, 110: 115, 111: 108, 112: 113, 113: 112, 114: 99, 115: 106,
              116: 107, 117: 98, 118: 104, 119: 105, 120: 97, 121: 102, 122: 103, 123: 96, 124: 101, 125: 28,
              126: 39, 127: 46, 128: 47, 129: 38, 130: 44, 131: 45, 132: 37, 133: 42, 134: 43, 135: 36, 136: 41,
              137: 40, 138: 27, 139: 34, 140: 35, 141: 26, 142: 32, 143: 33, 144: 25, 145: 30, 146: 31, 147: 24,
              148: 145, 149: 29}

        self.assertEqual(test, dd)

        m = "\t============== END   TEST_22 ================================"
        print(m) if self.log is None else self.log.info(m)

    ########################################################################
    def test_23_assignresidues_fromlist(self):

        m = "\n\t============== START TEST_23 ================================\n"
        m += "\t       Using an external file assign residues.   \n"
        m += "\t       Only one chain and one type of chain.   \n"
        print(m) if self.log is None else self.log.info(m)

        filenamepdb = "../data/P3HB_3mon_01c.pdb"
        headinfo_file = "./P3HB_3mon_01c_headtail.dat"
        residueinfo_file = "./P3HB_3mon_01c_residues.dat"
        pdb = ReadPdbFormat(filenamepdb)
        headlist, taillist = pdb.read_head_tail_info(headinfo_file)
        test = pdb.write_renumber_pdb(head_idx_atom=headlist, tail_idx_atom=taillist)
        dd = {0: 4, 1: 37, 2: 0, 3: 6, 4: 3, 5: 9, 6: 10, 7: 11, 8: 1, 9: 7, 10: 8, 11: 2, 12: 5, 13: 28, 14: 24,
              15: 30, 16: 27, 17: 33, 18: 34, 19: 35, 20: 25, 21: 31, 22: 32, 23: 26, 24: 29, 25: 16, 26: 12,
              27: 18, 28: 15, 29: 21, 30: 22, 31: 23, 32: 13, 33: 19, 34: 20, 35: 14, 36: 17, 37: 36, 38: 38}
        self.assertEqual(test, dd)

        filenamepdb = "./P3HB_3mon_01c_renumber.pdb"
        pdb_new = ReadPdbFormat(filenamepdb)
        pdb_new.assign_residues_chains(residueinfo_file)

        m = "\t============== END   TEST_23 ================================"
        print(m) if self.log is None else self.log.info(m)

    ########################################################################
    def test_24_assignresidues_fromlist(self):

        m = "\n\t============== START TEST_24 ================================\n"
        m += "\t       Using an external file assign residues.   \n"
        m += "\t       Several chains and one type of chain.   \n"
        print(m) if self.log is None else self.log.info(m)

        filenamepdb = "../data/P4HB_crystal_2x2x2.pdb"
        headinfo_file = "./P4HB_crystal_2x2x2_headtail.dat"
        residueinfo_file = "./P4HB_crystal_2x2x2_residues.dat"
        pdb = ReadPdbFormat(filenamepdb)
        headlist, taillist = pdb.read_head_tail_info(headinfo_file)
        test = pdb.write_renumber_pdb(head_idx_atom=headlist, tail_idx_atom=taillist)
        print(test)
        dd = {0: 364, 1: 399, 2: 375, 3: 382, 4: 383, 5: 374, 6: 380, 7: 381, 8: 373, 9: 378, 10: 379,
              11: 372, 12: 377, 13: 376, 14: 363, 15: 370, 16: 371, 17: 362, 18: 368, 19: 369, 20: 361,
              21: 366, 22: 367, 23: 360, 24: 365, 25: 172, 26: 183, 27: 190, 28: 191, 29: 182, 30: 188,
              31: 189, 32: 181, 33: 186, 34: 187, 35: 180, 36: 185, 37: 184, 38: 171, 39: 178, 40: 179,
              41: 170, 42: 176, 43: 177, 44: 169, 45: 174, 46: 175, 47: 168, 48: 391, 49: 173, 50: 268,
              51: 395, 52: 279, 53: 286, 54: 287, 55: 278, 56: 284, 57: 285, 58: 277, 59: 282, 60: 283,
              61: 276, 62: 281, 63: 280, 64: 267, 65: 274, 66: 275, 67: 266, 68: 272, 69: 273, 70: 265,
              71: 270, 72: 271, 73: 264, 74: 269, 75: 76, 76: 87, 77: 94, 78: 95, 79: 86, 80: 92, 81: 93,
              82: 85, 83: 90, 84: 91, 85: 84, 86: 89, 87: 88, 88: 75, 89: 82, 90: 83, 91: 74, 92: 80,
              93: 81, 94: 73, 95: 78, 96: 79, 97: 72, 98: 387, 99: 77, 100: 64, 101: 386, 102: 51, 103: 58,
              104: 59, 105: 50, 106: 56, 107: 57, 108: 49, 109: 54, 110: 55, 111: 48, 112: 53, 113: 52,
              114: 63, 115: 70, 116: 71, 117: 62, 118: 68, 119: 69, 120: 61, 121: 66, 122: 67, 123: 60,
              124: 65, 125: 256, 126: 243, 127: 250, 128: 251, 129: 242, 130: 248, 131: 249, 132: 241,
              133: 246, 134: 247, 135: 240, 136: 245, 137: 244, 138: 255, 139: 262, 140: 263, 141: 254,
              142: 260, 143: 261, 144: 253, 145: 258, 146: 259, 147: 252, 148: 394, 149: 257, 150: 160,
              151: 390, 152: 147, 153: 154, 154: 155, 155: 146, 156: 152, 157: 153, 158: 145, 159: 150,
              160: 151, 161: 144, 162: 149, 163: 148, 164: 159, 165: 166, 166: 167, 167: 158, 168: 164,
              169: 165, 170: 157, 171: 162, 172: 163, 173: 156, 174: 161, 175: 352, 176: 339, 177: 346,
              178: 347, 179: 338, 180: 344, 181: 345, 182: 337, 183: 342, 184: 343, 185: 336, 186: 341,
              187: 340, 188: 351, 189: 358, 190: 359, 191: 350, 192: 356, 193: 357, 194: 349, 195: 354,
              196: 355, 197: 348, 198: 398, 199: 353, 200: 220, 201: 393, 202: 231, 203: 238, 204: 239,
              205: 230, 206: 236, 207: 237, 208: 229, 209: 234, 210: 235, 211: 228, 212: 233, 213: 232,
              214: 219, 215: 226, 216: 227, 217: 218, 218: 224, 219: 225, 220: 217, 221: 222, 222: 223,
              223: 216, 224: 221, 225: 28, 226: 39, 227: 46, 228: 47, 229: 38, 230: 44, 231: 45, 232: 37,
              233: 42, 234: 43, 235: 36, 236: 41, 237: 40, 238: 27, 239: 34, 240: 35, 241: 26, 242: 32,
              243: 33, 244: 25, 245: 30, 246: 31, 247: 24, 248: 385, 249: 29, 250: 316, 251: 397, 252: 327,
              253: 334, 254: 335, 255: 326, 256: 332, 257: 333, 258: 325, 259: 330, 260: 331, 261: 324,
              262: 329, 263: 328, 264: 315, 265: 322, 266: 323, 267: 314, 268: 320, 269: 321, 270: 313,
              271: 318, 272: 319, 273: 312, 274: 317, 275: 124, 276: 135, 277: 142, 278: 143, 279: 134,
              280: 140, 281: 141, 282: 133, 283: 138, 284: 139, 285: 132, 286: 137, 287: 136, 288: 123,
              289: 130, 290: 131, 291: 122, 292: 128, 293: 129, 294: 121, 295: 126, 296: 127, 297: 120,
              298: 389, 299: 125, 300: 16, 301: 384, 302: 3, 303: 10, 304: 11, 305: 2, 306: 8, 307: 9, 308: 1,
              309: 6, 310: 7, 311: 0, 312: 5, 313: 4, 314: 15, 315: 22, 316: 23, 317: 14, 318: 20, 319: 21,
              320: 13, 321: 18, 322: 19, 323: 12, 324: 17, 325: 208, 326: 195, 327: 202, 328: 203, 329: 194,
              330: 200, 331: 201, 332: 193, 333: 198, 334: 199, 335: 192, 336: 197, 337: 196, 338: 207,
              339: 214, 340: 215, 341: 206, 342: 212, 343: 213, 344: 205, 345: 210, 346: 211, 347: 204,
              348: 392, 349: 209, 350: 112, 351: 388, 352: 99, 353: 106, 354: 107, 355: 98, 356: 104,
              357: 105, 358: 97, 359: 102, 360: 103, 361: 96, 362: 101, 363: 100, 364: 111, 365: 118,
              366: 119, 367: 110, 368: 116, 369: 117, 370: 109, 371: 114, 372: 115, 373: 108, 374: 113,
              375: 304, 376: 291, 377: 298, 378: 299, 379: 290, 380: 296, 381: 297, 382: 289, 383: 294,
              384: 295, 385: 288, 386: 293, 387: 292, 388: 303, 389: 310, 390: 311, 391: 302, 392: 308,
              393: 309, 394: 301, 395: 306, 396: 307, 397: 300, 398: 396, 399: 305}

        self.assertEqual(test, dd)

        filenamepdb = "./P4HB_crystal_2x2x2_renumber.pdb"
        pdb_new = ReadPdbFormat(filenamepdb)
        pdb_new.assign_residues_chains(residueinfo_file)

        m = "\t============== END   TEST_24 ================================"
        print(m) if self.log is None else self.log.info(m)

    ########################################################################
    def test_25_generate_topology_fromDATALammps_AA(self):

        m = "\n\t============== START TEST_25 ================================\n"
        m += "\t       Test Generate_topology fromDATA object (Lammps). All Atom    \n"
        m += "\t       Assign bond order FALSE"
        print(m) if self.log is None else self.log.info(m)

        filepsf = "../data/noctane_3_1_1.data"
        t = top.Topology(logger=self.log)
        t.get_bonds_topologyMDAnalysis(filepsf)

        try:
            import pygraphviz as pgv
            t.draw_graph_forest_pygraphviz(title="../data/25-Topology_lammps/graph01_test25_")

        except ModuleNotFoundError:
            try:
                import networkx as nx
                t.draw_graph_networkx(title="../data/25-Topology_lammps/graph01_test25_")
                t.draw_graph_forest_networkx(title="./data/25-Topology_lammps/graph01_test25_")
            except ModuleNotFoundError:
                m = "\ttest_25 is not run. The module pygraphviz and/or networkx is not installed."
                print(m) if self.log is None else self.log.error(m)
                pass

        ll = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
              [26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51],
              [52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77]]

        edges = [{0, 1}, {0, 2}, {0, 3}, {0, 4}, {4, 5}, {4, 6}, {4, 7}, {8, 7}, {9, 7},
                 {10, 7}, {10, 11}, {10, 12}, {10, 13}, {13, 14}, {13, 15}, {16, 13}, {16, 17},
                 {16, 18}, {16, 19}, {19, 20}, {19, 21}, {19, 22}, {22, 23}, {24, 22}, {25, 22},
                 {26, 27}, {26, 28}, {26, 29}, {26, 30}, {30, 31}, {32, 30}, {33, 30}, {33, 34},
                 {33, 35}, {33, 36}, {36, 37}, {36, 38}, {36, 39}, {40, 39}, {41, 39}, {42, 39},
                 {42, 43}, {42, 44}, {42, 45}, {45, 46}, {45, 47}, {48, 45}, {48, 49}, {48, 50},
                 {48, 51}, {52, 53}, {52, 54}, {52, 55}, {56, 52}, {56, 57}, {56, 58}, {56, 59},
                 {59, 60}, {59, 61}, {59, 62}, {62, 63}, {64, 62}, {65, 62}, {65, 66}, {65, 67},
                 {65, 68}, {68, 69}, {68, 70}, {68, 71}, {72, 71}, {73, 71}, {74, 71}, {74, 75},
                 {74, 76}, {74, 77}]

        self.assertEqual(t._nmols, ll)
        self.assertEqual(t.natoms, 78)
        self.assertEqual(t.get_edges(), edges)

        m = "\t============== END   TEST_25 ================================"
        print(m) if self.log is None else self.log.info(m)


if __name__ == '__main__':
    unittest.main()
