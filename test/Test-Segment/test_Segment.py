import unittest
import numpy as np
from copy import copy
from topology.Segment import Segment

class SegmentTests(unittest.TestCase):

    # ##################################################################################################################
    def setUp(self):

        print("=============== SEGMENT TEST =================")

   # ##################################################################################################################
    def test_00_empty_segment(self):

        """Test the creation of an empty Segment

        Returns
        -------
        None

        """
        print("Test 0 : Create an empty segment and minimun information segemnt.")
        print("\t"+self.id())

        s1 = Segment()

        self.assertEqual(s1._natoms, 0)
        self.assertEqual(s1._topology, None)

        # Exception if the file is in an unkwown format
        with self.assertRaises(Exception):
            s2 = Segment(natoms=10)

        with self.assertRaises(Exception):
            s3 = Segment(natoms=5, elementlist=['C', 'H', 'H', 'H', 'H'])

        # s4 = Segment(natoms=5, elementlist=['C', 'H', 'H', 'H', 'H'], xlist=[0.03205128,  0.32460314, 0.32462156,  0.32462156, -1.10205128],
        #                                                               ylist=[0.16025641, -0.84855359, 0.66465460,  0.66465460,  0.16026959],
        #                                                               zlist=[0.00000000,  0.00000000, 0.87365150, -0.87365150,  0.00000000])
        s4 = Segment(natoms=5, elementlist=['C', 'H', 'H', 'H', 'H'], xlist=[0.0320,  0.324, 0.324,  0.324, -1.102],
                                                                      ylist=[0.1602, -0.848, 0.664,  0.664,  0.160],
                                                                      zlist=[0.0000,  0.000, 0.874, -0.873,  0.000])


        self.assertEqual(s4._natoms, 5)
        self.assertNotEqual(s4._topology, None)

    # ##################################################################################################################
    def test_01_ethylene_UA_manual(self):

        """Test the creation of a Segment object manually

        Returns
        -------
        None

        """
        print("Test 1 : Create an ethylene segment UA model.")
        print("\t"+self.id())

        s1 = Segment(natoms=2, xlist=[0.0, -0.002], ylist=[0.765, -0.777], zlist=[0.0, 0.0], elementlist=['C','C'])

        self.assertEqual(s1._natoms, 2)

        result = [[0.0, 0.765, 0.0],
                  [-0.002, -0.777, 0.]]
        try:
            np.testing.assert_almost_equal(s1.get_coords(), result)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)

        self.assertEqual(s1._topology.elements,['C','C'])
        self.assertEqual(s1._topology._bonds, [{0,1}])

    # ##################################################################################################################
    def test_02_ethylene_UA_fromPDB(self):

        """Test the creation of a Segment object from PDB file

        Returns
        -------
        None

        """
        print("Test 2 : Create an ethylene segment UA model from a PDB file.")
        print("\t"+self.id())

        s1 = Segment(filecoord="../data/pe_1mon_ua.pdb")

        # Excpetion if the file does not exist
        with self.assertRaises(Exception):
            s1 = Segment(filecoord="../data/pe_1mon_uaa.pdb")

        # Exception if the file is in an unkwown format
        with self.assertRaises(Exception):
            s1 = Segment(filecoord="../data/pe_1mon_ua.psf")

        self.assertEqual(s1._natoms, 2)

        result = [[0.0, 0.765, 0.0],
                  [-0.002, -0.777, 0.]]
        try:
            np.testing.assert_almost_equal(s1.get_coords(), result)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)

        self.assertEqual(s1._topology.elements,['C','C'])
        self.assertEqual(s1._topology._bonds, [{0,1}])

    # ##################################################################################################################
    def test_03_ethylene_UA_fromXYZ(self):

        """Test the creation of a Segment object from XYZ file

        Returns
        -------
        None

        """
        print("Test 3 : Create an ethylene segment UA model from a XYZ file.")
        print("\t"+self.id())

        s1 = Segment(filecoord="../data/pe_1mon_ua.xyz")

        self.assertEqual(s1._natoms, 2)

        result = [[0.0, 0.765, 0.0],
                  [-0.002, -0.777, 0.]]
        try:
            np.testing.assert_almost_equal(s1.get_coords(), result)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)

        self.assertEqual(s1._topology.elements,['C','C'])
        self.assertEqual(s1._topology._bonds, [{0,1}])

    ################################################################################################################
    def test_04_ethylene_UA_fromGRO(self):

        """Test the creation of a Segment object from GRO file

        Returns
        -------
        None

        """
        print("Test 4 : Create an ethylene segment UA model from a GRO file.")
        print("\t"+self.id())

        s1 = Segment(filecoord="../data/pe_1mon_ua.gro")

        self.assertEqual(s1._natoms, 2)

        result = [[0.0, 0.76, 0.0],
                  [-0.00, -0.78, 0.]]
        try:
            np.testing.assert_almost_equal(s1.get_coords(), result)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)

        self.assertEqual(s1._topology.elements,['C','C'])
        self.assertEqual(s1._topology._bonds, [{0,1}])

    ################################################################################################################
    def test_05_ethylene_AA_fromPDB_withPDBtopology_guesstopology(self):

        """Test the creation of a Segment object from PDB file and PDB file.
        Guess the topology in base of distances between atoms

        Returns
        -------
        None

        """

        print("Test 5 : Create an ethylene segment AA model from a PDB file "
              "without \"CONECT\" and compare with a \"CONECT\".")
        print("\t"+self.id())

        # Without CONNECT
        s1 = Segment(filecoord="../data/pe_1mon_aa.pdb", filetop="../data/pe_1mon_aa_noconnect.pdb")
        # With CONNECT
        s2 = Segment(filecoord="../data/pe_1mon_aa.pdb", filetop="../data/pe_1mon_aa.pdb")

        self.assertEqual(s1._natoms, 8)
        self.assertEqual(s2._natoms, 8)

        result = [[0.0, 0.765, 0.0],
                  [0.513, 1.165, 0.887],
                  [0.513, 1.165, -0.887],
                  [-0.002, -0.777, 0.000],
                  [-0.515, -1.177, 0.887],
                  [1.021, -1.179, 0.000],
                  [-0.515, -1.177, -0.887],
                  [-1.024, 1.167, 0.]]

        try:
            np.testing.assert_almost_equal(s1.get_coords(), result)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)

        partial_result = [ [ 0.000,  0.765,  0.000],
                           [ 0.513,  1.165, -0.887]]

        try:
            np.testing.assert_almost_equal(s1.get_coords([0,2]), partial_result)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)

        try:
            np.testing.assert_almost_equal(s2.get_coords(), result)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)

        self.assertEqual(s1._topology.elements,['C','H','H','C','H','H','H','H'])
        self.assertEqual(s2._topology.elements,['C','H','H','C','H','H','H','H'])

        topo_res = {0:[1,2,3,7], 1:[0], 2:[0],
                    3:[0,4,5,6], 4:[3], 5:[3], 6:[3], 7:[0]}

        self.assertEqual(s1._topology._graphdict, topo_res)
        self.assertEqual(s2._topology._graphdict, topo_res)

    # ################################################################################################################
    def test_06_center_of_mass_and_translate_to_origin(self):

        print("Test 6 : Check the calculation of the center of mass and translation to the origin")
        print("\t"+self.id())

        s1 = Segment(filecoord="../data/propane_1mon_aa.pdb", filetop="../data/propane_1mon_aa.pdb")

        # Part 1 --> com
        com_result = [0.92694, -0.16618, 0.24637]
        com = s1.center_of_mass()
        try:
            np.testing.assert_almost_equal(com, com_result, 5)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)

        # Part 2 --> cog
        cog_result = [0.90400, -0.15582, 0.26927]
        cog = s1.center_of_geom()
        try:
            np.testing.assert_almost_equal(cog, cog_result, 5)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)

        # Part 3 --> Translate to the origin
        s1.translate_vector(-com)
        com_result = [0.0, 0.0, 0.0]
        try:
            np.testing.assert_almost_equal(s1.center_of_mass(), com_result, 10)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)

    # ################################################################################################################
    def test_07_translate_along_vector(self):

        print("Test 7 : Translate along the vector [2,2,2]")
        print("\t"+self.id())

        s1 = Segment(filecoord="../data/propane_1mon_aa.pdb", filetop="../data/propane_1mon_aa.pdb")
        n = [2.0, 2.0, 2.0]
        s1.translate_vector(n)
        d = s1.get_coords()

        result = [[1.774, 2.058, 2.022],
                  [1.659, 3.094, 2.372],
                  [1.259, 1.408, 2.743],
                  [1.244, 1.973, 1.063],
                  [3.257, 1.675, 1.88 ],
                  [3.372, 0.639, 1.53 ],
                  [3.773, 2.325, 1.159],
                  [3.773, 1.758, 2.814],
                  [4.8  , 1.485, 2.682],
                  [3.32 , 1.105, 3.532],
                  [3.713, 2.766, 3.165]]

        try:
            np.testing.assert_almost_equal(d, result)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)

    # ################################################################################################################
    def test_08_change_orientation(self):

        print("Test 8 : Orient Segment generating Euler angles")
        print("\t"+self.id())

        s1 = Segment(filecoord="../data/propane_1mon_aa.pdb", filetop="../data/propane_1mon_aa.pdb")
        # sim = Simulator(mol=(s1,), box=CubicBox(14.0, 14.0, 14.0))
        # sim.write_PDB_simulator(filenamePDB="simulator_append_orient.pdb")
        # cog = s1.center_of_geom()
        # s1.translate_vector(-cog)
        # sim.write_PDB_simulator(filenamePDB="simulator_append_orient.pdb", with_vmd_tcl=True, appendpdb=True)

        s1.euler_orientation(iseed=987655)
        # sim.write_PDB_simulator(filenamePDB="simulator_append_orient.pdb", with_vmd_tcl=True, appendpdb=True)

        r  = [[-0.06088824, -0.06608958,  0.21644581],
              [-1.0904413,   0.2388534,   0.45319734],
              [ 0.59720215,  0.380485,    0.97476908],
              [-0.00526952, -1.15758847,  0.33195658],
              [ 0.34111606,  0.3696992, -1.20292823],
              [ 1.37066912,  0.06475622, -1.43967976],
              [-0.31695929, -0.07648652, -1.96217268],
              [ 0.28672872,  1.43267833, -1.31525984],
              [ 0.57213627,  1.70326286, -2.31125609],
              [ 0.9486838 ,  1.90051976, -0.61520136],
              [-0.71398482,  1.75919019, -1.1286388 ]]

        try:
            np.testing.assert_almost_equal(s1.get_coords(), r)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)

    # ################################################################################################################
    def test_09_copy_segment(self):

        print("Test 9 : Copy Segment")
        print("\t"+self.id())

        s_original = Segment(filecoord="../data/propane_1mon_aa.pdb", filetop="../data/propane_1mon_aa.pdb")
        s_copy = copy(s_original)

        self.assertTrue(s_original == s_copy)

        # _natoms
        self.assertEqual(s_original._natoms, s_copy._natoms)
        s_copy._natoms += 1
        self.assertNotEqual(s_original._natoms, s_copy._natoms)

        # _topology
        self.assertTrue(s_copy._topology == s_original._topology)
        s_copy._topology.natoms += 1
        self.assertFalse(s_copy._topology == s_original._topology)

        # _coords
        try:
            np.testing.assert_equal(s_original._coords, s_copy._coords)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)
        s_copy._coords[0,0] = 10000
        res = np.array_equal(s_original._coords, s_copy._coords)
        self.assertFalse(res)

        # _elements
        try:
            np.testing.assert_equal(s_original._topology.elements, s_copy._topology.elements)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)
        s_copy._topology.elements[0] = 'Z'
        res = np.array_equal(s_original._topology.elements, s_copy._topology.elements)
        self.assertFalse(res)

    ########################################################################
    def test_10_assign_bonds(self):

        print("Test 10: Assign_bonds to the segment")
        print("\t"+self.id())

        # One aromatic ring --> Nitrobenzene ===========================================================
        filenamepdb = "../data/bond_orders/Conformer3D_CID_7416.pdb"
        s1 = Segment(filecoord=filenamepdb, filetop=filenamepdb)
        s1._topology.perception_rings()
        s1._topology.guess_nringsCauchy()

        try:
            import pygraphviz as pgv
            s1._topology.draw_graph_pygraphviz("graphs/test_10a_assignbonds")
        except ModuleNotFoundError:
            try:
                import networkx as nx
                s1._topology.draw_graph_networkx(title="graphs/test_10a_assignbonds_nx")
            except ModuleNotFoundError:
                print("test_10 is not run. The pygraphviz and/or networkx module is not installed.")
                pass

        s1._topology.assign_bond_orders()

        bo_test = {0: 1.5, 1: 1.5, 2: 1.0, 3: 1.5, 4: 1.5,
                   5: 1.5, 6: 1.0, 7: 1.5, 8: 1.0, 9: 1.5,
                   10: 1.0, 11: 1.5, 13: 1.0, 12: 1.0}

        try:
            np.testing.assert_equal(s1._topology._orderbonds, bo_test)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)

        # One non-aromatic ring --> D-glucose ===========================================================
        filenamepdb="../data/bond_orders/Conformer3D_CID_5793.pdb"
        s2 = Segment(filecoord=filenamepdb, filetop=filenamepdb)
        s2._topology.perception_rings()
        s2._topology.guess_nringsCauchy()

        try:
            import pygraphviz as pgv
            s2._topology.draw_graph_pygraphviz("graphs/test_10b_assignbonds")
        except ModuleNotFoundError:
            try:
                import networkx as nx
                s2._topology.draw_graph_networkx(title="graphs/test_10b_assignbonds_nx")
            except ModuleNotFoundError:
                print("test_10 is not run. The pygraphviz and/or networkx module is not installed.")
                pass

        s2._topology.assign_bond_orders()

        bo_test = {0: 1.0, 1: 1.0 , 13: 1.0, 17: 1.0, 18: 1.0,
                   19: 1.0, 8: 1.0, 21: 1.0, 2: 1.0, 12: 1.0,
                   14: 1.0, 10: 1.0, 22: 1.0, 23: 1.0, 3: 1.0,
                   4: 1.0, 15: 1.0, 16: 1.0, 5: 1.0, 6: 1.0,
                   20: 1.0, 7: 1.0, 9: 1.0, 11: 1.0}

        try:
            np.testing.assert_equal(s2._topology._orderbonds, bo_test)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)

        # Two condensed aromatic rings --> Naphthalene ===========================================================
        # This test is not enabled due to for
        filenamepdb="../data/bond_orders/Conformer3D_CID_931.pdb"
        s3 = Segment(filecoord=filenamepdb, filetop=filenamepdb)
        s3._topology.perception_rings()
        s3._topology.guess_nringsCauchy()

        try:
            import pygraphviz as pgv
            s3._topology.draw_graph_pygraphviz("graphs/test_10c_assignbonds")
        except ModuleNotFoundError:
            try:
                import networkx as nx
                s3._topology.draw_graph_networkx(title="graphs/test_10c_assignbonds_nx")
            except ModuleNotFoundError:
                print("test_10 is not run. The pygraphviz and/or networkx module is not installed.")
                pass

        s3._topology.assign_bond_orders()

        bo_test = {0: 1.5, 1: 1.5, 2: 1.5, 3: 1.5, 4: 1.5, 5: 1.5, 6: 1.0,
                   9: 1.5, 10: 1.0, 7: 1.5, 8: 1.0, 11: 1.5, 12: 1.0, 13: 1.5,
                   15: 1.0, 14: 1.0, 16: 1.5, 18: 1.0, 17: 1.0}

        try:
            np.testing.assert_equal(s3._topology._orderbonds, bo_test)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)

    # ################################################################################################################
    def test_11_volume_segment(self):

        print("Test 11: Volume Segment")
        print("\t"+self.id())

        s = Segment(filecoord="../data/n-hexane.pdb", filetop="../data/n-hexane.pdb")
        v_vdw, v_tsar = s.calc_vdw_volume_VABC()
        self.assertAlmostEqual(v_vdw, 112.36)
        self.assertAlmostEqual(v_tsar, 90.18036)

        s = Segment(filecoord="../data/nitrobenzene.pdb", filetop="../data/nitrobenzene.pdb")
        v_vdw, v_tsar = s.calc_vdw_volume_VABC()
        self.assertAlmostEqual(v_vdw, 107.12)
        self.assertAlmostEqual(v_tsar, 85.98312)

        s = Segment(filecoord="../data/bond_orders/Conformer3D_CID_12398.pdb", filetop="../data/bond_orders/Conformer3D_CID_12398.pdb")
        v_vdw, v_tsar = s.calc_vdw_volume_VABC()
        self.assertAlmostEqual(v_vdw, 302.66)
        self.assertAlmostEqual(v_tsar, 242.61066)

        s = Segment(filecoord="../data/bond_orders/Conformer3D_CID_7966.pdb", filetop="../data/bond_orders/Conformer3D_CID_7966.pdb")
        v_vdw, v_tsar = s.calc_vdw_volume_VABC()
        self.assertAlmostEqual(v_vdw, 108.79)
        self.assertAlmostEqual(v_tsar, 87.32079)

        s = Segment(filecoord="../data/bond_orders/Conformer3D_CID_7222.pdb", filetop="../data/bond_orders/Conformer3D_CID_7222.pdb")
        v_vdw, v_tsar = s.calc_vdw_volume_VABC()
        self.assertAlmostEqual(v_vdw, 102.09)
        self.assertAlmostEqual(v_tsar, 81.95409)

    # ################################################################################################################
    def test_12_assign_okuwaki(self):

        print("Test 12: Assign type atoms in Segment")
        print("\t"+self.id())

        s1 = Segment(filecoord="../data/n-hexane.pdb",
                     filetop="../data/n-hexane.pdb")
        self.assertEqual(s1._typeelements, None)

        s1.set_typeatoms("../data/n-hexane_types.dat")
        type_atoms_test = ['c3','c3', 'c3', 'hc', 'hc', 'hc', 'hc', 'hc',
                           'hc', 'hc', 'c3', 'c3', 'c3', 'hc',
                           'hc', 'hc', 'hc', 'hc', 'hc', 'hc']
        try:
            np.testing.assert_equal(s1._typeelements, type_atoms_test)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)

        s2 = Segment(filecoord="../data/n-hexane.pdb",
                     filetop="../data/n-hexane.pdb",
                     filetypeatoms="../data/n-hexane_types.dat")

        try:
            np.testing.assert_equal(s2._typeelements, type_atoms_test)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)

    # ################################################################################################################
    def test_13_nitrobenzene_AA_fromPDB_withPDBtopology_guesstopology(self):

        """
        Test the creation of a Segment object from PDB file with Conect
        and without conect records PDB file.
        Guess the topology in base of distances between atoms

        Returns
        -------
        None

        """

        print("Test 13 : Create an nitrobenzene segment AA model from a PDB file "
              "without \"CONECT\" and compare with a \"CONECT\".")
        print("\t"+self.id())

        # Without CONNECT
        s1 = Segment(filecoord="../data/nitrobenzene_noconnect.pdb", filetop="../data/nitrobenzene_noconnect.pdb")
        # With CONNECT
        s2 = Segment(filecoord="../data/nitrobenzene.pdb", filetop="../data/nitrobenzene.pdb")

        self.assertEqual(s1._natoms, 14)
        self.assertEqual(s2._natoms, 14)

        result = [[-0.152,   0.662,  -0.001],
                  [ 0.540,  -0.550,  -0.000],
                  [ 1.931,  -0.558,   0.001],
                  [ 2.602,   0.662,   0.002],
                  [ 1.931,   1.882,   0.001],
                  [ 0.540,   1.874,  -0.000],
                  [-1.241,   0.662,  -0.001],
                  [-0.004,  -1.493,  -0.000],
                  [ 2.503 , -1.482,   0.001],
                  [ 2.503,   2.806,   0.001],
                  [-0.004,   2.816,  -0.000],
                  [ 4.079,   0.662,   0.004],
                  [ 4.652,   1.759,  -0.003],
                  [ 4.652,  -0.435,  -0.003]]


        try:
            np.testing.assert_almost_equal(s1.get_coords(), result)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)

        try:
            np.testing.assert_almost_equal(s2.get_coords(), result)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)

        self.assertEqual(s1._topology.elements,['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'N', 'O', 'O'])
        self.assertEqual(s2._topology.elements,['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'N', 'O', 'O'])

        topo_res = {0:[1, 5, 6], 1:[0, 2, 7], 2:[1, 3, 8],
                    3:[2, 4, 11], 4:[3, 5, 9], 5:[0, 4, 10], 6:[0], 7:[1],
                    8:[2], 9:[4], 10:[5], 11:[3, 12, 13], 12:[11], 13:[11]}

        self.assertEqual(s1._topology._graphdict, topo_res)
        self.assertEqual(s2._topology._graphdict, topo_res)



    # ##################################################################################################################
    def tearDown(self):

        print("============= END SEGMENT TEST ===============")

if __name__ == '__main__':

    unittest.main(verbosity=2)



