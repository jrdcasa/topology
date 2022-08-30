import unittest
import numpy as np
from topology.MolecularGraph import MolecularGraph


class MolecularTests(unittest.TestCase):

    # ##################################################################################################################
    def setUp(self):

        """ g1 --------
           (0) -- (1) -- (2) -- (5) -- (6)
                   |
                  (3)
                   |
                  (4)

        """
        self.g1 = MolecularGraph(undirected=True)

        # Add Vertex
        for i in range(0, 7):
            self.g1.add_vertex(i)

        # Add get_edges
        self.g1.add_edge([0, 1])
        self.g1.add_edge([1, 2])
        self.g1.add_edge([1, 3])
        self.g1.add_edge([2, 5])
        self.g1.add_edge([3, 4])
        self.g1.add_edge([5, 6])

        """ g2 --------
        
                         (7) -- (8) -- (9)
                          |
           (0) -- (1) -- (2) -- (5) -- (6)
                   |             |
                  (3)           (10)
                   |             |
                  (4)           (11)

        """
        self.g2 = MolecularGraph(undirected=True)

        # Add Vertex
        for i in range(0, 12):
            self.g2.add_vertex(i)

        # Add get_edges
        self.g2.add_edge([0, 1])
        self.g2.add_edge([1, 2])
        self.g2.add_edge([1, 3])
        self.g2.add_edge([2, 5])
        self.g2.add_edge([3, 4])
        self.g2.add_edge([5, 6])
        self.g2.add_edge([2, 7])
        self.g2.add_edge([7, 8])
        self.g2.add_edge([8, 9])
        self.g2.add_edge([5, 10])
        self.g2.add_edge([10, 11])

        """ g3 --------
        
                         (7) -- (8) -- (9)--\
                          |                  (12)
           (0) -- (1) -- (2) -- (5) -- (6)--/
                   |             |
                  (3)           (10)
                   |             |
                  (4)           (11)

        """
        self.g3 = MolecularGraph(undirected=True)

        # Add Vertex
        for i in range(0, 13):
            self.g3.add_vertex(i)

        # Add get_edges
        self.g3.add_edge([0, 1])
        self.g3.add_edge([1, 2])
        self.g3.add_edge([1, 3])
        self.g3.add_edge([2, 5])
        self.g3.add_edge([3, 4])
        self.g3.add_edge([5, 6])
        self.g3.add_edge([5, 10])
        self.g3.add_edge([2, 7])
        self.g3.add_edge([7, 8])
        self.g3.add_edge([8, 9])
        self.g3.add_edge([10, 11])
        self.g3.add_edge([9, 12])
        self.g3.add_edge([6, 12])

        """ g4 --------
        
                         
                                           
           (0) -- (1) -- (2) -- (3) -- (4) -- (5)
                   
                  (6) -- (7) -- (8) -- (9)
                          |      |
                         (10)   (11)
            (12)
        """
        self.g4 = MolecularGraph(undirected=True)
        
        # Add Vertex
        for i in range(0, 13):
            self.g4.add_vertex(i)

        # Add get_edges
        self.g4.add_edge([0, 1])
        self.g4.add_edge([1, 2])
        self.g4.add_edge([2, 3])
        self.g4.add_edge([3, 4])
        self.g4.add_edge([4, 5])

        self.g4.add_edge([6, 7])
        self.g4.add_edge([7, 8])
        self.g4.add_edge([7, 10])
        self.g4.add_edge([8, 9])
        self.g4.add_edge([8, 11])

        """ g5 --------
        
                         
                                           
           (0) -- (1) -- (2) -- (3) -- (4) -- (5) -- (13)
                   
                  (6) -- (7) -- (8) -- (9) -- (14)
                          |      |
                         (10)   (11)
            (12)
        """
        self.g5 = MolecularGraph(undirected=True)

        # Add Vertex
        for i in range(0, 15):
            self.g5.add_vertex(i)

        # Add get_edges
        self.g5.add_edge([0, 1])
        self.g5.add_edge([1, 2])
        self.g5.add_edge([2, 3])
        self.g5.add_edge([3, 4])
        self.g5.add_edge([4, 5])
        self.g5.add_edge([5, 13])

        self.g5.add_edge([6, 7])
        self.g5.add_edge([7, 8])
        self.g5.add_edge([7, 10])
        self.g5.add_edge([8, 9])
        self.g5.add_edge([8, 11])
        self.g5.add_edge([9, 14])

        """ g6 --------
                                           
           (0) -- (1) -- (2) -- (3) -- (4) -- (5)
                   
                  (6) -- (7) -- (8) -- (9)
                          |      |
                         (10)   (11)
            (12)
        """

        self.g6 = MolecularGraph(nvert=13, undirected=True,
                                 listbonds=[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5],
                                            [6, 7], [7, 8], [7, 10], [8, 9], [8, 11]])
        """ g7 --------
        
                         (7) -- (8) -- (9)--\
                          |                  (12)
           (0) -- (1) -- (2) -- (5) -- (6)--/
                   |             |
                  (3)           (10)
                   |             |
                  (4)           (11)

            (13) -- (14) -- (15) -- (16)
                             |
                            (17)
        """
        self.g7 = MolecularGraph(undirected=True)

        # Add Vertex
        for i in range(0, 18):
            self.g7.add_vertex(i)

        # Add get_edges
        self.g7.add_edge([0, 1])
        self.g7.add_edge([1, 2])
        self.g7.add_edge([1, 3])
        self.g7.add_edge([2, 5])
        self.g7.add_edge([3, 4])
        self.g7.add_edge([5, 6])
        self.g7.add_edge([5, 10])
        self.g7.add_edge([2, 7])
        self.g7.add_edge([7, 8])
        self.g7.add_edge([8, 9])
        self.g7.add_edge([10, 11])
        self.g7.add_edge([9, 12])
        self.g7.add_edge([6, 12])
        self.g7.add_edge([13, 14])
        self.g7.add_edge([14, 15])
        self.g7.add_edge([15, 17])
        self.g7.add_edge([15, 16])

        """ g8 --------

                             (6) -- (7)
                              |        
        (0) -- (1) -- (2) -- (3)
                       |      |
                      (4) -- (5) -- (9) -- (10) -- (11)
                       |                     \     /
                      (8)                      (12)
        """
        self.g8 = MolecularGraph(undirected=True)

        # Add Vertex
        for i in range(0, 13):
            self.g8.add_vertex(i)

        # Add get_edges
        self.g8.add_edge([0, 1])
        self.g8.add_edge([1, 2])
        self.g8.add_edge([2, 3])
        self.g8.add_edge([2, 4])
        self.g8.add_edge([3, 5])
        self.g8.add_edge([3, 6])
        self.g8.add_edge([4, 5])
        self.g8.add_edge([4, 8])
        self.g8.add_edge([5, 9])
        self.g8.add_edge([6, 7])
        self.g8.add_edge([9, 10])
        self.g8.add_edge([10, 11])
        self.g8.add_edge([10, 12])
        self.g8.add_edge([11, 12])

        """ g9 --------
                 -------------------------------------
                 |                                   |
                 |                                   |
                (0) --- (1) -- (3) -- (7) -- (6) -- (2)
                 |       |
                (4)     (5)    Nafthalene
                 |       |
                (8) --- (9)


        """
        self.g9 = MolecularGraph(undirected=True)

        # Add Vertex
        for i in range(0, 10):
            self.g9.add_vertex(i)

        # Add get_edges
        self.g9.add_edge([0, 1])
        self.g9.add_edge([0, 2])
        self.g9.add_edge([0, 4])
        self.g9.add_edge([1, 3])
        self.g9.add_edge([1, 5])
        self.g9.add_edge([2, 6])
        self.g9.add_edge([3, 7])
        self.g9.add_edge([4, 8])
        self.g9.add_edge([5, 9])
        self.g9.add_edge([6, 7])
        self.g9.add_edge([8, 9])

        """ g10 --------

                            (0)     (2) 
                               \   /
                                (1)
                                 |
                        (10)    (3)
                       /    \  /   \
               (15)--(11)    (7)   (4)--(5)
              /       |       |      \   /
            (14)     (12)    (8)      (6)
                \   /    \ / 
                (13)     (9)  
                             
         Example: Figure 5
         Paper: A New Algorithm for Exhaustive Ring Perception in a Molecular Graph
                               
        """
        self.g10 = MolecularGraph(undirected=True)

        # Add Vertex
        for i in range(0, 16):
            self.g10.add_vertex(i)

        self.g10.add_edge([0, 1])
        self.g10.add_edge([1, 2])
        self.g10.add_edge([1, 3])
        self.g10.add_edge([3, 4])
        self.g10.add_edge([3, 7])
        self.g10.add_edge([4, 5])
        self.g10.add_edge([4, 6])
        self.g10.add_edge([5, 6])
        self.g10.add_edge([7, 8])
        self.g10.add_edge([7, 10])
        self.g10.add_edge([8, 9])
        self.g10.add_edge([9, 12])
        self.g10.add_edge([10, 11])
        self.g10.add_edge([11, 12])
        self.g10.add_edge([11, 15])
        self.g10.add_edge([12, 13])
        self.g10.add_edge([13, 14])
        self.g10.add_edge([14, 15])

        """ g11 --------
           (0) --> (1) --> (2) --> (5) --> (6)
                   |
                   \/
                  (3)
                   |
                   \/
                  (4)

        """
        self.g11 = MolecularGraph(undirected=False)

        # Add Vertex
        for i in range(0, 7):
            self.g11.add_vertex(i)

        # Add get_edges
        self.g11.add_edge([0, 1])
        self.g11.add_edge([1, 2])
        self.g11.add_edge([1, 3])
        self.g11.add_edge([2, 5])
        self.g11.add_edge([3, 4])
        self.g11.add_edge([5, 6])

        """ g12 --------
           (0) --> (1) --> (2) --> (5) --> (6)
                   |
                   \/
                  (3)
                   |
                   \/
                  (4)

            (7)
            (8)
        """
        self.g12 = MolecularGraph(undirected=False)

        # Add Vertex
        for i in range(0, 9):
            self.g12.add_vertex(i)

        # Add get_edges
        self.g12.add_edge([0, 1])
        self.g12.add_edge([1, 2])
        self.g12.add_edge([1, 3])
        self.g12.add_edge([2, 5])
        self.g12.add_edge([3, 4])
        self.g12.add_edge([5, 6])

    # ##################################################################################################################
    def test_01_emptyGraph(self):

        """Test an empty Molecular Graph

        Returns
        -------
        None

        """
        g = MolecularGraph(undirected=True)
        self.assertTrue(isinstance(g, MolecularGraph))

    # ##################################################################################################################
    def test_02_manualGraph(self):

        """Test the creation of a MolecularGraph by hand

           (0) -- (1) -- (2) -- (5) -- (6)
                   |
                  (3)
                   |
                  (4)

        Returns
        -------
        None

        """
        g = MolecularGraph()

        # Add Vertex
        for i in range(0, 7):
            g.add_vertex(i)

        # Check the vertex
        v = [0, 1, 2, 3, 4, 5, 6]
        vv = g.get_vertices()
        self.assertEqual(v, vv)

        # Add get_edges
        g.add_edge([0, 1])
        g.add_edge([1, 2])
        g.add_edge([1, 3])
        g.add_edge([2, 5])
        g.add_edge([3, 4])
        g.add_edge([5, 6])

        # Check number the get_edges
        lo = [{0, 1}, {1, 2}, {1, 3}, {2, 5}, {3, 4}, {5, 6}]
        ll = g.get_edges()
        self.assertEqual(lo, ll)

        # Draw the graph
        try:
            import pygraphviz as pgv
            g.draw_graph_pygraphviz("graphs/test_02")
        except ModuleNotFoundError:
            try:
                import networkx as nx
                g.draw_graph_networkx(title="graphs/test_02_nx")
            except ModuleNotFoundError:
                print("test_02 cannot create a graph. The pygraphviz and/or networkx modules are not installed.")
                pass

    # ##################################################################################################################
    def test_03_manualGraph(self):

        """Test the creation of a MolecularGraph by hand

           (0) -- (1) -- (2) -- (5) -- (6)
                   |
                  (3)
                   |
                  (4)

        The list of bonds is added in the initilization of the graph

        Returns
        -------
        None

        """
        g = MolecularGraph(nvert=7, listbonds=[[0, 1], [1, 2], [1, 3], [2, 5], [3, 4], [5, 6]])

        # # Add Vertex
        # for i in range(0,7):
        #     g.add_vertex(i)

        # Check the vertex
        v = [0, 1, 2, 3, 4, 5, 6]
        vv = g.get_vertices()
        self.assertEqual(v, vv)

        # # Add get_edges
        # g.add_edge([0,1])
        # g.add_edge([1,2])
        # g.add_edge([1,3])
        # g.add_edge([2,5])
        # g.add_edge([3,4])
        # g.add_edge([5,6])

        # Check number the get_edges
        lo = [{0, 1}, {1, 2}, {1, 3}, {2, 5}, {3, 4}, {5, 6}]
        ll = g.get_edges()
        self.assertEqual(lo, ll)

        # Draw the graph
        try:
            import pygraphviz as pgv
            g.draw_graph_pygraphviz("graphs/test_03")
        except ModuleNotFoundError:
            try:
                import networkx as nx
                g.draw_graph_networkx(title="graphs/test_03_nx")
            except ModuleNotFoundError:
                print("test_03 cannot create a graph. The pygraphviz and/or networkx modules are not installed.")
                pass

    # ##################################################################################################################
    def test_04_molecularGraph_withcycle(self):

        d = self.g3._graphdict

        dd = {0: [1],
              1: [0, 2, 3],
              2: [1, 5, 7],
              3: [1, 4],
              4: [3],
              5: [2, 6, 10],
              6: [5, 12],
              7: [2,  8],
              8: [7,  9],
              9: [8, 12],
              10: [5, 11],
              11: [10],
              12: [9, 6], }

        self.assertEqual(d, dd)

        iscyclic3 = self.g3.iscyclic()
        self.assertEqual(iscyclic3, [True])

        self.g3.perception_rings()
        self.assertTrue(self.g3._cycles == [[5, 2, 7, 8, 9, 12, 6, 5]])

    # # ###############################################################################################################
    @staticmethod
    def test_05_directed_undirected_graph():

        """ g2_undirected --------

                         (7) -- (8) -- (9)
                          |
           (0) -- (1) -- (2) -- (5) -- (6)
                   |             |
                  (3)           (10)
                   |             |
                  (4)           (11)

        """
        g2u = MolecularGraph(undirected=True)

        # Add Vertex
        for i in range(0, 12):
            g2u.add_vertex(i)

        # Add get_edges
        g2u.add_edge([0, 1])
        g2u.add_edge([1, 2])
        g2u.add_edge([1, 3])
        g2u.add_edge([2, 5])
        g2u.add_edge([3, 4])
        g2u.add_edge([5, 6])
        g2u.add_edge([2, 7])
        g2u.add_edge([7, 8])
        g2u.add_edge([8, 9])
        g2u.add_edge([5, 10])
        g2u.add_edge([10, 11])

        print("TEST 05 Undirected: {}".format(g2u._graphdict))

        g2d = MolecularGraph(undirected=False)

        # Add Vertex
        for i in range(0, 12):
            g2d.add_vertex(i)

        # Add get_edges
        g2d.add_edge([0, 1])
        g2d.add_edge([1, 2])
        g2d.add_edge([1, 3])
        g2d.add_edge([2, 5])
        g2d.add_edge([3, 4])
        g2d.add_edge([5, 6])
        g2d.add_edge([2, 7])
        g2d.add_edge([7, 8])
        g2d.add_edge([8, 9])
        g2d.add_edge([5, 10])
        g2d.add_edge([10, 11])

        print("TEST 05   Directed: {}".format(g2d._graphdict))

        try:
            import networkx as nx
            import matplotlib as plt
            g2d.draw_graph_networkx(title="graphs/g2dvsg2u")
        except:
            pass

    # ##################################################################################################################
    def test_07_molecularGraph_findallpaths_cycle(self):

        # Find 1
        path = self.g3.find_all_paths(2, 12)
        path_r = [[2, 5, 6, 12], [2, 7, 8, 9, 12]]
        self.assertEqual(path, path_r)

    # ##################################################################################################################
    def test_08_molecularGraph_findallpaths_length_cycle(self):

        # Find 1
        self.g3.draw_graph_networkx(title="graphs/g3_03_nx")
        path = self.g3.find_all_paths_length(2, 4)
        path_r = [[2, 5, 6, 12, 9], [2, 7, 8, 9, 12]]
        self.assertEqual(path, path_r)

        # Find 2
        path = self.g3.find_all_paths_length(2, 40)
        path_r = []
        self.assertEqual(path, path_r)

        # Find 3
        path = self.g3.find_all_paths_length(2, -1)
        path_r = []
        self.assertEqual(path, path_r)

        # Find 4
        path = self.g3.find_all_paths_length(2, 0)
        path_r = [[2]]
        self.assertEqual(path, path_r)

        # Find 5
        path1 = self.g3.find_all_paths_length(0, 3)
        path2 = self.g3.find_all_paths_length(5, 3)
        l1 = [0, 1, 2, 5]
        t1 = l1 in path1
        l2 = [5, 2, 1, 0]
        t2 = l2 in path2
        self.assertTrue(t1)
        self.assertTrue(t2)

    # #############################################################################################################
    def test_09_molecularGraph_drawgraph(self):

        # Draw the graph
        try:
            import pygraphviz as pgv
            self.g4.draw_graph_pygraphviz(title="graphs/graph04")
            print("TODO: Compare files graph04_test08")

            self.g7.draw_graph_pygraphviz(title="graphs/graph07")
            print("TODO: Compare files graph07_test08.dot")

        except ModuleNotFoundError:
            try:
                import networkx as nx
                self.g4.draw_graph_networkx(title="graphs/graph04_nx")
                self.g7.draw_graph_networkx(title="graphs/graph07_nx")

            except ModuleNotFoundError:
                print("test_09 is not run. The pygraphviz and/or networkx modules are not installed.")
                pass

    # ##################################################################################################################
    def test_10_molecularGraph_isolatedvertices(self):

        n1 = self.g3.find_isolated_vertices()
        n2 = self.g4.find_isolated_vertices()

        self.assertEqual(n1, [])
        self.assertEqual(n2, [12])

    # ##################################################################################################################
    def test_11_molecularGraph_isconnected(self):

        t1 = self.g1.is_connected()
        t3 = self.g1.is_connected()
        t4 = self.g4.is_connected()

        self.assertTrue(t1)
        self.assertTrue(t3)
        self.assertFalse(t4)

        t11 = self.g11.is_connected()
        print(self.g11.get_graph())
        v = set()
        t12 = self.g12.is_connected(vertices_encountered=v)
        vl = {0, 1, 2, 3, 4, 5, 6}
        self.assertEqual(v, vl)

        self.assertTrue(t11)
        self.assertFalse(t12)

        v = set()
        self.g7.is_connected(vertices_encountered=v)
        vl = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}
        self.assertEqual(v, vl)

    # ##################################################################################################################
    def test_12_molecularGraph_components(self):

        l1 = self.g4._nmols
        l2 = self.g4.get_forest()
        ll1 = ll2 = llr = []
        lr = [[0, 1, 2, 3, 4, 5], [6, 7, 8, 9, 11, 10], [12]]
        for item in l1:
            i = sorted(item)
            ll1.append(i)
        for item in l2:
            i = sorted(item)
            ll2.append(i)
        for item in lr:
            i = sorted(item)
            llr.append(i)
        self.assertEqual(ll1, llr)
        self.assertEqual(ll2, llr)

        l1 = self.g5._nmols
        lr = [[0, 1, 2, 3, 4, 5, 13], [6, 7, 8, 9, 14, 11, 10], [12]]
        ll1 = llr = []
        for item in l1:
            i = sorted(item)
            ll1.append(i)
        for item in lr:
            i = sorted(item)
            llr.append(i)

        self.assertEqual(ll1, llr)

        l1 = self.g7._nmols
        print(l1)

    # ##################################################################################################################
    def test_13_molecularGraph_autograph(self):

        try:
            import pygraphviz as pgv
            self.g6.draw_graph_pygraphviz(title="graphs/graph06_test12")
            print("TODO: Compare files graph06_test12.dot")

        except ModuleNotFoundError:
            try:
                import networkx as nx
                self.g6.draw_graph_networkx(title="graphs/graph06_test12_nx")
            except ModuleNotFoundError:
                print("test_13 is not run. The pygraphviz module is not installed.")
                pass

    # ##################################################################################################################
    def test_14_molecularGraph_bfs_dfs(self):

        """ g6 --------

           (0) -- (1) -- (2) -- (3) -- (4) -- (5)

                  (6) -- (7) -- (8) -- (9)
                          |      |
                         (10)   (11)
            (12)
        """

        p6_bfs = self.g6.bfs_iterative(start=6)
        p6_dfs = self.g6.dfs_iterative(start=6)
        self.assertEqual(p6_bfs, [6, 7, 8, 10,  9, 11])
        self.assertEqual(p6_dfs, [6, 7, 8,  9, 11, 10])

        """ g3 --------

                         (7) -- (8) -- (9)--\
                          |                  (12)
           (0) -- (1) -- (2) -- (5) -- (6)--/
                   |             |
                  (3)           (10)
                   |             |
                  (4)           (11)

        """
        p3_bfs = self.g3.bfs_iterative(start=0)
        p3_dfs = self.g3.dfs_iterative(start=0)
        self.assertEqual(p3_bfs, [0, 1, 2, 3, 5, 7, 4, 6, 10, 8, 12, 11, 9])
        self.assertEqual(p3_dfs, [0, 1, 2, 5, 6, 12, 9, 8, 7, 10, 11, 3, 4])

        """ g11 --------
           (0) --> (1) --> (2) --> (5) --> (6)
                   |
                   \/
                  (3)
                   |
                   \/
                  (4)

        """

        p11_bfs = self.g11.bfs_iterative(start=0)
        p11_dfs = self.g11.dfs_iterative()
        self.assertEqual(p11_bfs, [0, 1, 2, 3, 5, 4, 6])
        self.assertEqual(p11_dfs, [0, 1, 2, 5, 6, 3, 4])

    # ##################################################################################################################
    def test_15_molecularGraph_drawgraphforest(self):

        try:
            import pygraphviz as pgv
            self.g4.draw_graph_forest_pygraphviz(title="graphs/graph04_f")
            print("TODO: Compare files graph04_f_test14")
            self.g7.draw_graph_forest_pygraphviz(title="graphs/graph07_f")
            print("TODO: Compare files graph07_f_test14.dot")

        except ModuleNotFoundError:
            try:
                import networkx as nx
                self.g4.draw_graph_forest_networkx(title="graphs/graph04_f")
                self.g7.draw_graph_forest_networkx(title="graphs/graph07_f")
            except ModuleNotFoundError:
                print("test_15 is not run. The pygraphviz module is not installed.")
                pass

    # ##################################################################################################################
    def test_16_molecularGraph_iscyclic(self):

        iscyclic1 = self.g4.iscyclic()
        iscyclic2 = self.g7.iscyclic()

        self.assertEqual(iscyclic1, [False, False, False])
        self.assertEqual(iscyclic2, [True, False])

    # ##################################################################################################################
    def test_17_molecularGraph_allbonds(self):

        try:
            import networkx as nx
            import matplotlib as plt
        except:
            pass

        self.g7.draw_graph_networkx(title="graphs/graph07_nx")

        bl1 = self.g4.get_allbonds()
        self.assertEqual(bl1, [[0, 1], [1, 2], [2, 3], [3, 4],
                               [4, 5], [6, 7], [7, 8], [7, 10],
                               [8, 9], [8, 11]])

        bl1 = self.g7.get_allbonds()
        self.assertEqual(bl1, [[0, 1], [1, 2], [1, 3], [2, 5], [2, 7], [3, 4],
                               [5, 6], [5, 10], [6, 12], [7, 8], [8, 9], [9, 12],
                               [10, 11], [13, 14], [14, 15], [15, 17], [15, 16]])

    # ##################################################################################################################
    def test_18_molecularGraph_allbends(self):

        al1 = self.g4.get_allbends()
        self.assertEqual(al1, [[0, 1, 2], [1, 2, 3], [2, 3, 4], [3, 4, 5],
                               [6, 7, 8], [6, 7, 10], [7, 8, 9], [7, 8, 11],
                               [8, 7, 10], [9, 8, 11]])

        al1 = self.g7.get_allbends()
        self.assertEqual(al1, [[0, 1, 2], [0, 1, 3], [1, 2, 5], [1, 2, 7],
                               [1, 3, 4], [2, 1, 3], [2, 5, 6], [2, 5, 10],
                               [2, 7, 8], [5, 2, 7], [5, 6, 12], [5, 10, 11],
                               [6, 5, 10], [6, 12, 9], [7, 8, 9], [8, 9, 12],
                               [13, 14, 15], [14, 15, 17], [14, 15, 16], [16, 15, 17]])

    # ##################################################################################################################
    def test_19_molecularGraph_alldihedral(self):

        d1 = self.g4.get_alldihedrals()
        self.assertEqual(d1, [[3, 2, 1, 0], [4, 3, 2, 1], [5, 4, 3, 2],
                              [9, 8, 7, 6], [11, 8, 7, 6], [10, 7, 8, 9],
                              [11, 8, 7, 10]])
        d1 = self.g7.get_alldihedrals()
        self.assertEqual(d1, [[5, 2, 1, 0], [7, 2, 1, 0], [4, 3, 1, 0], [6, 5, 2, 1],
                              [10, 5, 2, 1], [8, 7, 2, 1], [4, 3, 1, 2], [12, 6, 5, 2],
                              [11, 10, 5, 2], [9, 8, 7, 2], [5, 2, 1, 3], [7, 2, 1, 3],
                              [8, 7, 2, 5], [9, 12, 6, 5], [7, 2, 5, 6], [11, 10, 5, 6],
                              [8, 9, 12, 6], [10, 5, 2, 7], [12, 9, 8, 7], [12, 6, 5, 10],
                              [17, 15, 14, 13], [16, 15, 14, 13]])

    # ##################################################################################################################
    def test_20_molecularGraph_allimproper(self):

        isbackbone_dict = {}
        for i in range(0, self.g4.natoms):
            isbackbone_dict[i] = True
        isbackbone_dict[10] = False
        isbackbone_dict[11] = False
        d1 = self.g4.get_allimpropers_ua_polymers(isbackbone_dict=isbackbone_dict)
        self.assertEqual(d1, [[7, 6, 8, 10], [8, 7, 9, 11]])

        isbackbone_dict = {}
        for i in range(0, self.g7.natoms):
            isbackbone_dict[i] = True
        isbackbone_dict[2] = False
        isbackbone_dict[5] = False
        isbackbone_dict[6] = False
        isbackbone_dict[7] = False
        isbackbone_dict[8] = False
        isbackbone_dict[9] = False
        isbackbone_dict[10] = False
        isbackbone_dict[11] = False
        isbackbone_dict[12] = False
        isbackbone_dict[16] = False
        d1 = self.g7.get_allimpropers_ua_polymers(isbackbone_dict=isbackbone_dict)
        self.assertEqual(d1, [[1, 0, 3, 2], [2, 5, 7, 1], [15, 14, 17, 16]])

    # ##################################################################################################################
    def test_21_molecularGraph_neigh(self):

        ln = self.g3.get_neighbours(2)
        lo = len(ln)
        self.assertEqual(lo, 3)
        self.assertEqual(ln, [1, 5, 7])

        # Node without neighbours
        ln = self.g6.get_neighbours(12)
        lo = len(ln)
        self.assertEqual(lo, 0)
        self.assertEqual(ln, [])

    # ##################################################################################################################
    def test_22_MolecularGraph_enumerate_fusedcycles(self):

        self.g9.perception_rings()
        self.assertTrue(self.g9._cycles == [[1, 0, 2, 6, 7, 3, 1],
                                            [1, 0, 4, 8, 9, 5, 1],
                                            [1, 3, 7, 6, 2, 0, 4, 8, 9, 5, 1]])

        self.g5.perception_rings()
        self.assertTrue(self.g5._cycles == [])

        self.g7.perception_rings()
        self.assertTrue(self.g7._cycles == [[5, 2, 7, 8, 9, 12, 6, 5]])

        self.g10.perception_rings()
        self.assertTrue(self.g10._cycles == [[4, 6, 5, 4],
                                             [12, 11, 15, 14, 13, 12],
                                             [12, 11, 10, 7, 8, 9, 12],
                                             [12, 13, 14, 15, 11, 10, 7, 8, 9, 12]])

    # ##################################################################################################################
    def test_23_MolecularGraph_addedgesbetweennotexistingvertex(self):

        """ g4 --------



           (0) -- (1) -- (2) -- (3) -- (4) -- (5)

                  (6) -- (7) -- (8) -- (9)
                          |      |
                         (10)   (11)
            (12)
        """

        e = self.g4._generate_edges()
        print(e)
        self.g4.find_isolated_vertices()

        self.g4.add_vertex(13)
        self.g4.add_edge([12, 13])
        self.g4.remove_edge([2, 3])

        try:
            import pygraphviz as pgv
            self.g4.draw_graph_pygraphviz(title="graphs/graph04_after_test23")
            print("TODO: Compare files graph06_test12.dot")

        except ModuleNotFoundError:
            try:
                import networkx as nx
                self.g4.draw_graph_networkx(title="graphs/graph04_after_test23_nx")
                self.g4.draw_graph_forest_networkx(title="graphs/graph04_after_test23_nx")
            except ModuleNotFoundError:
                print("test_23 is not run. The pygraphviz and/or networkx module is not installed.")
                pass

        try:
            import networkx as nx
            self.g4.draw_graph_forest_networkx(title="graphs/graph04_after_test23_nx")
        except ModuleNotFoundError:
            print("test_23 is not run. The pygraphviz and/or networkx module is not installed.")
            pass

    # ###############################################################################################################
    def test_24_MolecularGraph_directed(self):

        """ g11 --------
           (0) --> (1) --> (2) --> (5) --> (6)
                   |
                   \/
                  (3)
                   |
                   \/
                  (4)

        """

        v = self.g11.get_vertices()
        vl = [0, 1, 2, 3, 4, 5, 6]
        self.assertEqual(v, vl)

        try:
            import networkx as nx
            import matplotlib as plt
            self.g11.draw_graph_networkx(title="graphs/graph11_nx")
        except:
            pass

        e = self.g11.get_edges()
        el = [{0, 1}, {1, 2}, {1, 3}, {2, 5}, {3, 4}, {5, 6}]
        self.assertEqual(e, el)

        p1 = self.g11.find_all_paths(0, 3)
        p1l = [[0, 1, 3]]
        p2 = self.g11.find_all_paths(5, 0)
        p2l = []
        p3 = self.g11.find_all_paths_length(0, 3)
        p3l = [[0, 1, 2, 5], [0, 1, 3, 4]]
        p4 = self.g11.find_all_paths_length(2, 1)
        p4l = [[2, 5]]

        self.assertEqual(p1, p1l)
        self.assertEqual(p2, p2l)
        self.assertEqual(p3, p3l)
        self.assertEqual(p4, p4l)

        v1 = self.g11.find_isolated_vertices()
        v1l = []
        self.assertEqual(v1, v1l)

        v1 = self.g12.find_isolated_vertices()
        v1l = [7, 8]
        self.assertEqual(v1, v1l)

    # ##################################################################################################################
    def test_25_manualGraph(self):

        """Test the creation of a MolecularGraph by hand

           (0) -- (1) -- (2) -- (5) -- (6)
                   |
                  (3)
                   |
                  (4)

        The list of bonds is added in the initilization of the graph

        And after, three new vertex and edges are added.

        (7), (8), (9)

        Returns
        -------
        None

        """
        g = MolecularGraph(nvert=7, listbonds=[[0, 1], [1, 2], [1, 3], [2, 5], [3, 4], [5, 6]])

        # # Add Vertex
        # for i in range(0,7):
        #     g.add_vertex(i)

        # Check the vertex
        v = [0, 1, 2, 3, 4, 5, 6]
        vv = g.get_vertices()
        self.assertEqual(v, vv)

        # # Add get_edges
        # g.add_edge([0,1])
        # g.add_edge([1,2])
        # g.add_edge([1,3])
        # g.add_edge([2,5])
        # g.add_edge([3,4])
        # g.add_edge([5,6])

        # Check number the get_edges
        lo = [{0, 1}, {1, 2}, {1, 3}, {2, 5}, {3, 4}, {5, 6}]
        ll = g.get_edges()
        self.assertEqual(lo, ll)

        # Draw the graph
        try:
            import pygraphviz as pgv
            g.draw_graph_pygraphviz("graphs/test_25_before")
        except ModuleNotFoundError:
            try:
                import networkx as nx
                g.draw_graph_networkx(title="graphs/test_25_before_nx")
            except ModuleNotFoundError:
                print("test_03 cannot create a graph. The pygraphviz and/or networkx modules are not installed.")
                pass

        # Add vertex
        g.add_vertex(7)
        g.add_vertex(8)
        g.add_vertex(9)

        # Draw the graph
        try:
            import pygraphviz as pgv
            g.draw_graph_pygraphviz("graphs/test_25_2_before")
        except ModuleNotFoundError:
            try:
                import networkx as nx
                g.draw_graph_networkx(title="graphs/test_25_2_before_nx")
            except ModuleNotFoundError:
                print("test_03 cannot create a graph. The pygraphviz and/or networkx modules are not installed.")
                pass

        # Add edges
        g.add_edge((1, 7), setforest=False)
        g.add_edge((2, 8), setforest=False)
        g.add_edge((9, 8), setforest=True)

        # Draw the graph
        try:
            import pygraphviz as pgv
            g.draw_graph_pygraphviz("graphs/test_25_after")
        except ModuleNotFoundError:
            try:
                import networkx as nx
                g.draw_graph_networkx(title="graphs/test_25_after_nx")
            except ModuleNotFoundError:
                print("test_03 cannot create a graph. The pygraphviz and/or networkx modules are not installed.")
                pass

        print(g._nmols)

    # ##################################################################################################################
    def test_26_size_of_graphs(self):

        print("Graph  1: {} bytes".format(self.g1.__sizeof__()))
        print("Graph  2: {} bytes".format(self.g2.__sizeof__()))
        print("Graph  3: {} bytes".format(self.g3.__sizeof__()))
        print("Graph  7: {} bytes".format(self.g7.__sizeof__()))
        print("Graph 10: {} bytes".format(self.g10.__sizeof__()))
        print("Graph 12: {} bytes".format(self.g12.__sizeof__()))

    # ##################################################################################################################
    def test_27_iatch_graph(self):

        lo = np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2])
        res = np.array_equal(self.g4._iatch, lo)
        self.assertTrue(res)

        lo = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
        res = np.array_equal(self.g7._iatch, lo)
        self.assertTrue(res)


if __name__ == '__main__':

    unittest.main()
