#!/usr/bin/env python3

import math
import copy
import networkx as nx

from nested.utils import intervals

#COPIA
#[LTR, PBS, GAG, AP, INT, RT, RNaseH, PPT, LTR]
#GYPSY
#[LTR, PBS, GAG, AP, RT, RNaseH, INT, PPT, LTR]

class Graph(object):
    """Transposon path graph used to evaluate domain position and penalize incorrect order.
    It is constructed from one ltr pair and a list of domains found in gene.
    """
    def __init__(self, te, domain_list):
        self._graph = None
        self._build_graph(te, domain_list)

    def _build_graph(self, te, domain_list): #build graphs
        #TODO
        #check negative strands
        self._graph = nx.DiGraph()
        self._add_nodes(te, domain_list)
        self._add_edges()

    def _add_nodes(self, te, domain_list): #add all necessary nodes
        #LTR nodes
        self._graph.add_node(
            'ltr_left',
            location=te.ltr_left_location,
            score=1, #TODO introduce tsr into scoring
            node_class='ltr_left',
            strand='0',
            features={}
        )

        self._graph.add_node(
            'ltr_right',
            location=te.ltr_right_location,
            score=1, #TODO introduce tsr into scoring
            node_class='ltr_right',
            strand='0',
            features={}
        )

        #PBS/PPT nodes
        if not math.isnan(te.ppt[0]):
            location = copy.deepcopy(te.ppt)
            strand = '0'
            if location[0] > location[1]: 
                location = [location[1], location[0]]
                strand = '-'
            if (intervals.compare(te.ltr_left_location, location) != 1
              and intervals.compare(location, te.ltr_right_location) != 1):
                self._graph.add_node(
                    'ppt',
                    location=location,
                    score=1,
                    node_class='ppt',
                    strand=strand,
                    features={}
                )

        if not math.isnan(te.pbs[0]):
            location = copy.deepcopy(te.pbs)
            strand = '0'
            if location[0] > location[1]: 
                location = [location[1], location[0]]
                strand = '-'
            if (intervals.compare(te.ltr_left_location, location) != 1
              and intervals.compare(location, te.ltr_right_location) != 1):
                self._graph.add_node(
                    'pbs',
                    location=location,
                    score=1,
                    node_class='pbs',
                    strand=strand,
                    features={}
                )

        #DOMAINS
        #TODO negative strands
        i = 0
        for domain in domain_list:
            if domain.type not in ['GAG', 'AP', 'INT', 'RT', 'RNaseH']:
                continue
            if (intervals.compare(te.ltr_left_location, domain.location) != 1
              and intervals.compare(domain.location, te.ltr_right_location) != 1):
                self._graph.add_node(
                    'domain_{}'.format(i),
                    location=domain.location,
                    score=domain.score,
                    node_class=domain.type,
                    strand='+',
                    features={'domain': domain}
                )

                i += 1

    def _add_edges(self): #add all necessary edges
        for node1 in self._graph.nodes(data=True):
            for node2 in self._graph.nodes(data=True):
                self._add_edge(node1, node2)

    def _add_edge(self, node1, node2): #add edge between nodes
        if node1 == node2 or node1[1]['node_class'] == node2[1]['node_class']:
            return
        
        copia = ['ltr_left', 'pbs', 'GAG', 'AP', 'INT', 'RT', 'RNaseH', 'ppt', 'ltr_right']
        gypsy = ['ltr_left', 'pbs', 'GAG', 'AP', 'RT', 'RNaseH', 'INT', 'ppt', 'ltr_right']
        #node1 before node2
        if intervals.compare(node1[1]['location'], node2[1]['location']) != 1:
            copia_index1 = copia.index(node1[1]['node_class'])
            copia_index2 = copia.index(node2[1]['node_class'])
            gypsy_index1 = gypsy.index(node1[1]['node_class'])
            gypsy_index2 = gypsy.index(node2[1]['node_class'])

            copia_diff = copia_index2 - copia_index1
            if copia_diff == 1:
                self._graph.add_edge(node1[0], node2[0], weight=float(1)/node2[1]['score'])
            elif copia_diff > 1:
                self._graph.add_edge(node1[0], node2[0], weight=1000 * copia_diff)
            else:
                self._graph.add_edge(node1[0], node2[0], weight=1000 * 2 * (-copia_diff))

            """
            gypsy_diff = gypsy_index1 - gypsy_index2
            if gypsy_diff != copia_diff:                
                if gypsy_diff == 1:
                    self._graph.add_edge(node1[0], node2[0], weight=float(1)/node2[1]['score'])
                elif gypsy_diff > 1:
                    self._graph.add_edge(node1[0], node2[0], weight=1000 * gypsy_diff)
                else:
                    self._graph.add_edge(node1[0], node2[0], weight=1000 * 2 * (-gypsy_diff))
            """

    def get_score(self): 
        path_score = 0
        path_features = {
            'domains': [],
            'pbs': [float('nan'), float('nan')],
            'ppt': [float('nan'), float('nan')]
        }

        path = nx.dijkstra_path(self._graph, 'ltr_left', 'ltr_right', 'weight')[1:-1]
        scores = nx.get_node_attributes(self._graph, 'score')
        locations = nx.get_node_attributes(self._graph, 'location')
        features = nx.get_node_attributes(self._graph, 'features')

        for node in path:
            path_score += scores[node]
            if node == 'pbs' or node == 'ppt':
                path_features[node] = locations[node]
            else:
                path_features['domains'].append(features[node]['domain'])

        return path_score, path_features

