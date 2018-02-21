#!/usr/bin/env python3

import math
import copy
import networkx as nx

from python import intervals

#COPIA
#[LTR, PBS, GAG, AP, INT, RT, RNaseH, PPT, LTR]
#GYPSY
#[LTR, PBS, GAG, AP, RT, RNaseH, INT, PPT, LTR]

class TransposonGraph(object):
    """Transposon graph constructed from one ltr pair and list of domains of one gene.
    """
    def __init__(self, te, domainList):
        self._graph = None
        self._buildGraph(te, domainList)

    def _buildGraph(self, te, domainList): #build graphs
        #TODO
        #check negative strands
        self._graph = nx.DiGraph()
        self._addNodes(te, domainList)
        self._addEdges()

    def _addNodes(self, te, domainsList): #add all necessary nodes
        #LTR nodes
        self._graph.add_node(
            n='ltr_left',
            location=te.ltrLeftLocation,
            score=1,
            node_class='ltr_left',
            strand='0',
            features={}
        )

        self._graph.add_node(
            n='ltr_right',
            location=te.ltrRightLocation,
            score=1,
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
            if (intervals.compare(te.ltrLeftLocation, location) != 1
              and intervals.compare(location, te.ltrRightLocation) != 1):
                self._graph.add_node(
                    n='ppt',
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
            if (intervals.compare(te.ltrLeftLocation, location) != 1
              and intervals.compare(location, te.ltrRightLocation) != 1):
                self._graph.add_node(
                    n='pbs',
                    location=location,
                    score=1,
                    node_class='pbs',
                    strand=strand,
                    features={}
                )

        #DOMAINS
        #TODO negative strands
        i = 0
        for domain in domainsList:
            if domain.type not in ['GAG', 'AP', 'INT', 'RT', 'RNaseH']:
                continue
            if (intervals.compare(te.ltrLeftLocation, domain.location) != 1
              and intervals.compare(domain.location, te.ltrRightLocation) != 1):
                self._graph.add_node(
                    n='domain_{}'.format(i),
                    location=domain.location,
                    score=domain.score,
                    node_class=domain.type,
                    strand='+',
                    features={'domain': domain}
                )

                i += 1

    def _addEdges(self): #add all necessary edges
        for node1 in self._graph.nodes(data=True):
            for node2 in self._graph.nodes(data=True):
                self._addEdge(node1, node2)

    def _addEdge(self, node1, node2): #add edge between nodes
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

            gypsy_diff = gypsy_index1 - gypsy_index2
            if gypsy_diff != copia_diff:                
                if gypsy_diff == 1:
                    self._graph.add_edge(node1[0], node2[0], weight=float(1)/node2[1]['score'])
                elif gypsy_diff > 1:
                    self._graph.add_edge(node1[0], node2[0], weight=1000 * gypsy_diff)
                else:
                    self._graph.add_edge(node1[0], node2[0], weight=1000 * 2 * (-gypsy_diff))

    """Find best evaluated path and return its score

    Returns:
        float: score
    """
    def getScore(self): 
        pathScore = 0
        pathFeatures = {
            'domains': [],
            'pbs': [float('nan'), float('nan')],
            'ppt': [float('nan'), float('nan')]
        }

        path = nx.dijkstra_path(self._graph, 'ltr_left', 'ltr_right', 'weight')[1:-1]
        scores = nx.get_node_attributes(self._graph, 'score')
        locations = nx.get_node_attributes(self._graph, 'location')
        features = nx.get_node_attributes(self._graph, 'features')

        for node in path:
            pathScore += scores[node]
            if node == 'pbs' or node == 'ppt':
                pathFeatures[node] = locations[node]
            else:
                pathFeatures['domains'].append(features[node]['domain'])

        return pathScore, pathFeatures
