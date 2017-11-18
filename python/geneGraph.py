import math, copy
import itertools
import networkx as nx
import numpy as np
import pprint

import sys

#COPIA
#[LTR, PBS, GAG, AP, INT, RT, RNaseH, PPT, LTR]

#GYPSY
#[LTR, PBS, GAG, AP, RT, RNaseH, INT, PPT, LTR]

class GeneGraph:
    def __init__(self, ltr_pair, domains):
        self.__buildGraph(ltr_pair, domains)

    #compare intervals a, b
    def __compareIntervals(self, a, b):
        return -1 if a[1] < b[0] else (0 if a[1] == b[0] else 1)

    def __addNodes(self, ltr_pair, domains):
        ltr_left_location = [ltr_pair['location'][0], 
                                ltr_pair['location'][0] 
                                    + int(ltr_pair['ltr len'].split(',')[0])]

        ltr_right_location = [ltr_pair['location'][1], 
                                ltr_pair['location'][1] 
                                    + int(ltr_pair['ltr len'].split(',')[1])] 

        self._graph.add_node('ltr_left',
                location=ltr_left_location,
                score=1,
                node_class='ltr_left',
                features={})

        self._graph.add_node('ltr_right',
                location=ltr_right_location,
                score=1,
                node_class='ltr_right',
                features={})

        

        #pbs/ppt nodes
        if not math.isnan(ltr_pair['ppt'][0]):
            location = copy.deepcopy(ltr_pair['ppt'])
            #TODO negative strand, check positions
            if location[0] > location[1]: location = [location[1], location[0]] 
            if (self.__compareIntervals(ltr_left_location, location) != 1 
                    and self.__compareIntervals(location, ltr_right_location) != 1):
                #self._graph.add_node(node_ppt['id'], attr_dict=node_ppt)
                self._graph.add_node('ppt',
                        location=location,
                        score=1,
                        node_class='ppt',
                        features={})

        if not math.isnan(ltr_pair['pbs'][0]):
            location = copy.deepcopy(ltr_pair['pbs'])
            #TODO negative strand, check positions
            if location[0] > location[1]: location = [location[1], location[0]] 
            if (self.__compareIntervals(ltr_left_location, location) != 1 
                    and self.__compareIntervals(location, ltr_right_location) != 1):
                #self._graph.add_node(node_pbs['id'], attr_dict=node_pbs) 
                self._graph.add_node('pbs',
                        location=location,
                        score=1,
                        node_class='pbs',
                        features={})


        #add domain nodes
        #TODO domain orientation
        for domain_type in domains:
            if domain_type == 'XX':
                continue
            i = 0
            for domain in domains[domain_type]:                
                if (self.__compareIntervals(ltr_left_location, domain['location']) != 1 
                        and self.__compareIntervals(domain['location'], ltr_right_location) != 1):
                    self._graph.add_node('{}_{}'.format(domain_type, i), 
                            location=domain['location'],
                            score=domain['score'],
                            node_class=domain_type,
                            features={'domain':domain})
                    i += 1   
                      

    def _addEdge(self, node1, node2):
        if node1 == node2:
            return
        if node1[1]['node_class'] == node2[1]['node_class']:
            return

        copia = ['ltr_left', 'pbs', 'GAG', 'AP', 'INT', 'RT', 'RNaseH', 'ppt', 'ltr_right']

        #node1 is before node2
        if self.__compareIntervals(node1[1]['location'], node2[1]['location']) != 1:
            node1_index = copia.index(node1[1]['node_class'])
            node2_index = copia.index(node2[1]['node_class'])
            diff = node2_index - node1_index
            if diff == 1:
                self._graph.add_edge(node1[0], node2[0], weight=float(1)/node2[1]['score'])
            elif diff > 1:
                self._graph.add_edge(node1[0], node2[0], weight=1000 * diff)
            else:
                self._graph.add_edge(node1[0], node2[0], weight=1000 * 2 * (-diff))

    def __addEdges(self):
        for node1 in self._graph.nodes(data=True):
            for node2 in self._graph.nodes(data=True):
                self._addEdge(node1, node2)
                

    def getScore(self):
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

        for p in path:
            path_score += scores[p]
            if p == 'pbs':
                path_features['pbs'] = locations['pbs']
            elif p == 'ppt':
                path_features['ppt'] = locations['ppt']
            else:
                path_features['domains'].append(features[p]['domain'])

        return path_score, path_features


    def __buildGraph(self, ltr_pair, domains):
        #TODO
        #CHECK NEGATIVE STRAND
        self._graph = nx.DiGraph()
        
        self.__addNodes(ltr_pair, domains)
        self.__addEdges()
