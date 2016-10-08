#!/usr/bin/python3

import networkx as nx
import matplotlib.pyplot as plt

input_file = '/home/genome/joseph7e/SSRs/pairwise_stuff/ssr_matches_table2'

my_list = []
my_edges = []
my_nodes = []


with open(input_file, 'r') as i:
    for line in i.readlines():
        line = line.rstrip()
        count, group = line.split("\t")
        list = group.split("_")
        tl = set(list)
        my_list.append(tl)

        # dict( size=int(d[1]) )
        #             ]) for d in data]
        #
        tup = (group.replace("_", "+"), dict( size=int(count) ))
        my_nodes.append(tup)


for i in range(len(my_list)):
    for j in range(i + 1, len(my_list)):
        if my_list[i].issubset(my_list[j]) or my_list[j].issubset(my_list[i]):
            my_edges.append((my_list[i],my_list[j]))

my_edges = [('+'.join(sorted(e[0])), '+'.join(sorted(e[1]))) for e in my_edges]


if __name__ == '__main__':
    scale_factor = 1
    G = nx.Graph()
    nodes = my_nodes
    node_sizes = [ (n[1])['size']*scale_factor
                   for n in nodes ]
    edges = my_edges
    G.add_edges_from( edges )

    nx.draw_networkx(G,
                     pos=nx.spring_layout(G),
                     node_size = node_sizes,
                     font_size = 6,
                     label = None)
    plt.axis('off')
    plt.savefig("my_partitioned_venn_un.jpg")
    plt.show()