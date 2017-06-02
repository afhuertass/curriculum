

import networkx as nx

import pylab as p

G = nx.DiGraph()

G.add_path( [ 0 ,  1  , 2 ,3 ] )
G.add_edge ( 2 ,  0)
G.add_edge( 1 , 4 )

pred =  nx.dfs_predecessors( G ,  3  ) 
print  ( pred )
#print ( pred[2] )

suss =  nx.dfs_successors( G , 3  ) 

print suss


for keys, value in pred.iteritems():
    #key tiene por predecesor al value
    print value
    
for keys, value in suss.iteritems():
    # key tiene por sucesor al value 
    print value



nx.draw(G)
p.show()

neww = nx.reverse(G)

pred = nx.dfs_successors(neww, 2  )
print( pred )
for keys, value in pred.iteritems():
    #key tiene por predecesor al value
    print value

