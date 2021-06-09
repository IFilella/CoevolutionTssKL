import coevolution
import numpy as np

#print("TssGF---------------------------------------------------------------")
liGF = coevolution.CoevList('data/coevolution/TssGF_pairs_raptorGremlin.txt',366)
liGF.add_distances(pdb='data/structures/6n38.pdb',chainsA=['G'],chainsB=['H','I'],resomitA=list(range(1,122))+list(range(331,337)),resomitB=list(range(1,47))+list(range(393,407)))
interliGF = liGF.get_interpairs()
for pair in interliGF:
    #print(pair)
    pass

#print("TssKG---------------------------------------------------------------")
liKG = coevolution.CoevList('data/coevolution/TssKG_pairs_raptorGremlin.txt',445)
liKG.add_distances(pdb='data/structures/6n38.pdb',chainsA=['A','B','C','D','E','F'],chainsB=['G'],resomitA=list(range(1,2))+list(range(220,233))+list(range(312,446)),resomitB=list(range(1,122))+list(range(331,336)))
interliKG = liKG.get_interpairs()
for pair in interliKG:
    #print(pair)
    pass

distThreshold = 8

liGFKG = np.asarray(list(liGF.pairlist)+list(liKG.pairlist))
coevolution.plot_RocPlot(liGFKG,dist_T=distThreshold)

print("TssKctdL---------------------------------------------------------------")

liKctdL = coevolution.CoevList('data/coevolution/TssKctdL_pairs_raptorGremlin.txt',445,314)
interliKctdL = liKctdL.get_interpairs()
for pair in interliKctdL[0:10]:
    print(pair)
    pass
