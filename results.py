import coevolution
import numpy as np
import matplotlib.pyplot as plt

#print("TssGF---------------------------------------------------------------")
liGF = coevolution.CoevList('data/coevolution/TssGF_pairs_raptorGremlin.txt',366)
liGF.add_distances(pdb='data/structures/6n38.pdb',chainsA=['G'],chainsB=['H','I'],resomitA=list(range(1,122))+list(range(331,337)),resomitB=list(range(1,47))+list(range(393,407)))
interliGF = liGF.get_interpairs()
intraliGF = liGF.get_intrapairs()

#print("TssKG---------------------------------------------------------------")
liKG = coevolution.CoevList('data/coevolution/TssKG_pairs_raptorGremlin.txt',445)
liKG.add_distances(pdb='data/structures/6n38.pdb',chainsA=['A','B','C','D','E','F'],chainsB=['G'],resomitA=list(range(1,2))+list(range(220,233))+list(range(312,446)),resomitB=list(range(1,122))+list(range(331,336)))
interliKG = liKG.get_interpairs()
intraliKG = liKG.get_intrapairs()

distThreshold = 8

liinterGFKG = np.asarray(list(interliGF)+list(interliKG))
liintraGFKG = np.asarray(list(interliGF)+list(intraliKG))
liGFKG = np.asarray(list(liGF.pairlist)+list(liKG.pairlist))
coevolution.plot_contact_isotonicReg(liinterGFKG,distThreshold)
#coevolution.plot_contact_isotonicReg(liintraGFKG,distThreshold)
coevolution.plot_RocPlot(liinterGFKG,dist_T=distThreshold,intervals=50,reg=True)
#coevolution.plot_RocPlot(liintraGFKG,dist_T=distThreshold,intervals=50,reg=True)
#counts = coevolution.get_BinaryClass(liinterGFKG,dist_T=distThreshold,prob_T=0.88)
#print(counts)
#coevolution.plot_RocPlot(liGFKG,dist_T=distThreshold)
plt.show()
exit()

print("TssKctdL---------------------------------------------------------------")

liKctdL = coevolution.CoevList('data/coevolution/TssKctdL_pairs_raptorGremlin.txt',445,314)
interliKctdL = liKctdL.get_interpairs()
for pair in interliKctdL[0:10]:
    print(pair)
    pass
