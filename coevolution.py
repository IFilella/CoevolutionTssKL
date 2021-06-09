import numpy as np
from Bio.PDB import *
from Bio import *
import matplotlib.pyplot as plt

class CoevList(object):
    """
    Raptor (concatenated MSA) + Gremlin (coevolving pairs) list of coevolving pairs
    """
    def __init__(self,pairlist,Alen,Astart=1,Bstart=1):
        """
        pairlist -> .txt outcome from gremlin server with all pairs
        Alen -> integer lenght of the first gene
        Adisp -> integer starting residue of the first gene
        Bdisp -> integer startung residue of the second gene
        """
        self.Alen = Alen
        self.Astart = Astart
        self.Bstart = Bstart
        self.ABbreak = self.Alen-self.Astart+1
        self.pairlist = np.genfromtxt(pairlist,dtype=str,delimiter='\t',skip_header=1)
        self.parse_list()
     
    def parse_list(self):
        newpairlist=[]
        for pair in self.pairlist:
            newpair=[]
            #print(pair)
            if int(pair[0]) <= self.ABbreak:
                newpair.append(int(pair[0])+self.Astart-1)
                newpair.append('A')
                newpair.append(pair[2].split("_")[1])
            else:
                newpair.append(int(pair[0])-self.ABbreak+self.Bstart-1)
                newpair.append('B')
                newpair.append(pair[2].split("_")[1])
            if int(pair[1]) <= self.ABbreak:
                newpair.append(int(pair[1])+self.Astart-1)
                newpair.append('A')
                newpair.append(pair[3].split("_")[1])
            else:
                newpair.append(int(pair[1])-self.ABbreak+self.Bstart-1)
                newpair.append('B')
                newpair.append(pair[3].split("_")[1])
            #print(newpair)
            newpair.append(pair[-1])
            newpairlist.append(newpair)
        self.pairlist = np.asarray(newpairlist)

    def get_intrapairs(self):
        newpairlist=[]
        for pair in self.pairlist:
            if pair[1]==pair[4]: newpairlist.append(pair)
        return(np.asarray(newpairlist))
    
    def get_interpairs(self):
        newpairlist=[]
        for pair in self.pairlist:
            if pair[1]!=pair[4]: newpairlist.append(pair)
        return(np.asarray(newpairlist))

    def add_distances(self,pdb,chainsA,chainsB,resomitA,resomitB):
        newpairlist=[]
        parser = PDBParser()
        structure = parser.get_structure('structure',pdb)
        for i,pair in enumerate(self.pairlist):
            #print(self.pairlist[i])
            if pair[1]=='A' and pair[4]=='A':
                if pair[0]==pair[3]:
                    np.concatenate(self.pairlist[i],np.asarray(['-','-',float(0)]))
                else:
                    minDist = self.get_minmResDist_chains(structure,pair[0],pair[3],chainsA,chainsA,resomitA,resomitA)
            elif pair[1]=='B' and pair[4]=='B':
                if pair[0]==pair[3]:
                    self.pairlist[i].append('-','-',float(0))
                else:
                    minDist = self.get_minmResDist_chains(structure,pair[0],pair[3],chainsB,chainsB,resomitB,resomitB)
            elif pair[1]=='A' and pair[4]=='B':
                minDist = self.get_minmResDist_chains(structure,pair[0],pair[3],chainsA,chainsB,resomitA,resomitB)
            elif pair[1]=='B' and pair[4]=='A':
                minDist = self.get_minmResDist_chains(structure,pair[0],pair[3],chainsB,chainsA,resomitB,resomitA)
            minDist = np.array([minDist[0],minDist[1],minDist[2]])
            aux  = np.concatenate((self.pairlist[i],minDist))
            newpairlist.append(aux)
            #print(aux)
        self.pairlist = np.asarray(newpairlist)
            
    def get_minmResDist_chains(self,structure,res1,res2,chains1,chains2,resomit1,resomit2):
        if (int(res1) in resomit1) or (int(res2) in resomit2): return "-","-","-"
        mindistance = float("inf")
        for chain1 in chains1:
            for chain2 in chains2:
                aux_res1=[chain1,int(res1)]
                aux_res2=[chain2,int(res2)]
                dist=self.minimResDist(structure,aux_res1,aux_res2)
                if dist < mindistance:
                    mindistance = dist
                    ch1 = chain1
                    ch2 = chain2
        return str(mindistance),ch1,ch2

    #Minimum distance between all atoms of two residues
    def minimResDist(self,structure,res1,res2):
        mindistance = float("inf")
        res1_cont = structure[0][res1[0]][res1[1]]
        res2_cont = structure[0][res2[0]][res2[1]]
        for atom1 in res1_cont:
            for atom2 in res2_cont:
                distance = atom1-atom2
                if distance < mindistance:
                    mindistance = distance
        return mindistance

def get_BinaryClass(coevlist,dist_T,prob_T):
    """
    """
    keys =["TP","FN","FP","TN"]
    values = [0] * 4
    for pair in coevlist:
        if pair[7]!="-":
            if float(pair[7]) <= dist_T:
                if float(pair[6]) > prob_T:
                    values[0] += 1
                else:
                    values[1] += 1
            else:
                if float(pair[6]) > prob_T:
                    values[2] += 1
                else:
                    values[3] += 1
    counts = dict(zip(keys, values))
    return counts

def plot_RocPlot(coevlist,dist_T,intervals=30):
    opt_accuracy = 0
    opt_threshold_a = 0
    opt_precision = 0
    opt_threshold_p = 0
    probs = np.asarray(coevlist[:,6],dtype=float)
    thresholds = np.linspace(np.amin(probs),np.amax(probs),intervals)
    TPRs = []
    FPRs = []
    precisions = []
    for t in thresholds[:-1]:
        counts = get_BinaryClass(coevlist,dist_T,t)
        TPR = float(counts["TP"])/(float(counts["TP"])+float(counts["FN"]))
        FPR = float(counts["FP"])/(float(counts["FP"])+float(counts["TN"]))
        TPRs.append(TPR)
        FPRs.append(FPR)
        accuracy = (float(counts["TP"])+float(counts["TN"]))/(float(counts["TP"])+float(counts["TN"])+float(counts["FN"])+float(counts["FP"]))
        precision = (float(counts["TP"]))/(float(counts["TP"])+float(counts["FP"]))
        precisions.append(precision)
        if accuracy > opt_accuracy:
            opt_accuracy = accuracy
            opt_threshold_a = t
        if precision > opt_precision:
            opt_precision = precision
            opt_threshold_p = t
    values = zip(FPRs,TPRs)
    values = sorted(values,key=lambda x: x[0])
    values = np.asarray(values)
    FPRs = np.asarray(values[:,0])
    TPRs = np.asarray(values[:,1])
    auc = np.trapz(TPRs,x=FPRs)
    plt.plot(FPRs,TPRs)
    plt.xlabel("False positive rate",fontsize=13)
    plt.ylabel("True positive rate",fontsize=13)
    plt.plot([0, 1], [0, 1], 'k-')
    plt.axis('scaled')
    axes = plt.gca()
    plt.legend()
    axes.set_xlim([0, 1])
    axes.set_ylim([0, 1])
    plt.figure()
    precisions = np.asarray(precisions)
    thresholds = np.asarray(thresholds[:-1])
    plt.plot(thresholds,precisions)
    keys = ["Opt_Accuracy","Opt_Accuracy_T","Opt_Precision","Opt_Precision_T","AUC"]
    values = [opt_accuracy, opt_threshold_a, opt_precision, opt_threshold_p, auc]
    dic = dict(zip(keys, values))
    print(dic)
    plt.show()
    return dic
 
if __name__ == '__main__':
   
    #print("TssGF---------------------------------------------------------------") 
    liGF = CoevList('data/coevolution/TssGF_pairs_raptorGremlin.txt',366)
    liGF.add_distances(pdb='data/structures/6n38.pdb',chainsA=['G'],chainsB=['H','I'],resomitA=list(range(1,122))+list(range(331,337)),resomitB=list(range(1,47))+list(range(393,407)))
    interliGF = liGF.get_interpairs()
    for pair in interliGF:
        #print(pair)
        pass
   
    #print("TssKG---------------------------------------------------------------") 
    liKG = CoevList('data/coevolution/TssKG_pairs_raptorGremlin.txt',445)
    liKG.add_distances(pdb='data/structures/6n38.pdb',chainsA=['A','B','C','D','E','F'],chainsB=['G'],resomitA=list(range(1,2))+list(range(220,233))+list(range(312,446)),resomitB=list(range(1,122))+list(range(331,336)))    
    interliKG = liKG.get_interpairs()
    for pair in interliKG:
        #print(pair)
        pass

    distThreshold = 8
    
    liGFKG = np.asarray(list(liGF.pairlist)+list(liKG.pairlist))
    plot_RocPlot(liGFKG,dist_T=distThreshold)   
 
    print("TssKctdL---------------------------------------------------------------")

    liKctdL = CoevList('data/coevolution/TssKctdL_pairs_raptorGremlin.txt',445,314)
    interliKctdL = liKctdL.get_interpairs()
    for pair in interliKctdL[0:10]:
        print(pair)
        pass


