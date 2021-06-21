import numpy as np
from Bio.PDB import *
from Bio import *
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import sklearn
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LinearRegression

class CoevList(object):
    """
    Raptor (concatenated MSA) + Gremlin (coevolving pairs) list of coevolving pairs
    between gene A and gene B
    """
    def __init__(self,pairlist,Alen,Astart=1,Bstart=1):
        """
        pairlist -> .txt outcome from gremlin server with all pairs
        Alen -> integer lenght of the gene A
        Adisp -> integer starting residue of the gene A
        Bdisp -> integer startung residue of the gene B
        """
        self.Alen = Alen
        self.Astart = Astart
        self.Bstart = Bstart
        self.ABbreak = self.Alen-self.Astart+1
        self.pairlist = np.genfromtxt(pairlist,dtype=str,delimiter='\t',skip_header=1)
        self.parse_list()
     
    def parse_list(self):
        """
        Initialization function. It parses the list of coevolving residues according
        to the initialization variables Alen, Astart and Bstart
        """
        newpairlist=[]
        for pair in self.pairlist:
            newpair=[]
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
            newpair.append(pair[-1])
            newpairlist.append(newpair)
        self.pairlist = np.asarray(newpairlist)

    def get_intrapairs(self):
        """
        Return a list with the intra coevolving pairs (residues from the same
        gene A-A or B-B)
        """
        newpairlist=[]
        for pair in self.pairlist:
            if pair[1]==pair[4]: newpairlist.append(pair)
        return(np.asarray(newpairlist))
    
    def get_interpairs(self):
        """
        Return a list with the inter coevolving pairs (residues from
        different genes A-B)
        """
        newpairlist=[]
        for pair in self.pairlist:
            if pair[1]!=pair[4]: newpairlist.append(list(pair))
        return(np.asarray(newpairlist))
        #return(newpairlist)

    def add_distances(self,pdb,chainsA,chainsB,resomitA,resomitB):
        """
        Using a PDB structure of the proteins from which the coevolving pairs had been
        calculated add an extra column to the list with the residue pairs minimum
        distance among all their atoms
        pdb -> directory to the pdb file
        chainsA -> list of the pdb chains of the protein A
        chainsB -> list of the pdb chains of the protein B
        resomitA -> list of residues of the protein A not present in the PDB
        resomitB -> list of residues of the protein B not present in the PDB
        """
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
        """
        Compute the minimum distance between two residues (res1, res2) in 'structure'
        by selecting the minimum one between chains1 and chains2
        """
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
        """
        Get the minimum distance between two residues (res1, res2)  in 'structure'
        as the minimum distance of all their atoms
        """
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
    Given a list of coevolving residues with the minimum distance of its parts,
    perform a binary classification with a given distance and probability theshold
    coevlist -> CoevList object
    dist_T -> float distance threshold
    prob_T -> float probability threshold
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

def plot_contact_isotonicReg(coevlist,dist_T):
    problist=[]
    isincontact=[]
    for pair in coevlist:
        if pair[7] == '-': continue
        problist.append(float(pair[6]))
        isincontact.append(1.0 if float(pair[7]) < dist_T else 0.0)
    problist, isincontact = zip(*sorted(zip(problist, isincontact),reverse=False))
    ir = IsotonicRegression(increasing=False).fit(problist,isincontact)
    # Plot
    plt.figure()
    plt.yticks((0,1))
    plt.plot(problist,isincontact,'o',label='TssGF/TssKG pairs',color='blue')
    ir = IsotonicRegression(out_of_bounds="clip")
    isincontact_ = ir.fit_transform(problist, isincontact)
    plt.plot(problist,isincontact_, 'C1.-', markersize=5,color='orange',label='Isotonic Regression')
    plt.axvline(x=0.9,color='black',linestyle='--',label='Threshold (0.9)')
    plt.xlabel("Contact Probability",fontsize=13)
    plt.legend()
    plt.savefig("plots/isincontact.png")

def plot_RocPlot(coevlist,dist_T,intervals=30,reg=False):
    """
    Given a list of coevolving residues with the minimum distance of its parts and
    a distance threshold, calculate its receaving operating curve (ROC curve)
    coevlist -> CoevList object
    dist_T -> float distance threshold
    intervals -> integer
    """
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
    precisions = np.asarray(precisions)
    thresholds = np.asarray(thresholds[:-1])
    values = zip(FPRs,TPRs)
    values = sorted(values,key=lambda x: x[0])
    values = np.asarray(values)
    FPRs = np.asarray(values[:,0])
    TPRs = np.asarray(values[:,1])
    auc = np.trapz(TPRs,x=FPRs)
    # Plot ROC
    plt.figure()
    plt.plot(FPRs,TPRs)
    plt.xlabel("False positive rate",fontsize=13)
    plt.ylabel("True positive rate",fontsize=13)
    plt.plot([0, 1], [0, 1], 'k-')
    plt.axis('scaled')
    axes = plt.gca()
    plt.legend()
    axes.set_xlim([0, 1])
    axes.set_ylim([0, 1])
    # Get isotonic and linear regression for precision
    if reg:
        ir = IsotonicRegression(out_of_bounds="clip")
        precisions_ = ir.fit_transform(thresholds, precisions)
        #lr = LinearRegression()
        #lr.fit(thresholds[:, np.newaxis], precisions)
    # Plot evolution of precision
    plt.figure()
    if reg:
        plt.plot(thresholds,precisions_, 'C1.-', markersize=5,color='orange',label='Isotonic Regresion')
        #plt.plot(thresholds, lr.predict(thresholds[:, np.newaxis]), 'C2-',color='green')
    plt.plot(thresholds,precisions,'o',markersize=5,color='blue')
    plt.xlabel("Contact Probability",fontsize=13)
    plt.ylabel("Precision",fontsize=13)
    plt.axvline(x=0.9,color='black',linestyle='--',label='Threshold (0.9)')
    axes.set_xlim([0, 1])
    axes.set_ylim([0, 1])
    plt.legend()
    plt.savefig("plots/precision.png")
    # Return
    keys = ["Opt_Accuracy","Opt_Accuracy_T","Opt_Precision","Opt_Precision_T","AUC"]
    values = [opt_accuracy, opt_threshold_a, opt_precision, opt_threshold_p, auc]
    dic = dict(zip(keys, values))
    print(dic)
    return dic 
