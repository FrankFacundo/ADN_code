import lit_maf
import copy
m=lit_maf.lit_maf("ENm010.maf")
align=m[0]#un alignement

def lettre(align,n):
    dico={}
    esp=list(align.keys())
    esp.remove('score')
    for e in esp:
        d=align[str(e)]
        seq=d['seq']
        dico[str(e)]=seq[n]
    return dico
    
def sequence(align,n,m):
    dico={}
    esp=list(align.keys())
    esp.remove('score')
    for e in esp:
        d=align[str(e)]
        seq=d['seq']
        dico[str(e)]=seq[n:m]
    return dico

esp=list(align.keys())
esp.remove('score')
for e in esp:
    d=align[str(e)]
    seq=d['seq']
mot='ATG'

def find (seq, mot):
    pos=[]
    indel=-1
    seqcopie=seq[:]
    seqcopie=seqcopie.replace('-','')
    if len('mot')<len(seqcopie):
        for i in range(len(seqcopie)-len(mot)+1) :
            k=0
            indel+=1
            while seq[indel]=='-' :
                indel+=1
            while k<len(mot) and seqcopie[i+k]==mot[k]:
                    if k==len(mot)-1:
                        pos+=[indel]
                    k+=1
    return pos
    
#appliquer sur les alignements m la fonction find
mot='GATTGGCA'
mot='AAGTGAGT'

def faire_tourner(m,mot):
    liste=[]
    for align in m:
        esp=list(align.keys())
        esp.remove('score')
        e=0
        while e in [i for i in range (len(esp))]:
            if find (align[str(esp[e])]['seq'], mot)!=[] :
                liste+=find(align[str(esp[e])]['seq'], mot) #espèce,sequence,position
            e+=1
    return liste


def alignement(m,mot):
    liste=[]
    for align in m:
        esp=list(align.keys())
        esp.remove('score')
        e=0
        deja_fait=False
        while e in [i for i in range (len(esp))] and deja_fait==False :
            if find (align[str(esp[e])]['seq'], mot)!=[]:
                for pos in find (align[str(esp[e])]['seq'], mot):
                    if int(pos)-10<0:
                        min=0
                    else:
                        min=int(pos)-10
                    if int(pos)+10+len(mot)>len(align[str(esp[e])]['seq'])-1:
                        max=len(align[str(esp[e])]['seq'])-1
                    else:
                        max=int(pos)+10+len(mot)
                    dico=sequence(align,min,max)
                    clés=dico.keys()
                    for c in clés:
                        liste+=[[c,dico[c],pos]] #espèce,sequence,position
                deja_fait=True
            e+=1
    return liste
    
def lg_sr(seq):
    compt=0
    for i in range (len(seq)):
        if seq[i].islower():
            compt+=1
    return compt
        
def elt_rep(m):
    tot={'panTro':0,'hg':0,'rheMac':0} #tot des lettres
    elt={'panTro':0,'hg':0,'rheMac':0} #elt répété pour les trois espèces, panTro, hg, rheMac
    for align in m:
        esp=list(align.keys())
        esp.remove('score')
        for e in esp:
            if e=='panTro':
                elt[e]+=lg_sr(align[e]['seq'])
                tot[e]+=len(align[e]['seq'])
            elif e=='hg':
                elt[e]+=lg_sr(align[e]['seq'])
                tot[e]+=len(align[e]['seq'])
            elif e=='rheMac':
                elt[e]+=lg_sr(align[e]['seq'])
                tot[e]+=len(align[e]['seq'])
            else :
                 elt.append(lg_sr(align[e]['seq']))
                 tot[e]+=len(align[e]['seq'])
    return elt,tot
    
align=m[-1]    
def correction(m):
    for align in m:
        esp=list(align.keys())
        esp.remove('score')
        dico={}
        for es in esp:
            dico[es]=list(align[es]['seq'])
        for i in range (len(align[es]['seq'])):
            for es in esp:
                indice_copy=esp[:]
                indice_copy.remove(es)
                if dico[es][i].isupper() and dico[es][i]!='-' :
                    for f in indice_copy:
                        if  dico[f][i].islower() and dico[f][i]!='-' :
                            dico[es][i]=dico[es][i].lower()
        for es in esp :
            dico[es]=''.join(dico[es])
        for es in esp :
            align[es]['seq']=dico[es]
    return m

def proportion (m):
    elt_m,tot_m=elt_rep(m)
    elt_corr,tot_corr=elt_rep(correction(m))
    avant={}
    apres={}
    for e in elt_m.keys():
        avant[e]=elt_m[e]/tot_m[e]
        apres[e]=elt_corr[e]/tot_corr[e]
    return avant, apres
    
def masque(m):
    for align in m:
        esp=list(align.keys())
        esp.remove('score')
        dico={}
        for es in esp:
            dico[es]=list(align[es]['seq'])
        for i in range (len(align[es]['seq'])):
            for es in esp:
                if dico[es][i].islower() and dico[es][i]!='-' :
                    dico[es][i]='N'
        for es in esp :
            dico[es]=''.join(dico[es])
        for es in esp :
            align[es]['seq']=dico[es]
    return m
    
m_masque=masque(correction(m))

def distance_genome(m_masque):
    dist={}
    dico={}
    esp=list(m_masque[0].keys())
    esp.remove('score')
    for es in esp:
        indice_copy=esp[:]
        indice_copy.remove(es)
        for f in indice_copy:
            dist[(es,f)]=0         #construction de dist
    for align in m_masque:
        esp=list(align.keys())
        esp.remove('score')
        dico={}
        for es in esp:
            dico[es]=list(align[es]['seq']) #reconstruire le dico pour chaque align
        for i in range (len(dico[es])):
            for es in esp:
                indice_copy=esp[:]
                indice_copy.remove(es)
                for f in indice_copy:
                    if dico[es][i]!=dico[f][i] :
                        dist[(es,f)]+=1
    return dist
    
def parcisite(homme,chimp,mac):
    if homme==chimp:
        return homme
    else:
        if homme==mac or chimp==mac:
            return mac
        else:
            return '-'
            
align=m[0]
def parciseq_Hs_Pt(align):
    seq=[]
    esp=list(align.keys())
    esp.remove('score')
    dico={}
    for es in esp:
        dico[es]=list(align[es]['seq']) #reconstruire le dico
    for i in range (len(dico[es])):
        seq+=[parcisite(dico['hg'][i],dico['panTro'][i],dico['rheMac'][i])]
    return ''.join(seq)
    
def matrice_substi(seq_anc,homme):
    mat={}  #subst A,T,C,G
    for lett in ['A','T','G','C','N','-']: #construction de la matrice
        for lettre in ['A','T','G','C','N','-']:
            mat[(lett,lettre)]=0
    for i in range (len(seq_anc)):
        mat[(seq_anc[i],homme[i])]+=1
    return mat
    
def evo_CN(seq_anc,homme):
    dico={}; mat={}
    for lett in ['A','T','G','C','N','-']: #construction de la matrice
        for lettre in ['A','T','G','C','N','-']:
            dico[(lett,lettre)]=0
    for lett in ['A','T','G','C','N','-']: #construction de la matrice
        for lettre in ['A','T','G','C','N','-']:
            mat[(lett,lettre)]=dico.copy()
    for i in range (len(seq_anc)-1):
        if seq_anc[i]=='C':
            mat[(seq_anc[i],seq_anc[i+1])][(homme[i],homme[i+1])]+=1
    return mat
        