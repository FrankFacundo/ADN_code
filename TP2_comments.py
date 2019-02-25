import lit_maf                                                                  #agrega la libreria lit_maf.py
import copy
m=lit_maf.lit_maf("ENm010.maf")                                                 #Hace el parsage del .maf
align=m[0]#un alignemen                                                         #align es un diccionario

"""la funcion "lettre" permite a partir de un diccionario con un triplete de especies (align) y un entero (n), obtener la letra n del triplete de secuencias (dico)"""
"""ejemplo
lettre(align,5)
Out[81]: {'hg': 'C', 'panTro': 'C', 'rheMac': 'C'}"""
def lettre(align,n):                                                            #lettre toma como parametros a el diccionatio "align" y a un entero "n" 
    dico={}                                                                     # creacion del diccionario "dico" 
    esp=list(align.keys())                                                      #"esp" es una lista de todos los campos de un "align" osea "hg", "panTro", "rheMac", "score"
    esp.remove('score')                                                         #quita el index 0 osea 'score'
    for e in esp:                                                               # 'e' va a tomar todos los valores de la lista 'esp' excepto score, osea hg panTro rheMac en cada ciclo de la funcion 'for'
        d=align[str(e)]                                                         # "str(e)" convierte en string a 'e' // 'd' toma el diccionario de 'hg' 'panTro' y 'rheMac'
        seq=d['seq']                                                            # "seq" toma las secuencias de hg panTro y rheMac
        dico[str(e)]=seq[n]                                                     # guarda en el diccionario dico con keys hg, panTro y rheMac las letras de la position 'n'
    return dico                                                                 #retorna el diccionario "dico"
    

"""la funcion "sequence" a partir de un diccionario de un triplete de especies (align), un inicio (n) y un final (m), 
permite obtener las palabras del triplete de secuencias (dico), 
si solo utilizas dos parametros funciona como la funcion 'lettre' osea esta funcion es mas general o generica.
La funcion sequence la he modificado"""
"""ejemplo1
sequence(align,5)
Out[82]: {'hg': 'C', 'panTro': 'C', 'rheMac': 'C'}"""
"""ejemplo2
sequence(align,5,8)
Out[83]: {'hg': 'CATA', 'panTro': 'CATA', 'rheMac': 'CACA'}"""
def sequence(align,n,m=0):                                                        #es casi el mismo codigo que la funcion anterior
    dico={}
    if (m==0):                                                                  #si no utilzamos el tercer parametro la letra final sera "n"
        m = n
    esp=list(align.keys())
    esp.remove('score')
    for e in esp:
        d=align[str(e)]
        seq=d['seq']
        dico[str(e)]=seq[n:m+1]                                                   #la diferencia esta en que dico toma los valores desde n hasta m (es m+1 porque la lista comienza en 0 por lo tanto se le tiene que agregar una unidad)
    return dico

mot='ATG'

"""la funcion "find" a partir de una secuencia (string seq) y una palabra (string mot)
permite obtener la posicion de mot en seq, tienes que tener cuidado con el hecho de que si se repiten varias veces mot, la funcion solo te retornara la ultima repeticion"""

"""ejemplo
align["hg"]["seq"]
Out[92]: 'CTTTACATAACAAGGTTGGAGATCTGAATATAAAAAT'

find(align["hg"]["seq"],"GGT")
Out[93]: [13]
"""

def find (seq, mot):
    pos=[]                                                                      #crea una lista "pos"
    indel=-1                                        
    seqcopie=seq[:]                                                             #copia todo el string "seq"
    seqcopie=seqcopie.replace('-','')                                           #borra todos los '-'
    if len('mot')<len(seqcopie):                                                #nos aseguramos que la palabra que buscamos es mas pequeña que la secuencia misma
        for i in range(len(seqcopie)-len(mot)+1) :                              #hacemos un bucle de "len(seqcopie)-len(mot)+1" repeticiones
            k=0
            indel+=1
            while seq[indel]=='-' :                                             #si existen '-' las obviamos
                indel+=1
            while k<len(mot) and seqcopie[i+k]==mot[k]:                         #verifica si "mot" es igual a "seq" a partir de la casilla 'n'
                    if k==len(mot)-1:                                           #si toda la palabra 'mot' se verifica entonces:
                        pos+=[indel]                                            #asignamos la position "pos" actual como resultado para retornar
                    k+=1
    return pos
    
#appliquer sur les alignements m la fonction find
mot='GATTGGCA'
mot='AAGTGAGT'

"""la funcion "faire tourner" a partir de todo la base de datos (m) y una palabra (string mot)
permite obtener la posicion de mot en m """

"""ejemplo
faire_tourner(m,'AAGTGAGT')
Out[99]: 
[10,
 10,
 10,
 47,
 47,
 47,
 36,
 36,
 98,
 98,
 98,
 23,
 23,
 23,
 125,
 125,
 13,
 13,
 13,
 26,
 20,
 20,
 20,
 103,
 12,
 12,
 12]
"""

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
        