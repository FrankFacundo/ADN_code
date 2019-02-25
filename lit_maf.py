# -*- coding: utf-8 -*-

def lit_seq(ligne):
  """retourne un dictionnaire qui contient chacun des éléments de la
    ligne.
    """
  ret_dico = {}     # ou ret_dico = dict()
  seq = ligne.split()
  desc=seq[1]
  pos=0
  while pos<len(desc) and not desc[pos].isdigit():
    pos+=1
  pos2=desc.find(".")

  ret_dico['esp'] = desc[:pos]
  ret_dico['vers'] = int(desc[pos:pos2])
  ret_dico['chrom'] = desc[pos2+1:]
  ret_dico['pos'] = int(seq[2])
  ret_dico['nbl'] = int(seq[3])
  ret_dico['sens'] = seq[4]
  ret_dico['lgchr'] = int(seq[5])
  ret_dico['seq'] = seq[6]
  return(ret_dico)

def lit_ali(fichier, lig_score):
  """ficher est un file ouvert en lecture au début de la première
    ligne de séquence d'un alignement. Cette fonction retourne un
    dictionnaire dont les clefs sont les noms d'espèces et les valeurs les
    dictionnaires construits par lit_seq. Ce dictionnaire a aussi la clef
    score pour le score de l'alignement."""

  dico = {}

  dico["score"]=float(lig_score[lig_score.find("=")+1:])

  ligne=fichier.readline()
  while len(ligne)>1:
    if ligne[0]=="s":
      ret_dico = lit_seq(ligne)
      dico[ret_dico['esp']] = ret_dico
    ligne = fichier.readline()
  return(dico)

def lit_maf(nom_fichier):
  """Sur un nom de fichier .maf, cette fonction retourne une liste de
  dictionnaires construits par lit_ali."""

  liste_align = []
  try:                     # Pour gerer un fichier inconnu
    fichier = open(nom_fichier, 'r')
  except IOError:       #  gestion des exceptions
    print("Fichier inconnu: ", nom_fichier)
    return

  lig_score=fichier.readline()
  while lig_score:
    if lig_score.find("a score")!=-1:
      dico = lit_ali(fichier, lig_score)
      liste_align.append(dico)
    lig_score=fichier.readline()

  fichier.close()
  return(liste_align)
