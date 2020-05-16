def pfas_brics(smiles, patt):
  from rdkit import Chem 
  from rdkit.Chem import BRICS
  import numpy as np
  import pandas as pd
  mol = Chem.MolFromSmiles(smiles)
  patt = Chem.MolFromSmarts(patt)
  if mol.HasSubstructMatch(patt):
    bonds = mol.GetSubstructMatches(patt)
    bonds = [((x,y),(0,0)) for x,y in bonds]
    p = BRICS.BreakBRICSBonds(mol, bonds = bonds)
    mols = Chem.GetMolFrags(p,asMols=True)
    molwt = [Chem.rdMolDescriptors.CalcExactMolWt(m) for m in mols]
    flct = [len(m.GetSubstructMatches(Chem.MolFromSmarts('[F]'))) for m in mols]
    d = pd.DataFrame({"Molwt":molwt, "F_ct":flct}).sort_values(by = ['Molwt', 'F_ct'])
    mols = [mols[i] for i in d.index.values]
    smis = [Chem.MolToSmiles(x,True) for x in mols]
  else:
    smis = []
  return(smis)
