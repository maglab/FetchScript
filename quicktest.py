import fetch

f = fetch.FetchReference()
#print f.fetchPubmed(22737090)
#print f.fetchPubmed(20177970)
# has terms, from our lab
#print f.fetchPubmed(23039964, withTerms=True)
# has no terms
#print f.fetchPubmed(22998224, withTerms=True)

f = fetch.FetchDetails()
#f.fetchDetailsFromNucleotide('NM_000163.4')
#up = f.entrezToUniProtID(71)
#uni = f.translateID(2690, 'P_ENTREZGENEID', 'ID')
#nuc = f.translateID(uni, 'ID', 'REFSEQ_NT_ID')
#print f.fetchDetailsFromNucleotide('NM_000163.4')
#print f.fetchDetailsFromUniProt(uni)
#print f.fetchDetailsFromEntrez(2690)
#print f.symbolToEntrezGeneID('HLA-DRB1')
print f.fetchDetailsFromdbSNP('rs10736086')
print f.fetchDetailsFromdbSNP('rs1205035')
print f.fetchDetailsFromdbSNP('rs11568820')
