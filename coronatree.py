import os
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import pylab
from Bio import Phylo
Entrez.email = "juliagolebiowskaa@gmail.com"

#analiza dla korona wirusów, tworzę listę nazw
name1="2019-nCov"
name2="SARS"
name3="bat coronavirus"
name4="MERS"
name5="influenza A"
name6="hepatitis A"

listofnames=[name1, name2, name3, name4, name5, name6]
#interesują nas białka "spike"
proteinname="spike"

def read(listofnames, proteinname, proteinfilename, genomefilename):
    listofproteins=[]
    listofgenomes=[]
    for name in listofnames:
        #pobierz genomy
        #obsługa wyjątków
        recordgenome = SeqRecord(seq="")
        try:
            search_term = name
            handle = Entrez.esearch(db="nuccore", term=search_term)
            rec = Entrez.read(handle)
            rec_handle = Entrez.efetch(db="nucleotide", id=rec["IdList"][2], rettype="fasta")
            genome = SeqIO.read(rec_handle, "fasta")
            recordgenome=SeqRecord(seq=genome.seq, id=name)
            listofgenomes+=[recordgenome]

        except:
            try:         #zamist kategorii "organizm" użyj kategorii "virus"

                search_term = name+"[virus]"
                handle = Entrez.esearch(db="nuccore", term=search_term)
                rec = Entrez.read(handle)
                rec_handle = Entrez.efetch(db="nucleotide", id=rec["IdList"][2], rettype="fasta")
                genome = SeqIO.read(rec_handle, "fasta")
                recordgenome=SeqRecord(seq=genome.seq, id=name)
                listofgenomes+=[recordgenome]

            except:
                try:
                    # spróbuj użyć pierwszego członu z nazwy
                    name=name.split(" ")
                    search_term = name[0]
                    handle = Entrez.esearch(db="nuccore", term=search_term)
                    rec = Entrez.read(handle)
                    rec_handle = Entrez.efetch(db="nucleotide", id=rec["IdList"][2], rettype="fasta")
                    genome = SeqIO.read(rec_handle, "fasta")
                    recordgenome=SeqRecord(seq=genome.seq, id=name)
                    listofgenomes+=[recordgenome]

                except:

                    pass

        #pobierz białka
        #obsługa wyjątków

        proteinrecord = SeqRecord(seq="", id="")

        try:
            search_term=name+" "+proteinname
            handle = Entrez.esearch(db="protein", term=search_term)
            protein_rec = Entrez.read(handle)
            protein_handle = Entrez.efetch(db="protein", id=protein_rec["IdList"][0], rettype="fasta")
            protein = SeqIO.read(protein_handle, "fasta")
            proteinrecord = SeqRecord(seq=protein.seq, id=name)
            listofproteins += [proteinrecord]

        except:
            #zamist kategorii organizm użyj kategorii virus
            try:
                search_term = "("+name+"[virus])" +" AND " + proteinname+"[Protein Name]"
                handle = Entrez.esearch(db="protein", term=search_term)
                protein_rec = Entrez.read(handle)
                #ogarnąc co jest nie tak dla hepatitis A
                protein_handle = Entrez.efetch(db="protein", id=protein_rec["IdList"][0], rettype="fasta")
                protein = SeqIO.read(protein_handle, "fasta")
                proteinrecord = SeqRecord(seq=protein.seq, id=name)
                listofproteins += [proteinrecord]

            except:
                #spróbuj użyć pierwszego członu z nazwy
                try:
                    name=name.split(" ")
                    search_term = "(" + name[0] + "[virus])" + " AND " + proteinname + "[Protein Name]"
                    handle = Entrez.esearch(db="protein", term=search_term)
                    protein_rec = Entrez.read(handle)
                    # ogarnąc co jest nie tak dla hepatitis A
                    protein_handle = Entrez.efetch(db="protein", id=protein_rec["IdList"][0], rettype="fasta")
                    protein = SeqIO.read(protein_handle, "fasta")
                    proteinrecord = SeqRecord(seq=protein.seq, id=name)
                    listofproteins += [proteinrecord]
                except:
                    pass


        handle.close()
    if len(listofproteins)!=0 and len(listofgenomes)!=0:
        #wrzycamy nasze sekwencje do pliku
        SeqIO.write(listofgenomes, genomefilename, "fasta")
        SeqIO.write(listofproteins, proteinfilename, "fasta")
    else:
        print('ERROR:Data not aviable')
        raise ValueError('Data not aviable')

filegenome ="genome"
fileprotein="spiekeprotein"

try:
    #pobieram  i tworzę pliki fasta dla wybranych genomów i białek

    read(listofnames,proteinname, fileprotein, filegenome)
    #muliple alignment

    #znajdź clustalw
    def find_all(name, path):
        result = []
        for root, dirs, files in os.walk(path):
            if name in files:
                result.append(os.path.join(root, name))
        return result

    clustalpath=find_all("clustalw2","/home")

    #multiple alignment genomów
    clustalw_cline = ClustalwCommandline(clustalpath[0], infile=filegenome)
    clustalw_cline()
    genomealigname=filegenome+".aln"
    aligngenome = AlignIO.read(filegenome+".aln", "clustal")

    #multiple alignment białek
    clustalw_cline = ClustalwCommandline(clustalpath[0], infile=fileprotein)
    clustalw_cline()
    proteinalignname=fileprotein+".aln"
    alignproteine = AlignIO.read(proteinalignname, "clustal")


    #liczę macierz odległości dla genomów
    calculator = DistanceCalculator('identity')
    genomedist= calculator.get_distance(aligngenome)
    constructor = DistanceTreeConstructor()

    #liczę macierz odległości dla białek
    calculator = DistanceCalculator('blosum62')
    proteindist= calculator.get_distance(alignproteine)
    constructor = DistanceTreeConstructor()

    #drzewa genom
    genomedisttreenj=constructor.nj(genomedist)
    genomedisttreeupgma=constructor.upgma(genomedist)

    #drzewa bialka
    proteindisttreenj=constructor.nj(proteindist)
    proteindisttreeupgma=constructor.upgma(proteindist)
    #usuwam etykiety węzłów wewnętrznych - funkcja która zostanie użyta do wizualizacji
    def get_label(leaf):
        if "Inner" in leaf.name:
            leaf.name=""
        return leaf.name
    #wizualizacja
    Phylo.draw(genomedisttreenj, label_func=get_label)
    pylab.savefig('genometreenj')
    Phylo.draw(genomedisttreeupgma, label_func=get_label)
    pylab.savefig('genometreeupgma')
    Phylo.draw(proteindisttreenj, label_func=get_label)
    pylab.savefig('proteinetreenj')
    Phylo.draw(proteindisttreeupgma, label_func=get_label)
    pylab.savefig('proteinetreeupgma')
except:
    print("analysis failed")