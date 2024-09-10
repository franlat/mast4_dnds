### Modules ###

import sys
import argparse
import itertools
import numpy as np

parser = argparse.ArgumentParser()
parser._action_groups.pop()

required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')

required.add_argument("-gff",
                                        "--gff_file",
                                        dest="gff_file",
                                        action="store",
                                        type=str,
                                        required=True,
                                        help="GFF file input")

required.add_argument("-vcf",
                                        "--vcf_file",
                                        dest="vcf_file",
                                        action="store",
                                        type=str,
                                        required=True,
                                        help="VCF file from SnpEff output")

required.add_argument("-o",
                                        "--output",
                                        dest="outfix",
                                        type=str,
                                        action="store",
                                        required=True,
                                        help="Name of the output file to store all the data")

required.add_argument("-f",
                                        "--fasta",
                                        dest="infasta",
                                        action="store",
                                        type=str,
                                        required=True,
                                        help="Fasta file with contigs matching the gff3 file")




arguments = parser.parse_args()

def read_fasta_file(infile):
        """
        Parses a fasta file with one or more sequences and generates a
        dictionary, in which each sequence_name has a sequence string 
        associated.

        i.e.
        > final_list = read_fasta_file(sequences.file)
        > final_list[contig1]
        "ATGTTCCGCGAACATGCGAGC"

        """
        ### Input ###

        fh1 = open(infile, 'r') ## reading the fasta file

        ### Variables ###

        h = 0 ## counter for number of lines
        all_sequences = {}  ## list of sequences

        ### Body of the function ###

        for line in fh1:
                h += 1
                if line.startswith(">"):                                       ## if it's a header
                        if h > 1:                                                  ## if it's not the first header
                                all_sequences[sequence_name] = sequence        ## append the info to the final list

                        sequence_name = line.lstrip(">").rstrip()                  ## remove end of line characters and the initial '>'
                        sequence = ''                                              ## generate an empty sequence string

                else:                                                          ## if it's a sequence line
                        sequence += line.rstrip()                                  ## concatenate them

        all_sequences[sequence_name] = sequence                    ## because it's not going to to find any header, we need to store the last sequence info into the array

        ### Closing input ###

        fh1.close()

        ### Return ###

        return all_sequences

def get_CDS_from_gff3(infile):
        """
        Parses a GFF file to get the three dictionaries: 
        1) Contig position for the CDS region for each gene,
        2) In which strand the gene is coded for each gene.
        3) To which contig each gene belongs to

        i.e.
        > gene_cds, gene_strand, gene_contig = get_CDS_from_gff3(file.gff3)
        > gene_cds[gene1]
        "[100,101,102,103,500,501,503]"
        > gene_strand[gene1]
        "+"

        """
        ### Input ###

        fh1 = open(infile, 'r') ## reading the gff file

        ### Variables ###

        h = 0                   ## counter for number of CDS lines
        gene_list = []          ## list of genes
        gene_cds = {}           ## Dictionary for gene CDS positions
        gene_strand = {}        ## Dictionary for gene strand
        gene_contig = {}        ## Dictionary for gene and contig they belong

        ### Body of the function ###

        # Regular GFF3 format for my genes:

        #scaffold124_size45503  AUGUSTUS    gene    1   1283    0.32    +   .   ID=g1
        #scaffold124_size45503  AUGUSTUS    transcript  1   1283    0.32    +   .   ID=g1.t1;Parent=g1
        #scaffold124_size45503  AUGUSTUS    intron  1   523 0.74    +   .   Parent=g1.t1
        #scaffold124_size45503  AUGUSTUS    intron  535 591 0.52    +   .   Parent=g1.t1
        #scaffold124_size45503  AUGUSTUS    CDS 524 534 0.64    +   1   ID=g1.t1.cds;Parent=g1.t1
        #scaffold124_size45503  AUGUSTUS    CDS 592 1283    0.38    +   2   ID=g1.t1.cds;Parent=g1.t1
        #scaffold124_size45503  AUGUSTUS    stop_codon  1281    1283    .   +   0   Parent=g1.t1

        for line in fh1:

                if line.startswith("#"):                                       ## if it's a comment, skip it
                    continue
                else:
                    record = line.split("\t")
                    feature = record[2]
                    if feature == "CDS":                                        ## get only the CDSs lines
                        h += 1
                        gene = record[8].split("=")[1].split(".")[0]  # Parse the gene name >>>> WARNING: This may have to change depending on the name of your genes!!! <<<<<

                        if gene in gene_list:  
                            spos = int(record[3])
                            fpos = int(record[4])
                            strand = record[6]
                            contig = record[0]
                            positions = positions + list(range(spos,fpos+1,1))
                        else:
                            gene_list.append(gene)
                            if h > 1:
                                gene_cds[gene_old] = positions
                                gene_strand[gene_old] = strand
                                gene_contig[gene_old] = contig
                            

                            spos = int(record[3])
                            fpos = int(record[4])
                            strand = record[6]
                            contig = record[0]
                            positions = list(range(spos, fpos+1, 1))

                        gene_old = gene

        gene_cds[gene_old] = positions      ## Add the last gene, which won't meet the criteria inside the loop.
        gene_strand[gene_old] = strand
        gene_contig[gene_old] = contig

    ### Closing inoput

        fh1.close()

    ### Return ###

        return (gene_cds, gene_strand, gene_contig)

def read_vcf_file(infile):
    '''
    Reads a VCF file from the SnpEff output and returns 
    a a list for each position including contig_name, 
    syn_genes, miss_genes and counts 

    > pos_info = read_vcf_file(vcf_file)
    > pos_info[1]
    "contig_num1"
    > pos_info[0]
    "176"
    '''
    h = 0 # counter for lines

    ## typical line of a VCF file: 
    ## Contig   Pos  ID   Ref  Allele  Quality  SMTH   INFO FORMAT  [SAMPLES]
    ## In INFO you can find the ANN parameter, which includes:
        ## ##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
    
    ## scaffold68_size42219 795 .   A   G   3.349   .   AB=0;ABP=0;AC=1;AF=0.111111;AN=9;AO=4;CIGAR=1X;DP=24;DPB=24;DPRA=2.66667;EPP=11.6962;EPPR=4.74748;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=60;NS=9;NUMALT=1;ODDS=0.150253;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=144;QR=741;RO=20;RPL=2;RPP=3.0103;RPPR=3.44459;RPR=2;RUN=1;SAF=2;SAP=3.0103;SAR=2;SRF=11;SRP=3.44459;SRR=9;TYPE=snp;technology.illumina=1;ANN=G|missense_variant|MODERATE|g1|g1|transcript|g1.t1|protein_coding|1/1|c.794A>G|p.Asp265Gly|795/2764|794/2763|265/920||WARNING_TRANSCRIPT_NO_START_CODON,G|downstream_gene_variant|MODIFIER|g2|g2|transcript|g2.t1|protein_coding||c.*2480T>C|||||2480|   GT:DP:AD:RO:QR:AO:QA:GL .   .   ..  .   .   .   .   .   0:2:2,0:2:71:0:0:0,-6.74038 .   .   .   .   .   0:1:1,0:1:41:0:0:0,-4.09457 0:1:1,0:1:41:0:0:0,-4.09457 0:4:4,0:4:134:0:0:0,-12.3908    .   .   .   .   .   .   .   .   .   ..  .   .   .   .   .   .   .   .   1:6:2,4:2:81:4:144:-5.62771,0   .   .   .   .   .   ..  .   0:2:2,0:2:75:0:0:0,-7.11959 .   .   0:5:5,0:5:193:0:0:0,-17.7396    0:1:1,0:1:38:0:0:0,-3.79727 .   0:2:2,0:2:67:0:0:0,-6.36304 .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   ..  .   .   .   .   .   .   .   .   .   .   .   . 
    ## For FORMAT we have GT:DP:AD:RO:QR:AO:QA:GL for each sample: 0:2:2,0:2:71:0:0:0,-6.74038
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype"> 
        ##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
        ##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count">
        ##FORMAT=<ID=QR,Number=1,Type=Integer,Description="Sum of quality of the reference observations">
        ##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">
        ##FORMAT=<ID=QA,Number=A,Type=Integer,Description="Sum of quality of the alternate observations">

    ## What we want is: Pos Contig ListGenesAfecting Syn/MissSNP [DP RO AO]persample
        
    fh1 = open(infile, 'r')

    pos_info = []

    for line in fh1:
        if line.startswith("#"):
            continue
        else:
            record = line.rstrip().lstrip().split("\t")  ## split lines by tab
            contig = record[0]   # get the contig name
            pos = str(record[1].rstrip().lstrip())    # get the pos in the contig
            ANN = record[7].split(";")[42].split(",")  # get the ANN field
            vtype = record[3].rstrip().lstrip()

            genes_syn = [gene.split("|")[3] for gene in ANN if gene.split("|")[1] == "synonymous_variant" and len(vtype) == 1 ]  ## split the ANN field to get the gene names if they have synonymous variants, length of this list can vary
            genes_miss = [gene.split("|")[3] for gene in ANN if gene.split("|")[1] == "missense_variant" and len(vtype) == 1 ]   ## split the ANN field to get the gene names if they have missense variants, length of this list can vary
            
            if len(genes_syn) != 0 or len(genes_miss) !=0:  ## If a position doesn't have any of these variants, go to next line
                SAMPLES = record[9:]    ## get all the samples counts
                counts = []
                for sp in SAMPLES:
                    if sp != ".":
                        dp = int(sp.split(":")[1])
                        ro = int(sp.split(":")[3])
                        ao = dp - ro
                        counts.append(ao/float(dp))
                    else:
                        counts.append(float(0))
                #counts = [[int(sp.split(":")[1]), int(sp.split(":")[3]), int(sp.split(":")[5])] for sp in SAMPLES]   ## get only the the Depth, Reference coverge, Allele coverage, depth = RO + AO
                #counts = [[0,0,0,0,0,0,0,0] for sp in SAMPLES if sp == "."]   ## change the missing values to 0s

                pos_info.append([pos, contig, genes_syn, genes_miss, counts])

    fh1.close()

    return pos_info

def get_cds_seq(gene_cds, gene_strand, gene_contig, sequences):

    '''
    Reads the output from read_fasta_file() and get_CDS_from_gff3() to
    generate a dictionary with the CDS sequence for each gene.
    
    i.e.
    > cds_seq = get_cds_seq(gene_cds, gene_strand, gene_contig, sequences)
    > cds_seq["g1"]
    'AAATTGATGATTGAT'
    '''
    
    gene_cds_seq = {}

    for gene in gene_cds.keys():

        ## extract CDS sequence from scaffold
        cds_pos = gene_cds[gene]
        cds_strand = gene_strand[gene]
        cds_contig = gene_contig[gene]
        contig_seq = sequences[cds_contig]

        cds_seq = ''

        if cds_strand == "-":
            cds_pos = list(reversed(cds_pos))             ## if the gene is read in the neg direction, reverse the cds sequence
        
        if len(cds_pos) > len(contig_seq):
            print(len(cds_pos), len(contig_seq))

        for pos in cds_pos:
            
            cds_seq += contig_seq[int(pos)-1]

        gene_cds_seq[gene] = [cds_seq, cds_pos]

    return gene_cds_seq

def get_codon_pos(gene_seqs):
    """
    Get the output from get_cds_seq and return
    a new dictionary with the codon belonging 
    to each position in the gene.
    
    i.e
    > codon_pos = get_codon_pos(gene_seqs)
    > codon_pos["g1"][0][0]
    "ATG"
    > codon_pos["g1"][1][0]
    345
    """

    codon_pos = {} 
    for g in gene_seqs.keys():
        s = gene_seqs[g][0]
        codon = list(itertools.chain.from_iterable(itertools.repeat(s[i:i+3], 3) for i in range(0, len(s), 3)))
        codon_pos[g] = [codon, gene_seqs[g][1]]

    return codon_pos

def get_genes_counts(pos_info):
    """
    Parse the input from read_vcf_file and returns 
    two dictionaries of genes, with counts and pos
    and contig info for both synonymous and missense
    variants.

    i.e.

    counts_genes_syn, counts_genes_miss = get_genes_counts(pos_info)
    > counts_genes_syn["g1"]
    """

    counts_genes_syn = {}
    counts_genes_miss = {}
    list_syn = []
    list_miss = []

    for record in pos_info:
        pos = record[0]
        contig = record[1]
        genes_syn = record[2]
        genes_miss = record[3]
        counts = record[4]

        for g in genes_syn:
            if g in list_syn:
                counts_genes_syn[g].append([contig, pos, counts])
            else:
                list_syn.append(g)
                counts_genes_syn[g] = [[contig, pos, counts]]

        for g in genes_miss:
            if g in list_miss:
                counts_genes_miss[g].append([contig, pos, counts])
            else:
                list_miss.append(g)
                counts_genes_miss[g] = [[contig, pos, counts]]

    return (counts_genes_syn, counts_genes_miss)

 #counts = [map(sum, zip(*t)) for t in zip(X, Y)]

def column_sum(lst):

    return [sum(i) for i in zip(*lst)]

## Dictionarites 

miss_codon = { 'TTT':'2.666666667','TTC':'2.666666667','TTA':'2.333333333','TTG':'2.333333333','CTT':'2','CTC':'2','CTA':'1.666666667','CTG':'1.666666667','ATT':'2.333333333','ATC':'2.333333333','ATA':'2.333333333','ATG':'3','GTT':'2','GTC':'2','GTA':'2','GTG':'2','TCT':'2','TCC':'2','TCA':'2','TCG':'2','CCT':'2','CCC':'2','CCA':'2','CCG':'2','ACT':'2','ACC':'2','ACA':'2','ACG':'2','GCT':'2','GCC':'2','GCA':'2','GCG':'2','TAT':'2.666666667','TAC':'2.666666667','TAA':'2.333333333','TAG':'2.666666667','CAT':'2.666666667','CAC':'2.666666667','CAA':'2.666666667','CAG':'2.666666667','AAT':'2.333333333','AAC':'2.333333333','AAA':'2.666666667','AAG':'2.666666667','GAT':'2.333333333','GAC':'2.333333333','GAA':'2.666666667','GAG':'2.666666667','TGT':'2.666666667','TGC':'2.666666667','TGA':'2.666666667','TGG':'3','CGT':'2','CGC':'2','CGA':'1.666666667','CGG':'1.666666667','AGT':'2.666666667','AGC':'2.666666667','AGA':'2.333333333','AGG':'2.333333333','GGT':'2','GGC':'2','GGA':'2','GGG':'2' } 
syn_codon = { 'TTT':'0.333333333','TTC':'0.333333333','TTA':'0.666666667','TTG':'0.666666667','CTT':'1','CTC':'1','CTA':'1.333333333','CTG':'1.333333333','ATT':'0.666666667','ATC':'0.666666667','ATA':'0.666666667','ATG':'0','GTT':'1','GTC':'1','GTA':'1','GTG':'1','TCT':'1','TCC':'1','TCA':'1','TCG':'1','CCT':'1','CCC':'1','CCA':'1','CCG':'1','ACT':'1','ACC':'1','ACA':'1','ACG':'1','GCT':'1','GCC':'1','GCA':'1','GCG':'1','TAT':'0.333333333','TAC':'0.333333333','TAA':'0.666666667','TAG':'0.333333333','CAT':'0.333333333','CAC':'0.333333333','CAA':'0.333333333','CAG':'0.333333333','AAT':'0.666666667','AAC':'0.666666667','AAA':'0.333333333','AAG':'0.333333333','GAT':'0.666666667','GAC':'0.666666667','GAA':'0.333333333','GAG':'0.333333333','TGT':'0.333333333','TGC':'0.333333333','TGA':'0.333333333','TGG':'0','CGT':'1','CGC':'1','CGA':'1.333333333','CGG':'1.333333333','AGT':'0.333333333','AGC':'0.333333333','AGA':'0.666666667','AGG':'0.666666667','GGT':'1','GGC':'1','GGA':'1','GGG':'1' }

if __name__ == "__main__":


    sequences = read_fasta_file(arguments.infasta)  ## contig sequences 
    gene_cds, gene_strand, gene_contig = get_CDS_from_gff3(arguments.gff_file)  ## gene info, cds positions in contig, strand and contig to which they belong
    pos_info = read_vcf_file(arguments.vcf_file)  ## vcf info from SNPs syn and missense mutations for each position per contig
    gene_seqs = get_cds_seq(gene_cds, gene_strand, gene_contig, sequences) ## get the CDS sequences to read CODONs per gene
    counts_genes_syn, counts_genes_miss = get_genes_counts(pos_info)  ## get for each gene, all the positions that have syn and miss variations
    codon_pos = get_codon_pos(gene_seqs)  ## per gene, get all the positions and its codon associated

    fh2 = open(arguments.outfix,'w')

    header="Gene\tTARA100\tTARA102\tTARA106\tTARA109\tTARA110\tTARA111\tTARA112\tTARA113\tTARA11\tTARA122\tTARA123\tTARA124\tTARA125\tTARA128\tTARA129\tTARA130\tTARA131\tTARA132\tTARA135\tTARA136\tTARA137\tTARA138\tTARA139\tTARA142\tTARA143\tTARA144\tTARA145\tTARA146\tTARA147\tTARA148\tTARA149\tTARA150\tTARA151\tTARA152\tTARA16\tTARA18\tTARA20\tTARA22\tTARA23\tTARA25\tTARA30\tTARA32\tTARA34\tTARA36\tTARA38\tTARA39\tTARA41\tTARA42\tTARA43\tTARA45\tTARA46\tTARA4\tTARA51\tTARA52\tTARA58\tTARA5\tTARA64\tTARA65\tTARA66\tTARA67\tTARA68\tTARA6\tTARA70\tTARA72\tTARA76\tTARA78\tTARA7\tTARA80\tTARA81\tTARA82\tTARA83\tTARA84\tTARA85\tTARA86\tTARA89\tTARA92\tTARA93\tTARA95\tTARA96\tTARA97\tTARA98\tTARA9"

    ######fh2.write("\t".joint())(len(codon_pos.keys()))

    fh2.write(header+"\n")

    for g in codon_pos.keys():
        r = len(codon_pos[g][0])/3  # number of codons in gene
        #print(g,r)

        if g in counts_genes_syn.keys():
            norm_syn_counts = []
            for i in range(0,len(counts_genes_syn[g])):
                pos = counts_genes_syn[g][i][1]
                counts = counts_genes_syn[g][i][2]

                codon = codon_pos[g][0][codon_pos[g][1].index(int(pos))]

                if codon not in syn_codon.keys():
                    codon = 'ATG'
                
                n = syn_codon[codon]    ## get the number of possible syn mutations for a specific position
                if float(n) == 0:
                    n = 1
                counts = [float(nd)/float(n) for nd in counts]                            ## divide the normalized counts by coverage by the possible number of mutations
                norm_syn_counts.append(counts)                                     

            total_syn_counts = column_sum(norm_syn_counts)                             ## sum all the positions of a gene per sample
            total_syn_counts = [float(ts)/r for ts in total_syn_counts]
            ds_values = [(-3/4)*np.log(1-(4*float(ps)/3)) for ps in total_syn_counts]
        else:
            ds_values = [0]*82          ## this can change depending on the number of samples 

        if g in counts_genes_miss.keys():
            norm_miss_counts = []
            for i in range(0,len(counts_genes_miss[g])):
                pos = counts_genes_miss[g][i][1]
                counts = counts_genes_miss[g][i][2]
                
                codon = codon_pos[g][0][codon_pos[g][1].index(int(pos))]

                if codon not in miss_codon.keys():
                    codon = 'ATG'

                n = miss_codon[codon]    ## get the number of possible syn mutations for a specific position
                if float(n) == 0:
                    n = 1
                counts = [float(nd)/float(n) for nd in counts]                            ## divide the normalized counts by coverage by the possible number of mutations
                norm_miss_counts.append(counts)                                     

            total_miss_counts = column_sum(norm_miss_counts)                             ## sum all the positions of a gene per sample
            total_miss_counts = [float(ts)/r for ts in total_miss_counts]
            dn_values = [(-3/4)*np.log(1-(4*float(pn)/3)) for pn in total_miss_counts]   
        else:
            dn_values = [0]*82            ##this can change depending on the number of samples       

        dnds_values = [abs(float(dn)/(ds+1))+0 if ds == 0 else abs(float(dn)/ds)+0 for dn,ds in zip(dn_values, ds_values)]

        res = "\t".join([format(x, "10.5f") for x in dnds_values])

        fh2.write(g.rstrip().lstrip()+"\t"+res+"\n")

    fh2.close()
         







        


        
     
    

        






    
    


