# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 12:29:04 2016

@author: SpiffKringle
"""
from itertools import combinations
#need to factor in number of sequences represented at a position! more = better
#peanlizing number of gaps does this (kind of?)
def getLen(alignment):
    """
    Returns the length of an alignment for use in the GenFrequencies Function
    """
    AlignmentLength = 0
    with open(alignment,'rU') as F:
        F.readline()
        while True:
            Data = F.readline().strip()
            if Data[0] != '>':
                AlignmentLength += len(Data)
            else:
                return AlignmentLength

def GenPositionConsensus(position, cutoff = 0.65,MinPercentBetter = 0.1,GapCutoff = 0.5):
    """
    >Generates the consensus of a position. Takes a dictionary of base frequencies as input
    >e.g: {'A':0.15,'C':0.1,'G':O.4,'T':0.35.'-':0.0}, cutoff 0.7,MinPercentBetter = 0.1
    >scans through till checking combinations of 2, at whcih point G+T = 0.75, is > cutoff
    >Checks against MinPercentBetter, next better 2 combinations is G+A = 0.55, 
    >0.75 - 0.55 = 0.2 > MinPercentBetter, so passes. Otherwise would have gone
     on to to check 3 combinations then 4 combinations.
    >Gaps are subtracted from cutoff so that a combination ofthe bases will always
     be greater than the cutoff. Always returns somethin
    >output format for bases (('A','C'),0.8,0.1) where 
     0.8 = summed frequencies for A + C. 0.1 = gap frequency.
    >used later in scoring primers
    """
    if   GapCutoff - position['-'] < 0.001:#Check to see if it has too many gaps
        return (('-',),position['-'])
    cutoff -= position['-']#otherwise ignore gaps and adjust the cutoff proportionally, basically 1.0 now = 1.0 - %gaps
    threshold = (1-cutoff)*MinPercentBetter#set the threshold, this is the minimum amount better than the next best that a score must be to get accepted, means that ACG = 0.81, AGT = 0.79 returns AGCT not AGC 
    for i in [1,2,3,4]:
        permutations = combinations(['A','C','G','T'],i)
        scores = []
        for permutation in permutations:
            total = 0            
            for base in permutation:
                total += position[base]            
            scores.append((total,permutation))        
        scores.sort(reverse=True)
        if i == 4:
            return  (scores[0][1],scores[0][0],position['-'])#default to returning ACGT (N)

        if cutoff - scores[0][0]< 0.001 and \
        scores[0][0] - scores[1][0] > threshold:
            return  (scores[0][1],scores[0][0],position['-'])
    
    return ('*',)
            
    
def GenFrequencies(alignment):
    """
    >Takes a plain text alignment from Muscle (or other)
    >Returns the frequencies of bases at each position
    >Output is a list of dictionaries
    """
    bases = {'A':0,'C':0,'G':0,'T':0,'-':0}
    FreqArray = []
    SeqLen = getLen(alignment)
    for i in range(SeqLen):
        FreqArray.append(bases.copy())
    count = 0
    SeqNum = 0
    with open(alignment,'rU') as F:
        data = 'placeHolder'
        while data:
            data =  F.readline().strip()
            if data and not data[0] == '>':
                for char in data:
                    FreqArray[count][char] += 1
                    count +=1
            elif data:
                count = 0
                SeqNum += 1
            else:
                break
        for position in FreqArray:
            for base in position:
                position[base] /= float(SeqNum)
        return FreqArray

IUPACDEGEN = {
('A',):'A',
('C',):'C',
('G',):'G',
('T',):'T',
('A','G'):'R',
('C','T'):'Y',
('G','C'):'S',
('A','T'):'W',
('G','T'):'K',
('A','C'):'M',
('A','G','T'):'D',
('A','C','T'):'H',
('A','C','G'):'V',
('C','G','T'):'B',
('A','C','G','T'):'N',
('-',):''
}#make a dict with tuples as keys and IUPAC degenerate bases as values                

def GenConsensus(alignment):
    """
    >Integrates the functions GenFrequencies and GenPositionConsensus to generate
     A consensus sequence
    >Input is an alignment file
    >Output is a list of tuples, where each tuple is generated according to
     GenPositionConsensus, format --> (('A','C'),0.8,0.1)
     """
    count = 0
    consensus = []
    FreqArray = GenFrequencies(alignment)
    for position in FreqArray:
        count +=1
        consensus.append(GenPositionConsensus(position))
    return consensus

consensus = GenConsensus('alignment1.txt')



def DeGap(consensus,cutoff = 0.95):
    """
    Removes gaps from the consensus where greater than cutoff seqs have a gap
    at a particular position. Returns a new consensus in the same format
    """
    newCon = []
    for i in consensus:
        if i[0][0] == '-' and i[1] > cutoff:
            print i
            continue
        newCon.append(i)
    return newCon


DG_consensus = DeGap(consensus)

GapCount = 0
for i in DG_consensus:
    if i[0] == ('-',):
        GapCount +=1


           
for i in range(len(DG_consensus)):
    print ''.join(DG_consensus[i][0]) +'*' +str(round(DG_consensus[i][1],2)), 
print '***'
print len(consensus), len(DG_consensus), GapCount

#print consensus

#(kmer,Start,Stop)
#sort into potential fwd vs revs based on 3' and 5' characteristics

def GenKmers(consensus,MinLen=18,MaxLen=22):
    """
    Generates all possible kmers between MinLen and MaxLen
    These are then scored for primer 'goodnesss' according to the function
    ScoreKmer
    """
    lengths = [i+MinLen for i in range(MaxLen+1-MinLen)]
    kmers = []
    for length in lengths:
        for i in range(len(consensus)+1 - length):
            kmer = consensus[i:i+length]
            kmers.append((i,kmer))
    return kmers
        
kmers = GenKmers(DG_consensus) 

#make dictionary of iupac bases where keys are tuples

#print x           

def ScoreKmer(kmer, n = 1.0):#increasing n will decrease the -ve impact of degeneracy on score
    score = 0
    gaps = 0
   # print kmer[0], kmer[1]
    #print
    #print
    for base in kmer[1]:
     #   print base, base[0][0]
        if base[0][0] == '-':
      #      print 'gappy'
            score -= (1-base[1])
            gaps += 1
        else:
       #     print 'not gappy'
            fraction = base[1]
            degen = len(base[0])
            FractionGaps = base[2]
        #    print fraction, FractionGaps, degen
            score += (n - degen + fraction - FractionGaps)
    
    if len(kmer[1]) - gaps > 17: 
        return (score,kmer)
        
        
def ScoreAndRankKmers(kmers):
    ScoredKmers = []
    for kmer in kmers:
        x = ScoreKmer(kmer)
        if x:
            ScoredKmers.append(x)
    return ScoredKmers
        
ToRank = ScoreAndRankKmers(kmers)
ToRank.sort(reverse=True)

def ParseAndPair(kmers,MinDist = 150, MaxDist = 400):
    FwKmers = []
    RevKmers = []
    pairs = []
    for kmer in kmers:
        seq = kmer[1][1]
        #print seq
        #print seq[0][0], seq[-1][0]
        if len(seq[0][0]) == 1:
            #print "Found Rev"
            RevKmers.append(kmer)
        if len(seq[-1][0]) == 1:
            #print "Found Fw"
            FwKmers.append((kmer))
    for FwKmer in FwKmers:
        FwScore = FwKmer[0]
        FwStart = FwKmer[1][0]
        #print FwScore, FwStart
        for RevKmer in RevKmers:
           RevScore = RevKmer[0]
           RevStart = RevKmer[1][0] +-len(RevKmer[1][1])
           ProductLength = RevStart - FwStart
           if MinDist <= ProductLength <= MaxDist:
               pairs.append((FwScore+RevScore,(FwKmer,RevKmer)))
               #print FwStart, RevStart, FwScore+RevScore 
    pairs.sort(reverse = True)
    return (pairs)
        
            
z = ParseAndPair(ToRank)


#need to penalize degeneracy more heavily   


#How are there still gaps in this!!

# consensus[i][0] = fraction that IUPAC base represents
#len(consensus[i][1]) = degeneracy, len(consensus[i][1]) = the tuple of bases at that position
#consensus[i][2] = fraction gaps
#treat gap positions differently. 
#The closer to 100% gaps, the closer to no -ve impact
#Don't subtract gap from overall primer len
#score by consensus[i][0]/len(consensus[i][1]) or consensus[i][0]/len(consensus[i][1])**n 
#could make n lower or higher depending on how much you value lack of degeneracy
#do we even need to make varible length? Can we just set default to 22?
#also factor % '-'s. This is consensus[i][2]
#also also, 3' nt must not be degenerate and must pass a stringency cutoff for fraction
#apply extra weight to additional 3' end nts?
#divide overall score by degeneracy to further penalize?
#Postition score = n*consensus[i][0] - (len(consensus[i][1])+consensus[i][2])
# or  n + consensus[i][0] -   (len(consensus[i][1])+consensus[i][2])  
#if position = gap, ignore and score -= (1-consensus[i][2]) --> tends to zero as more gaps
#or + n-len, set n to 2 so two fold degeneracy is neutral, one fold positve, 3 and four are -ve
#multiply the whole lot by 1-numgaps
#score = n - len + fraction - gaps, default n = 2      
#write a function that checks the number of sequences to which a primer successully anneals     
