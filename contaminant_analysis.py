import os
import math
import sys
import statistics
import pandas as pd
import numpy as np
import seaborn as sns
import pysam
import json
import matplotlib.pyplot as plt

def calculate_positional_depths(bam):
    """
    Return actual sequence matrix.
    """
    print("Calculating positional depths")
    query_names = []
    samfile = pysam.AlignmentFile(bam, "rb")
    #print(thing.query_alignment_sequence)
    #print(thing.get_reference_sequence())
    #print(thing.query_alignment_qualities)
    #print(thing.get_reference_positions())
       
    unique_headers = []
    position_dict = {}    
    for count,thing in enumerate(samfile): 
        if count % 100000 == 0:
            print(count)
        #loop over the aligned sequence parts
        for (pos, letter, ref, qual,test) in zip(thing.get_reference_positions(), thing.query_alignment_sequence, thing.get_reference_sequence(), thing.query_alignment_qualities, thing.get_aligned_pairs(matches_only=True)):        
            if qual < 20:
                continue
            #sometimes we have misalignments
            #if pos == 6057:    
            #    print(pos, letter, ref, qual, test)
            
            #if we've seen this position before
            if pos-1 in position_dict:
                #check if we've seen this nuc before
                if letter in position_dict[pos-1]["allele"]:
                    position_dict[pos-1]['allele'][letter] += 1
                    position_dict[pos-1]['total_depth'] += 1 
                #else we will add the letter            
                else:
                    position_dict[pos-1]['allele'][letter] = 1
                    position_dict[pos-1]['total_depth'] += 1
                if ref in position_dict[pos-1]['ref']:
                    position_dict[pos-1]['ref'][ref] += 1
                else:
                    position_dict[pos-1]['ref'][ref] = 1
                    
            #else we'll add this position
            else:
                position_dict[pos-1] = {}
                position_dict[pos-1]['ref'] = {ref:1}
                position_dict[pos-1]['allele'] = {letter: 1}
                position_dict[pos-1]['total_depth'] = 1
            
        query_names.append(thing.query_name)
    return(position_dict)

def remove_low_depth_positions(position_dict):
    temp_dict = {}
    for key,value in position_dict.items():
        if value['total_depth'] < 10:
            continue
        else:
            temp_dict[key] = value
    return(temp_dict)

def calculate_read_probability(position_depths, bam, usable_pos):
    """
    Calculates the probabilites per read and then creates
    a curve to represent the distribution of read probabilities.
    
    Parameters
    ----------
    position_depths : dict
        A dictionary representing position based allele depths.
    
    bam : str
        The path to the bam file we're looking at.
    """
    samfile = pysam.AlignmentFile(bam, "rb")
    all_read_probs = []
    all_read_arrays = []
    all_nt_arrays = []
    all_read_pos = []
    
    for count, thing in enumerate(samfile):
        if count % 100000 == 0:
            print(count)
            #if count >= 100000:
            #    break       
        
        temp_read_prob = 1
        temp_read_arrays = []
        #print(dir(thing))
        #print(thing.query_alignment_sequence)
        #print(thing.get_reference_sequence())
        #print(thing.query_alignment_qualities)
        #print(thing.get_reference_positions())
        #sys.exit(0)       
        for (pos, nuc, ref, qual) in zip(thing.get_reference_positions(), thing.query_alignment_sequence, thing.get_reference_sequence(), thing.query_alignment_qualities): 
            pos = str(pos-1)
            try:
                pos_prob = position_depths[pos] 
            except:
                print("key not found")
                continue
            pos_ref = max(pos_prob["ref"], key=pos_prob["ref"].get)
            #if a nuc quality is low, we ignore it and don't even bother to calculate the probability
            if qual < 20:
                continue

            #if the query nuc is same as the actual nuc we pass
            if nuc == pos_ref:
                continue
            if nuc.lower() == pos_ref.lower():
                continue
            if nuc.lower() == 'n':
                continue
            
            total_depth = pos_prob['total_depth']
            nuc_depth = pos_prob['allele'][nuc]
            probability = nuc_depth/total_depth
            temp_read_arrays.append(probability)
            temp_read_prob += math.log(probability)
            
            if int(pos) == 5984:
                #print("false", probability, nuc, pos_prob["ref"], ref)
                pass
            if int(pos) == 6056:
                #print(ref)
                #print("true", probability, nuc, ref, pos_ref)
                pass
        all_nt_arrays.append(thing.query_sequence)
        all_read_arrays.append(temp_read_arrays)
        all_read_pos.append(thing.positions)
        all_read_probs.append(temp_read_prob/len(thing.positions))
    
    
    #print(position_depths.keys())
    print("Max: ", max(all_read_probs))
    print("Min: ", min(all_read_probs))
    q1 = np.percentile(all_read_probs, 25)
    print("Q1: ", q1)
    median = statistics.median(all_read_probs)
    print("Median: ", median)
    mean = statistics.mean(all_read_probs)
    print("Mean: ", mean)
    ten_lower = np.percentile(all_read_probs, 5)
    print("Lower 5: ", ten_lower)

    #pca on the probabilities
    #all_read_arrays
    """
    from sklearn.decomposition import PCA
    pca = PCA(n_components=2)
    np_all_read_arrays = np.array(all_read_arrays)
    print(np_all_read_arrays.shape)
    principalComponents = pca.fit_transform(np_all_read_arrays)
    df_pca = pd.DataFrame(principalComponents)
    print(df_pca)
    """
    
    #sns.boxplot(all_read_probs)
    #plt.show()
    flagged_nucs = []
    flagged_read = []
    flagged_pos = []
    for count, (prob,nt,pos) in enumerate(zip(all_read_probs, all_nt_arrays, all_read_pos)):
        #the condition for determining which reads are flagged
        if prob < ten_lower:
            flagged_read.append(count)   
            flagged_nucs.append(nt)
            flagged_pos.append(pos)
    return(flagged_read, flagged_pos, flagged_nucs, all_read_probs)

def condense_flagged_read_NT(flagged_read_NT, flagged_read_pos):
    """
    Make a summary dict by position of all flagged read NTs.
    
    Parameters
    ----------
    flagged_read_NT : list
    flaged_read_pos : list    

    Returns
    -------
    flagged_nucs : dict
    
    """
    flagged_nucs = {}
    for positions,sequences in zip(flagged_read_pos, flagged_read_NT):
        for pos, letter in zip(positions, sequences):
            if pos in flagged_nucs:
                #check if we've seen this nuc before
                if letter in flagged_nucs[pos]["allele"]:
                    flagged_nucs[pos]['allele'][letter] += 1
                    flagged_nucs[pos]['total_depth'] += 1 
                
                #else we will add the letter            
                else:
                    flagged_nucs[pos]['allele'][letter] = 1
                    flagged_nucs[pos]['total_depth'] += 1
                    
            #else we'll add this position
            else:
                flagged_nucs[pos] = {}
                flagged_nucs[pos]['allele'] = {letter: 1}
                flagged_nucs[pos]['total_depth'] = 1
    
    return(flagged_nucs)

def snv(position_depths, flagged_positions, bam, usable_pos, ground_truth):
    """
    Determines whether or not a nuc is: 
    (1) Matches the ref, ignore it!
    (2) Low quality, ignore it!
    (3) SNV
    (4) Potential contamination

    Parameters
    ----------
    position_depths : dict
        Depths at every position.
    flagged_positions : dict 

    Returns
    -------
    
    """
    #make ground truth dict
    zip_iterator = zip(usable_pos, ground_truth) 
    gt_dict = dict(zip_iterator)

    #for normal reads
    frequency_cutoff = 0.03
    
    #for flagged reads
    frequency_cutoff_harsh = 0.50
    
    samfile = pysam.AlignmentFile(bam, "rb")
   
    fake_consensus_threshold = 0.03

    #this is where we store LOW DEPTH, SNV, CONTAM
    choice = []

    #now we want to iterate every position, and say what we think is the likely SNV vs. contaminating NT, vs. matches ref
    for key, value in position_depths.items():
        if int(key) == 6055:
            #print(value) 
            pass
            #print(value['allele'])
            #sys.exit(0)
        if int(key)+2 in gt_dict:
            pass
        else:
            continue
        total_depth = value["total_depth"]
        possible_alleles = value["allele"]
        reference_allele = max(value["ref"], key=value["ref"].get)
        if total_depth < 10:
            continue
        
        #print(possible_alleles)
        
        #let's try and call consensus
        consensus_alleles = {k: v / total_depth for k, v in possible_alleles.items()}
        #consensus_alleles = {k:v for (k,v) in consensus_alleles.items() if v > fake_consensus_threshold}       
        
        #iterate all indels that make our cuttoff
        for ckey, cvalue in consensus_alleles.items():
            if int(key)+2 == 6056:
                #print(consensus_alleles)
                #print(ckey, cvalue)
                #sys.exit(0)
                pass
            #matches the reference
            if ckey == reference_allele:
                continue
            if ckey.lower() == reference_allele.lower():
                continue
            if ckey == "N":
                continue
            if ckey.lower() == "n":
                continue
            
            #we think this is a contaminant
            if cvalue <= 0.03:
                if gt_dict[int(key)+2] == "TRUE":
                    pass
                    #print("1: false negative: ", ckey, cvalue, consensus_alleles, reference_allele, int(key)+2)
                temp_dict = {"pos":int(key)+2, "ref": reference_allele, "allele":ckey, \
                    "status": "low frequency", "freq":cvalue}
                choice.append(temp_dict)
                
            #it's not consensus but it's not low enough to be flagged as contamination
            else:
                #we go look for it's occurence in the pos in contaminated samples
                if int(key) in flagged_positions:
                    #this is how often it occurs in the contamination
                    temp_allele = flagged_positions[int(key)]
                    #normalize it
                    temp_allele = {k: v / temp_allele['total_depth'] for k, v in temp_allele['allele'].items()}
                    if int(key)+2 == 10376:
                        print(temp_allele)
                        print(ckey, cvalue)
                        print(temp_allele[ckey])
           
                    #the occurence in flagged reads as a whole is is under 30%
                    if ckey not in temp_allele: 
                        print("MISTAKE??")
                        continue
                    
                    #we have contamination
                    if (temp_allele[ckey] < 0.30):
                        if int(key)+2 == 10376:
                            print("HERE", ckey, cvalue)
                        if gt_dict[int(key)+2] == "FALSE":
                            pass
                            #print("false positive: ", gt_dict[int(key)+2], reference_allele, temp_allele)
                        temp_dict = {"pos":int(key)+2, "ref": reference_allele, "allele":ckey, \
                            "status": "contaminant", "freq":cvalue}
                        choice.append(temp_dict)
                    
                    #we have an SNV? 
                    else:               
                        if gt_dict[int(key)+2] == "TRUE":
                            pass
                            #print("false negative: ",gt_dict[int(key)+2], reference_allele, temp_allele, cvalue, ckey)          
                        temp_dict = {"pos":int(key)+2, "ref":reference_allele,"allele":ckey, \
                            "status": "normal base", "freq":cvalue}
                        choice.append(temp_dict)
                else:
                    print("No occurences in a flagged read.")
               
    return(choice)
     
def manipulate_spreadsheet(filename, ground_truth_filename):
    df = pd.read_csv(filename)
    ground_truth_df = pd.read_csv(ground_truth_filename)
    print(df.columns)
    print(ground_truth_df)
    
    merge_df = df.merge(ground_truth_df, right_on="Position", left_on="pos", how='left')
    merge_df = merge_df[merge_df['True/False/Remove'] != "Remove"]
    
    false_pos = 0
    false_neg = 0
    true_pos = 0
    true_neg = 0
    new_row = []
    for index,row in merge_df.iterrows():
        prediction = row["status"]
        gt = row["True/False/Remove"]
        freq = row['freq']
        
        if (prediction == "normal base" or prediction == "low frequency") and gt == "FALSE":
            true_neg += 1
            new_row.append("true neg")
        elif prediction == "contaminant" and gt == "TRUE":
            #condition relevant to the dataset 
            if freq > 0.03 and freq < 0.39:
                true_pos += 1
                new_row.append("true pos")
            #it's the wrong base, false pos
            else:
                new_row.append("false pos")
                false_pos += 1
        elif (prediction == "normal base" or prediction == "low frequency") and gt == "TRUE":
            #this is our known contam range
            if freq < 0.03 or freq > 0.39:
                new_row.append("true neg")
                true_neg += 1
            else:
                new_row.append("false neg")
                false_neg += 1
        elif prediction == "contaminant" and gt == "FALSE":
            new_row.append("false pos")
            false_pos += 1

    print("True Pos: ", true_pos)
    print("True Neg: ", true_neg)
    print("False Pos: ", false_pos)
    print("False Neg: ", false_neg)
    merge_df['TFP'] = new_row
    merge_df.to_csv("test.csv")

def look_for_viable_pos(ground_truth_filename):
    ground_truth_df = pd.read_csv(ground_truth_filename)
    df = ground_truth_df[ground_truth_df['True/False/Remove'] != "Remove"]
    return(df["Position"].tolist(), df['True/False/Remove'].tolist())
 
def main():
    """
	Proof of concept code for identifying contaminants vs SNV's in a bam file
	"""
     
    filename = "./test.csv"
    ground_truth_filename = "./ZIKV-intrahost_True-False-Remove.csv"
    """
    manipulate_spreadsheet(filename, ground_truth_filename)
    sys.exit(0)
    """
    ground_truth_filename = "./ZIKV-intrahost_True-False-Remove.csv"
    usable_pos,ground_truth = look_for_viable_pos(ground_truth_filename)
    
    bam = "./bam/primalseq_mix10_repA.sorted.bam"
    #bam = "./bam/spike_in_metagenomics_mix10_repA.sorted.bam"
 
    """ 
    position_depths = calculate_positional_depths(bam) 
    with open('pos_depths.json','w') as jfile:
        json.dump(position_depths, jfile)
    sys.exit(0) 
    """
    
    with open('pos_depths.json', 'r') as jfile:
        position_depths = json.load(jfile)
    """ 
    #using probabilties flag reads as potential contaminants
    flagged_reads, flagged_read_pos, flagged_read_NT, all_read_probs= calculate_read_probability(position_depths, bam, usable_pos)
    dump_dict = {"flagged_reads": flagged_reads, "flagged_read_pos":flagged_read_pos, "flagged_read_NT":flagged_read_NT, \
    "all_read_probs": all_read_probs}      
    
    with open("flagged_depths.json", "w") as jfile:
        json.dump(dump_dict, jfile)
    sys.exit(0)
    """
    with open("flagged_depths.json", "r") as jfile:
        dump_dict = json.load(jfile)
    
    flagged_reads = dump_dict["flagged_reads"]
    flagged_read_NT = dump_dict["flagged_read_NT"]
    flagged_read_pos = dump_dict["flagged_read_pos"]
    all_read_probs = dump_dict["all_read_probs"]
   
    #visualize_distribution(all_read_probs, 
    #sns.distplot(all_read_probs)
    #plt.show()
    print("Flagged Reads: ", len(flagged_reads))
     
    flagged_positions = condense_flagged_read_NT(flagged_read_NT, flagged_read_pos)   
    #go through and try and use frequency to determine SNV, using flagged reads
    temp_list = snv(position_depths, flagged_positions, bam, usable_pos, ground_truth)
    #print(temp_list)
    df = pd.DataFrame(temp_list)
    df.to_csv("test.csv")
    #print(df)
    manipulate_spreadsheet(filename, ground_truth_filename)
    

if __name__ == "__main__":
    main()
