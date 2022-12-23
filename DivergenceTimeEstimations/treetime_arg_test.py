#!/usr/bin/python

import sys, getopt, os
import random
import numpy as np
import pandas as pd
import json
from treetime.utils import parse_dates
from treetime.arg import parse_arg, setup_arg
from treetime import TreeTime
from Bio import Phylo
import json

def main(argv):
    """
    Read in input and output directory name as well as round number (will be added to output name)
    """
    try:
        opts, args = getopt.getopt(argv,"hi:o:r:k:n:",["inputdir=","outputdir=","key=","no_simulations="])
    except getopt.GetoptError:
        print("treetime_arg_test.py -i <inputdir> -o <outputdir> -k <key_folder_name> -n <no_simulations>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("treetime_arg_test.py -i <inputdir> -o <outputdir> -k <key_folder_name> -n <no_simulations>")
            sys.exit()
        elif opt=="-i":
            inputdir = arg
        elif opt=="-o":
            outputdir = arg
        elif opt=="-k":
            key = int(arg)
        elif opt=="-n":
            no_simulations = int(arg)
    return inputdir, outputdir, key, no_simulations

def find_which_tree_subset(MCC_file):
    f = open(MCC_file)
    data = json.load(f)
    trees = set()
    for key in data["MCC_dict"]:
        for s in data["MCC_dict"][key]["trees"]:
            trees.add(s)
    return list(trees)

def run_treetime_with_arg(tree_name, trees, alignments, MCC_file, dates, method_anc = "probabilistic", branch_length_mode="marginal"):
    """
    Read in trees, alignments and MCC files using the setup_arg function and then run treetime on the resulting tree,
    this must be done twice as setup_arg will always use the first tree as the main tree and add information from the second tree 
    """
    
    arg_params = parse_arg(trees, alignments, MCC_file)

    tt = setup_arg(tree_name, arg_params['trees_dict'], arg_params['alignment'], dates, arg_params['MCCs_dict'], arg_params['masks_dict'],
                    gtr='JC69', verbose=False, fixed_clock_rate = 1, alphabet="nuc_nogap")
    
    tt.optimize_tree(infer_gtr=False, max_iter=10, method_anc = method_anc, branch_length_mode=branch_length_mode,
            resolve_polytomies=False, prune_short=False)

    return tt.tree

def run_inference_on_1_folder(output_dir, pair=False):
    dates = parse_dates(output_dir+'/metadata.csv')
    dist_dict_list = []
    if pair:
        k_vals = range(2,3)
        with open(output_dir+'/shared_branches.json') as json_file:
            shared_dict = json.load(json_file)
    else:
        k_vals = range(2,9)
        shared_dict = None
    for k in k_vals:
        MCC_file = output_dir+ '/' + str(k) + '/MCCs.json'
        tree_subset = find_which_tree_subset(MCC_file)
        if pair:
            sample_tree_name = tree_subset[0]
        else:
            sample_tree_name = random.sample(tree_subset, 1)[0]
        tree_nwk = [output_dir+"/" + t+".nwk" for t in tree_subset]
        aln = [output_dir+"/" +t+".fasta" for t in tree_subset]
        tt_arg_tree = run_treetime_with_arg(sample_tree_name, tree_nwk, aln, MCC_file, dates, method_anc = "probabilistic")
        tt_old = TreeTime(dates=dates, tree=output_dir+"/" + sample_tree_name+".nwk", aln=output_dir+"/" +sample_tree_name+".fasta",
            gtr='JC69', verbose=False, keep_node_order=True, compress=False, fixed_clock_rate = 1, alphabet="nuc_nogap")
        tt_old.optimize_tree(infer_gtr=False, max_iter=10, method_anc = "probabilistic", branch_length_mode="marginal",
            resolve_polytomies=False, prune_short=False)
        tt_old_tree = tt_old.tree
        true_tree = Phylo.read(output_dir+"/" +"true_segment_"+sample_tree_name+".nwk", "newick")
        dist_dict = compute_average_dist(true_tree, tt_old_tree, tt_arg_tree, pair, shared_dict)
        dist_dict_list.append(dist_dict)
    return dist_dict_list

def compute_average_dist_dict(true_tree, tt_old, tt_arg, val=None, shared_dict=None):
    branch_lengths_true = {}
    branch_lengths_old = {}
    branch_lengths_arg = {}

    for node in true_tree.find_clades():
        if node != true_tree.root and node not in true_tree.root.clades:
            branch_lengths_true[node.name] = node.branch_length
    for node in tt_old.find_clades():
        if node != tt_old.root and node not in tt_old.root.clades:
            branch_lengths_old[node.name] = node.branch_length
    for node in tt_arg.find_clades():
        if node != tt_arg.root and node not in tt_arg.root.clades:
            branch_lengths_arg[node.name] = node.branch_length
    node_names = branch_lengths_true.keys()
    if shared_dict is not None:
        node_names = [n for n in node_names if shared_dict.get(n, 0) == val]
    if len(node_names) == 0:
        return {"median_old": np.nan, "mean_old": np.nan, "var_old": np.nan, "median_arg": np.nan, "mean_arg": np.nan, "var_arg": np.nan}

    epsilon = 1e-12
    old_dist_vector = [(branch_lengths_old.get(n,0) - branch_lengths_true.get(n,0))/np.sqrt(branch_lengths_true.get(n,0)+ epsilon) for n in node_names]
    arg_dist_vector = [(branch_lengths_arg.get(n,0) - branch_lengths_true.get(n,0))/np.sqrt(branch_lengths_true.get(n,0)+ epsilon) for n in node_names]
    
    old_dist_vector_full_norm = [(branch_lengths_old.get(n,0) - branch_lengths_true.get(n,0))/(branch_lengths_true.get(n,0)+ epsilon) for n in node_names]
    arg_dist_vector_full_norm = [(branch_lengths_arg.get(n,0) - branch_lengths_true.get(n,0))/(branch_lengths_true.get(n,0)+ epsilon) for n in node_names]
    
    return {"median_old": np.median(old_dist_vector_full_norm), "mean_old": np.mean(old_dist_vector_full_norm), "var_old": np.var(old_dist_vector), "median_arg": np.median(arg_dist_vector_full_norm), "mean_arg": np.mean(arg_dist_vector_full_norm), "var_arg": np.var(arg_dist_vector)}


def compute_average_dist(true_tree, tt_old, tt_arg, pair=False, shared_dict=None):
    if not pair:
        return compute_average_dist_dict(true_tree, tt_old, tt_arg, val=None, shared_dict=None)
    else:
        dist_dict = {}
        for val in ["TP", "TN", "FP", "FN"]:
            dist_dict[val] = compute_average_dist_dict(true_tree, tt_old, tt_arg, val=val, shared_dict=shared_dict)
        return dist_dict

def append_val(dict_, key, vector, name):
    if key in dict_:
        dict_[key].append(vector[name])
    else:
        dict_[key] = [vector[name]]

def run_inference_on_all_rec_folders(out_dir, key, no_sim):
    all_rec_folders = [out_dir+"/sim_"+str(key)+"_"+str(n) for n in range(no_sim)]
    folder = out_dir+"/sim_"+str(key)+"_0"
    if os.path.isdir(folder + "/2") and not(os.path.isdir(folder + "/8")):
        pair = True
        k_vals = range(2,3)
    else:
        pair = False
        k_vals = range(2,9)
    median_old = {}
    mean_old = {}
    var_old = {}
    median_arg = {}
    mean_arg = {}
    var_arg = {}
    if pair:
        dict_ = {}
        for val in ["TP", "TN", "FP", "FN"]:
            median_old[val] = {}
            mean_old[val] = {}
            var_old[val] = {}
            median_arg[val] = {}
            mean_arg[val] = {}
            var_arg[val] = {}
        for folder in all_rec_folders:
            diff_list = run_inference_on_1_folder(folder, pair)
            for val in ["TP", "TN", "FP", "FN"]:
                for k in k_vals:
                    vector = diff_list[k-2][val]
                    append_val(median_old[val], k, vector, "median_old")
                    append_val(mean_old[val], k, vector, "mean_old")
                    append_val(var_old[val], k, vector, "var_old")
                    append_val(median_arg[val], k, vector, "median_arg")
                    append_val(mean_arg[val], k, vector, "mean_arg")
                    append_val(var_arg[val], k, vector, "var_arg")
        for val in ["TP", "TN", "FP", "FN"]:
            dict_[val] = {"median_old": median_old[val], "mean_old": mean_old[val], "var_old": var_old[val], "median_arg": median_arg[val], "mean_arg": mean_arg[val], "var_arg": var_arg[val]}
        return dict_
    else:
        for folder in all_rec_folders:
            diff_list = run_inference_on_1_folder(folder, pair)
            for k in k_vals:
                vector = diff_list[k-2]
                append_val(median_old, k, vector, "median_old")
                append_val(mean_old, k, vector, "mean_old")
                append_val(var_old, k, vector, "var_old")
                append_val(median_arg, k, vector, "median_arg")
                append_val(mean_arg, k, vector, "mean_arg")
                append_val(var_arg, k, vector, "var_arg")
        return {"median_old": median_old, "mean_old": mean_old, "var_old": var_old, "median_arg": median_arg, "mean_arg": mean_arg, "var_arg": var_arg}

def write_results(data_full, results_dir, key):
    if len(data_full.keys()) == 4:
        for k in range(2,3):
            data = {}
            for val in ["TP", "TN", "FP", "FN"]:
                for name in ["median_old", "mean_old", "var_old", "median_arg", "mean_arg", "var_arg"]:
                    data[val+"_"+name] = data_full[val][name][k]
            df = pd.DataFrame(data)
            if not os.path.isdir(results_dir):
                os.mkdir(results_dir)
            df.to_csv(results_dir+"/tt_summary_"+ str(key)+"_"+str(k)+".csv")
    else:
        for k in range(2,9):
            data = {}
            for name in ["median_old", "mean_old", "var_old", "median_arg", "mean_arg", "var_arg"]:
                data[name] = data_full[name][k]
            df = pd.DataFrame(data)
            if not os.path.isdir(results_dir):
                os.mkdir(results_dir)
            df.to_csv(results_dir+"/tt_summary_"+ str(key)+"_"+str(k)+".csv")
    return None


if __name__ == '__main__':

    inputdir, outputdir, key, no_simulations = main(sys.argv[1:])
    # inputdir = "DivergenceTimeEstimations/sim_flu_0.3_true"
    # no_simulations = 1
    # key = 4
    # outputdir = "results_test"
    dict_ = run_inference_on_all_rec_folders(inputdir, key, no_simulations)
    write_results(dict_, outputdir, key)