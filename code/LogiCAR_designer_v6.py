#! usr/bin/env python

"""
Author: Tiangen Chang
date: 2023-08-03
E-mail: changtiangen@gmail.com
Description: Design triplet-gene CAR-T using genetic algorithm
    This script is used to find the (near-)optimal triplet-gene combination for CAR-T design
    usage: python Triplet_CART_designer_v5.py [args]
    Input (files):
        -i1: the TME scRNAseq single cell gene expression matrix (e.g. HNSC_GSE103322_matrix.txt)
        -i2: the normal tissue scRNAseq single cell gene expression matrix (e.g. Tabula_Sapiens_subsampled_matrix.csv)
        -g: the candidate CAR target genes encoding surface proteins (e.g. intersection_and_clinical.txt)
    Input (parameters):
        -N: output topN triplet-gene combinations (e.g. 200)
        -c: logic-gate combination type (e.g. AA)
        -r: random seed
        -ps: population size of each generation
        -sg: stop after certain generations without new solution
        -mp: genome mutation probability of patients
        -cp: genome crossover probability between two patients
    Output:
        -o: the identified (near-)optimal triplet-gene combination for CAR-T (e.g. ./output/HNSC_GSE103322_optimal_triplets.txt)
# improvement:
1. extend from 3-gene combination to N-gene combination
2. move pred_power_calculator() and evalLogicGate*() functions to utils_v2.py for better readability
3. set max_generation = 400 to force quit evolution after 400 generations

How to run this program? e.g.
# 2-gene
python LogiCAR_designer_v6.py -i1 ../../02.Input/15BRCAdatasets_tumorCells_LogiCAR_wgene_and_clinical.txt -i2 ../../02.Input/HPA/version24_20241022/rna/HPA_LogiCAR_cell_surface_genes_wgene_and_clinical.txt -g ../../02.Input/Constants/cell_surface_genes_wgene_and_clinical_uniqueStandardNames.txt -o ./output/BRCA_AllDatasets_all_optimal_doublets_A.txt -N 200 -c A -mp 0.2 -cp 0.5 -r 1
# 3-gene
python LogiCAR_designer_v6.py -i1 ../../02.Input/15BRCAdatasets_tumorCells_LogiCAR_wgene_and_clinical.txt -i2 ../../02.Input/HPA/version24_20241022/rna/HPA_LogiCAR_cell_surface_genes_wgene_and_clinical.txt -g ../../02.Input/Constants/cell_surface_genes_wgene_and_clinical_uniqueStandardNames.txt -o ./output/BRCA_AllDatasets_all_optimal_triplets_AA.txt -N 200 -c AA -mp 0.2 -cp 0.5 -r 1
# 5-gene
python LogiCAR_designer_v6.py -i1 ../../02.Input/15BRCAdatasets_tumorCells_LogiCAR_wgene_and_clinical.txt -i2 ../../02.Input/HPA/version24_20241022/rna/HPA_LogiCAR_cell_surface_genes_wgene_and_clinical.txt -g ../../02.Input/Constants/cell_surface_genes_wgene_and_clinical_uniqueStandardNames.txt -o ./output/BRCA_AllDatasets_all_optimal_quintuplets_OA_O_O.txt -N 200 -c OA_O_O -mp 0.2 -cp 0.5 -r 1


# 2-gene logic combination types (2 unique combinations):
A: G1 & G2
O: G1 | G2
# 3-gene gate code (4 unique combinations):
AA: (G1 & G2) & G3
AO: (G1 & G2) | G3
OA: (G1 | G2) & G3
OO: (G1 | G2) | G3
# 4-gene logic combination types (12 unique combinations):
AA_A: ((G1 & G2) & G3) & G4
AO_A: ((G1 & G2) | G3) & G4
OA_A: ((G1 | G2) & G3) & G4
OO_A: ((G1 | G2) | G3) & G4
AA_O: ((G1 & G2) & G3) | G4
AO_O: ((G1 & G2) | G3) | G4
OA_O: ((G1 | G2) & G3) | G4
OO_O: ((G1 | G2) | G3) | G4
A_A_O: (G1 & G2) & (G3 | G4)
A_O_A: (G1 & G2) | (G3 & G4)
O_A_O: (G1 | G2) & (G3 | G4)
O_O_A: (G1 | G2) | (G3 & G4)

# 5-gene logic combination types (40 unique combinations):
AA_A_A: (((G1 & G2) & G3) & G4) & G5
AO_A_A: (((G1 & G2) | G3) & G4) & G5
OA_A_A: (((G1 | G2) & G3) & G4) & G5
OO_A_A: (((G1 | G2) | G3) & G4) & G5
AA_O_A: (((G1 & G2) & G3) | G4) & G5
AO_O_A: (((G1 & G2) | G3) | G4) & G5
OA_O_A: (((G1 | G2) & G3) | G4) & G5
OO_O_A: (((G1 | G2) | G3) | G4) & G5
A_A_O_A: ((G1 & G2) & (G3 | G4)) & G5
A_O_A_A: ((G1 & G2) | (G3 & G4)) & G5
O_A_O_A: ((G1 | G2) & (G3 | G4)) & G5
O_O_A_A: ((G1 | G2) | (G3 & G4)) & G5

AA_A_O: (((G1 & G2) & G3) & G4) | G5
AO_A_O: (((G1 & G2) | G3) & G4) | G5
OA_A_O: (((G1 | G2) & G3) & G4) | G5
OO_A_O: (((G1 | G2) | G3) & G4) | G5
AA_O_O: (((G1 & G2) & G3) | G4) | G5
AO_O_O: (((G1 & G2) | G3) | G4) | G5
OA_O_O: (((G1 | G2) & G3) | G4) | G5
OO_O_O: (((G1 | G2) | G3) | G4) | G5
A_A_O_O: ((G1 & G2) & (G3 | G4)) | G5
A_O_A_O: ((G1 & G2) | (G3 & G4)) | G5
O_A_O_O: ((G1 | G2) & (G3 | G4)) | G5
O_O_A_O: ((G1 | G2) | (G3 & G4)) | G5

AA_a_A: ((G1 & G2) & G3) & (G4 & G5)
AO_a_A: ((G1 & G2) | G3) & (G4 & G5)
OA_a_A: ((G1 | G2) & G3) & (G4 & G5)
OO_a_A: ((G1 | G2) | G3) & (G4 & G5)

AA_o_A: ((G1 & G2) & G3) | (G4 & G5)
AO_o_A: ((G1 & G2) | G3) | (G4 & G5)
OA_o_A: ((G1 | G2) & G3) | (G4 & G5)
OO_o_A: ((G1 | G2) | G3) | (G4 & G5)

AA_a_O: ((G1 & G2) & G3) & (G4 | G5)
AO_a_O: ((G1 & G2) | G3) & (G4 | G5)
OA_a_O: ((G1 | G2) & G3) & (G4 | G5)
OO_a_O: ((G1 | G2) | G3) & (G4 | G5)

AA_o_O: ((G1 & G2) & G3) | (G4 | G5)
AO_o_O: ((G1 & G2) | G3) | (G4 | G5)
OA_o_O: ((G1 | G2) & G3) | (G4 | G5)
OO_o_O: ((G1 | G2) | G3) | (G4 | G5)

"""

import argparse
import random
from deap import creator, base, tools, algorithms
from sklearn.metrics import confusion_matrix
import numpy as np
import pandas as pd
import time

def createHelp():
    """
    Create the command line interface of the program.
    """
    epilog_string="Any bug is welcome reported to changtiangen@gmail.com"
    description_string='The program is going to find the near-optimal triplet-gene combination for CAR-T design'
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i1', '--input1', dest='fnIn1', default='../../02.Data/HNSC_GSE103322_matrix_train.txt', type=str,help='input raw TME scRNAseq gene expression file')
    parser.add_argument('-i2', '--input2', dest='fnIn2', default='../../02.Data/Tabula_Sapiens_subsampled_matrix.csv', type=str,help='input raw normal tissue scRNAseq gene expression file')
    parser.add_argument('-g', '--gene', dest='fnGene', default='../intersection_and_clinical.txt', type=str,
                        help='input file for candidate CAR target genes encoding surface proteins')
    parser.add_argument('-o', '--output', dest='fnOut', default='./output/HNSC_GSE103322_optimal_triplets.txt', type=str,help='output file')
    parser.add_argument('-N', '--topN', dest='topN', default=200, type=int,help='output topN triplet-gene combinations')
    parser.add_argument('-c', '--comb_type', dest='ct', default='AA', type=str, help='combination type')
    parser.add_argument('-r', '--random_seed', dest='rs', default=0, type=int, help='random seed')
    parser.add_argument('-ps', '--population_size', dest='ps', default=1000, type=int, help='population size of each generation')
    parser.add_argument('-sg', '--stop_after_generation', dest='sg', default=200, type=int, help='stop after N generations without new solution')
    parser.add_argument('-mp', '--mutate_prob', dest='mp', default=0.2, type=float, help='genome mutation probability of patients')
    parser.add_argument('-cp', '--crossover_prob', dest='cp', default=0.5, type=float, help='genome crossover probability between two patients')
    op=parser.parse_args()
    return op


#exec(open('./utils_v2.py').read())
exec(open('./utils_v3.py').read())

def setup_toolbox(df, combination_type):
    gene_num = len(combination_type.replace("_", "")) + 1
    # Define the fitness and individual
    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax)
    toolbox = base.Toolbox()
    # Redefine the individual and population creation operations to generate gene_num unique genes.
    toolbox.register("indices", random.sample, range(len(df.columns) - 3), gene_num)
    toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.indices)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    command = 'toolbox.register("evaluate", evalLogicGate' + combination_type + ')'
    exec(command)
    toolbox.register("mate", tools.cxTwoPoint)
    toolbox.register("mutate", tools.mutUniformInt, low=0, up=len(df.columns) - 4, indpb=0.33)
    toolbox.register("select", tools.selTournament, tournsize=3)  # select the highest score one from 3 random pick
    #toolbox.register("select", tools.selRoulette) # another patient selection atrategy
    return toolbox


def main(toolbox, random_seed, population_size, elite_size, stop_gen, crossover_prob, mutate_prob,
         addNewInd_freq, write2File_freq, fnOut, max_generation):
    random.seed(random_seed)
    fh_log = open(fnOut[0:-4]+'_GAlog.txt', 'w')
    fh_log.write('\t'.join(['Generation', 'GeneCombo_index', 'GeneCombo_geneNA', 'Score'])+'\n')

    # Initialization
    pop = toolbox.population(n=population_size)
    hof = tools.HallOfFame(elite_size)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("min", np.min)
    stats.register("max", np.max)

    # Evolution
    continue_evolution = 1
    no_update_gen = 0
    best_score = -1
    current_generation = 0
    while continue_evolution and (current_generation < max_generation):
        current_generation += 1
        print('  Generation %d ...'%current_generation)
        # At every N generation, inject 10% new individuals
        if current_generation % addNewInd_freq == 0 and current_generation > 0:
            num_new_individuals = len(pop) // 10
            new_individuals = [toolbox.individual() for _ in range(num_new_individuals)]
            pop.sort(key=lambda ind: ind.fitness)  # sort fitness in ascending order
            pop[0:num_new_individuals] = new_individuals # replace the lowest scored ind with new ind
        # Apply crossover (mating) and mutation
        offspring = algorithms.varAnd(pop, toolbox, cxpb=crossover_prob, mutpb=mutate_prob)
        # Evaluate the individuals
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit
        # Select the next generation individuals
        offspring = toolbox.select(offspring, k=population_size)
        ### Implement elitism: ensure the best individuals remain in the next generation
        # offspring.sort(key=lambda ind: ind.fitness)  # sorts in ascending order
        # offspring[:elite_size] = hof.items
        pop = offspring

        # Update the hall of fame with the generated individuals
        hof.update(pop)
        # Gather all the fitnesses in one list and print the stats
        logbook = tools.Logbook()
        logbook.record(gen=current_generation, evals=len(invalid_ind), **stats.compile(pop))
        current_best_score = hof[0].fitness.values[0]
        top1_info = "%s\t%s\t%.5f" % (hof[0], df.columns[hof[0]].tolist(), current_best_score)
        print("  Current best individual (indices, genes, fitness): %s" %top1_info)
        fh_log.write('\t'.join([str(current_generation), str(hof[0]), ','.join(df.columns[hof[0]].tolist()), str(current_best_score)])+'\n')
        #print(logbook.stream)
        if current_best_score - best_score > 0.0001:
            no_update_gen = 0
            best_score = current_best_score
        else:
            no_update_gen += 1
        if no_update_gen > stop_gen:
            continue_evolution = 0
        # Write current global elite individuals into file every 100 generations
        if current_generation % write2File_freq == 0 and current_generation > 0:
            fhOut = open(fnOut, 'w')
            fhOut.write('Generation: %d'%current_generation + '\n')
            # Access the best individuals in the hall of fame
            best_individuals = hof.items
            # Retrieve the fitness values of the best individuals
            best_fitness_values = [ind.fitness.values for ind in best_individuals]
            for i in range(elite_size):
                topN_info = "%s\t%s\t%.5f" % (best_individuals[i], df.columns[best_individuals[i]].tolist(), best_fitness_values[i][0])
                fhOut.write(topN_info + '\n')
            fhOut.close()
    fh_log.close()
    return pop, logbook, hof


if __name__ == "__main__":
    time_start = time.time()
    op = createHelp()
    fnIn1 = op.fnIn1
    fnIn2 = op.fnIn2
    fnOut = op.fnOut
    combination_type = op.ct
    elite_size = op.topN
    fnGene = op.fnGene
    random_seed = op.rs
    population_size = op.ps
    stop_gen = op.sg
    crossover_prob = op.cp # 0.8
    mutate_prob = op.mp # 0.2

    #max_generation = 400 # force quit evolution after max_generation
    max_generation = 1000 # force quit evolution after max_generation


    addNewInd_freq = 10 # every 10 generations, adding 10% new individuals to avoid local optimum
    write2File_freq = 10 # every 10 generations, writing the current optimal solutions to file

    print('Step 1: Preprocessing ...')
    df = pre_processing(fnIn1,fnIn2, fnGene)
    true_labels = df['cell_type']
    not_true_labels = (df['cell_type'] == False)
    sampleNUM = df.shape[0]
    print('  Gene number: %d. Sample number: %d.'%((df.shape[1]-3)//2,df.shape[0]))
    print('Step 2: Setup_toolbox ...')
    toolbox = setup_toolbox(df, combination_type)
    print('Step 3: Run GA ...')
    pop, log, hof = main(toolbox, random_seed, population_size, elite_size, stop_gen, crossover_prob, mutate_prob,
                         addNewInd_freq, write2File_freq, fnOut, max_generation)
    print("Best individual (indices, genes, fitness): %s\t%s\t%s." % (hof[0], df.columns[hof[0]].tolist(), hof[0].fitness.values[0]))

    print('All done! Time used: %d s'%(time.time() - time_start))
