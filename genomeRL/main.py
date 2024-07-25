from overlap import OverlapResolver
from node import Node
from itertools import permutations 
import random
import sys
import numpy as np
import multiprocessing
import csv
import os

### GA ###

class GA:
    def __init__(self):
        self.crossover_prob = 0.7
        self.mutation_prob = 0.2
        self.ring_size = 3
        self.buffer = {}

# calculates each instance of a recurrent function to calculate sequences overlaps
    def _getOverlapValue(self, i, j, matrix, s1, s2, match, mismatch, gap):
        score = match if s1[i-1] == s2[j-1] else mismatch
        aux = max(
            matrix[i-1][j-1] + score,
            matrix[i-1][j] + gap,
            matrix[i][j-1] + gap,
            0
        )
        return aux

    def compute_overlap(self, left_read, right_read):
        for i in range(len(left_read)):
            l = left_read[i:]
            size = len(l)
            r = right_read[:size]
            if l == r:
                return l, size
        return "", 0

    def _getOverlap(self, s1, s2, match, mismatch, gap):
        l = len(s1)+1
        c = len(s2)+1
        matrix = np.array([0.0 for _ in range(l * c)]).reshape(l, c)
        for i in range(1, l):
            for j in range(1, c):
                matrix[i][j] = self._getOverlapValue(i, j, matrix, s1, s2, match, mismatch, gap)
        return np.max(matrix)

    def _getSuffixPrefixOverlap(self, left, right):
        return self.compute_overlap(left, right)[1]

    def _findOverlap(self, reads, id1, id2, match = 1.0, mismatch = -0.33, gap = -1.33):
        left = reads[id1]
        right = reads[id2]
        if left in self.buffer and right in self.buffer[left]:
            overlap = self.buffer[left][right]
        else:
            overlap = self._getSuffixPrefixOverlap(left, right)
            if left not in self.buffer:
                self.buffer[left] = {}
            self.buffer[left][right] = overlap
        return overlap

    # calculates the score of overlap between two strings
    def _findOverlap2(self, reads, id1, id2, match = 1.0, mismatch = -0.33, gap = -1.33):
        if id1 < id2:
            minId = id1
            maxId = id2
        else:
            minId = id2
            maxId = id1
        if minId in self.buffer and maxId in self.buffer[minId]:
            overlap = self.buffer[minId][maxId]
        else:
            overlap = self._getOverlap(reads[id1], reads[id2], match, mismatch, gap)
            if minId not in self.buffer:
                self.buffer[minId] = {}
            self.buffer[minId][maxId] = overlap
        return overlap

    #combina due cromosomi
    def _crossover(self, cromossome1, cromossome2):
        if type(cromossome1) == list:
            cromossome1 = np.array(cromossome1)
        if type(cromossome2) == list:
            cromossome2 = np.array(cromossome2)
        genes = np.random.choice(len(cromossome1), size=2, replace=False)
        genes.sort()

        aux1 = cromossome1[genes[0]:genes[1]+1]
        aux2 = cromossome2[genes[0]:genes[1]+1]

        diff2 = cromossome2[~np.in1d(cromossome2, aux1, assume_unique=True)]
        diff1 = cromossome1[~np.in1d(cromossome1, aux2, assume_unique=True)]

        child1 = np.append(aux1, diff2).copy()
        child2 = np.append(aux2, diff1).copy()

        return child1, child2
    
    #checks if there are no repeated reads
    def _checkCross(self, part_size, A, C, A1, C1):
        for i in range(part_size):
            for j in range(part_size):
                if A[i] == C1[j]:
                    return 0
                if A1[i] == C[j]:
                    return 0
        return 1
        


    def _RLmutation(self, reads, cromossome):
        if type(cromossome) == list:
            cromossome = np.array(cromossome)


        #divide cromosome in three chunks
        lenght = len(cromossome)
        part_size = lenght // 3
        remainder = lenght % 3

        part1 = cromossome[:part_size]
        part2 = cromossome[part_size:2 * part_size + (1 if remainder == 1 else 2 if remainder == 2 else 0)]
        part3 = cromossome[2 * part_size + (1 if remainder == 1 else 2 if remainder == 2 else 0):]

        reads1 = []
        picked_index = []
 
        #print('cromossome1', cromossome1, flush=True)
        #print('part1', part1, flush=True)
        #print('part6', part6, flush=True)
        #print('part4', part4, flush=True)
        #print('part3', part3, flush=True)

        #picks training reads
        for i in part2:
            reads1.append(reads[i])
            picked_index.append(i)
            #print('reads1', reads1, flush=True)

        
        rl = np.array(qlearning(reads1, 40, dm_rm_metrics=False))
 
            
        #rl1 = np.array(multiprocessing.Process(target=qlearning, args=(reads1, 40, dm_rm_metrics=False)))
        #rl2 = np.array(multiprocessing.Process(target=qlearning, args=(reads2, 40, dm_rm_metrics=False)))
        #rl1.start()
        #rl1.join()
        #rl2.start()
        #rl2.join()
        #print('Rl1 e rl2', rl1, rl2, flush=True)

        new = []
        #new = np.empty(len(reads1))

        #maps rl reads to cromosome indexes
        seen_elements = set()
        for k in rl:
            for l in range(len(reads)):
                if reads1[int(rl[k])] == reads[l]:
                    if picked_index[int(rl[k])] == l:
                        if l in seen_elements:
                            continue
                        new.append(l)
                        seen_elements.add(l)
                        break

        #child
        cromossome = np.concatenate((part1, np.array(new), part3))
 
        #print('child1 e child2', child1, child2, flush=True)
        return cromossome
    

    def _RLcrossover(self, reads, cromossome1, cromossome2, counter):
        if type(cromossome1) == list:
            cromossome1 = np.array(cromossome1)
        if type(cromossome2) == list:
            cromossome2 = np.array(cromossome2)

        #divide cromosome in three chunks
        lenght = len(cromossome1)
        part_size = lenght // 3
        remainder = lenght % 3

        part1 = cromossome1[:part_size]
        part2 = cromossome1[part_size:2 * part_size + (1 if remainder == 1 else 2 if remainder == 2 else 0)]
        part3 = cromossome1[2 * part_size + (1 if remainder == 1 else 2 if remainder == 2 else 0):]

        part4 = cromossome2[:part_size]
        part5 = cromossome2[part_size:2 * part_size + (1 if remainder == 1 else 2 if remainder == 2 else 0)]
        part6 = cromossome2[2 * part_size + (1 if remainder == 1 else 2 if remainder == 2 else 0):]

        reads1 = []
        reads2 = []
        picked_index1 = []
        picked_index2 = []
        
        if self._checkCross(part_size, part1, part3, part4, part6) == 1:
            #print('cromossome1', cromossome1, flush=True)
            #print('part1', part1, flush=True)
            #print('part6', part6, flush=True)
            #print('part4', part4, flush=True)
            #print('part3', part3, flush=True)

            #picks training reads
            for i in cromossome1:
                if i not in part1:
                    if i not in part6:
                        reads1.append(reads[i])
                        picked_index1.append(i)
                        #print('reads1', reads1, flush=True)
                if i not in part4:
                    if i not in part3:
                        reads2.append(reads[i])
                        picked_index2.append(i)
                        #print('reads2', reads2, flush=True)
        
                
            rl1 = np.array(qlearning(reads1, 40, test_each_episode = False, dm_rm_metrics=False))
            rl2 = np.array(qlearning(reads2, 40, test_each_episode = False, dm_rm_metrics=False))
            #print('Rl1 e rl2', rl1, rl2, flush=True)

            #maps rl reads to cromosome indexes
            seen_elements1 = set()
            seen_elements2 = set()   
            new2 = []
            new5 = []
            for k in rl1:
                for l in range(len(reads)):
                    if reads1[int(rl1[k])] == reads[l]:
                        if picked_index1[int(rl1[k])] == l:
                            if l in seen_elements1:
                                continue
                            new2.append(l)
                            seen_elements1.add(l)
                            break
                        #print('new2', new2, flush=True)
            for k in rl2:
                for l in range(len(reads)):
                    if reads2[int(rl2[k])] == reads[l]:
                        if picked_index2[int(rl2[k])] == l:
                            if l in seen_elements2:
                                continue
                            new5.append(l)
                            seen_elements2.add(l)
                            break
                        #print('new5', new5, flush=True)
        
            child1 = np.concatenate((part1, np.array(new2), part6))
            child2 = np.concatenate((part4, np.array(new5), part3))

            counter +=1
            #print('child1 e child2', child1, child2, flush=True)
            return child1, child2, counter 
        return cromossome1, cromossome2, counter
        

    #scambia 2 posizioni delle read all'interno del cromosoma
    def _mutation(self, cromossome): 
        if type(cromossome) == list:
            cromossome = np.array(cromossome)
        genes = np.random.choice(len(cromossome), size=2, replace=False)
        
        cromossome[genes[[0,1]]] = cromossome[genes[[1,0]]]
        return cromossome

    #calcola punteggio overlap tra tutte le read che compongono il cromosoma
    def _fitness(self, cromossome, reads):
        score = 0
        for i in range(1, len(cromossome)):
            score += self._findOverlap(reads, cromossome[i-1], cromossome[i])
        return score


    #calcola fit di tutti i cromosomi
    def _evaluatePopulation(self, population, reads):
        scores = np.zeros(len(population))
        for i in range(len(population)):
            scores[i] = self._fitness(population[i], reads)
        return scores

    #poppo 3 elementi a caso della popolazione e prendo elemento con punteggio maggiore, selezione da torneo
    def _ring(self, pop_fitness, ring_size):
        fighters = np.random.choice(len(pop_fitness), size=self.ring_size, replace=False)
        fit = pop_fitness[fighters]
        winner = fit.argmax()
        return fighters[winner]



    def run_ga(self, env, memory, reads, generations):
        pop_size = len(memory)
        cromo_size = len(reads)
        counter = 0

        population = np.array(memory)
        pop_fitness = self._evaluatePopulation(population, reads) #punteggio overlap complessivo di ogni cromosoma
        best_ind = population[pop_fitness.argmax()] #cromosoma con punteggio maggiore
        best_fit = self._fitness(best_ind, reads) #punteggio di fit piu alto

        fitness_evolution = np.zeros(generations)

        for generation in range(generations):
            #print('Generation {}'.format(generation+1))

            # Tournament selection
            selected = []
            for i in range(pop_size):
                winner = self._ring(pop_fitness, self.ring_size) #ring size 3
                selected.append(population[winner].copy())


            # Crossover
            for i in range(0, pop_size, 2):
                population[i], population[i+1], counter = self._RLcrossover(reads, selected[i], selected[i+1], counter) #sostituisco crossover ai cromosomi originali altrimenti lascio uguale
            print('Effettuati ' + str(counter) + ' crossing-over', flush=True)

            #RL mutation 
            for i in range(pop_size):
                population[i]= self._RLmutation(reads, selected[i])


            pop_fitness = self._evaluatePopulation(population, reads)
            # Elitism
            fitness_evolution[generation] = pop_fitness.max() #segna valore punteggio piu alto generazione
            if fitness_evolution[generation] < best_fit: #se fit di migliore individuo < di generazione precedente lo inietto nel primo slot della popolazione
                population[0] = best_ind.copy()
                fitness_evolution[generation] = best_fit 
            else: #se fit di migliore individuo > fit generazione precedente aggiorno il migliore
                best_ind = population[pop_fitness.argmax()].copy()
                best_fit = pop_fitness.max()

        return best_ind

def printTree(cur, target, verbose, offset = ""):
    if verbose:
        tag = "*" if cur == target else ""
        print(offset + "|" + cur.get_read_content() + "(" + str(cur.acc_overlap) + ")" + tag, flush=True)
    fully = True
    for _, v in cur.children.items():
        fully = fully and printTree(v, target, verbose, offset + "\t")
    return cur.is_fully_explored() and fully

def manual(reads, verbose = True):
    ovr = OverlapResolver(reads)
    root = Node.createRootNode(ovr)
    cur = root
    while True:
        printTree(root, cur, verbose)
        code = input()
        if cur.is_leaf():
            cur = root
        else:
            aux = cur.get_child(int(code))
            if aux is not None:
                cur = aux
            elif cur.parent_node is None:
                cur = root

def gt():
    reads = ['ACATTAGG', 'TTAGGCCCTT', 'GCCCTTA', 'CCTTACA']
    ovr = OverlapResolver(reads)
    for p in permutations([0,1,2,3]):
        ov = None
        for r in p:
            if ov is None:
                ov = 0
            else:
                aux = ovr.get_overlap_by_read_ids(last,r)
                if aux == 0:
                    ov = -1
                    break
                ov += aux
            last = r
        if ov > 0:
            print(p, ov)

def auto(reads, step_by_step = False, verbose = False, auto_interrupt = False):
    ovr = OverlapResolver(reads)
    root = Node.createRootNode(ovr)
    cur = root
    actions_taken, leafs_reached = 0, 0
    while True:
        codes = cur.candidate_children[:]
        codes.extend(list(cur.children.keys()))
        if cur.is_leaf() or len(codes) == 0:
            leafs_reached += 1
            if verbose or auto_interrupt:
                if printTree(root, cur, verbose) and auto_interrupt:
                    break
            if step_by_step:
                input()
            cur = root
        else:
            actions_taken += 1
            code = random.sample(codes, 1)[0]
            aux = cur.get_child(code)
            if aux is not None:
                cur = aux
            elif cur.parent_node is None:
                cur = root

    print("actions_taken", actions_taken)
    print("leafs_reached", leafs_reached)

def _get_maximum_shifted_path(root_node, leaf_node):
    path = []
    cur_node = leaf_node
    while cur_node.parent_node is not None:
        path.insert(0, cur_node.read_id)
        cur_node = cur_node.parent_node
    cur_node = root_node.get_child(leaf_node.read_id)
    if cur_node is None:
        return None
    cur_node = cur_node.get_child(path[0])
    if cur_node is None:
        return None

    i = len(path) - 1
    max_value = leaf_node.acc_overlap
    max_i = 0
    j = 1
    while i > 0:
        while True:
            cur_node = cur_node.get_child(path[j])
            j += 1
            if j >= len(path):
                j = 0
            if cur_node is None or cur_node.is_leaf():
                if cur_node is not None and cur_node.acc_overlap >= max_value:
                    max_value = cur_node.acc_overlap
                    max_i = i
                break
        i -= 1
        j = i
        cur_node = root_node
    if max_i > 0 and max_value == root_node.overlap_resolver.max_acc:
        return [path[(i + max_i) % len(path) if i + max_i >= len(path) else i + max_i] for i in range(len(path))]
    return None

def qlearning(reads, episodes, genome = None, test_each_episode = True, dm_rm_metrics=True):
    generations = 15
    num_ind = 60
    ovr = OverlapResolver(reads)
    root = Node.createRootNode(ovr)
    actions_taken, leafs_reached = 0, 0
    epsilon = 1.0
    gamma = 0.9
    alpha = 0.8
    epsilon_decay = 1.0 / episodes
    factor = 1.0 / (len(reads) * max([len(read) for read in reads]))
    forced_path = None
    succ_run_counter_dm = 0
    succ_run_counter_rm = 0

    
    #debug = [4, 15, 13, 2, 5, 3, 18, 9, 7, 16, 11, 19, 0, 8, 17, 12, 14, 1, 6, 10]
    #cur_node = root
    #for read in debug:
    #    cur_node = cur_node.get_child(read)
    #forced_path = _get_maximum_shifted_path(root, cur_node)
    #print(forced_path)

    ga = GA() #empty GA population and unset evolution flag
    memory = []
    ind_evolved = []
    pointer = -1
    for episode in range(episodes):
        cur_node = root
        total_reward = 0.0
        actions_train = []
        ind = []
        #printTree(root, ovr, verbose=True)
        while True:
            if len(ind_evolved) > 0: #reset del puntatore se individuo e evoluto e ho preso gene
                if pointer >= len(ind_evolved[0]):
                    del ind_evolved[0] #unset evolution flag
                    pointer = 0

            #pruning
            candidates = cur_node.get_outputs() #candidate children and their keys
            if len(candidates) == 0: #foglia se non ho azioni
                leafs_reached += 1
                if forced_path is None:
                    forced_path = _get_maximum_shifted_path(root, cur_node) #se non c'e percorso calcolo percorso forzato da foglia alla radice considerando gli overlap
                else:
                    forced_path = None
                break


            if len(ind_evolved) == 0: #no evolved
                if forced_path is None or len(forced_path) == 0: #seleziono casualmente una azione dai candidati se un forced path non e specificato
                    action = random.sample(candidates, 1)[0]
                else:
                    action = forced_path[0] 
                    del forced_path[0]
                rand = random.random()
                if rand > epsilon: #and cur_node != root: #epsilon greedy approach, se rand maggiore di epsilon vado greedy
                    a = cur_node.get_max_action()
                    if a is not None:
                        action = a 
                #print('Ind ', ind, flush=True)
                #print('Action ', int(action), flush=True)
                if int(action) not in ind: #risoluzione problema read ripetute
                    ind.append(int(action)) #insert action as a gene after the last gene that individual has 
                #print('ind after action appended', ind, flush=True)
            else: #evolved, ottengo prima azione del genoma evoluto
                action = ind_evolved[0][pointer] #remove the first gene from individual and set such gene as action
                pointer += 1
            next_node = cur_node.get_child(action) 
            if next_node is None:
                break

            #equazione 1.4, aggiorno qvalues e stato per azione intrapresa
            reward = 0.1 if cur_node == root else next_node.pairwise_overlap * factor 
            reward += 1.0 if next_node.is_leaf() else 0.0
            total_reward += reward
            cur_node.update_q(action, reward + gamma * next_node.get_max_qvalue(), alpha)
            actions_taken += 1
            actions_train.append(action)
            cur_node = next_node

        
        #absorbing state
        if len(ind_evolved) == 0: #evolutionary flag not set
            if len(ind) < len(reads):
                remaining = list(set(range(len(reads))) - set(ind))
                random.shuffle(remaining)
                #print('ind prima estensione', ind, flush=True)
                ind.extend(remaining)
                #print('ind dopo estensione', ind, flush=True)
            if len(ind) == len(reads) and len(set(ind)) == len(reads):
                memory.append(ind) #add individual to GA population
        if len(memory) == num_ind: #set evolution flag when GA population is completely filled, evolve this population and set its most fitted individual as evolved
            ind_evolved = [ga.run_ga(None, memory, reads, generations)]
            pointer = 0
            memory = []#empty GA population

        if test_each_episode or episode + 1 == episodes:
            test = test_qlearning(root, factor, genome)
        else:
            test = (None, 0.0, None)
        if dm_rm_metrics == True:
            print("ep.:", episode+1, "max_acc:", ovr.max_acc, "train_rw:", "%.5f" % total_reward, "test_rw:", "%.5f" % test[1], "test:", test[0], "train", actions_train, "dist:", test[2])
            if(test[2] == 0):
                succ_run_counter_dm += 1
            if(test[1] >= total_reward):
                succ_run_counter_rm += 1
        epsilon -= epsilon_decay
        '''
        if episode +1 == episodes and dm_rm_metrics:
            filename = 'results.csv'
            file_exists = os.path.isfile(filename)
            with open(filename, 'a', newline='') as csvfile:
                csvwriter = csv.writer(csvfile)
                
                if not file_exists:
                    csvwriter.writerow(['Successful RM runs', 'Episodes', 'Average RM'])
                csvwriter.writerow([succ_run_counter_rm, episodes])
            total_succ_runs, total_episodes = calc_average_rm(filename)
            if total_episodes != 0:
                average_rm = round((total_succ_runs / total_episodes) * 100, 2)
            else:
                average_rm = 0
            with open(filename, 'a', newline='') as csvfile:
                csvwriter = csv.writer(csvfile)
                csvwriter.writerow(['average_rm', '', average_rm])    

            print("Successful DM runs: ", succ_run_counter_dm, "/", episodes, flush=True)
            print("Successful RM runs: ", succ_run_counter_rm, "/", episodes, flush=True)
            print("Media DM: ", "%.2f" % (succ_run_counter_dm / episodes*100), "%", flush=True)
            print("Media RM: ", "%.2f" % (succ_run_counter_rm / episodes*100), "%", flush=True)
            '''

    return ind
    #print("actions_taken", actions_taken)
    #print("leafs_reached", leafs_reached)

def calc_average_rm(filename):
    total_succ_runs = 0
    total_episodes = 0
    if os.path.isfile(filename):
        with open(filename, 'r', newline='') as csvfile:
            csvreader = csv.reader(csvfile)
            next(csvreader)  # Salta l'intestazione
            for row in csvreader:
                if row and row[0] != 'average_rm':  # Assicurati che la riga non sia vuota e non sia la riga della media
                    total_succ_runs += float(row[0])
                    total_episodes += float(row[1])
    return total_succ_runs, total_episodes


def test_qlearning(root_node, factor, genome):
    def levenshtein(s, t, costs=(1, 1, 1)):
        """
            iterative_levenshtein(s, t) -> ldist
            ldist is the Levenshtein distance between the strings
            s and t.
            For all i and j, dist[i,j] will contain the Levenshtein
            distance between the first i characters of s and the
            first j characters of t

            costs: a tuple or a list with three integers (d, i, s)
                   where d defines the costs for a deletion
                         i defines the costs for an insertion and
                         s defines the costs for a substitution
        """
        rows = len(s)+1
        cols = len(t)+1
        deletes, inserts, substitutes = costs

        dist = [[0 for x in range(cols)] for x in range(rows)]
        # source prefixes can be transformed into empty strings
        # by deletions:
        for row in range(1, rows):
            dist[row][0] = row * deletes
        # target prefixes can be created from an empty source string
        # by inserting the characters
        for col in range(1, cols):
            dist[0][col] = col * inserts

        for col in range(1, cols):
            for row in range(1, rows):
                if s[row-1] == t[col-1]:
                    cost = 0
                else:
                    cost = substitutes
                dist[row][col] = min(dist[row-1][col] + deletes,
                                     dist[row][col-1] + inserts,
                                     dist[row-1][col-1] + cost) # substitution
        return dist[rows-1][cols-1]    
    cur_node = root_node
    actions = []
    total_reward = 0.0
    while True:
        a = cur_node.get_max_action()
        if a is None:
            break
        actions.append(a)
        aux = cur_node.get_child(a)
        if aux is None:
            break
        cur_node = aux
        reward = 0.1 if cur_node.parent_node == root_node else cur_node.pairwise_overlap * factor
        reward += 1.0 if cur_node.is_leaf() else 0.0
        total_reward += reward #reward system 1.4
    dist = None
    if genome is not None: 
        dist = levenshtein(cur_node.get_consensus(), genome) #numero di edits per arrivare a genoma originale
    return actions, total_reward, dist

if __name__ == "__main__":
    if len(sys.argv) > 3:
        seed = int(sys.argv[3])
    else:
        seed = random.randrange(sys.maxsize)
    random.seed(seed)
    print(seed, file=sys.stderr)
    
    dataset = {}
    genome25 = "TACTAGCAATACGCTTGCGTTCGGT"
    genome50 = "CCTAACCATTTTAACAGCAACATAACAGGCTAAGAGGGGCCGGACACCCA"
    genome381 = "ATGGCAATATTAGGTTTAGGCACGGATATTGTGGAGATCGCTCGCATCGAAGCGGTGATCGCCCGATCCGGTGATCGCCTGGCACGCCGCGTATTAAGCGATAACGAATGGGCTATCTGGAAAACGCACCACCAGCCGGTGCGTTTTCTGGCGAAGCGTTTTGCTGTGAAAGAAGCCGCAGCAAAAGCGTTTGGCACCGGGATCCGCAATGGTCTGGCGTTTAATCAATTTGAAGTATTCAATGATGAGCTCGGCAAACCACGGCTACGGCTATGGGGCGAGGCATTAAAACTGGCGGAAAAGCTGGGCGTTGCAAATATGCATGTAACGCTGGCAGATGAGCGGCACTATGCTTGTGCCACGGTAATTATTGAAAGTTAA"
    genome567 = "ATGAGCAAAGCAGGTGCGTCGCTTGCGACCTGTTACGGCCCTGTCAGCGCCGACGTTATAGCAAAAGCAGAGAACATTCGTCTGCTGATCCTCGATGTCGATGGCGTACTGTCAGATGGCCTGATTTATATGGGCAATAATGGCGAAGAGCTGAAAGCGTTCAATGTTCGTGACGGTTATGGCATTCGTTGTGCGCTCACCTCTGATATTGAAGTCGCTATCATTACCGGGCGAAAGGCTAAACTGGTAGAAGATCGTTGTGCCACATTGGGGATCACTCACTTGTATCAGGGGCAGTCAAACAAACTGATCGCCTTTAGCGATCTGCTGGAAAAACTGGCGATTGCCCCGGAAAATGTGGCTTATGTCGGCGATGATCTCATCGACTGGCCGGTAATGGAAAAAGTGGGTTTAAGCGTCGCCGTGGCCGATGCGCATCCACTGTTGATCCCGCGCGCCGATTACGTGACGCGCATTGCTGGCGGTCGTGGCGCAGTGCGCGAAGTTTGCGACTTATTACTCCTGGCGCAGGGCAAACTGGATGAAGCCAAAGGGCAATCGATATGA"
    genome726 = 'ATGGCAACTGTTTCCATGCGCGACATGCTCAAGGCTGGTGTTCACTTCGGTCACCAGACCCGTTACTGGAACCCGAAAATGAAGCCGTTCATCTTCGGTGCGCGTAACAAAGTTCACATCATCAACCTTGAGAAAACTGTACCGATGTTCAACGAAGCTCTGGCTGAACTGAACAAGATTGCTTCTCGCAAAGGTAAAATCCTTTTCGTTGGTACTAAACGCGCTGCAAGCGAAGCGGTGAAAGACGCTGCTCTGAGCTGCGACCAGTTCTTCGTGAACCATCGCTGGCTGGGCGGTATGCTGACTAACTGGAAAACCGTTCGTCAGTCCATCAAACGTCTGAAAGACCTGGAAACTCAGTCTCAGGACGGTACTTTCGACAAGCTGACCAAGAAAGAAGCGCTGATGCGCACTCGTGAGCTGGAGAAACTGGAAAACAGCCTGGGCGGTATCAAAGACATGGGCGGTCTGCCGGACGCTCTGTTTGTAATCGATGCTGACCACGAACACATTGCTATCAAAGAAGCAAACAACCTGGGTATTCCGGTATTTGCTATCGTTGATACCAACTCTGATCCGGACGGTGTTGACTTCGTTATCCCGGGTAACGACGACGCAATCCGTGCTGTGACCCTGTACCTGGGCGCTGTTGCTGCAACCGTACGTGAAGGCCGTTCTCAGGATCTGGCTTCCCAGGCGGAAGAAAGCTTCGTAGAAGCTGAGTAA'
    genome930 = 'ATGACGCAATTTGCATTTGTGTTCCCTGGACAGGGTTCTCAAACCGTTGGAATGCTGGCTGATATGGCGGCGAGCTATCCAATTGTCGAAGAAACGTTTGCTGAAGCTTCTGCGGCGCTGGGCTACGACCTGTGGGCGCTGACCCAGCAGGGGCCAGCTGAAGAACTGAATAAAACCTGGCAAACTCAGCCTGCGCTGTTGACTGCATCTGTTGCGCTGTATCGCGTATGGCAGCAGCAGGGCGGTAAAGCACCGGCAATGATGGCCGGTCACAGCCTGGGGGAATACTCCGCGCTGGTTTGCGCTGGTGTGATTGATTTCGCTGATGCGGTGCGTCTGGTTGAGATGCGCGGCAAGTTCATGCAAGAAGCCGTACCGGAAGGCACGGGCGCTATGGCGGCAATCATCGGTCTGGATGATGCGTCTATTGCGAAAGCGTGTGAAGAAGCTGCAGAAGGTCAGGTCGTTTCTCCGGTAAACTTTAACTCTCCGGGACAGGTGGTTATTGCCGGTCATAAAGAAGCGGTTGAGCGTGCTGGCGCTGCCTGTAAAGCGGCGGGCGCAAAACGCGCGCTGCCGTTACCAGTGAGCGTACCGTCTCACTGTGCGCTGATGAAACCAGCAGCCGACAAACTGGCAGTAGAATTAGCGAAAATCACCTTTAACGCACCAACAGTTCCTGTTGTGAATAACGTTGATGTGAAATGCGAAACCAATGGTGATGCCATCCGTGACGCACTGGTACGTCAGTTGTATAACCCGGTTCAGTGGACGAAGTCTGTTGAGTACATGGCAGCGCAAGGCGTAGAACATCTCTATGAAGTCGGCCCGGGCAAAGTGCTTACTGGCCTGACGAAACGCATTGTCGACACCCTGACCGCCTCGGCGCTGAACGAACCTTCAGCGATGGCAGCGGCGCTCGAGCTTTAA'
    genome4224 = 'GTGAAAGATTTATTAAAGTTTCTGAAAGCGCAGACTAAAACCGAAGAGTTTGATGCGATCAAAATTGCTCTGGCTTCGCCAGACATGATCCGTTCATGGTCTTTCGGTGAAGTTAAAAAGCCGGAAACCATCAACTACCGTACGTTCAAACCAGAACGTGACGGCCTTTTCTGCGCCCGTATCTTTGGGCCGGTAAAAGATTACGAGTGCCTGTGCGGTAAGTACAAGCGCCTGAAACACCGTGGCGTCATCTGTGAGAAGTGCGGCGTTGAAGTGACCCAGACTAAAGTACGCCGTGAGCGTATGGGCCACATCGAACTGGCTTCCCCGACTGCGCACATCTGGTTCCTGAAATCGCTGCCGTCCCGTATCGGTCTGCTGCTCGATATGCCGCTGCGCGATATCGAACGCGTACTGTACTTTGAATCCTATGTGGTTATCGAAGGCGGTATGACCAACCTGGAACGTCAGCAGATCCTGACTGAAGAGCAGTATCTGGACGCGCTGGAAGAGTTCGGTGACGAATTCGACGCGAAGATGGGGGCGGAAGCAATCCAGGCTCTGCTGAAGAGCATGGATCTGGAGCAAGAGTGCGAACAGCTGCGTGAAGAGCTGAACGAAACCAACTCCGAAACCAAGCGTAAAAAGCTGACCAAGCGTATCAAACTGCTGGAAGCGTTCGTTCAGTCTGGTAACAAACCAGAGTGGATGATCCTGACCGTTCTGCCGGTACTGCCGCCAGATCTGCGTCCGCTGGTTCCGCTGGATGGTGGTCGTTTCGCGACTTCTGACCTGAACGATCTGTATCGTCGCGTCATTAACCGTAACAACCGTCTGAAACGTCTGCTGGATCTGGCTGCGCCGGACATCATCGTACGTAACGAAAAACGTATGCTGCAGGAAGCGGTAGACGCCCTGCTGGATAACGGTCGTCGCGGTCGTGCGATCACCGGTTCTAACAAGCGTCCTCTGAAATCTTTGGCCGACATGATCAAAGGTAAACAGGGTCGTTTCCGTCAGAACCTGCTCGGTAAGCGTGTTGACTACTCCGGTCGTTCTGTAATCACCGTAGGTCCATACCTGCGTCTGCATCAGTGCGGTCTGCCGAAGAAAATGGCACTGGAGCTGTTCAAACCGTTCATCTACGGCAAGCTGGAACTGCGTGGTCTTGCTACCACCATTAAAGCTGCGAAGAAAATGGTTGAGCGCGAAGAAGCTGTCGTTTGGGATATCCTGGACGAAGTTATCCGCGAACACCCGGTACTGCTGAACCGTGCACCGACTCTGCACCGTCTGGGTATCCAGGCATTTGAACCGGTACTGATCGAAGGTAAAGCTATCCAGCTGCACCCGCTGGTTTGTGCGGCATATAACGCCGACTTCGATGGTGACCAGATGGCTGTTCACGTACCGCTGACGCTGGAAGCCCAGCTGGAAGCGCGTGCGCTGATGATGTCTACCAACAACATCCTGTCCCCGGCGAACGGCGAACCAATCATCGTTCCGTCTCAGGACGTTGTACTGGGTCTGTACTACATGACCCGTGACTGTGTTAACGCCAAAGGCGAAGGCATGGTGCTGACTGGCCCGAAAGAAGCAGAACGTCTGTATCGCTCTGGTCTGGCTTCTCTGCATGCGCGCGTTAAAGTGCGTATCACCGAGTATGAAAAAGATGCTAACGGTGAATTAGTAGCGAAAACCAGCCTGAAAGACACGACTGTTGGCCGTGCCATTCTGTGGATGATTGTACCGAAAGGTCTGCCTTACTCCATCGTCAACCAGGCGCTGGGTAAAAAAGCAATCTCCAAAATGCTGAACACCTGCTACCGCATTCTCGGTCTGAAACCGACCGTTATTTTTGCGGACCAGATCATGTACACCGGCTTCGCCTATGCAGCGCGTTCTGGTGCATCTGTTGGTATCGATGACATGGTCATCCCGGAGAAGAAACACGAAATCATCTCCGAGGCAGAAGCAGAAGTTGCTGAAATTCAGGAGCAGTTCCAGTCTGGTCTGGTAACTGCGGGCGAACGCTACAACAAAGTTATCGATATCTGGGCTGCGGCGAACGATCGTGTATCCAAAGCGATGATGGATAACCTGCAAACTGAAACCGTGATTAACCGTGACGGTCAGGAAGAGAAGCAGGTTTCCTTCAACAGCATCTACATGATGGCCGACTCCGGTGCGCGTGGTTCTGCGGCACAGATTCGTCAGCTTGCTGGTATGCGTGGTCTGATGGCGAAGCCGGATGGCTCCATCATCGAAACGCCAATCACCGCGAACTTCCGTGAAGGTCTGAACGTACTCCAGTACTTCATCTCCACCCACGGTGCTCGTAAAGGTCTGGCGGATACCGCACTGAAAACTGCGAACTCCGGTTACCTGACTCGTCGTCTGGTTGACGTGGCGCAGGACCTGGTGGTTACCGAAGACGATTGTGGTACCCATGAAGGTATCATGATGACTCCGGTTATCGAGGGTGGTGACGTTAAAGAGCCGCTGCGCGATCGCGTACTGGGTCGTGTAACTGCTGAAGACGTTCTGAAGCCGGGTACTGCTGATATCCTCGTTCCGCGCAACACGCTGCTGCACGAACAGTGGTGTGACCTGCTGGAAGAGAACTCTGTCGACGCGGTTAAAGTACGTTCTGTTGTATCTTGTGACACCGACTTTGGTGTATGTGCGCACTGCTACGGTCGTGACCTGGCGCGTGGCCACATCATCAACAAGGGTGAAGCAATCGGTGTTATCGCGGCACAGTCCATCGGTGAACCGGGTACACAGCTGACCATGCGTACGTTCCACATCGGTGGTGCGGCATCTCGTGCGGCTGCTGAATCCAGCATCCAAGTGAAAAACAAAGGTAGCATCAAGCTCAGCAACGTGAAGTCGGTTGTGAACTCCAGCGGTAAACTGGTTATCACTTCCCGTAATACTGAACTGAAACTGATCGACGAATTCGGTCGTACTAAAGAAAGCTACAAAGTACCTTACGGTGCGGTACTGGCGAAAGGCGATGGCGAACAGGTTGCTGGCGGCGAAACCGTTGCAAACTGGGACCCGCACACCATGCCGGTTATCACCGAAGTAAGCGGTTTTGTACGCTTTACTGACATGATCGACGGCCAGACCATTACGCGTCAGACCGACGAACTGACCGGTCTGTCTTCGCTGGTGGTTCTGGATTCCGCAGAACGTACCGCAGGTGGTAAAGATCTGCGTCCGGCACTGAAAATCGTTGATGCTCAGGGTAACGACGTTCTGATCCCAGGTACCGATATGCCAGCGCAGTACTTCCTGCCGGGTAAAGCGATTGTTCAGCTGGAAGATGGCGTACAGATCAGCTCTGGTGACACCCTGGCGCGTATTCCGCAGGAATCCGGCGGTACCAAGGACATCACCGGTGGTCTGCCGCGCGTTGCGGACCTGTTCGAAGCACGTCGTCCGAAAGAGCCGGCAATCCTGGCTGAAATCAGCGGTATCGTTTCCTTCGGTAAAGAAACCAAAGGTAAACGTCGTCTGGTTATCACCCCGGTAGACGGTAGCGATCCGTACGAAGAGATGATTCCGAAATGGCGTCAGCTCAACGTGTTCGAAGGTGAACGTGTAGAACGTGGTGACGTAATTTCCGACGGTCCGGAAGCGCCGCACGACATTCTGCGTCTGCGTGGTGTTCATGCTGTTACTCGTTACATCGTTAACGAAGTACAGGACGTATACCGTCTGCAGGGCGTTAAGATTAACGATAAACACATCGAAGTTATCGTTCGTCAGATGCTGCGTAAAGCTACCATCGTTAACGCGGGTAGCTCCGACTTCCTGGAAGGCGAACAGGTTGAATACTCTCGCGTCAAGATCGCAAACCGCGAACTGGAAGCGAACGGCAAAGTGGGTGCAACTTACTCCCGCGATCTGCTGGGTATCACCAAAGCGTCTCTGGCAACCGAGTCCTTCATCTCCGCGGCATCGTTCCAGGAGACCACTCGCGTGCTGACCGAAGCAGCCGTTGCGGGCAAACGCGACGAACTGCGCGGCCTGAAAGAGAACGTTATCGTGGGTCGTCTGATCCCGGCAGGTACCGGTTACGCGTACCACCAGGATCGTATGCGTCGCCGTGCTGCGGGTGAAGCTCCGGCTGCACCGCAGGTGACTGCAGAAGACGCATCTGCCAGCCTGGCAGAACTGCTGAACGCAGGTCTGGGCGGTTCTGATAACGAGTAA' 
    
    # reads = ['ACATTAGG', 'TTAGGCCCTT', 'GCCCTTA', 'CCTTACA']
    #reads = ['GCTTGCGT','GCGTTCGG','CTAGCAAT','GCAATACG','CTAGCAAT','TAGCAATA','CGTTCGGT','TACGCTTG','CTTGCGTT','TACGCTTG','TACTAGCA','ACGCTTGC','TGCGTTCG','CTAGCAAT','TGCGTTCG','CTAGCAAT','ACGCTTGC','TTGCGTTC','ATACGCTT','CGCTTGCG']
    #reads = ['CAACATAA','TTTAACAG','AGGCTAAG','GGGCCGGA','GGCTAAGA','AACATAAC','TAACAGCA','ACATAACA','AGGCTAAG','CAACATAA','GGGCCGGA','CCATTTTA','AACATAAC','CCATTTTA','CCTAACCA','GAGGGGCC','GACACCCA','TAACAGGC','TAACAGCA','GCCGGACA']
    #reads = ['GGGGCCGG','ACAGCAAC','AGAGGGGC','CCTAACCA','CCATTTTA','AAGAGGGG','TTAACAGC','AGGCTAAG','AGGGGCCG','GGGCCGGA','CCGGACAC','TAAGAGGG','CATAACAG','TTAACAGC','GAGGGGCC','TTTAACAG','CAGCAACA','GGGGCCGG','GGCTAAGA','CGGACACC','CCTAACCA','AACCATTT','AGGGGCCG','GACACCCA','AGCAACAT','CCGGACAC','ACCATTTT','TAAGAGGG','GGCCGGAC','GCTAAGAG']
    #gt()
    #manual(reads)
    #auto(reads)
    reads_25_10_8 = ['CGTTCGGT','TTGCGTTC','CTTGCGTT','ACGCTTGC','ATACGCTT','AATACGCT','AGCAATAC','CTAGCAAT','ACTAGCAA','TACTAGCA']
    reads_25_10_10 = ['TGCGTTCGGT','CTAGCAATAC','ACTAGCAATA','CAATACGCTT','GCTTGCGTTC','CTTGCGTTCG','GCTTGCGTTC','ACTAGCAATA','TACTAGCAAT','CTAGCAATAC']
    reads_25_10_15 = ['CAATACGCTTGCGTT','AGCAATACGCTTGCG','GCAATACGCTTGCGT','TACTAGCAATACGCT','ACGCTTGCGTTCGGT','ATACGCTTGCGTTCG','TAGCAATACGCTTGC','CTAGCAATACGCTTG','AATACGCTTGCGTTC','AGCAATACGCTTGCG']
    reads_50_10_8 = ['GGGGCCGG','GGACACCC','ATAACAGG','GACACCCA','AAGAGGGG','TTAACAGC','AACATAAC','ATTTTAAC','CCTAACCA','ACAGGCTA']
    reads_50_10_10 = ['CGGACACCCA','CCTAACCATT','TAACCATTTT','CTAAGAGGGG','ACATAACAGG','GGGCCGGACA','GAGGGGCCGG','TTAACAGCAA','CATAACAGGC','ACAGCAACAT']
    reads_50_10_15 = ['CAACATAACAGGCTA','CCTAACCATTTTAAC','ATTTTAACAGCAACA','GGGGCCGGACACCCA','CCATTTTAACAGCAA','TAACCATTTTAACAG','TAAGAGGGGCCGGAC','AACAGCAACATAACA','ACAGGCTAAGAGGGG','TTTTAACAGCAACAT']
    reads_25_20_8 = ['GCTTGCGT','GCGTTCGG','CTAGCAAT','GCAATACG','CTAGCAAT','TAGCAATA','CGTTCGGT','TACGCTTG','CTTGCGTT','TACGCTTG','TACTAGCA','ACGCTTGC','TGCGTTCG','CTAGCAAT','TGCGTTCG','CTAGCAAT','ACGCTTGC','TTGCGTTC','ATACGCTT','CGCTTGCG']
    reads_25_20_10 = ['CAATACGCTT','GCTTGCGTTC','AATACGCTTG','CAATACGCTT','TACTAGCAAT','ATACGCTTGC','TGCGTTCGGT','TAGCAATACG','AGCAATACGC','TACGCTTGCG','TACGCTTGCG','CTAGCAATAC','CAATACGCTT','ACGCTTGCGT','ATACGCTTGC','TAGCAATACG','ACGCTTGCGT','GCAATACGCT','CAATACGCTT','AGCAATACGC']
    reads_25_20_15 = ['ACGCTTGCGTTCGGT','AATACGCTTGCGTTC','ACGCTTGCGTTCGGT','AGCAATACGCTTGCG','TAGCAATACGCTTGC','AATACGCTTGCGTTC','ATACGCTTGCGTTCG','CTAGCAATACGCTTG','ACGCTTGCGTTCGGT','ACGCTTGCGTTCGGT','ATACGCTTGCGTTCG','ACGCTTGCGTTCGGT','ACGCTTGCGTTCGGT','TACTAGCAATACGCT','CTAGCAATACGCTTG','ACGCTTGCGTTCGGT','GCAATACGCTTGCGT','TACGCTTGCGTTCGG','CAATACGCTTGCGTT','TACGCTTGCGTTCGG']
    reads_50_20_8 = ['CAACATAA','TTTAACAG','AGGCTAAG','GGGCCGGA','GGCTAAGA','AACATAAC','TAACAGCA','ACATAACA','AGGCTAAG','CAACATAA','GGGCCGGA','CCATTTTA','AACATAAC','CCATTTTA','CCTAACCA','GAGGGGCC','GACACCCA','TAACAGGC','TAACAGCA','GCCGGACA']
    reads_50_20_10 = ['GAGGGGCCGG','TAACCATTTT','TAACAGGCTA','TTTTAACAGC','CATTTTAACA','ACCATTTTAA','GGCTAAGAGG','CAGCAACATA','AGAGGGGCCG','TAACCATTTT','CCATTTTAAC','AACAGGCTAA','CTAACCATTT','GGGCCGGACA','CCTAACCATT','CCATTTTAAC','TTTTAACAGC','CAACATAACA','ACATAACAGG','CGGACACCCA']
    reads_50_20_15 = ['ATAACAGGCTAAGAG','TTTTAACAGCAACAT','CAGGCTAAGAGGGGC','GGCTAAGAGGGGCCG','ATAACAGGCTAAGAG','CATAACAGGCTAAGA','GCTAAGAGGGGCCGG','TTTTAACAGCAACAT','CAACATAACAGGCTA','ACAGGCTAAGAGGGG','GGCTAAGAGGGGCCG','TAAGAGGGGCCGGAC','AACAGGCTAAGAGGG','TAACAGGCTAAGAGG','TAACAGCAACATAAC','CATAACAGGCTAAGA','CCTAACCATTTTAAC','GGGGCCGGACACCCA','ACAGCAACATAACAG','TTTTAACAGCAACAT']
    reads_25_30_8 = ['TACTAGCA','AATACGCT','CGTTCGGT','CGCTTGCG','AGCAATAC','TGCGTTCG','TAGCAATA','CTTGCGTT','TAGCAATA','CTAGCAAT','GCGTTCGG','TTGCGTTC','TTGCGTTC','TGCGTTCG','GCGTTCGG','TACGCTTG','CAATACGC','ACGCTTGC','GCTTGCGT','CGCTTGCG','ATACGCTT','CGTTCGGT','CGCTTGCG','GCAATACG','GCTTGCGT','ACGCTTGC','CTTGCGTT','TTGCGTTC','GCTTGCGT','TACTAGCA'] 
    reads_25_30_10 = ['TAGCAATACG','CGCTTGCGTT','TGCGTTCGGT','ACGCTTGCGT','CTAGCAATAC','ACTAGCAATA','TTGCGTTCGG','TGCGTTCGGT','ATACGCTTGC','CGCTTGCGTT','CTAGCAATAC','CTAGCAATAC','TAGCAATACG','GCAATACGCT','TACGCTTGCG','AGCAATACGC','GCTTGCGTTC','ACGCTTGCGT','ATACGCTTGC','CAATACGCTT','AATACGCTTG','TAGCAATACG','GCAATACGCT','TACGCTTGCG','AGCAATACGC','TAGCAATACG','CGCTTGCGTT','TAGCAATACG','TACTAGCAAT','GCTTGCGTTC']
    reads_25_30_15 = ['ACTAGCAATACGCTT','ACTAGCAATACGCTT','AATACGCTTGCGTTC','ACTAGCAATACGCTT','ACGCTTGCGTTCGGT','TACGCTTGCGTTCGG','AATACGCTTGCGTTC','AATACGCTTGCGTTC','CTAGCAATACGCTTG','TACTAGCAATACGCT','CTAGCAATACGCTTG','CAATACGCTTGCGTT','TACGCTTGCGTTCGG','AATACGCTTGCGTTC','AATACGCTTGCGTTC','CAATACGCTTGCGTT','TACGCTTGCGTTCGG','CAATACGCTTGCGTT','AATACGCTTGCGTTC','ACTAGCAATACGCTT','ACTAGCAATACGCTT','ACTAGCAATACGCTT','AGCAATACGCTTGCG','CTAGCAATACGCTTG','ATACGCTTGCGTTCG','GCAATACGCTTGCGT','ATACGCTTGCGTTCG','ACTAGCAATACGCTT','ATACGCTTGCGTTCG','TAGCAATACGCTTGC']
    reads_50_30_8 = ['GGGGCCGG','ACAGCAAC','AGAGGGGC','CCTAACCA','CCATTTTA','AAGAGGGG','TTAACAGC','AGGCTAAG','AGGGGCCG','GGGCCGGA','CCGGACAC','TAAGAGGG','CATAACAG','TTAACAGC','GAGGGGCC','TTTAACAG','CAGCAACA','GGGGCCGG','GGCTAAGA','CGGACACC','CCTAACCA','AACCATTT','AGGGGCCG','GACACCCA','AGCAACAT','CCGGACAC','ACCATTTT','TAAGAGGG','GGCCGGAC','GCTAAGAG'] 
    reads_50_30_10 = ['ACATAACAGG','CTAACCATTT','CAGGCTAAGA','AGCAACATAA','CTAAGAGGGG','AGAGGGGCCG','CAGGCTAAGA','CGGACACCCA','AACAGGCTAA','CAGCAACATA','GGCCGGACAC','GGCCGGACAC','TAACAGCAAC','CCTAACCATT','GGGGCCGGAC','CAGGCTAAGA','GGCTAAGAGG','TAAGAGGGGC','AACATAACAG','CAGCAACATA','TAACAGGCTA','TTTTAACAGC','ACCATTTTAA','AACATAACAG','AGAGGGGCCG','GCAACATAAC','TAACAGGCTA','GGCTAAGAGG','TAACCATTTT','CAGGCTAAGA']
    reads_50_30_15 = ['CTAAGAGGGGCCGGA','AACAGCAACATAACA','GGGGCCGGACACCCA','AGCAACATAACAGGC','AACAGGCTAAGAGGG','ATAACAGGCTAAGAG','TAACAGCAACATAAC','GAGGGGCCGGACACC','CTAAGAGGGGCCGGA','ATAACAGGCTAAGAG','TTTTAACAGCAACAT','CATTTTAACAGCAAC','CCTAACCATTTTAAC','TTTTAACAGCAACAT','GCAACATAACAGGCT','GAGGGGCCGGACACC','GGCTAAGAGGGGCCG','ATTTTAACAGCAACA','ACAGGCTAAGAGGGG','TTTAACAGCAACATA','CAACATAACAGGCTA','CAACATAACAGGCTA','TTTTAACAGCAACAT','AACATAACAGGCTAA','CATTTTAACAGCAAC','TAACCATTTTAACAG','AACATAACAGGCTAA','CTAAGAGGGGCCGGA','AGGCTAAGAGGGGCC','CTAAGAGGGGCCGGA']
    reads_381_20_75 = ['AAAGAAGCCGCAGCAAAAGCGTTTGGCACCGGGATCCGCAATGGTCTGGCGTTTAATCAATTTGAAGTATTCAAT',
    'GGCGTTGCAAATATGCATGTAACGCTGGCAGATGAGCGGCACTATGCTTGTGCCACGGTAATTATTGAAAGTTAA',
    'AAGCGTTTGGCACCGGGATCCGCAATGGTCTGGCGTTTAATCAATTTGAAGTATTCAATGATGAGCTCGGCAAAC',
    'GCAGCAAAAGCGTTTGGCACCGGGATCCGCAATGGTCTGGCGTTTAATCAATTTGAAGTATTCAATGATGAGCTC',
    'ATGATGAGCTCGGCAAACCACGGCTACGGCTATGGGGCGAGGCATTAAAACTGGCGGAAAAGCTGGGCGTTGCAA',
    'CTGGCACGCCGCGTATTAAGCGATAACGAATGGGCTATCTGGAAAACGCACCACCAGCCGGTGCGTTTTCTGGCG',
    'TGAAGTATTCAATGATGAGCTCGGCAAACCACGGCTACGGCTATGGGGCGAGGCATTAAAACTGGCGGAAAAGCT',
    'GCACCGGGATCCGCAATGGTCTGGCGTTTAATCAATTTGAAGTATTCAATGATGAGCTCGGCAAACCACGGCTAC',
    'GGCGAGGCATTAAAACTGGCGGAAAAGCTGGGCGTTGCAAATATGCATGTAACGCTGGCAGATGAGCGGCACTAT',
    'GGCAATATTAGGTTTAGGCACGGATATTGTGGAGATCGCTCGCATCGAAGCGGTGATCGCCCGATCCGGTGATCG',
    'TTAGGCACGGATATTGTGGAGATCGCTCGCATCGAAGCGGTGATCGCCCGATCCGGTGATCGCCTGGCACGCCGC',
    'TATCTGGAAAACGCACCACCAGCCGGTGCGTTTTCTGGCGAAGCGTTTTGCTGTGAAAGAAGCCGCAGCAAAAGC',
    'ATGGCAATATTAGGTTTAGGCACGGATATTGTGGAGATCGCTCGCATCGAAGCGGTGATCGCCCGATCCGGTGAT',
    'GCGAGGCATTAAAACTGGCGGAAAAGCTGGGCGTTGCAAATATGCATGTAACGCTGGCAGATGAGCGGCACTATG',
    'GGTGCGTTTTCTGGCGAAGCGTTTTGCTGTGAAAGAAGCCGCAGCAAAAGCGTTTGGCACCGGGATCCGCAATGG',
    'TTAATCAATTTGAAGTATTCAATGATGAGCTCGGCAAACCACGGCTACGGCTATGGGGCGAGGCATTAAAACTGG',
    'GCGTTTTCTGGCGAAGCGTTTTGCTGTGAAAGAAGCCGCAGCAAAAGCGTTTGGCACCGGGATCCGCAATGGTCT',
    'TTTTCTGGCGAAGCGTTTTGCTGTGAAAGAAGCCGCAGCAAAAGCGTTTGGCACCGGGATCCGCAATGGTCTGGC',
    'CCGATCCGGTGATCGCCTGGCACGCCGCGTATTAAGCGATAACGAATGGGCTATCTGGAAAACGCACCACCAGCC',
    'GGGCGTTGCAAATATGCATGTAACGCTGGCAGATGAGCGGCACTATGCTTGTGCCACGGTAATTATTGAAAGTTA']
    reads_567_30_75 = ['GTCGCCGTGGCCGATGCGCATCCACTGTTGATCCCGCGCGCCGATTACGTGACGCGCATTGCTGGCGGTCGTGGC',
    'GTATCAGGGGCAGTCAAACAAACTGATCGCCTTTAGCGATCTGCTGGAAAAACTGGCGATTGCCCCGGAAAATGT',
    'GATCCCGCGCGCCGATTACGTGACGCGCATTGCTGGCGGTCGTGGCGCAGTGCGCGAAGTTTGCGACTTATTACT',
    'ATGAGCAAAGCAGGTGCGTCGCTTGCGACCTGTTACGGCCCTGTCAGCGCCGACGTTATAGCAAAAGCAGAGAAC',
    'TGTGCGCTCACCTCTGATATTGAAGTCGCTATCATTACCGGGCGAAAGGCTAAACTGGTAGAAGATCGTTGTGCC',
    'TCCACTGTTGATCCCGCGCGCCGATTACGTGACGCGCATTGCTGGCGGTCGTGGCGCAGTGCGCGAAGTTTGCGA',
    'GCGTACTGTCAGATGGCCTGATTTATATGGGCAATAATGGCGAAGAGCTGAAAGCGTTCAATGTTCGTGACGGTT',
    'GTGACGCGCATTGCTGGCGGTCGTGGCGCAGTGCGCGAAGTTTGCGACTTATTACTCCTGGCGCAGGGCAAACTG',
    'TGGCGAAGAGCTGAAAGCGTTCAATGTTCGTGACGGTTATGGCATTCGTTGTGCGCTCACCTCTGATATTGAAGT',
    'ACGCGCATTGCTGGCGGTCGTGGCGCAGTGCGCGAAGTTTGCGACTTATTACTCCTGGCGCAGGGCAAACTGGAT',
    'CGCCGTGGCCGATGCGCATCCACTGTTGATCCCGCGCGCCGATTACGTGACGCGCATTGCTGGCGGTCGTGGCGC',
    'CTCACTTGTATCAGGGGCAGTCAAACAAACTGATCGCCTTTAGCGATCTGCTGGAAAAACTGGCGATTGCCCCGG',
    'TAAGCGTCGCCGTGGCCGATGCGCATCCACTGTTGATCCCGCGCGCCGATTACGTGACGCGCATTGCTGGCGGTC',
    'GCGCATCCACTGTTGATCCCGCGCGCCGATTACGTGACGCGCATTGCTGGCGGTCGTGGCGCAGTGCGCGAAGTT',
    'TTTAGCGATCTGCTGGAAAAACTGGCGATTGCCCCGGAAAATGTGGCTTATGTCGGCGATGATCTCATCGACTGG',
    'GGCTTATGTCGGCGATGATCTCATCGACTGGCCGGTAATGGAAAAAGTGGGTTTAAGCGTCGCCGTGGCCGATGC',
    'GCGTTCAATGTTCGTGACGGTTATGGCATTCGTTGTGCGCTCACCTCTGATATTGAAGTCGCTATCATTACCGGG',
    'TCCTCGATGTCGATGGCGTACTGTCAGATGGCCTGATTTATATGGGCAATAATGGCGAAGAGCTGAAAGCGTTCA',
    'TCGTTGTGCGCTCACCTCTGATATTGAAGTCGCTATCATTACCGGGCGAAAGGCTAAACTGGTAGAAGATCGTTG',
    'GCAGTGCGCGAAGTTTGCGACTTATTACTCCTGGCGCAGGGCAAACTGGATGAAGCCAAAGGGCAATCGATATGA',
    'AGGTGCGTCGCTTGCGACCTGTTACGGCCCTGTCAGCGCCGACGTTATAGCAAAAGCAGAGAACATTCGTCTGCT',
    'GCAGTGCGCGAAGTTTGCGACTTATTACTCCTGGCGCAGGGCAAACTGGATGAAGCCAAAGGGCAATCGATATGA',
    'TCGTTGTGCCACATTGGGGATCACTCACTTGTATCAGGGGCAGTCAAACAAACTGATCGCCTTTAGCGATCTGCT',
    'GCCACATTGGGGATCACTCACTTGTATCAGGGGCAGTCAAACAAACTGATCGCCTTTAGCGATCTGCTGGAAAAA',
    'CGAAGAGCTGAAAGCGTTCAATGTTCGTGACGGTTATGGCATTCGTTGTGCGCTCACCTCTGATATTGAAGTCGC',
    'ATGATCTCATCGACTGGCCGGTAATGGAAAAAGTGGGTTTAAGCGTCGCCGTGGCCGATGCGCATCCACTGTTGA',
    'TATCAGGGGCAGTCAAACAAACTGATCGCCTTTAGCGATCTGCTGGAAAAACTGGCGATTGCCCCGGAAAATGTG',
    'CGTTCAATGTTCGTGACGGTTATGGCATTCGTTGTGCGCTCACCTCTGATATTGAAGTCGCTATCATTACCGGGC',
    'TCAGCGCCGACGTTATAGCAAAAGCAGAGAACATTCGTCTGCTGATCCTCGATGTCGATGGCGTACTGTCAGATG',
    'AAACAAACTGATCGCCTTTAGCGATCTGCTGGAAAAACTGGCGATTGCCCCGGAAAATGTGGCTTATGTCGGCGA']
    reads_726_40_75 = ['AAAATGAAGCCGTTCATCTTCGGTGCGCGTAACAAAGTTCACATCATCAACCTTGAGAAAACTGTACCGATGTTC',
    'CCCTGTACCTGGGCGCTGTTGCTGCAACCGTACGTGAAGGCCGTTCTCAGGATCTGGCTTCCCAGGCGGAAGAAA',
    'AGCTGCGACCAGTTCTTCGTGAACCATCGCTGGCTGGGCGGTATGCTGACTAACTGGAAAACCGTTCGTCAGTCC',
    'TCTCGCAAAGGTAAAATCCTTTTCGTTGGTACTAAACGCGCTGCAAGCGAAGCGGTGAAAGACGCTGCTCTGAGC',
    'CTAACTGGAAAACCGTTCGTCAGTCCATCAAACGTCTGAAAGACCTGGAAACTCAGTCTCAGGACGGTACTTTCG',
    'TCAACGAAGCTCTGGCTGAACTGAACAAGATTGCTTCTCGCAAAGGTAAAATCCTTTTCGTTGGTACTAAACGCG',
    'AAAGACGCTGCTCTGAGCTGCGACCAGTTCTTCGTGAACCATCGCTGGCTGGGCGGTATGCTGACTAACTGGAAA',
    'CAAAGACATGGGCGGTCTGCCGGACGCTCTGTTTGTAATCGATGCTGACCACGAACACATTGCTATCAAAGAAGC',
    'CGAACACATTGCTATCAAAGAAGCAAACAACCTGGGTATTCCGGTATTTGCTATCGTTGATACCAACTCTGATCC',
    'GACTTCGTTATCCCGGGTAACGACGACGCAATCCGTGCTGTGACCCTGTACCTGGGCGCTGTTGCTGCAACCGTA',
    'CCGAAAATGAAGCCGTTCATCTTCGGTGCGCGTAACAAAGTTCACATCATCAACCTTGAGAAAACTGTACCGATG',
    'TTTCGACAAGCTGACCAAGAAAGAAGCGCTGATGCGCACTCGTGAGCTGGAGAAACTGGAAAACAGCCTGGGCGG',
    'CGCTGCAAGCGAAGCGGTGAAAGACGCTGCTCTGAGCTGCGACCAGTTCTTCGTGAACCATCGCTGGCTGGGCGG',
    'GATACCAACTCTGATCCGGACGGTGTTGACTTCGTTATCCCGGGTAACGACGACGCAATCCGTGCTGTGACCCTG',
    'GAAGCTCTGGCTGAACTGAACAAGATTGCTTCTCGCAAAGGTAAAATCCTTTTCGTTGGTACTAAACGCGCTGCA',
    'CGGTGAAAGACGCTGCTCTGAGCTGCGACCAGTTCTTCGTGAACCATCGCTGGCTGGGCGGTATGCTGACTAACT',
    'CCTGGAAACTCAGTCTCAGGACGGTACTTTCGACAAGCTGACCAAGAAAGAAGCGCTGATGCGCACTCGTGAGCT',
    'CAGGACGGTACTTTCGACAAGCTGACCAAGAAAGAAGCGCTGATGCGCACTCGTGAGCTGGAGAAACTGGAAAAC',
    'AATGAAGCCGTTCATCTTCGGTGCGCGTAACAAAGTTCACATCATCAACCTTGAGAAAACTGTACCGATGTTCAA',
    'ATGGCAACTGTTTCCATGCGCGACATGCTCAAGGCTGGTGTTCACTTCGGTCACCAGACCCGTTACTGGAACCCG',
    'TCGTGAACCATCGCTGGCTGGGCGGTATGCTGACTAACTGGAAAACCGTTCGTCAGTCCATCAAACGTCTGAAAG',
    'GCTGCAACCGTACGTGAAGGCCGTTCTCAGGATCTGGCTTCCCAGGCGGAAGAAAGCTTCGTAGAAGCTGAGTAA',
    'TACTGGAACCCGAAAATGAAGCCGTTCATCTTCGGTGCGCGTAACAAAGTTCACATCATCAACCTTGAGAAAACT',
    'ATCAAAGAAGCAAACAACCTGGGTATTCCGGTATTTGCTATCGTTGATACCAACTCTGATCCGGACGGTGTTGAC',
    'GACTTCGTTATCCCGGGTAACGACGACGCAATCCGTGCTGTGACCCTGTACCTGGGCGCTGTTGCTGCAACCGTA',
    'GACCACGAACACATTGCTATCAAAGAAGCAAACAACCTGGGTATTCCGGTATTTGCTATCGTTGATACCAACTCT',
    'GGTAAAATCCTTTTCGTTGGTACTAAACGCGCTGCAAGCGAAGCGGTGAAAGACGCTGCTCTGAGCTGCGACCAG',
    'GCTGCAACCGTACGTGAAGGCCGTTCTCAGGATCTGGCTTCCCAGGCGGAAGAAAGCTTCGTAGAAGCTGAGTAA',
    'TCCTTTTCGTTGGTACTAAACGCGCTGCAAGCGAAGCGGTGAAAGACGCTGCTCTGAGCTGCGACCAGTTCTTCG',
    'TTCCATGCGCGACATGCTCAAGGCTGGTGTTCACTTCGGTCACCAGACCCGTTACTGGAACCCGAAAATGAAGCC',
    'TATTTGCTATCGTTGATACCAACTCTGATCCGGACGGTGTTGACTTCGTTATCCCGGGTAACGACGACGCAATCC',
    'CGTTGATACCAACTCTGATCCGGACGGTGTTGACTTCGTTATCCCGGGTAACGACGACGCAATCCGTGCTGTGAC',
    'CGTCAGTCCATCAAACGTCTGAAAGACCTGGAAACTCAGTCTCAGGACGGTACTTTCGACAAGCTGACCAAGAAA',
    'CTTTCGACAAGCTGACCAAGAAAGAAGCGCTGATGCGCACTCGTGAGCTGGAGAAACTGGAAAACAGCCTGGGCG',
    'TGACTTCGTTATCCCGGGTAACGACGACGCAATCCGTGCTGTGACCCTGTACCTGGGCGCTGTTGCTGCAACCGT',
    'TCAGTCTCAGGACGGTACTTTCGACAAGCTGACCAAGAAAGAAGCGCTGATGCGCACTCGTGAGCTGGAGAAACT',
    'ATCCTTTTCGTTGGTACTAAACGCGCTGCAAGCGAAGCGGTGAAAGACGCTGCTCTGAGCTGCGACCAGTTCTTC',
    'CGTTACTGGAACCCGAAAATGAAGCCGTTCATCTTCGGTGCGCGTAACAAAGTTCACATCATCAACCTTGAGAAA',
    'CACTCGTGAGCTGGAGAAACTGGAAAACAGCCTGGGCGGTATCAAAGACATGGGCGGTCTGCCGGACGCTCTGTT',
    'CGCAAAGGTAAAATCCTTTTCGTTGGTACTAAACGCGCTGCAAGCGAAGCGGTGAAAGACGCTGCTCTGAGCTGC']
    reads_930_50_75 = ['TTAGCGAAAATCACCTTTAACGCACCAACAGTTCCTGTTGTGAATAACGTTGATGTGAAATGCGAAACCAATGGT',
    'TACGTCAGTTGTATAACCCGGTTCAGTGGACGAAGTCTGTTGAGTACATGGCAGCGCAAGGCGTAGAACATCTCT',
    'TGGCCTGACGAAACGCATTGTCGACACCCTGACCGCCTCGGCGCTGAACGAACCTTCAGCGATGGCAGCGGCGCT',
    'CAAGTTCATGCAAGAAGCCGTACCGGAAGGCACGGGCGCTATGGCGGCAATCATCGGTCTGGATGATGCGTCTAT',
    'ATGACGCAATTTGCATTTGTGTTCCCTGGACAGGGTTCTCAAACCGTTGGAATGCTGGCTGATATGGCGGCGAGC',
    'CGGTTCAGTGGACGAAGTCTGTTGAGTACATGGCAGCGCAAGGCGTAGAACATCTCTATGAAGTCGGCCCGGGCA',
    'TTCGCTGATGCGGTGCGTCTGGTTGAGATGCGCGGCAAGTTCATGCAAGAAGCCGTACCGGAAGGCACGGGCGCT',
    'TTAGCGAAAATCACCTTTAACGCACCAACAGTTCCTGTTGTGAATAACGTTGATGTGAAATGCGAAACCAATGGT',
    'GGAATGCTGGCTGATATGGCGGCGAGCTATCCAATTGTCGAAGAAACGTTTGCTGAAGCTTCTGCGGCGCTGGGC',
    'AAACGCATTGTCGACACCCTGACCGCCTCGGCGCTGAACGAACCTTCAGCGATGGCAGCGGCGCTCGAGCTTTAA',
    'TGATTTCGCTGATGCGGTGCGTCTGGTTGAGATGCGCGGCAAGTTCATGCAAGAAGCCGTACCGGAAGGCACGGG',
    'TGGCAAACTCAGCCTGCGCTGTTGACTGCATCTGTTGCGCTGTATCGCGTATGGCAGCAGCAGGGCGGTAAAGCA',
    'TTCCCTGGACAGGGTTCTCAAACCGTTGGAATGCTGGCTGATATGGCGGCGAGCTATCCAATTGTCGAAGAAACG',
    'TCGGCCCGGGCAAAGTGCTTACTGGCCTGACGAAACGCATTGTCGACACCCTGACCGCCTCGGCGCTGAACGAAC',
    'CTGAAGCTTCTGCGGCGCTGGGCTACGACCTGTGGGCGCTGACCCAGCAGGGGCCAGCTGAAGAACTGAATAAAA',
    'AAGCTGCAGAAGGTCAGGTCGTTTCTCCGGTAAACTTTAACTCTCCGGGACAGGTGGTTATTGCCGGTCATAAAG',
    'ACCAGTGAGCGTACCGTCTCACTGTGCGCTGATGAAACCAGCAGCCGACAAACTGGCAGTAGAATTAGCGAAAAT',
    'TGAAACCAGCAGCCGACAAACTGGCAGTAGAATTAGCGAAAATCACCTTTAACGCACCAACAGTTCCTGTTGTGA',
    'TATCCAATTGTCGAAGAAACGTTTGCTGAAGCTTCTGCGGCGCTGGGCTACGACCTGTGGGCGCTGACCCAGCAG',
    'TTTCTCCGGTAAACTTTAACTCTCCGGGACAGGTGGTTATTGCCGGTCATAAAGAAGCGGTTGAGCGTGCTGGCG',
    'GTGCTGGCGCTGCCTGTAAAGCGGCGGGCGCAAAACGCGCGCTGCCGTTACCAGTGAGCGTACCGTCTCACTGTG',
    'CCGGAAGGCACGGGCGCTATGGCGGCAATCATCGGTCTGGATGATGCGTCTATTGCGAAAGCGTGTGAAGAAGCT',
    'AGCGGTTGAGCGTGCTGGCGCTGCCTGTAAAGCGGCGGGCGCAAAACGCGCGCTGCCGTTACCAGTGAGCGTACC',
    'AGCAGCCGACAAACTGGCAGTAGAATTAGCGAAAATCACCTTTAACGCACCAACAGTTCCTGTTGTGAATAACGT',
    'CATGGCAGCGCAAGGCGTAGAACATCTCTATGAAGTCGGCCCGGGCAAAGTGCTTACTGGCCTGACGAAACGCAT',
    'TTACCAGTGAGCGTACCGTCTCACTGTGCGCTGATGAAACCAGCAGCCGACAAACTGGCAGTAGAATTAGCGAAA',
    'AAACTTTAACTCTCCGGGACAGGTGGTTATTGCCGGTCATAAAGAAGCGGTTGAGCGTGCTGGCGCTGCCTGTAA',
    'GCATTTGTGTTCCCTGGACAGGGTTCTCAAACCGTTGGAATGCTGGCTGATATGGCGGCGAGCTATCCAATTGTC',
    'GCGAAACCAATGGTGATGCCATCCGTGACGCACTGGTACGTCAGTTGTATAACCCGGTTCAGTGGACGAAGTCTG',
    'GGCTGATATGGCGGCGAGCTATCCAATTGTCGAAGAAACGTTTGCTGAAGCTTCTGCGGCGCTGGGCTACGACCT',
    'CGGGCGCTATGGCGGCAATCATCGGTCTGGATGATGCGTCTATTGCGAAAGCGTGTGAAGAAGCTGCAGAAGGTC',
    'ACCGTTGGAATGCTGGCTGATATGGCGGCGAGCTATCCAATTGTCGAAGAAACGTTTGCTGAAGCTTCTGCGGCG',
    'ATCCAATTGTCGAAGAAACGTTTGCTGAAGCTTCTGCGGCGCTGGGCTACGACCTGTGGGCGCTGACCCAGCAGG',
    'TGACGCAATTTGCATTTGTGTTCCCTGGACAGGGTTCTCAAACCGTTGGAATGCTGGCTGATATGGCGGCGAGCT',
    'GTGCGTCTGGTTGAGATGCGCGGCAAGTTCATGCAAGAAGCCGTACCGGAAGGCACGGGCGCTATGGCGGCAATC',
    'AAGGCACGGGCGCTATGGCGGCAATCATCGGTCTGGATGATGCGTCTATTGCGAAAGCGTGTGAAGAAGCTGCAG',
    'TGTGAATAACGTTGATGTGAAATGCGAAACCAATGGTGATGCCATCCGTGACGCACTGGTACGTCAGTTGTATAA',
    'TGGCAGTAGAATTAGCGAAAATCACCTTTAACGCACCAACAGTTCCTGTTGTGAATAACGTTGATGTGAAATGCG',
    'TGGCAGCGCAAGGCGTAGAACATCTCTATGAAGTCGGCCCGGGCAAAGTGCTTACTGGCCTGACGAAACGCATTG',
    'AACCAGCAGCCGACAAACTGGCAGTAGAATTAGCGAAAATCACCTTTAACGCACCAACAGTTCCTGTTGTGAATA',
    'ATGCTGGCTGATATGGCGGCGAGCTATCCAATTGTCGAAGAAACGTTTGCTGAAGCTTCTGCGGCGCTGGGCTAC',
    'GAAATGCGAAACCAATGGTGATGCCATCCGTGACGCACTGGTACGTCAGTTGTATAACCCGGTTCAGTGGACGAA',
    'GTTCCCTGGACAGGGTTCTCAAACCGTTGGAATGCTGGCTGATATGGCGGCGAGCTATCCAATTGTCGAAGAAAC',
    'CTCCGGTAAACTTTAACTCTCCGGGACAGGTGGTTATTGCCGGTCATAAAGAAGCGGTTGAGCGTGCTGGCGCTG',
    'TGTGAAGAAGCTGCAGAAGGTCAGGTCGTTTCTCCGGTAAACTTTAACTCTCCGGGACAGGTGGTTATTGCCGGT',
    'CTGTATCGCGTATGGCAGCAGCAGGGCGGTAAAGCACCGGCAATGATGGCCGGTCACAGCCTGGGGGAATACTCC',
    'ATGGCCGGTCACAGCCTGGGGGAATACTCCGCGCTGGTTTGCGCTGGTGTGATTGATTTCGCTGATGCGGTGCGT',
    'TTGCTGAAGCTTCTGCGGCGCTGGGCTACGACCTGTGGGCGCTGACCCAGCAGGGGCCAGCTGAAGAACTGAATA',
    'CGCTGGGCTACGACCTGTGGGCGCTGACCCAGCAGGGGCCAGCTGAAGAACTGAATAAAACCTGGCAAACTCAGC',
    'GAAACGTTTGCTGAAGCTTCTGCGGCGCTGGGCTACGACCTGTGGGCGCTGACCCAGCAGGGGCCAGCTGAAGAA']
    reads_4224_230_75 = ['TTATCGATATCTGGGCTGCGGCGAACGATCGTGTATCCAAAGCGATGATGGATAACCTGCAAACTGAAACCGTGA',
    'TCTGCGTCCGGCACTGAAAATCGTTGATGCTCAGGGTAACGACGTTCTGATCCCAGGTACCGATATGCCAGCGCA',
    'CGTGCACCGACTCTGCACCGTCTGGGTATCCAGGCATTTGAACCGGTACTGATCGAAGGTAAAGCTATCCAGCTG',
    'GTTCTGTTGTATCTTGTGACACCGACTTTGGTGTATGTGCGCACTGCTACGGTCGTGACCTGGCGCGTGGCCACA',
    'CTGCTGGATAACGGTCGTCGCGGTCGTGCGATCACCGGTTCTAACAAGCGTCCTCTGAAATCTTTGGCCGACATG',
    'ATCACCGCGAACTTCCGTGAAGGTCTGAACGTACTCCAGTACTTCATCTCCACCCACGGTGCTCGTAAAGGTCTG',
    'TCGGTTGTGAACTCCAGCGGTAAACTGGTTATCACTTCCCGTAATACTGAACTGAAACTGATCGACGAATTCGGT',
    'AACACGCTGCTGCACGAACAGTGGTGTGACCTGCTGGAAGAGAACTCTGTCGACGCGGTTAAAGTACGTTCTGTT',
    'TAACAAGCGTCCTCTGAAATCTTTGGCCGACATGATCAAAGGTAAACAGGGTCGTTTCCGTCAGAACCTGCTCGG',
    'GTACGTAACGAAAAACGTATGCTGCAGGAAGCGGTAGACGCCCTGCTGGATAACGGTCGTCGCGGTCGTGCGATC',
    'AACGTACCGCAGGTGGTAAAGATCTGCGTCCGGCACTGAAAATCGTTGATGCTCAGGGTAACGACGTTCTGATCC',
    'GCTACAAAGTACCTTACGGTGCGGTACTGGCGAAAGGCGATGGCGAACAGGTTGCTGGCGGCGAAACCGTTGCAA',
    'AGAGCAGTATCTGGACGCGCTGGAAGAGTTCGGTGACGAATTCGACGCGAAGATGGGGGCGGAAGCAATCCAGGC',
    'AGAGCCGGCAATCCTGGCTGAAATCAGCGGTATCGTTTCCTTCGGTAAAGAAACCAAAGGTAAACGTCGTCTGGT',
    'AAACTGGTTATCACTTCCCGTAATACTGAACTGAAACTGATCGACGAATTCGGTCGTACTAAAGAAAGCTACAAA',
    'AGAACGTTATCGTGGGTCGTCTGATCCCGGCAGGTACCGGTTACGCGTACCACCAGGATCGTATGCGTCGCCGTG',
    'AACCTGCTCGGTAAGCGTGTTGACTACTCCGGTCGTTCTGTAATCACCGTAGGTCCATACCTGCGTCTGCATCAG',
    'TCGGTCTGAAACCGACCGTTATTTTTGCGGACCAGATCATGTACACCGGCTTCGCCTATGCAGCGCGTTCTGGTG',
    'AACGATCGTGTATCCAAAGCGATGATGGATAACCTGCAAACTGAAACCGTGATTAACCGTGACGGTCAGGAAGAG',
    'CGACTTCGATGGTGACCAGATGGCTGTTCACGTACCGCTGACGCTGGAAGCCCAGCTGGAAGCGCGTGCGCTGAT',
    'CGCTGCGCGATCGCGTACTGGGTCGTGTAACTGCTGAAGACGTTCTGAAGCCGGGTACTGCTGATATCCTCGTTC',
    'ACTCCGAAACCAAGCGTAAAAAGCTGACCAAGCGTATCAAACTGCTGGAAGCGTTCGTTCAGTCTGGTAACAAAC',
    'ACAGGACGTATACCGTCTGCAGGGCGTTAAGATTAACGATAAACACATCGAAGTTATCGTTCGTCAGATGCTGCG',
    'GTGAAAGATTTATTAAAGTTTCTGAAAGCGCAGACTAAAACCGAAGAGTTTGATGCGATCAAAATTGCTCTGGCT',
    'GTACCGAAAGGTCTGCCTTACTCCATCGTCAACCAGGCGCTGGGTAAAAAAGCAATCTCCAAAATGCTGAACACC',
    'AAAAGATTACGAGTGCCTGTGCGGTAAGTACAAGCGCCTGAAACACCGTGGCGTCATCTGTGAGAAGTGCGGCGT',
    'CGTGGTCTGATGGCGAAGCCGGATGGCTCCATCATCGAAACGCCAATCACCGCGAACTTCCGTGAAGGTCTGAAC',
    'CGGTGCTCGTAAAGGTCTGGCGGATACCGCACTGAAAACTGCGAACTCCGGTTACCTGACTCGTCGTCTGGTTGA',
    'TGGATCTGGAGCAAGAGTGCGAACAGCTGCGTGAAGAGCTGAACGAAACCAACTCCGAAACCAAGCGTAAAAAGC',
    'CTCCAGCGGTAAACTGGTTATCACTTCCCGTAATACTGAACTGAAACTGATCGACGAATTCGGTCGTACTAAAGA',
    'GTTGCGGGCAAACGCGACGAACTGCGCGGCCTGAAAGAGAACGTTATCGTGGGTCGTCTGATCCCGGCAGGTACC',
    'GTGACTGCAGAAGACGCATCTGCCAGCCTGGCAGAACTGCTGAACGCAGGTCTGGGCGGTTCTGATAACGAGTAA',
    'CGTCTGCTGGATCTGGCTGCGCCGGACATCATCGTACGTAACGAAAAACGTATGCTGCAGGAAGCGGTAGACGCC',
    'AAAGCTGCGAAGAAAATGGTTGAGCGCGAAGAAGCTGTCGTTTGGGATATCCTGGACGAAGTTATCCGCGAACAC',
    'AAAGATCTGCGTCCGGCACTGAAAATCGTTGATGCTCAGGGTAACGACGTTCTGATCCCAGGTACCGATATGCCA',
    'TCGTGGGTCGTCTGATCCCGGCAGGTACCGGTTACGCGTACCACCAGGATCGTATGCGTCGCCGTGCTGCGGGTG',
    'TGGGTCGTCTGATCCCGGCAGGTACCGGTTACGCGTACCACCAGGATCGTATGCGTCGCCGTGCTGCGGGTGAAG',
    'AGAGCATGGATCTGGAGCAAGAGTGCGAACAGCTGCGTGAAGAGCTGAACGAAACCAACTCCGAAACCAAGCGTA',
    'AAGGACATCACCGGTGGTCTGCCGCGCGTTGCGGACCTGTTCGAAGCACGTCGTCCGAAAGAGCCGGCAATCCTG',
    'GACCCGCACACCATGCCGGTTATCACCGAAGTAAGCGGTTTTGTACGCTTTACTGACATGATCGACGGCCAGACC',
    'CTGCTACGGTCGTGACCTGGCGCGTGGCCACATCATCAACAAGGGTGAAGCAATCGGTGTTATCGCGGCACAGTC',
    'CCTGGTGGTTACCGAAGACGATTGTGGTACCCATGAAGGTATCATGATGACTCCGGTTATCGAGGGTGGTGACGT',
    'GTTGCTGAAATTCAGGAGCAGTTCCAGTCTGGTCTGGTAACTGCGGGCGAACGCTACAACAAAGTTATCGATATC',
    'TCCGGTTATCGAGGGTGGTGACGTTAAAGAGCCGCTGCGCGATCGCGTACTGGGTCGTGTAACTGCTGAAGACGT',
    'GACATCATCGTACGTAACGAAAAACGTATGCTGCAGGAAGCGGTAGACGCCCTGCTGGATAACGGTCGTCGCGGT',
    'AGCGCAGACTAAAACCGAAGAGTTTGATGCGATCAAAATTGCTCTGGCTTCGCCAGACATGATCCGTTCATGGTC',
    'AGCCCAGCTGGAAGCGCGTGCGCTGATGATGTCTACCAACAACATCCTGTCCCCGGCGAACGGCGAACCAATCAT',
    'CAAGCTGGAACTGCGTGGTCTTGCTACCACCATTAAAGCTGCGAAGAAAATGGTTGAGCGCGAAGAAGCTGTCGT',
    'ACCTGCTGGAAGAGAACTCTGTCGACGCGGTTAAAGTACGTTCTGTTGTATCTTGTGACACCGACTTTGGTGTAT',
    'GTTGACGTGGCGCAGGACCTGGTGGTTACCGAAGACGATTGTGGTACCCATGAAGGTATCATGATGACTCCGGTT',
    'CAACCTGGAACGTCAGCAGATCCTGACTGAAGAGCAGTATCTGGACGCGCTGGAAGAGTTCGGTGACGAATTCGA',
    'TGGATGATTGTACCGAAAGGTCTGCCTTACTCCATCGTCAACCAGGCGCTGGGTAAAAAAGCAATCTCCAAAATG',
    'GGTAAACGTCGTCTGGTTATCACCCCGGTAGACGGTAGCGATCCGTACGAAGAGATGATTCCGAAATGGCGTCAG',
    'AGGTTTCCTTCAACAGCATCTACATGATGGCCGACTCCGGTGCGCGTGGTTCTGCGGCACAGATTCGTCAGCTTG',
    'TCTGGACGCGCTGGAAGAGTTCGGTGACGAATTCGACGCGAAGATGGGGGCGGAAGCAATCCAGGCTCTGCTGAA',
    'TCTGCCGGTACTGCCGCCAGATCTGCGTCCGCTGGTTCCGCTGGATGGTGGTCGTTTCGCGACTTCTGACCTGAA',
    'CACGCTGCTGCACGAACAGTGGTGTGACCTGCTGGAAGAGAACTCTGTCGACGCGGTTAAAGTACGTTCTGTTGT',
    'AAAGTACGTTCTGTTGTATCTTGTGACACCGACTTTGGTGTATGTGCGCACTGCTACGGTCGTGACCTGGCGCGT',
    'GTACTGGGTCTGTACTACATGACCCGTGACTGTGTTAACGCCAAAGGCGAAGGCATGGTGCTGACTGGCCCGAAA',
    'AAAGCTATCCAGCTGCACCCGCTGGTTTGTGCGGCATATAACGCCGACTTCGATGGTGACCAGATGGCTGTTCAC',
    'AAACCAGAGTGGATGATCCTGACCGTTCTGCCGGTACTGCCGCCAGATCTGCGTCCGCTGGTTCCGCTGGATGGT',
    'GTTATCACCGAAGTAAGCGGTTTTGTACGCTTTACTGACATGATCGACGGCCAGACCATTACGCGTCAGACCGAC',
    'CGCTGATGATGTCTACCAACAACATCCTGTCCCCGGCGAACGGCGAACCAATCATCGTTCCGTCTCAGGACGTTG',
    'CGAAAACCAGCCTGAAAGACACGACTGTTGGCCGTGCCATTCTGTGGATGATTGTACCGAAAGGTCTGCCTTACT',
    'TCCGGTTATCGAGGGTGGTGACGTTAAAGAGCCGCTGCGCGATCGCGTACTGGGTCGTGTAACTGCTGAAGACGT',
    'GAAAGGCGATGGCGAACAGGTTGCTGGCGGCGAAACCGTTGCAAACTGGGACCCGCACACCATGCCGGTTATCAC',
    'ATGGCGAAGCCGGATGGCTCCATCATCGAAACGCCAATCACCGCGAACTTCCGTGAAGGTCTGAACGTACTCCAG',
    'CTGGGTCGTGTAACTGCTGAAGACGTTCTGAAGCCGGGTACTGCTGATATCCTCGTTCCGCGCAACACGCTGCTG',
    'GGAAGCGCGTGCGCTGATGATGTCTACCAACAACATCCTGTCCCCGGCGAACGGCGAACCAATCATCGTTCCGTC',
    'AACCTGCTCGGTAAGCGTGTTGACTACTCCGGTCGTTCTGTAATCACCGTAGGTCCATACCTGCGTCTGCATCAG',
    'CATCGGTGAACCGGGTACACAGCTGACCATGCGTACGTTCCACATCGGTGGTGCGGCATCTCGTGCGGCTGCTGA',
    'CGTATGCGTCGCCGTGCTGCGGGTGAAGCTCCGGCTGCACCGCAGGTGACTGCAGAAGACGCATCTGCCAGCCTG',
    'AGCCGGCAATCCTGGCTGAAATCAGCGGTATCGTTTCCTTCGGTAAAGAAACCAAAGGTAAACGTCGTCTGGTTA',
    'TATCCAGCTGCACCCGCTGGTTTGTGCGGCATATAACGCCGACTTCGATGGTGACCAGATGGCTGTTCACGTACC',
    'CGGTCGTGACCTGGCGCGTGGCCACATCATCAACAAGGGTGAAGCAATCGGTGTTATCGCGGCACAGTCCATCGG',
    'ACGTGTAGAACGTGGTGACGTAATTTCCGACGGTCCGGAAGCGCCGCACGACATTCTGCGTCTGCGTGGTGTTCA',
    'CGGTACTGGCGAAAGGCGATGGCGAACAGGTTGCTGGCGGCGAAACCGTTGCAAACTGGGACCCGCACACCATGC',
    'ACATCATCGTACGTAACGAAAAACGTATGCTGCAGGAAGCGGTAGACGCCCTGCTGGATAACGGTCGTCGCGGTC',
    'GCGAACAGGTTGAATACTCTCGCGTCAAGATCGCAAACCGCGAACTGGAAGCGAACGGCAAAGTGGGTGCAACTT',
    'CGGTCGTGACCTGGCGCGTGGCCACATCATCAACAAGGGTGAAGCAATCGGTGTTATCGCGGCACAGTCCATCGG',
    'GACTGCGCACATCTGGTTCCTGAAATCGCTGCCGTCCCGTATCGGTCTGCTGCTCGATATGCCGCTGCGCGATAT',
    'AACACCGTGGCGTCATCTGTGAGAAGTGCGGCGTTGAAGTGACCCAGACTAAAGTACGCCGTGAGCGTATGGGCC',
    'CAGCTGACCATGCGTACGTTCCACATCGGTGGTGCGGCATCTCGTGCGGCTGCTGAATCCAGCATCCAAGTGAAA',
    'GCACGAACAGTGGTGTGACCTGCTGGAAGAGAACTCTGTCGACGCGGTTAAAGTACGTTCTGTTGTATCTTGTGA',
    'ACCATCAACTACCGTACGTTCAAACCAGAACGTGACGGCCTTTTCTGCGCCCGTATCTTTGGGCCGGTAAAAGAT',
    'TATCGATATCTGGGCTGCGGCGAACGATCGTGTATCCAAAGCGATGATGGATAACCTGCAAACTGAAACCGTGAT',
    'CTGGTTCCTGAAATCGCTGCCGTCCCGTATCGGTCTGCTGCTCGATATGCCGCTGCGCGATATCGAACGCGTACT',
    'TTTCCGACGGTCCGGAAGCGCCGCACGACATTCTGCGTCTGCGTGGTGTTCATGCTGTTACTCGTTACATCGTTA',
    'ACAACAAAGTTATCGATATCTGGGCTGCGGCGAACGATCGTGTATCCAAAGCGATGATGGATAACCTGCAAACTG',
    'GTAAAGCTATCCAGCTGCACCCGCTGGTTTGTGCGGCATATAACGCCGACTTCGATGGTGACCAGATGGCTGTTC',
    'GCGGCCTGAAAGAGAACGTTATCGTGGGTCGTCTGATCCCGGCAGGTACCGGTTACGCGTACCACCAGGATCGTA',
    'ACGGTCCGGAAGCGCCGCACGACATTCTGCGTCTGCGTGGTGTTCATGCTGTTACTCGTTACATCGTTAACGAAG',
    'AAGCAATCCAGGCTCTGCTGAAGAGCATGGATCTGGAGCAAGAGTGCGAACAGCTGCGTGAAGAGCTGAACGAAA',
    'ACGAAAAACGTATGCTGCAGGAAGCGGTAGACGCCCTGCTGGATAACGGTCGTCGCGGTCGTGCGATCACCGGTT',
    'GTGGTGCGGCATCTCGTGCGGCTGCTGAATCCAGCATCCAAGTGAAAAACAAAGGTAGCATCAAGCTCAGCAACG',
    'GCCTGAAACACCGTGGCGTCATCTGTGAGAAGTGCGGCGTTGAAGTGACCCAGACTAAAGTACGCCGTGAGCGTA',
    'CTGAAACCGACCGTTATTTTTGCGGACCAGATCATGTACACCGGCTTCGCCTATGCAGCGCGTTCTGGTGCATCT',
    'ACTCCGGTTATCGAGGGTGGTGACGTTAAAGAGCCGCTGCGCGATCGCGTACTGGGTCGTGTAACTGCTGAAGAC',
    'GACACCCTGGCGCGTATTCCGCAGGAATCCGGCGGTACCAAGGACATCACCGGTGGTCTGCCGCGCGTTGCGGAC',
    'GGCTCTGCTGAAGAGCATGGATCTGGAGCAAGAGTGCGAACAGCTGCGTGAAGAGCTGAACGAAACCAACTCCGA',
    'CAGTCCATCGGTGAACCGGGTACACAGCTGACCATGCGTACGTTCCACATCGGTGGTGCGGCATCTCGTGCGGCT',
    'AGGAAGCGGTAGACGCCCTGCTGGATAACGGTCGTCGCGGTCGTGCGATCACCGGTTCTAACAAGCGTCCTCTGA',
    'CACCCGCTGGTTTGTGCGGCATATAACGCCGACTTCGATGGTGACCAGATGGCTGTTCACGTACCGCTGACGCTG',
    'ACATGATCAAAGGTAAACAGGGTCGTTTCCGTCAGAACCTGCTCGGTAAGCGTGTTGACTACTCCGGTCGTTCTG',
    'AGATCAGCTCTGGTGACACCCTGGCGCGTATTCCGCAGGAATCCGGCGGTACCAAGGACATCACCGGTGGTCTGC',
    'ACCTGCAAACTGAAACCGTGATTAACCGTGACGGTCAGGAAGAGAAGCAGGTTTCCTTCAACAGCATCTACATGA',
    'AAGAGTTCGGTGACGAATTCGACGCGAAGATGGGGGCGGAAGCAATCCAGGCTCTGCTGAAGAGCATGGATCTGG',
    'TCAACAGCATCTACATGATGGCCGACTCCGGTGCGCGTGGTTCTGCGGCACAGATTCGTCAGCTTGCTGGTATGC',
    'TGCTGCTCGATATGCCGCTGCGCGATATCGAACGCGTACTGTACTTTGAATCCTATGTGGTTATCGAAGGCGGTA',
    'CTGCGTGGTGTTCATGCTGTTACTCGTTACATCGTTAACGAAGTACAGGACGTATACCGTCTGCAGGGCGTTAAG',
    'CGATCAAAATTGCTCTGGCTTCGCCAGACATGATCCGTTCATGGTCTTTCGGTGAAGTTAAAAAGCCGGAAACCA',
    'GTCTGCCTTACTCCATCGTCAACCAGGCGCTGGGTAAAAAAGCAATCTCCAAAATGCTGAACACCTGCTACCGCA',
    'CTCCATCGTCAACCAGGCGCTGGGTAAAAAAGCAATCTCCAAAATGCTGAACACCTGCTACCGCATTCTCGGTCT',
    'TGGGACCCGCACACCATGCCGGTTATCACCGAAGTAAGCGGTTTTGTACGCTTTACTGACATGATCGACGGCCAG',
    'CCGTCTGGGTATCCAGGCATTTGAACCGGTACTGATCGAAGGTAAAGCTATCCAGCTGCACCCGCTGGTTTGTGC',
    'TCTGGTCTGGTAACTGCGGGCGAACGCTACAACAAAGTTATCGATATCTGGGCTGCGGCGAACGATCGTGTATCC',
    'CAGCTCTGGTGACACCCTGGCGCGTATTCCGCAGGAATCCGGCGGTACCAAGGACATCACCGGTGGTCTGCCGCG',
    'CTGCCTTACTCCATCGTCAACCAGGCGCTGGGTAAAAAAGCAATCTCCAAAATGCTGAACACCTGCTACCGCATT',
    'GGCGAACCAATCATCGTTCCGTCTCAGGACGTTGTACTGGGTCTGTACTACATGACCCGTGACTGTGTTAACGCC',
    'CACTCGCGTGCTGACCGAAGCAGCCGTTGCGGGCAAACGCGACGAACTGCGCGGCCTGAAAGAGAACGTTATCGT',
    'GAACTTCCGTGAAGGTCTGAACGTACTCCAGTACTTCATCTCCACCCACGGTGCTCGTAAAGGTCTGGCGGATAC',
    'TACCGTACGTTCAAACCAGAACGTGACGGCCTTTTCTGCGCCCGTATCTTTGGGCCGGTAAAAGATTACGAGTGC',
    'TAAAAGATTACGAGTGCCTGTGCGGTAAGTACAAGCGCCTGAAACACCGTGGCGTCATCTGTGAGAAGTGCGGCG',
    'ATGCTCAGGGTAACGACGTTCTGATCCCAGGTACCGATATGCCAGCGCAGTACTTCCTGCCGGGTAAAGCGATTG',
    'AATCCGGCGGTACCAAGGACATCACCGGTGGTCTGCCGCGCGTTGCGGACCTGTTCGAAGCACGTCGTCCGAAAG',
    'GAGATGATTCCGAAATGGCGTCAGCTCAACGTGTTCGAAGGTGAACGTGTAGAACGTGGTGACGTAATTTCCGAC',
    'GGCGAAACCGTTGCAAACTGGGACCCGCACACCATGCCGGTTATCACCGAAGTAAGCGGTTTTGTACGCTTTACT',
    'GCACGACATTCTGCGTCTGCGTGGTGTTCATGCTGTTACTCGTTACATCGTTAACGAAGTACAGGACGTATACCG',
    'GGTCTTTCGGTGAAGTTAAAAAGCCGGAAACCATCAACTACCGTACGTTCAAACCAGAACGTGACGGCCTTTTCT',
    'TTCAACAGCATCTACATGATGGCCGACTCCGGTGCGCGTGGTTCTGCGGCACAGATTCGTCAGCTTGCTGGTATG',
    'AACTGGAAGCGAACGGCAAAGTGGGTGCAACTTACTCCCGCGATCTGCTGGGTATCACCAAAGCGTCTCTGGCAA',
    'CCTATGCAGCGCGTTCTGGTGCATCTGTTGGTATCGATGACATGGTCATCCCGGAGAAGAAACACGAAATCATCT',
    'AAAGGTAAACGTCGTCTGGTTATCACCCCGGTAGACGGTAGCGATCCGTACGAAGAGATGATTCCGAAATGGCGT',
    'CAGCTCAACGTGTTCGAAGGTGAACGTGTAGAACGTGGTGACGTAATTTCCGACGGTCCGGAAGCGCCGCACGAC',
    'TGTCCCCGGCGAACGGCGAACCAATCATCGTTCCGTCTCAGGACGTTGTACTGGGTCTGTACTACATGACCCGTG',
    'CCTTACGGTGCGGTACTGGCGAAAGGCGATGGCGAACAGGTTGCTGGCGGCGAAACCGTTGCAAACTGGGACCCG',
    'CCTGCGTCTGCATCAGTGCGGTCTGCCGAAGAAAATGGCACTGGAGCTGTTCAAACCGTTCATCTACGGCAAGCT',
    'AAGCGGTTTTGTACGCTTTACTGACATGATCGACGGCCAGACCATTACGCGTCAGACCGACGAACTGACCGGTCT',
    'TCTGTATCGCTCTGGTCTGGCTTCTCTGCATGCGCGCGTTAAAGTGCGTATCACCGAGTATGAAAAAGATGCTAA',
    'AGTCCTTCATCTCCGCGGCATCGTTCCAGGAGACCACTCGCGTGCTGACCGAAGCAGCCGTTGCGGGCAAACGCG',
    'GATCAAAATTGCTCTGGCTTCGCCAGACATGATCCGTTCATGGTCTTTCGGTGAAGTTAAAAAGCCGGAAACCAT',
    'CTGTCGTTTGGGATATCCTGGACGAAGTTATCCGCGAACACCCGGTACTGCTGAACCGTGCACCGACTCTGCACC',
    'TATCGTCGCGTCATTAACCGTAACAACCGTCTGAAACGTCTGCTGGATCTGGCTGCGCCGGACATCATCGTACGT',
    'TGCTAACGGTGAATTAGTAGCGAAAACCAGCCTGAAAGACACGACTGTTGGCCGTGCCATTCTGTGGATGATTGT',
    'TCGTTCCAGGAGACCACTCGCGTGCTGACCGAAGCAGCCGTTGCGGGCAAACGCGACGAACTGCGCGGCCTGAAA',
    'TAGTAGCGAAAACCAGCCTGAAAGACACGACTGTTGGCCGTGCCATTCTGTGGATGATTGTACCGAAAGGTCTGC',
    'GTACTGCCGCCAGATCTGCGTCCGCTGGTTCCGCTGGATGGTGGTCGTTTCGCGACTTCTGACCTGAACGATCTG',
    'ACTGTTGGCCGTGCCATTCTGTGGATGATTGTACCGAAAGGTCTGCCTTACTCCATCGTCAACCAGGCGCTGGGT',
    'AGCCGCTGCGCGATCGCGTACTGGGTCGTGTAACTGCTGAAGACGTTCTGAAGCCGGGTACTGCTGATATCCTCG',
    'GGCCTGAAAGAGAACGTTATCGTGGGTCGTCTGATCCCGGCAGGTACCGGTTACGCGTACCACCAGGATCGTATG',
    'ACGCCGTGAGCGTATGGGCCACATCGAACTGGCTTCCCCGACTGCGCACATCTGGTTCCTGAAATCGCTGCCGTC',
    'TACGCGTACCACCAGGATCGTATGCGTCGCCGTGCTGCGGGTGAAGCTCCGGCTGCACCGCAGGTGACTGCAGAA',
    'CCGCTGGTTTGTGCGGCATATAACGCCGACTTCGATGGTGACCAGATGGCTGTTCACGTACCGCTGACGCTGGAA',
    'GTGTTGACTACTCCGGTCGTTCTGTAATCACCGTAGGTCCATACCTGCGTCTGCATCAGTGCGGTCTGCCGAAGA',
    'CTGACCGAAGCAGCCGTTGCGGGCAAACGCGACGAACTGCGCGGCCTGAAAGAGAACGTTATCGTGGGTCGTCTG',
    'GGGTGAAGCTCCGGCTGCACCGCAGGTGACTGCAGAAGACGCATCTGCCAGCCTGGCAGAACTGCTGAACGCAGG',
    'TTCTGTGGATGATTGTACCGAAAGGTCTGCCTTACTCCATCGTCAACCAGGCGCTGGGTAAAAAAGCAATCTCCA',
    'TACTGTACTTTGAATCCTATGTGGTTATCGAAGGCGGTATGACCAACCTGGAACGTCAGCAGATCCTGACTGAAG',
    'GTTGGCCGTGCCATTCTGTGGATGATTGTACCGAAAGGTCTGCCTTACTCCATCGTCAACCAGGCGCTGGGTAAA',
    'GGGTAGCTCCGACTTCCTGGAAGGCGAACAGGTTGAATACTCTCGCGTCAAGATCGCAAACCGCGAACTGGAAGC',
    'CTGGTTATCACCCCGGTAGACGGTAGCGATCCGTACGAAGAGATGATTCCGAAATGGCGTCAGCTCAACGTGTTC',
    'GACCATTACGCGTCAGACCGACGAACTGACCGGTCTGTCTTCGCTGGTGGTTCTGGATTCCGCAGAACGTACCGC',
    'GTACTTTGAATCCTATGTGGTTATCGAAGGCGGTATGACCAACCTGGAACGTCAGCAGATCCTGACTGAAGAGCA',
    'GCCGCACGACATTCTGCGTCTGCGTGGTGTTCATGCTGTTACTCGTTACATCGTTAACGAAGTACAGGACGTATA',
    'CTGAAGAGCAGTATCTGGACGCGCTGGAAGAGTTCGGTGACGAATTCGACGCGAAGATGGGGGCGGAAGCAATCC',
    'TGGATCTGGCTGCGCCGGACATCATCGTACGTAACGAAAAACGTATGCTGCAGGAAGCGGTAGACGCCCTGCTGG',
    'ACCATGCGTACGTTCCACATCGGTGGTGCGGCATCTCGTGCGGCTGCTGAATCCAGCATCCAAGTGAAAAACAAA',
    'ACAAAGTACCTTACGGTGCGGTACTGGCGAAAGGCGATGGCGAACAGGTTGCTGGCGGCGAAACCGTTGCAAACT',
    'AGTGAAAAACAAAGGTAGCATCAAGCTCAGCAACGTGAAGTCGGTTGTGAACTCCAGCGGTAAACTGGTTATCAC',
    'TGAACGAAACCAACTCCGAAACCAAGCGTAAAAAGCTGACCAAGCGTATCAAACTGCTGGAAGCGTTCGTTCAGT',
    'AGAAGTGCGGCGTTGAAGTGACCCAGACTAAAGTACGCCGTGAGCGTATGGGCCACATCGAACTGGCTTCCCCGA',
    'AACGTCGTCTGGTTATCACCCCGGTAGACGGTAGCGATCCGTACGAAGAGATGATTCCGAAATGGCGTCAGCTCA',
    'GCTGACGCTGGAAGCCCAGCTGGAAGCGCGTGCGCTGATGATGTCTACCAACAACATCCTGTCCCCGGCGAACGG',
    'ATGGTGCTGACTGGCCCGAAAGAAGCAGAACGTCTGTATCGCTCTGGTCTGGCTTCTCTGCATGCGCGCGTTAAA',
    'TATCGATATCTGGGCTGCGGCGAACGATCGTGTATCCAAAGCGATGATGGATAACCTGCAAACTGAAACCGTGAT',
    'CGTGCACCGACTCTGCACCGTCTGGGTATCCAGGCATTTGAACCGGTACTGATCGAAGGTAAAGCTATCCAGCTG',
    'ACTTTGAATCCTATGTGGTTATCGAAGGCGGTATGACCAACCTGGAACGTCAGCAGATCCTGACTGAAGAGCAGT',
    'GTGACACCCTGGCGCGTATTCCGCAGGAATCCGGCGGTACCAAGGACATCACCGGTGGTCTGCCGCGCGTTGCGG',
    'GACTAAAGTACGCCGTGAGCGTATGGGCCACATCGAACTGGCTTCCCCGACTGCGCACATCTGGTTCCTGAAATC',
    'GTGGTGTGACCTGCTGGAAGAGAACTCTGTCGACGCGGTTAAAGTACGTTCTGTTGTATCTTGTGACACCGACTT',
    'TGCATGCGCGCGTTAAAGTGCGTATCACCGAGTATGAAAAAGATGCTAACGGTGAATTAGTAGCGAAAACCAGCC',
    'GGTTCTGGATTCCGCAGAACGTACCGCAGGTGGTAAAGATCTGCGTCCGGCACTGAAAATCGTTGATGCTCAGGG',
    'ATATCCTCGTTCCGCGCAACACGCTGCTGCACGAACAGTGGTGTGACCTGCTGGAAGAGAACTCTGTCGACGCGG',
    'TAAACGTCGTCTGGTTATCACCCCGGTAGACGGTAGCGATCCGTACGAAGAGATGATTCCGAAATGGCGTCAGCT',
    'AAACCAGAGTGGATGATCCTGACCGTTCTGCCGGTACTGCCGCCAGATCTGCGTCCGCTGGTTCCGCTGGATGGT',
    'GAACTGGCTTCCCCGACTGCGCACATCTGGTTCCTGAAATCGCTGCCGTCCCGTATCGGTCTGCTGCTCGATATG',
    'AAGCGATTGTTCAGCTGGAAGATGGCGTACAGATCAGCTCTGGTGACACCCTGGCGCGTATTCCGCAGGAATCCG',
    'AATCCGGCGGTACCAAGGACATCACCGGTGGTCTGCCGCGCGTTGCGGACCTGTTCGAAGCACGTCGTCCGAAAG',
    'AACTGCGTGGTCTTGCTACCACCATTAAAGCTGCGAAGAAAATGGTTGAGCGCGAAGAAGCTGTCGTTTGGGATA',
    'TCGCGTGCTGACCGAAGCAGCCGTTGCGGGCAAACGCGACGAACTGCGCGGCCTGAAAGAGAACGTTATCGTGGG',
    'ACAAGCGTCCTCTGAAATCTTTGGCCGACATGATCAAAGGTAAACAGGGTCGTTTCCGTCAGAACCTGCTCGGTA',
    'CATGGTCATCCCGGAGAAGAAACACGAAATCATCTCCGAGGCAGAAGCAGAAGTTGCTGAAATTCAGGAGCAGTT',
    'CATCGTTCCAGGAGACCACTCGCGTGCTGACCGAAGCAGCCGTTGCGGGCAAACGCGACGAACTGCGCGGCCTGA',
    'AACCGGGTACACAGCTGACCATGCGTACGTTCCACATCGGTGGTGCGGCATCTCGTGCGGCTGCTGAATCCAGCA',
    'TGGCGCGTGGCCACATCATCAACAAGGGTGAAGCAATCGGTGTTATCGCGGCACAGTCCATCGGTGAACCGGGTA',
    'TCCTGGACGAAGTTATCCGCGAACACCCGGTACTGCTGAACCGTGCACCGACTCTGCACCGTCTGGGTATCCAGG',
    'TGGAAGCGCGTGCGCTGATGATGTCTACCAACAACATCCTGTCCCCGGCGAACGGCGAACCAATCATCGTTCCGT',
    'ACAAAGTACCTTACGGTGCGGTACTGGCGAAAGGCGATGGCGAACAGGTTGCTGGCGGCGAAACCGTTGCAAACT',
    'ATATCTGGGCTGCGGCGAACGATCGTGTATCCAAAGCGATGATGGATAACCTGCAAACTGAAACCGTGATTAACC',
    'TCGCTCTGGTCTGGCTTCTCTGCATGCGCGCGTTAAAGTGCGTATCACCGAGTATGAAAAAGATGCTAACGGTGA',
    'GTAAAGATCTGCGTCCGGCACTGAAAATCGTTGATGCTCAGGGTAACGACGTTCTGATCCCAGGTACCGATATGC',
    'AAAGCGTCTCTGGCAACCGAGTCCTTCATCTCCGCGGCATCGTTCCAGGAGACCACTCGCGTGCTGACCGAAGCA',
    'GTGACCTGCTGGAAGAGAACTCTGTCGACGCGGTTAAAGTACGTTCTGTTGTATCTTGTGACACCGACTTTGGTG',
    'TGATTCCGAAATGGCGTCAGCTCAACGTGTTCGAAGGTGAACGTGTAGAACGTGGTGACGTAATTTCCGACGGTC',
    'GCCTGAAAGACACGACTGTTGGCCGTGCCATTCTGTGGATGATTGTACCGAAAGGTCTGCCTTACTCCATCGTCA',
    'CAGCTGGAAGATGGCGTACAGATCAGCTCTGGTGACACCCTGGCGCGTATTCCGCAGGAATCCGGCGGTACCAAG',
    'TCGTTCAGTCTGGTAACAAACCAGAGTGGATGATCCTGACCGTTCTGCCGGTACTGCCGCCAGATCTGCGTCCGC',
    'CTGCAGGAAGCGGTAGACGCCCTGCTGGATAACGGTCGTCGCGGTCGTGCGATCACCGGTTCTAACAAGCGTCCT',
    'CTCCGGTTATCGAGGGTGGTGACGTTAAAGAGCCGCTGCGCGATCGCGTACTGGGTCGTGTAACTGCTGAAGACG',
    'CCCGACTGCGCACATCTGGTTCCTGAAATCGCTGCCGTCCCGTATCGGTCTGCTGCTCGATATGCCGCTGCGCGA',
    'ATCCTGACTGAAGAGCAGTATCTGGACGCGCTGGAAGAGTTCGGTGACGAATTCGACGCGAAGATGGGGGCGGAA',
    'CCAGAGTGGATGATCCTGACCGTTCTGCCGGTACTGCCGCCAGATCTGCGTCCGCTGGTTCCGCTGGATGGTGGT',
    'TATGAAAAAGATGCTAACGGTGAATTAGTAGCGAAAACCAGCCTGAAAGACACGACTGTTGGCCGTGCCATTCTG',
    'GCCCTGCTGGATAACGGTCGTCGCGGTCGTGCGATCACCGGTTCTAACAAGCGTCCTCTGAAATCTTTGGCCGAC',
    'CGCGTTAAAGTGCGTATCACCGAGTATGAAAAAGATGCTAACGGTGAATTAGTAGCGAAAACCAGCCTGAAAGAC',
    'TCTGTTGGTATCGATGACATGGTCATCCCGGAGAAGAAACACGAAATCATCTCCGAGGCAGAAGCAGAAGTTGCT',
    'ACCCAGACTAAAGTACGCCGTGAGCGTATGGGCCACATCGAACTGGCTTCCCCGACTGCGCACATCTGGTTCCTG',
    'TGACACCGACTTTGGTGTATGTGCGCACTGCTACGGTCGTGACCTGGCGCGTGGCCACATCATCAACAAGGGTGA',
    'GACATGATCGACGGCCAGACCATTACGCGTCAGACCGACGAACTGACCGGTCTGTCTTCGCTGGTGGTTCTGGAT',
    'CTGGCGCGTGGCCACATCATCAACAAGGGTGAAGCAATCGGTGTTATCGCGGCACAGTCCATCGGTGAACCGGGT',
    'GTGCCTGTGCGGTAAGTACAAGCGCCTGAAACACCGTGGCGTCATCTGTGAGAAGTGCGGCGTTGAAGTGACCCA',
    'TACCGAAGACGATTGTGGTACCCATGAAGGTATCATGATGACTCCGGTTATCGAGGGTGGTGACGTTAAAGAGCC',
    'AAGCGTGTTGACTACTCCGGTCGTTCTGTAATCACCGTAGGTCCATACCTGCGTCTGCATCAGTGCGGTCTGCCG',
    'TCGTCTGGTTATCACCCCGGTAGACGGTAGCGATCCGTACGAAGAGATGATTCCGAAATGGCGTCAGCTCAACGT',
    'ACACATCGAAGTTATCGTTCGTCAGATGCTGCGTAAAGCTACCATCGTTAACGCGGGTAGCTCCGACTTCCTGGA',
    'GGTAAAGATCTGCGTCCGGCACTGAAAATCGTTGATGCTCAGGGTAACGACGTTCTGATCCCAGGTACCGATATG',
    'GTTATCGCGGCACAGTCCATCGGTGAACCGGGTACACAGCTGACCATGCGTACGTTCCACATCGGTGGTGCGGCA',
    'GGCGGAAGCAATCCAGGCTCTGCTGAAGAGCATGGATCTGGAGCAAGAGTGCGAACAGCTGCGTGAAGAGCTGAA',
    'GCTGGGTAAAAAAGCAATCTCCAAAATGCTGAACACCTGCTACCGCATTCTCGGTCTGAAACCGACCGTTATTTT',
    'GAGTGCCTGTGCGGTAAGTACAAGCGCCTGAAACACCGTGGCGTCATCTGTGAGAAGTGCGGCGTTGAAGTGACC']

    dataset[1] = (genome25, reads_25_10_8) #va
    dataset[2] = (genome25, reads_25_10_10)
    dataset[3] = (genome25, reads_25_10_15)
    dataset[4] = (genome50, reads_50_10_8) #va
    dataset[5] = (genome50, reads_50_10_10)
    dataset[6] = (genome50, reads_50_10_15)
    dataset[7] = (genome25, reads_25_20_8)
    dataset[8] = (genome25, reads_25_20_10)
    dataset[9] = (genome25, reads_25_20_15)
    dataset[10] = (genome50, reads_50_20_8)
    dataset[11] = (genome50, reads_50_20_10)
    dataset[12] = (genome50, reads_50_20_15)
    dataset[13] = (genome25, reads_25_30_8)
    dataset[14] = (genome25, reads_25_30_10)
    dataset[15] = (genome25, reads_25_30_15)
    dataset[16] = (genome50, reads_50_30_8)
    dataset[17] = (genome50, reads_50_30_10)
    dataset[18] = (genome50, reads_50_30_15)
    dataset[19] = (genome381, reads_381_20_75) #va
    dataset[20] = (genome567, reads_567_30_75)
    dataset[21] = (genome726, reads_726_40_75)
    dataset[22] = (genome930, reads_930_50_75)
    dataset[23] = (genome4224, reads_4224_230_75)

    genome, reads = dataset[int(sys.argv[2])]
    qlearning(reads, int(sys.argv[1]), genome)
    


