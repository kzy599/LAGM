import numpy as np
import json
import argparse
import random
from deap import base, creator, tools, algorithms
import math

def load_input_data(file_path):
    """
    Loads and validates input data from a JSON file.

    Parameters:
        file_path (str): Path to the input JSON file.

    Returns:
        tuple: A tuple containing individual IDs, EBV vector, kinship matrix,
               female IDs, and male IDs.
    """
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except Exception as e:
        raise IOError(f"Error reading input file: {e}")

    # Extract data
    individual_ids = data.get("individual_ids")
    ebv_vector = data.get("ebv_vector")
    kinship_matrix = data.get("kinship_matrix")
    female_ids = data.get("female_ids")
    male_ids = data.get("male_ids")
    geno_matrix = data.get("geno_matrix")
    gen = data.get("gen")
    He = data.get("He")
    A_matrix = data.get("A_matrix")
    gen_pre = data.get("gen_pre")
    # Validate data
    if not (individual_ids and ebv_vector and kinship_matrix and female_ids and male_ids):
        raise ValueError("Input JSON file is missing required fields.")

    total_individuals = len(individual_ids)

    if len(ebv_vector) != total_individuals:
        raise ValueError("Length of EBV vector does not match number of individual IDs.")

    if len(kinship_matrix) != total_individuals:
        raise ValueError("Kinship matrix size does not match number of individuals.")

    for row in kinship_matrix:
        if len(row) != total_individuals:
            raise ValueError("Each row of the kinship matrix must have a length equal to the number of individuals.")

    if isinstance(ebv_vector, dict):
        ebv_vector = list(ebv_vector.values())
    # Convert to numpy arrays
    ebv_vector = np.array(ebv_vector, dtype=np.float32)
    kinship_matrix = np.array(kinship_matrix, dtype=np.float32)
    geno_matrix = np.array(geno_matrix, dtype=np.float32)
    if A_matrix is not None:
        A_matrix = np.array(A_matrix, dtype=np.float32)
    # print("Type of ebv_vector after loading:", type(ebv_vector))
    # print("Content of ebv_vector:", ebv_vector)

    # Check if kinship matrix is symmetric
    if not np.allclose(kinship_matrix, kinship_matrix.T, atol=1e-5):
        raise ValueError("Kinship matrix is not symmetric.")

    # Check if female and male IDs are subsets of individual IDs
    individual_set = set(individual_ids)
    if not set(female_ids).issubset(individual_set):
        raise ValueError("Some female IDs are not present in individual IDs.")
    if not set(male_ids).issubset(individual_set):
        raise ValueError("Some male IDs are not present in individual IDs.")

    return individual_ids, ebv_vector, kinship_matrix, female_ids, male_ids, geno_matrix,gen,He,A_matrix,gen_pre


def calHe(geno_arr: np.ndarray) -> float:
    # 1. 输入验证
    if not isinstance(geno_arr, np.ndarray) or geno_arr.ndim != 2:
        raise ValueError("输入必须是二维NumPy数组")
    n_samples, n_loci = geno_arr.shape
    if n_samples == 0:
        raise ValueError("样本数不能为0")
    
    # 2. 计算每个位点的杂合子数量
    hetero_counts = np.sum(geno_arr == 1, axis=0)  # axis=0正确：按列统计
    
    # 3. 计算各点位杂合率
    he_per_locus = hetero_counts / n_samples
    
    # 4. 返回全局均值
    return np.mean(he_per_locus)

def repair_individual(individual, female_ids, male_ids, female_num, male_num, female_mates, male_mates):
    """
    Repairs an individual to ensure it satisfies all constraints.
    """
    # Count the number of times each female and male is used
    female_counts = {female: 0 for female in female_ids}
    male_counts = {male: 0 for male in male_ids}

    # Filter valid pairs
    valid_pairs = []
    for female, male in individual:
        if female_counts[female] < female_mates and male_counts[male] < male_mates:
            valid_pairs.append((female, male))
            female_counts[female] += 1
            male_counts[male] += 1

    # Ensure we have exactly female_num females and male_num males
    selected_females = list(set(pair[0] for pair in valid_pairs))
    selected_males = list(set(pair[1] for pair in valid_pairs))

    while True:
        
        # 终止条件：数量精确匹配
        if len(selected_females) == female_num and len(selected_males) == male_num:
            break
            
         # Add more females if needed
        if len(selected_females) < female_num:
            available_females = [f for f in female_ids if female_counts[f] < female_mates]
            needed_females = female_num - len(selected_females)
            selected_females.extend(np.random.choice(available_females, size=needed_females, replace=False))

        elif len(selected_females) > female_num:
            remove_num = len(selected_females) - female_num
            remove_females = np.random.choice(selected_females, remove_num, replace=False)
            
            # 找出所有涉及被删雌性的配对，并更新雄性计数
            removed_pairs = [ (f,m) for f,m in valid_pairs if f in remove_females ]
            for f, m in removed_pairs:
                male_counts[m] -= 1  # 关键修改：减少相关雄性的计数
            
            # 更新valid_pairs和雌性列表
            valid_pairs = [ (f,m) for f,m in valid_pairs if f not in remove_females ]
            selected_females = [f for f in selected_females if f not in remove_females]
            # 重置被删雌性的计数
            for f in remove_females:
                female_counts[f] = 0

        # Add more males if needed
        if len(selected_males) < male_num:
            available_males = [m for m in male_ids if male_counts[m] < male_mates]
            needed_males = male_num - len(selected_males)
            selected_males.extend(np.random.choice(available_males, size=needed_males, replace=False))

        elif len(selected_males) > male_num:
            remove_num = len(selected_males) - male_num
            remove_males = np.random.choice(selected_males, remove_num, replace=False)
            
            # 找出所有涉及被删雄性的配对，并更新雌性计数
            removed_pairs = [ (f,m) for f,m in valid_pairs if m in remove_males ]
            for f, m in removed_pairs:
                female_counts[f] -= 1  # 关键修改：减少相关雌性的计数
            
            # 更新valid_pairs和雄性列表
            valid_pairs = [ (f,m) for f,m in valid_pairs if m not in remove_males ]
            selected_males = [m for m in selected_males if m not in remove_males]
            
            # 重置被删雄性的计数
            for m in remove_males:
                male_counts[m] = 0
    
    # 阶段3：生成最终配对
    new_pairs = []
    f_counts = {f:0 for f in selected_females}
    m_counts = {m:0 for m in selected_males}
    
    # 优先保留原始有效配对
    for f, m in valid_pairs:
        if f in selected_females and m in selected_males:
            new_pairs.append((f, m))
            f_counts[f] += 1
            m_counts[m] += 1
    
    # 补充新配对
    for f in selected_females:
        for m in selected_males:
            if f_counts[f] < female_mates and m_counts[m] < male_mates:
                if (f, m) not in new_pairs:
                    new_pairs.append((f, m))
                    f_counts[f] += 1
                    m_counts[m] += 1

    # Create an individual with fitness attribute
    individual = creator.Individual(new_pairs)
    individual.fitness = creator.FitnessMax()

    # Debugging: Check individual type and fitness
    # print("Type of repaired individual:", type(individual))
    # print("Fitness of repaired individual:", individual.fitness)

    return individual

def normalize(value, min_value, max_value):
    return (value - min_value) / (max_value - min_value + 1e-8)


def genetic_gain(breeding_pairs, ebv_vector, id_to_index, geno_matrix, gen, He,A_matrix,gen_pre):
    # 1. 输入验证
    if not breeding_pairs:
        return 0.0
    if He <= 0:
        raise ValueError("Initial heterozygosity (He) must be positive")
    if geno_matrix.shape[0] != len(ebv_vector):
        raise ValueError("geno_matrix rows must match ebv_vector length")
    
    total_gain = 0.0
    total_deltaC = 0.0
    
    # 2. 遍历所有配对
    for female, male in breeding_pairs:
        try:
            female_idx = id_to_index[female]
            male_idx = id_to_index[male]
        except KeyError as e:
            raise ValueError(f"ID {e} not found in id_to_index") from e
        
        # 3. 计算遗传增益分量
        pair_gain = (ebv_vector[female_idx] + ebv_vector[male_idx]) / 2
        total_gain += pair_gain
        
        # 4. 计算杂合率分量
        Cfm = A_matrix[female_idx,male_idx] / 2
        
        total_C += Cfm
    
    # 5. 计算平均增益和杂合率
    avg_gain = total_gain / len(breeding_pairs)
    avg_C = total_C / len(breeding_pairs)
    
    # 6. 计算多样性衰减
    deltaC = 1 - (1 - avg_C)**(1 / gen_pre)
    # decay_factor = (1-deltaC)**gen

    eps = 1e-12
    gain_pos = max(avg_gain, eps)
    ratio = max((1-deltaC),eps)
    
    return np.log(gain_pos) + gen * np.log(ratio)
    

def setup_deap_gain(female_ids, male_ids,id_to_index,ebv_vector,geno_matrix, gen, He,A_matrix,gen_pre,female_num=50, male_num=25, female_mates=1, male_mates=2):
    """
    Sets up DEAP for the genetic algorithm.
    """
    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax)

    toolbox = base.Toolbox()

    def create_individual():
        selected_females = np.random.choice(female_ids, size=female_num, replace=False)
        selected_males = np.random.choice(male_ids, size=male_num, replace=False)

        breeding_pairs = []
        female_counts = {female: 0 for female in selected_females}
        male_counts = {male: 0 for male in selected_males}

        for female in selected_females:
            for male in selected_males:
                if female_counts[female] < female_mates and male_counts[male] < male_mates:
                    breeding_pairs.append((female, male))
                    female_counts[female] += 1
                    male_counts[male] += 1

        # Create an individual with fitness attribute
        individual = creator.Individual(breeding_pairs)
        individual.fitness = creator.FitnessMax()
        return individual

    toolbox.register("individual", tools.initIterate, creator.Individual, create_individual)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    # Fitness evaluation
    def evaluate(individual):
        """
        Evaluate an individual based on genetic gain and genetic diversity.
        """
        gain = genetic_gain(individual,ebv_vector, id_to_index,geno_matrix, gen, He,A_matrix,gen_pre)
        return gain,

    toolbox.register("evaluate", evaluate)

    def crossover(ind1, ind2):
        # Combine breeding pairs from both parents
        child1 = ind1[:len(ind1) // 2] + ind2[len(ind2) // 2:]
        child2 = ind2[:len(ind2) // 2] + ind1[len(ind1) // 2:]

        # Repair the children to ensure constraints are satisfied
        child1 = repair_individual(child1, female_ids, male_ids, female_num, male_num, female_mates, male_mates)
        child2 = repair_individual(child2, female_ids, male_ids, female_num, male_num, female_mates, male_mates)

        ind1[:], ind2[:] = child1[:], child2[:]

        return ind1, ind2

    toolbox.register("mate", crossover)

    # Mutation
    def mutate(individual):
        idx = random.randint(0, len(individual) - 1)
        new_female = np.random.choice(female_ids)
        new_male = np.random.choice(male_ids)
        individual[idx] = (new_female, new_male)
        # Repair the individual to ensure constraints are satisfied
        individual_repair = repair_individual(individual, female_ids, male_ids, female_num, male_num, female_mates, male_mates)
        individual[:] = individual_repair[:]
        return individual,

    toolbox.register("mutate", mutate)

    # Selection
    toolbox.register("select", tools.selTournament, tournsize=3)

    return toolbox

def custom_ga(population, toolbox, cxpb, mutpb, ngen, convergence_window,stats=None,convergence_eps=1e-5):
    logbook = tools.Logbook()
    if stats is not None:
        logbook.header = ["gen"] + stats.fields

    # 初始评估
    invalid_ind = [ind for ind in population if not ind.fitness.valid]
    fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit

    # 记录初始统计信息
    if stats is not None:
        record = stats.compile(population)
        logbook.record(gen=0,**record)

    for gen in range(1, ngen + 1):
        # 选择后代（从当前种群中选择，但不直接修改当前种群）
        offspring = toolbox.select(population, len(population))

        # 克隆后代，避免修改原个体
        offspring = [toolbox.clone(ind) for ind in offspring]

        # 交叉：作用于后代种群
        for child1, child2 in zip(offspring[::2], offspring[1::2]):
            if random.random() < cxpb:
                toolbox.mate(child1, child2)
                del child1.fitness.values  # 重置适应度
                del child2.fitness.values

        # 变异：作用于后代种群
        for mutant in offspring:
            if random.random() < mutpb:
                toolbox.mutate(mutant)
                del mutant.fitness.values  # 重置适应度

        # 评估新生成的后代
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        # 替换当前种群（完全替换为后代）
        population[:] = offspring

        # 记录统计信息
        if stats is not None:
            record = stats.compile(population)
            logbook.record(gen=gen,**record)
            print(logbook.stream)

        # 检查收敛条件（例如适应度变化小于阈值）
        if gen > convergence_window:
            # 提取最近10代的最大适应度值
            max_values = [entry["max"] for entry in logbook[-convergence_window:]]
            # 计算相邻代的适应度差异
            diffs = np.abs(np.diff(max_values))
            mean_diff = np.mean(diffs)
            max_diff = np.max(diffs)
            if mean_diff < convergence_eps and max_diff < 2 * convergence_eps:
                print(f"在第 {gen} 代收敛（适应度变化 < {convergence_eps}）")
                break

    return population, logbook

def save_breeding_pairs(breeding_pairs, output_path):
    """
    Saves the breeding pairs to a JSON file.

    Parameters:
        breeding_pairs (list): List of tuples representing breeding pairs (female_id, male_id).
        output_path (str): Path to the output JSON file.
    """
    breeding_pairs_list = [{"female_id": pair[0], "male_id": pair[1]} for pair in breeding_pairs]
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump({"breeding_pairs": breeding_pairs_list}, f, ensure_ascii=False, indent=4)
    print(f"Breeding pairs have been saved to {output_path}")

def train_actor_critic(
    individual_ids,
    ebv_vector,
    female_ids,
    male_ids,
    geno_matrix,
    gen,
    He,
    A_matrix,
    gen_pre,
    output_breeding_pairs_path="breeding_pairs.json"
):
    """
    Trains the actor-critic networks to select optimal breeding pairs.

    Parameters:
        individual_ids (list): List of individual IDs.
        ebv_vector (np.ndarray): EBV vector.
        kinship_matrix (np.ndarray): Kinship matrix.
        female_ids (list): List of female individual IDs.
        male_ids (list): List of male individual IDs.
        epochs (int): Number of training epochs.
        early_stop_patience (int): Patience for early stopping.
        weights (tuple): Weights for genetic gain, diversity, and penalty.
        kinship_threshold (float): Threshold for applying penalties on kinship.
        output_breeding_pairs_path (str): Path to save the final breeding pairs JSON file.

    Returns:
        tuple: Trained actor and critic networks.
    """

    # Mapping from ID to index
    id_to_index = {id_: idx for idx, id_ in enumerate(individual_ids)}

    stats = tools.Statistics(lambda ind: ind.fitness.values[0])
    stats.register("avg", np.mean)
    stats.register("min", np.min)
    stats.register("max", np.max)

    toolbox_gain = setup_deap_gain(female_ids, male_ids, id_to_index,ebv_vector,geno_matrix, gen, He,A_matrix,gen_pre)
    
    population = toolbox_gain.population(n=100)
    
    population, log = custom_ga(
    population, toolbox_gain,
    cxpb=0.5, mutpb=0.2, ngen=1000,
    convergence_window=10,stats=stats, convergence_eps=1e-5
    )
    # results,log = algorithms.eaSimple(population, toolbox_gain, cxpb=0.5, mutpb=0.2, ngen=50, stats = stats,verbose=True)

    final_breeding_pairs = max(population, key=lambda ind: ind.fitness.values[0])
    # Save the final breeding pairs to a JSON file
    save_breeding_pairs(final_breeding_pairs, output_breeding_pairs_path)

    # Print the breeding pairs
    print(f"reward: {final_breeding_pairs.fitness.values[0]}")
    print("Final Breeding Pairs:")
    for pair in final_breeding_pairs:
        print(f"Female ID: {pair[0]}, Male ID: {pair[1]}")
    print(f"Number of Pairs: {len(final_breeding_pairs)}")
    print(f"Number of Females: {len({pair[0] for pair in final_breeding_pairs})}")  # 正确获取唯一雌性数量
    print(f"Number of Males: {len({pair[1] for pair in final_breeding_pairs})}")  
    return final_breeding_pairs

def main():
    """
    Main function to execute the breeding pair selection and policy network training.
    """
    parser = argparse.ArgumentParser(description="Breeding Pair Selection using Actor-Critic Reinforcement Learning with GNN")
    parser.add_argument('input_file', type=str, help='Path to the input JSON file.')
    parser.add_argument('--output_pairs', type=str, default='breeding_pairs.json', help='Path to save the breeding pairs JSON file.')
    args = parser.parse_args()

    # Load input data
    individual_ids, ebv_vector, kinship_matrix, female_ids, male_ids, geno_matrix,gen,He,A_matrix,gen_pre = load_input_data(args.input_file)
    if He is None:
        He = calHe(geno_arr=geno_matrix)
    # Train the actor-critic networks and get the final breeding pairs
    final_breeding_pairs = train_actor_critic(
        individual_ids,
        ebv_vector,
        female_ids,
        male_ids,
        geno_matrix,
        gen,
        He,
        A_matrix,
        gen_pre,
        output_breeding_pairs_path=args.output_pairs
    )

    # After training, you can save the models or perform further evaluation
    # Example: 
    # trained_actor.save('actor_network.h5')
    # trained_critic.save('critic_network.h5')

if __name__ == "__main__":
    main()
