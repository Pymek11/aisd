mport math

def load_cities(filename):
    cities = []
    with open(filename, 'r') as file:
        for line in file:
            data = line.strip().split()
            city_id = int(data[0])
            x = float(data[1])
            y = float(data[2])
            cities.append((city_id, x, y))
    return cities

def euclidean_distance(city1, city2):
    return math.sqrt((city1[1] - city2[1])**2 + (city1[2] - city2[2])**2)

def path_length(path, cities):
    length = 0.0
    for i in range(len(path) - 1):
        length += euclidean_distance(cities[path[i] - 1], cities[path[i + 1] - 1])
    length += euclidean_distance(cities[path[-1] - 1], cities[path[0] - 1])  # return to the starting city
    return length

# Example usage
filename = 'TSP.txt'
cities = load_cities(filename)
path = [i+1 for i in range(len(cities))]  # Example path (1, 2, 3, ..., 100)
print(f"Długość ścieżki: {path_length(path, cities)}")

path_in_file_order = [i+1 for i in range(len(cities))]
length_in_file_order = path_length(path_in_file_order, cities)
print(f"Długość ścieżki w kolejności z pliku: {length_in_file_order}")

def greedy_tsp(cities):
    n = len(cities)
    unvisited = set(range(1, n + 1))
    current_city = 1
    path = [current_city]
    unvisited.remove(current_city)
    
    while unvisited:
        next_city = min(unvisited, key=lambda city: euclidean_distance(cities[current_city - 1], cities[city - 1]))
        path.append(next_city)
        unvisited.remove(next_city)
        current_city = next_city
    
    return path

greedy_path = greedy_tsp(cities)
length_greedy = path_length(greedy_path, cities)
print(f"Długość ścieżki algorytmu zachłannego: {length_greedy}")

import random

def simulated_annealing(cities, initial_temp, cooling_rate, max_iterations):
    def swap_random(path):
        new_path = path[:]
        i, j = random.sample(range(len(path)), 2)
        new_path[i], new_path[j] = new_path[j], new_path[i]
        return new_path
    
    current_path = [i + 1 for i in range(len(cities))]
    random.shuffle(current_path)
    current_length = path_length(current_path, cities)
    
    best_path = current_path[:]
    best_length = current_length
    
    temp = initial_temp
    
    for _ in range(max_iterations):
        new_path = swap_random(current_path)
        new_length = path_length(new_path, cities)
        
        if new_length < current_length or random.random() < math.exp((current_length - new_length) / temp):
            current_path = new_path
            current_length = new_length
        
        if current_length < best_length:
            best_path = current_path
            best_length = current_length
        
        temp *= cooling_rate
    
    return best_path, best_length

initial_temp = 10000
cooling_rate = 0.995
max_iterations = 100000

sa_path, length_sa = simulated_annealing(cities, initial_temp, cooling_rate, max_iterations)
print(f"Długość ścieżki algorytmu heurystycznego (Simulated Annealing): {length_sa}")

# Results
print("Podsumowanie wyników:")
print(f"Kolejność z pliku: Długość = {length_in_file_order}")
print(f"Algorytm zachłanny: Długość = {length_greedy}")
print(f"Algorytm heurystyczny (Simulated Annealing): Długość = {length_sa}")

# Improvement comparison
improvement_greedy = (length_in_file_order - length_greedy) / length_in_file_order * 100
improvement_sa = (length_in_file_order - length_sa) / length_in_file_order * 100

print(f"Poprawa zachłannego względem kolejności z pliku: {improvement_greedy:.2f}%")
print(f"Poprawa heurystycznego względem kolejności z pliku: {improvement_sa:.2f}%")
