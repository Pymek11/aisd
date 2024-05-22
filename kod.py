import math

def load_cities(filename):
    cities = {}
    with open(filename, 'r') as file:
        for line in file:
            parts = line.split()
            city_id = int(parts[0])
            x = float(parts[1])
            y = float(parts[2])
            cities[city_id] = (x, y)
    return cities

def euclidean_distance(point1, point2):
    return math.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)

def calculate_path_length(cities, path):
    length = 0.0
    num_cities = len(path)
    for i in range(num_cities):
        city1 = cities[path[i]]
        city2 = cities[path[(i + 1) % num_cities]]  # wrap around to the first city
        length += euclidean_distance(city1, city2)
    return length

# Example usage
filename = 'TSP.txt'
cities = load_cities(filename)
path = list(cities.keys())  # Example path, can be any permutation
length = calculate_path_length(cities, path)
print(f"Długość ścieżki: {length}")

#2
import time

# Ścieżka zgodna z kolejnością w pliku
start_time = time.time()
default_path = list(cities.keys())
default_length = calculate_path_length(cities, default_path)
default_time = time.time() - start_time

print(f"Długość ścieżki według kolejności w pliku: {default_length}, Czas: {default_time}s")
def greedy_algorithm(cities):
    unvisited = set(cities.keys())
    current_city = default_path[0]
    path = [current_city]
    unvisited.remove(current_city)
    
    while unvisited:
        next_city = min(unvisited, key=lambda city: euclidean_distance(cities[current_city], cities[city]))
        path.append(next_city)
        unvisited.remove(next_city)
        current_city = next_city
    
    return path

# Algorytm zachłanny
start_time = time.time()
greedy_path = greedy_algorithm(cities)
greedy_length = calculate_path_length(cities, greedy_path)
greedy_time = time.time() - start_time

print(f"Długość ścieżki algorytmu zachłannego: {greedy_length}, Czas: {greedy_time}s")

import random

def alg_heurystyczny(cities, initial_temp, cooling_rate):
    def swap(path):
        new_path = path[:]
        i, j = random.sample(range(len(path)), 2)
        new_path[i], new_path[j] = new_path[j], new_path[i]
        return new_path

    current_path = list(cities.keys())
    current_length = calculate_path_length(cities, current_path)
    temp = initial_temp
    
    while temp > 1:
        new_path = swap(current_path)
        new_length = calculate_path_length(cities, new_path)
        
        if new_length < current_length or random.uniform(0, 1) < math.exp((current_length - new_length) / temp):
            current_path = new_path
            current_length = new_length
        
        temp *= cooling_rate
    
    return current_path, current_length

# Algorytm heurystyczny
initial_temp = 10000
cooling_rate = 0.995
start_time = time.time()
heuristic_path, heuristic_length = alg_heurystyczny(cities, initial_temp, cooling_rate)
heuristic_time = time.time() - start_time

print(f"Długość ścieżki algorytmu heurystycznego: {heuristic_length}, Czas: {heuristic_time}s")

print(f"Kolejność z pliku: {default_length}, Czas: {default_time}s")
print(f"Algorytm zachłanny: {greedy_length}, Czas: {greedy_time}s")
print(f"Algorytm heurystyczny: {heuristic_length}, Czas: {heuristic_time}s")

improvement_greedy = ((default_length - greedy_length) / default_length) * 100
improvement_heuristic = ((default_length - heuristic_length) / default_length) * 100

print(f"Poprawa zachłanny (%): {improvement_greedy:.2f}%")
print(f"Poprawa heurystyczny (%): {improvement_heuristic:.2f}%")
