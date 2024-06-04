def read_data_from_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    items_data = []
    for line in lines[1:]:  
        item = tuple(map(int, line.strip().split(',')[1:]))  # Skipping the first element (id)
        items_data.append(item)

    return items_data

def greedy_knapsack(items, capacity):
    items.sort(key=lambda x: x[2]/(x[0]*x[1]), reverse=True)  # Sort items by value per unit area
    knapsack = []
    total_value = 0
    total_width = 0
    total_height = 0

    for item in items:
        if total_width + item[0] <= capacity and total_height + item[1] <= capacity:
            knapsack.append(item)
            total_value += item[2]
            total_width += item[0]
            total_height += item[1]

    return knapsack, total_value

def dynamic_knapsack(items, capacity):
    n = len(items)
    dp = [[0] * (capacity + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        for j in range(1, capacity + 1):
            if items[i - 1][0] <= j and items[i - 1][1] <= j:
                dp[i][j] = max(dp[i - 1][j], dp[i - 1][j - items[i - 1][0]] + items[i - 1][2], dp[i - 1][j - items[i - 1][1]] + items[i - 1][2])
            elif items[i - 1][0] <= j:
                dp[i][j] = max(dp[i - 1][j], dp[i - 1][j - items[i - 1][0]] + items[i - 1][2])
            elif items[i - 1][1] <= j:
                dp[i][j] = max(dp[i - 1][j], dp[i - 1][j - items[i - 1][1]] + items[i - 1][2])
            else:
                dp[i][j] = dp[i - 1][j]

    knapsack = []
    i, j = n, capacity
    while i > 0 and j > 0:
        if dp[i][j] != dp[i - 1][j]:
            knapsack.append(items[i - 1])
            j -= items[i - 1][0]
        i -= 1

    return knapsack, dp[n][capacity]

filename = "packages500.txt"

items_data = read_data_from_file(filename)

# Example capacity
capacity = 100

greedy_result = greedy_knapsack(items_data, capacity)
print("Greedy Algorithm:")
print("Knapsack:", greedy_result[0])
print("Total Value:", greedy_result[1])

dynamic_result = dynamic_knapsack(items_data, capacity)
print("\nDynamic Programming Algorithm:")
print("Knapsack:", dynamic_result[0])
print("Total Value:", dynamic_result[1])
