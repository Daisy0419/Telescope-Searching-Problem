# One fixed start node (with net outflow = 1).
# One dynamically chosen end node (enforced by ∑i is_end[i]=1).
# At most one incoming/outgoing edge per node (so no node is used more than once in the path).
# Net flow balance for each node (except start) to ensure that if a node is not the end, incoming = outgoing; if it is the end, incoming - outgoing = 1.
# A budget constraint and an objective to maximize the collected prize.
# MTZ constraints to eliminate subtours (so the visited nodes form a single connected path).

set NODES;

param Prize{NODES};
param Cost{NODES, NODES};
param Budget;

# Start node (fixed)
param start;

# is_end[i] = 1 if node i is the chosen end node, 0 otherwise
var is_end{NODES} binary;

# x[i,j] = 1 if traveling from node i to node j, 0 otherwise
var x{NODES, NODES} binary;

# For MTZ subtour elimination
var u{NODES} >= 2 <= card(NODES);


maximize Total_Prize:
    sum{i in NODES, j in NODES: i != j} Prize[j] * x[i, j];


subject to One_End_Node:
    sum{i in NODES} is_end[i] = 1;

# At most one incoming edge for each node
subject to At_Most_One_Incoming{i in NODES}:
    sum{j in NODES: j != i} x[j,i] <= 1;

# At most one outgoing edge for each node
subject to At_Most_One_Outgoing{i in NODES}:
    sum{j in NODES: j != i} x[i,j] <= 1;


# Flow balance (excluding start):
subject to Flow_Balance{i in NODES: i != start}:
    sum{j in NODES: j != i} x[i,j] - sum{j in NODES: j != i} x[j,i] 
        = - is_end[i];

# For the start node: outflow - inflow = 1
subject to Start_Flow:
    sum{j in NODES: j != start} x[start,j] - sum{j in NODES: j != start} x[j,start] = 1;


# Budget constraint
subject to Budget_Constraint:
    sum{i in NODES, j in NODES: i != j} Cost[i, j] * x[i, j] <= Budget;

# Subtour elimination constraints (Miller-Tucker-Zemlin)
subject to Subtour_Elimination{i in NODES, j in NODES: i != j}:
    u[i] - u[j] + card(NODES) * x[i, j] <= card(NODES) - 1;

# display Total_Prize;
# display sum{i in NODES, j in NODES: i != j} Cost[i,j] * x[i,j];

# # 1. Compute next[i] = the node j for which x[i,j] = 1
# var next{NODES} integer;

# subject to Next_Definition{i in NODES}:
#     next[i] = sum{j in NODES: j != i} j * x[i,j];

# # 2. Build a route array to list the sequence of visited nodes.
# param route_len integer := card(NODES);  # an upper bound on path length
# param route{1..route_len} default 0;     # store the sequence of nodes

# # We can do an iterative assignment in the AMPL "command" section:
# let route[1] := start;
# for {k in 1..(route_len-1)} {
#   let route[k+1] := next[route[k]];
#   if is_end[route[k+1]] = 1 then {  # If we've reached the end node
#     # Could exit the loop if you want to stop here
#     break;
#   }
# }

# display route;



# # Fixed start node constraints
# subject to Fixed_Start:
#     sum{j in NODES: j != start} x[start, j] = 1;

# subject to No_Incoming_To_Start:
#     sum{j in NODES: j != start} x[j, start] = 0;


# # Enforce exactly one end node
# subject to One_End_Node:
#     sum{i in NODES} is_end[i] = 1;

# # End node constraints
# subject to End_Node_Incoming{i in NODES}:
#     sum{j in NODES: j != i} x[j, i] is_end[i] = 1;

# subject to No_Outgoing_From_End{i in NODES}:
#     sum{j in NODES: j != i} x[i, j] <= 1 - is_end[i];


# # Flow balance (excluding start):
# subject to Flow_Balance{i in NODES: i != start}:
#     sum{j in NODES: j != i} x[i,j] - sum{j in NODES: j != i} x[j,i] 
#         = - is_end[i];

# # For the start node: outflow - inflow = 1
# subject to Start_Flow:
#     sum{j in NODES: j != start} x[start,j] - sum{j in NODES: j != start} x[j,start] = 1;




# # Objective: Maximize the total prize collected
# maximize Total_Prize:
#     sum{i in NODES} Prize[i] * sum{j in NODES: j != i} x[i, j];

# # Budget constraint: Total cost should not exceed the budget
# subject to Budget_Constraint:
#     sum{i in NODES, j in NODES: i != j} Cost[i, j] * x[i, j] <= Budget;

# # Flow constraint: All nodes except start and end must have one incoming and outgoing connection

# subject to Edge_Incoming{i in NODES}:
#     sum{j in NODES: j != i} x[j, i] = (if i = start then 0 else 1);

# subject to Edge_Outgoing{i in NODES}:
#     sum{j in NODES: j != i} x[i, j] = (if i = end then 0 else 1);

# # Subtour elimination constraints (Miller-Tucker-Zemlin formulation)
# subject to Subtour_Elimination{i in NODES, j in NODES: i != j}:
#     u[i] - u[j] + card(NODES) * x[i, j] <= card(NODES) - 1;


# display Total_Prize;
# display sum{i in NODES, j in NODES: i != j} Cost[i, j] * x[i, j];
# param start := min{n in NODES} n;  # Start at the first node (or define a custom start)
# var next{NODES} integer;           # Store the next node for each node in the tour

# # Calculate the route
# subject to Next_Node_Calculation{i in NODES}:
#     next[i] = sum{j in NODES: j != i} j * x[i, j];

# # Output the route in sequence
# param current integer := start;
# param route{1..card(NODES)} integer;  # Record the route

# let {i in 1..card(NODES)} route[i] := current;
# let {i in 1..card(NODES)} current := next[route[i]];

# display route;




# set NODES;  # Set of nodes (cities)

# param Prize{NODES};      # Prize associated with each node
# param Cost{NODES, NODES}; # Cost to travel from node i to node j
# param Budget;             # Maximum budget constraint

# var x{NODES, NODES} binary;  # x[i, j] = 1 if traveling from i to j, 0 otherwise
# var u{NODES} >= 2 <= card(NODES);  # For subtour elimination (Miller-Tucker-Zemlin)

# # Objective: Maximize the total prize collected
# maximize Total_Prize:
#     sum{i in NODES} Prize[i] * sum{j in NODES: j != i} x[i, j];

# # Budget constraint: Total cost should not exceed the budget
# subject to Budget_Constraint:
#     sum{i in NODES, j in NODES: i != j} Cost[i, j] * x[i, j] <= Budget;

# # Flow constraint: Each node has exactly one incoming and one outgoing connection
# subject to Flow_Incoming{i in NODES}:
#     sum{j in NODES: j != i} x[j, i] = 1;

# subject to Flow_Outgoing{i in NODES}:
#     sum{j in NODES: j != i} x[i, j] = 1;

# # Subtour elimination constraints (Miller-Tucker-Zemlin formulation)
# subject to Subtour_Elimination{i in NODES, j in NODES: i != j}:
#     u[i] - u[j] + card(NODES) * x[i, j] <= card(NODES) - 1;


# printf "Optimal Route:\n";
# for {(i, j) in EDGES} {
#     if x[i, j] > 0.5 then
#         printf "Node %d -> Node %d\n", i, j;
# }

# param Total_Cost := sum {(i, j) in EDGES} Cost[i, j] * x[i, j];
# display Total_Cost;
# display Total_Profit;


# data;
# set NODES := 1 2 3;
# set EDGES := (1, 2) (1, 3) (2, 1) (2, 3) (3, 1) (3, 2);

# param Prize :=
# 1 0.800000
# 2 0.500000
# 3 0.700000
# ;
# param Cost :=
# [1, 2] 10.500000
# [1, 3] 5.200000
# [2, 1] 10.500000
# [2, 3] 8.300000
# [3, 1] 5.200000
# [3, 2] 8.300000
# ;
# param Budget := 20;

# set NODES;
# set EDGES within {NODES, NODES};

# param Prize{NODES};
# param Cost{EDGES};
# param Budget;

# var x{EDGES} binary;
# var u{NODES} >= 2 <= card(NODES);  # u[i] defines the order of visiting nodes

# maximize Total_Profit:
#     sum {(i, j) in EDGES} Prize[j] * x[i, j];

# s.t. Start:
#     sum {(1, j) in EDGES} x[1, j] = 1;

# s.t. End:
#     sum {(i, card(NODES)) in EDGES} x[i, card(NODES)] = 1;

# s.t. Flow_Conservation {k in NODES diff {1, card(NODES)}}:
#     sum {(i, k) in EDGES} x[i, k] = sum {(k, j) in EDGES} x[k, j];

# s.t. Travel_Budget:
#     sum {(i, j) in EDGES} Cost[i, j] * x[i, j] <= Budget;

# s.t. Subtour_Elimination {(i, j) in EDGES: i != 1 && j != card(NODES)}:
#     u[i] - u[j] + card(NODES) * x[i, j] <= card(NODES) - 1;

# printf "Optimal Route:\n";
# for {(i, j) in EDGES: x[i, j] > 0.5} {
#     printf "Node %d -> Node %d\n", i, j;
# }

# param Total_Cost := sum {(i, j) in EDGES} Cost[i, j] * x[i, j];
# display Total_Cost;
# display Total_Profit;



# # Sets and Indices
# set NODES;        # Set of nodes (1 to |N|)
# set EDGES within {NODES, NODES};  # Directed EDGES between nodes

# # Parameters
# param Prize{NODES};       # Profit at each node
# param Cost{EDGES};        # Travel cost between nodes
# param Budget;             # Maximum allowed travel cost

# # Variables
# var x{EDGES} binary;      # x[i, j] = 1 if arc from i to j is used
# var u{NODES} >= 2 <= card(NODES);  # u[i] defines order of visit for subtour elimination

# # Objective: Maximize total collected profit
# maximize Total_Profit:
#     sum {(i, j) in EDGES} Prize[j] * x[i, j];

# # Constraints

# # Start and end at specified nodes
# s.t. Start:
#     sum {(1, j) in EDGES} x[1, j] = 1;

# s.t. End:
#     sum {(i, card(NODES)) in EDGES} x[i, card(NODES)] = 1;

# # Ensure each node is visited at most once
# s.t. Flow_Conservation {k in NODES diff {1, card(NODES)}}:
#     sum {(i, k) in EDGES} x[i, k] = sum {(k, j) in EDGES} x[k, j];

# # Travel cost constraint
# s.t. Travel_Budget:
#     sum {(i, j) in EDGES} Cost[i, j] * x[i, j] <= Budget;

# # Subtour elimination
# s.t. Subtour_Elimination {(i, j) in EDGES: i != 1 && j != card(NODES)}:
#     u[i] - u[j] + card(NODES) * x[i, j] <= card(NODES) - 1;

# # Display the route (only selected EDGES)
# printf "Optimal Route:\n";
# for {(i, j) in EDGES: x[i, j] > 0.5} {
#     printf "Node %d -> Node %d\n", i, j;
# }

# # Compute the total travel cost
# param Total_Cost := sum {(i, j) in EDGES} Cost[i, j] * x[i, j];
# display Total_Cost;

# # Display the total profit (objective value)
# display Total_Profit;
