set NODES;

param Prize{NODES};
param Cost{NODES, NODES};
param Budget;

param start;

var is_end{NODES} binary;
var x{NODES, NODES} binary;
var u{NODES} >= 2 <= card(NODES);

maximize Total_Prize:
    sum{i in NODES, j in NODES: i != j} Prize[j] * x[i, j];

subject to One_End_Node:
    sum{i in NODES} is_end[i] = 1;

subject to At_Most_One_Incoming{i in NODES}:
    sum{j in NODES: j != i} x[j, i] <= 1;

subject to At_Most_One_Outgoing{i in NODES}:
    sum{j in NODES: j != i} x[i, j] <= 1;

subject to Flow_Balance{i in NODES: i != start}:
    sum{j in NODES: j != i} x[i,j] - sum{j in NODES: j != i} x[j,i] 
        = - is_end[i];

subject to Start_Flow:
    sum{j in NODES: j != start} x[start,j] 
      - sum{j in NODES: j != start} x[j,start] 
    = 1;

subject to Budget_Constraint:
    sum{i in NODES, j in NODES: i != j} Cost[i, j] * x[i, j] <= Budget;

# Subtour elimination constraints (Miller-Tucker-Zemlin)
subject to Subtour_Elimination{i in NODES, j in NODES: i != j}:
# subject to Subtour_Elimination {i in NODES, j in NODES: i != j and i != start and j != start}:
    u[i] - u[j] + card(NODES) * x[i, j] <= card(NODES) - 1;



