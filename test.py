import pickle
from gurobipy import Model, GRB, quicksum, gurobi
import os
import numpy as np

# Specify the full path to your Gurobi license file
#gurobi_license_path = "/Users/mariannapiperigou/Documents/gurobi.lic"  # marianna
#gurobi_license_path = "C:/Users/Jacob/OneDrive - Delft University of Technology/Documents/gurobi.lic"  # job
os.environ['GRB_LICENSE_FILE'] = '/Users/julia/Documents/OperationsOptimisation/gurobi.lic'
# Set the environment variable
#os.environ["GRB_LICENSE_FILE"] = gurobi_license_path

#Data Extraction
with open("B.pickle", "rb") as file1:
    bins = pickle.load(file1)

with open("I.pickle", "rb") as file2:
    items = pickle.load(file2)
print(items)
print(bins)

'''
Parameter Definition
'''

mbins = int(sum(entry[1][2] for entry in bins.values())/2)  # number of bins -- should be halved i think
<<<<<<< Updated upstream
nitems = 10                                      # number of items --> why???
=======
nitems = 8                               # number of items --> why???
>>>>>>> Stashed changes
li = [values[0] for values in items.values()]        # length of item
hi = [values[1] for values in items.values()]        # height of item
ai = [li[i] * hi[i] for i in range(len(li))]         # area of iteam
Lj = [values[1][0] for values in bins.values()]        # Length of binf
L = max(Lj)
Hj = [values[1][1] for values in bins.values()]        # height of bin
H = max(Hj)
Aj = [Lj[i] * Hj[i] for i in range(len(Lj))]         # area of bin
Cj = [values[1][3] for values in bins.values()]         # cost of bin
a = [values[1][4] for values in bins.values()]                  # corner shape of bin
b = [values[1][5] for values in bins.values()]                 # corner shape of bin

# Orientation parameters
lip = [values[2] for values in items.values()]       # item can be rotated by pi/2 or no
#dont get points of this??:
lplus = [1 if lip[i] == 1 else 0 for i in range(nitems)]  # 1 if can rotate along length, 0 otherwise
hplus = [1 if lip[i] == 1 else 0 for i in range(nitems)]  # 1 if can rotate along height, 0 otherwise

'''
How to tackle the cut???
'''
# Indices for bins without a cut (0-3) and bins with a cut (4-7)
indices_no_cut = [k for k, v in bins.items() if v[1][-2:] == [-1, -1]]#range(2)  # 0, 1, 2, 3
indices_with_cut = [k for k, v in bins.items() if v[1][-2:] != [-1, -1]]#range(2, 4)  # 4, 5, 6, 7

# Divide into subsets
bins_no_cut = {
    "length": [Lj[i] for i in indices_no_cut],
    "height": [Hj[i] for i in indices_no_cut],
}

bins_with_cut = {
    "length": [Lj[i] for i in indices_with_cut],
    "height": [Hj[i] for i in indices_with_cut],
    "a": [a[i] for i in indices_with_cut],
    "b": [b[i] for i in indices_with_cut],
}

Lnc = bins_no_cut['length']
Hnc = bins_no_cut['height']
Lcut = bins_with_cut['length']
Hcut = bins_with_cut['height']
acut = bins_with_cut['a']
bcut = bins_with_cut['b']

'''
Model Definition
'''
model = Model("2DBPP")
model.setParam('TimeLimit', 60*5)
model.setParam('Method', 2)
'''
Variables Definition
'''
p_ij = model.addVars(nitems, mbins, vtype=GRB.BINARY, name="p_ij")      # if box i in container j
u_j = model.addVars(mbins, vtype=GRB.BINARY, name='u_j')                # if container j is used

xp = model.addVars(nitems, nitems, vtype=GRB.BINARY, name="x_p")        # if box i is to the right of box k
zp = model.addVars(nitems, nitems, vtype=GRB.BINARY, name="z_p")        # if box i is above box k

# Define variables (coordinates)
x = model.addVars(nitems, vtype=GRB.CONTINUOUS, name="xi")              # Bottom-left x-coordinate
z = model.addVars(nitems, vtype=GRB.CONTINUOUS, name="zi")              # Bottom-left z-coordinate
xprime = model.addVars(nitems, vtype=GRB.CONTINUOUS, name="x_i_prime")  # Top-right x-coordinate
zprime = model.addVars(nitems, vtype=GRB.CONTINUOUS, name="z_i_prime")  # Top-right z-coordinate

r = model.addVars(nitems, 2, 2, vtype=GRB.BINARY, name="r")             # if
rho = model.addVars(nitems, vtype=GRB.BINARY, name='rho')               # if item is rotated?

g = model.addVars(nitems, vtype=GRB.BINARY, name='g')                   # 1 if item i lies on the ground of the bin

beta1 = model.addVars(nitems, nitems, vtype=GRB.BINARY, name='beta1')   # 1 if vertex 1 of item i is supported by item j
beta2 = model.addVars(nitems, nitems, vtype=GRB.BINARY, name='beta2')   # 1 if vertex 2 of item i is supported by item j
gamma = model.addVars(nitems, vtype=GRB.BINARY, name='gamma')           # 1 if vertex 1 of item i is supported by the cut of the bin where it is placed

<<<<<<< Updated upstream

=======
v=model.addVars(nitems,nitems, vtype=GRB.CONTINUOUS, lb=0, name="v_i_k")
h=model.addVars(nitems,nitems, vtype=GRB.BINARY, name="h_i_k")
n1=model.addVars(nitems,nitems, vtype=GRB.BINARY,  name="n1_i_k")
n2=model.addVars(nitems,nitems, vtype=GRB.BINARY,  name="n2_i_k")
model.update()
>>>>>>> Stashed changes
'''
Constraints Definition
'''

#constraint 3: area of items not larger than area of bin
for i in range(nitems):
    for j in range(mbins):
        model.addConstr(
            ai[i] * p_ij[i, j] <= Aj[j] * u_j[j],
            name=f"AreaConstraint_{i}_{j}"
        )

# Constraint 4: each item i is assigned to one bin j
for i in range(nitems):
    model.addConstr(quicksum(p_ij[i, j] for j in range(mbins)) == 1, name=f"OneItemOneBin_{i}")

# Constraint 5: item should not exceed container size
for i in range(nitems):
    model.addConstr(xprime[i] <= quicksum(Lj[j] * p_ij[i, j] for j in range(mbins)),
                    name=f"ItemFitToBin_{i}")

# Constraint 7: item should not exceed container size
for i in range(nitems):
    model.addConstr(zprime[i] <= quicksum(Hj[j] * p_ij[i,j] for j in range(mbins)),
                    name=f"ItemFitToBinZ_{i}")

#Constraint 8 and 10 : orthogonal rotation
# why is rho not used here?
for i in range(nitems):
    model.addConstr(
        xprime[i] - x[i] == sum(r[i, 0, b] * [li[i], hi[i]][b] for b in range(2)),
        name=f"TransformX_{i}"
    )
    model.addConstr(
        zprime[i] - z[i] == sum(r[i, 1, b] * [li[i], hi[i]][b] for b in range(2)),
        name=f"TransformZ_{i}"
    )


#constraint 11 and 12:
for i in range(nitems):
    # Each side aligns with exactly one axis
    model.addConstr(r[i, 0, 0] + r[i, 1, 0] == 1, name=f"LengthAlign_{i}")
    model.addConstr(r[i, 0, 1] + r[i, 1, 1] == 1, name=f"HeightAlign_{i}")

    # Each axis has exactly one side aligned
    model.addConstr(r[i, 0, 0] + r[i, 0, 1] == 1, name=f"XAxisAlign_{i}")
    model.addConstr(r[i, 1, 0] + r[i, 1, 1] == 1, name=f"ZAxisAlign_{i}")

# Constraint 13: Overlap occurs only within the same bin j
for i in range(nitems):
    for k in range(nitems):
        for j in range(mbins):
            if i != k:
                model.addConstr(xp[i, k] + xp[k, i] + zp[i, k] + zp[k, i] >= (p_ij[i, j] + p_ij[k, j]) - 1,
                                name=f"OverlapInBin_{i}_{k}_toBin_{j}")


# Constraint for relative positioning (this ensures item relative positions) (added by chatgpt, so far they work )
for i in range(nitems):
    for k in range(nitems):
        if i != k:  # Avoid self-comparison
            model.addConstr(xp[i, k] + xp[k, i] + zp[i, k] + zp[k, i] >= 1,
                            name=f"Overlap_{i}_{k}")

            model.addConstr(xp[i, k] <= (1 - xp[k, i]), name=f"RightLeft_{i}_{k}")
            model.addConstr(xp[k, i] <= (1 - xp[i, k]), name=f"LeftRight_{i}_{k}")

            model.addConstr(zp[i, k] <= (1 - zp[k, i]), name=f"AboveBelow_{i}_{k}")
            model.addConstr(zp[k, i] <= (1 - zp[i, k]), name=f"BelowAbove_{i}_{k}")

# Constraint 14:
for i in range(nitems):
    for k in range(nitems):
        if i != k:  # Avoid self-comparison
            model.addConstr(xprime[k] <= x[i] + (1 - xp[i, k]) * L,
                            name=f"Constraint14_{i}_{k}")

# Constraint 15:
for i in range(nitems):
    for k in range(nitems):
        if i != k:  # Avoid self-comparison
            model.addConstr(x[i] + 1 <= xprime[k] + xp[i, k] * L,
                            name=f"Constraint15_{i}_{k}")

# Constraint 18:
for i in range(nitems):
    for k in range(nitems):
        if i != k:  # Avoid self-comparison
            model.addConstr(zprime[k] <= z[i] + (1 - zp[i, k]) * H,
                            name=f"Constraint16_{i}_{k}")


# Orientation constraints (19–21)
for i in range(nitems):
<<<<<<< Updated upstream
    model.addConstr(r[i, 1, 0] <= lplus[i], name=f"OrientationLength_{i}") #if length in vertical position
    model.addConstr(r[i, 1, 1] <= hplus[i], name=f"OrientationHeight_{i}")  #if heigth
=======
    model.addConstr(r[i, 0, 1] <= lplus[i], name=f"OrientationLength_{i}")
    model.addConstr(r[i, 1, 0] <= hplus[i], name=f"OrientationHeight_{i}")
>>>>>>> Stashed changes

# #Constraint 22
for i in range(nitems):
    for j in indices_with_cut:
        model.addConstr(
            z[i] + b[j]/a[j] * x[i] >= b[j] - 50000*(1 - p_ij[i, j]), name=f'Constraint_{i}_{j}'
        )
        model.addConstr(
            z[i] + b[j] / a[j] * x[i] >= b[j] - 50000 * (1 - p_ij[i, j]) + 50000* (1 - gamma[i]), name=f'Constraint_{i}_{j}'
        )

'''
Constraints from lecture, mostly for vertical stability and cut
'''

## Stability extra

#constraint 16 from lecture (flagging constraint)
for i in range(nitems):
    for j in range(mbins):
        model.addConstr(
            u_j[j] >= p_ij[i, j]
        )

#Constraint for stability
for i in range(nitems):
    model.addConstr(
        gamma[i] + quicksum(beta1[i, j] for j in range(mbins)) + quicksum(beta2[i, j] for j in range(mbins)) + 2*g[i]
        >= 2
    )



#Constraints for cut of box

# #constraint for when theres no cut
# for i in range(nitems):
#     for j in range(len(Lnc)):
#         model.addConstr(
#             gamma[i] <= 1 - p_ij[i, j]
#         )
#
# for i in range(nitems):
#     for j in range(len(Lcut)):
#         model.addConstr(
#             u_j[j] >= - bcut[j]/acut[j] * x[i] + bcut[j] - (1 - p_ij[i,j])
#         )
#
# for i in range(nitems):
#     for j in range(len(Lcut)):
#         model.addConstr(
#             u_j[j] >= - bcut[j]/acut[j] * x[i] + bcut[j] - (1 - p_ij[i,j]) + (1 - gamma[i])
#         )

'''
Objective Function
'''
objective = quicksum(Cj[j] * u_j[j] for j in range(len(Cj))) #sum of Cj[i] * u_j for each i
model.setObjective(objective, GRB.MINIMIZE)


'''
Print constraints
'''
'''def my_callback(model, where):
    if where == GRB.Callback.MIPSOL:  # At integer feasible solutions
        print("\n--- Constraint Values at Current Iteration ---")
        for constr in model.getConstrs():
            lhs = sum(model.getVarByName(var.varName).X * coeff
                      for var, coeff in zip(constr.getVars(), constr.getCoeff()))
            rhs = constr.RHS
            print(f"{constr.ConstrName}: LHS = {lhs}, RHS = {rhs}, Residual = {rhs - lhs}")'''


model.optimize()

if model.status == GRB.INFEASIBLE:
    print("The model is infeasible. Computing IIS...")
    model.computeIIS()
    model.write("infeasible.ilp")

import matplotlib.pyplot as plt
import numpy as np

def visualize_with_overlap(items, nitems, mbins, Lj, Hj, xi, zi, x_i_prime, z_i_prime, p_ij, bins_with_cut, a, b):
    fig, axs = plt.subplots(1, mbins, figsize=(15, 5))

    # Ensure the model is optimized before visualization
    if model.status in [GRB.OPTIMAL, GRB.SUBOPTIMAL]:
        for j in range(mbins):
            axs[j].set_xlim(0, Lj[j])
            axs[j].set_ylim(0, Hj[j])
            axs[j].set_title(f"Bin {j}")
            axs[j].set_aspect('equal')
            # Set x and y ticks every 25 units
            axs[j].set_xticks(np.linspace(0, Lj[j], 5))
            axs[j].set_yticks(np.linspace(0, Hj[j], 5))

            bin_items = []
            for i in range(nitems):
                if p_ij[i, j].X > 0.5:  # Only visualize items assigned to bin j
                    x_start = xi[i].X
                    z_start = zi[i].X
                    width = x_i_prime[i].X - xi[i].X
                    height = z_i_prime[i].X - zi[i].X
                    bin_items.append((x_start, z_start, width, height, i))

            # Draw items and check overlaps
            for i, (x, z, w, h, item) in enumerate(bin_items):
                rect_color = "green"  # Default color for non-overlapping items
                for x2, z2, w2, h2, other_item in bin_items:
                    if item != other_item:  # Don't compare an item with itself
                        if not (x + w <= x2 or x2 + w2 <= x or z + h <= z2 or z2 + h2 <= z):
                            rect_color = "red"  # Overlapping items are marked in red
                            break
                axs[j].add_patch(plt.Rectangle((x, z), w, h, color=rect_color, alpha=0.5))

                axs[j].text(x + w / 2, z + h / 2, f"{item}\n {items[i][-3:]}", ha='center', va='center')

            # **Draw the ULD outline**
            if j in bins_with_cut:  # Checking if bin has a cut
                cut_a = a[j]
                cut_b = b[j]
                if cut_a != -1 and cut_b != -1:  # Only plot if the bin has a defined cut
                    # Calculate the cut line
                    x_cut_vals = np.array([0, Lj[j]])  # X starts at 0, ends at bin width
                    z_cut_vals = - (cut_b / cut_a) * x_cut_vals + cut_b  # Compute cut line equation

                    # Find intersection of cut with bin edge
                    x_cut_intersect = cut_b / (cut_b / cut_a)  # Solves - (b/a) * x + b = 0
                    z_cut_intersect = 0  # At the bottom

                    # **Draw ULD outline**
                    outline_x = [1, 1, Lj[j], Lj[j], x_cut_intersect, 0]  # x-coordinates
                    outline_z = [cut_b, Hj[j], Hj[j], 0, z_cut_intersect, cut_b]  # z-coordinates

                    axs[j].plot(outline_x, outline_z, 'k-', linewidth=2, label="ULD Outline")  # Draw outline

            else:  # **Regular bins (without a cut)**
                outline_x = [0, 0, Lj[j], Lj[j], 0]  # Full rectangle
                outline_z = [0, Hj[j], Hj[j], 0, 0]  # Full rectangle

                axs[j].plot(outline_x, outline_z, 'k-', linewidth=2, label="Bin Outline")

        plt.tight_layout()
        plt.show()
    else:
        print("Model didn't find a solution within the time limit.")

# Call the function
if model.status == GRB.OPTIMAL or model.status == GRB.SUBOPTIMAL:
    visualize_with_overlap(items, nitems, mbins, Lj, Hj, x, z, xprime, zprime, p_ij, indices_with_cut, a, b)
else:
    print("Model didn't find a solution within the time limit.")

if model.status in [GRB.OPTIMAL, GRB.SUBOPTIMAL]:
    print("\n--- Final Constraint Values ---")
    for constr in model.getConstrs():
        expr = model.getRow(constr)  # Get the constraint's left-hand side expression
        lhs_value = sum(expr.getVar(i).X * expr.getCoeff(i) for i in range(expr.size()))
        rhs_value = constr.RHS  # Right-hand side of the constraint
        residual = rhs_value - lhs_value  # Difference between RHS and LHS

        print(f"{constr.ConstrName}: LHS = {lhs_value}, RHS = {rhs_value}, Residual = {residual}")
<<<<<<< Updated upstream
=======

    for i in range(nitems):
        for j in range(mbins):
                if p_ij[i, j].X == 1:
                    print(f"p_ij[{i},{j}]: {p_ij[i, j].X}")

    for j in range(mbins):
        if u_j[j].X == 1:
            print(f"u_j[{j}]: {u_j[j].X}")

    for i in range(nitems):
        for j in range(nitems):
            if i != j:
                if xp[i, j].X == 1:
                    print(f"x_p[{i},{j}]: {xp[i, j].X}")

    for i in range(nitems):
        for j in range(nitems):
            if i!=j:
                if zp[i, j].X == 1:
                    print(f"z_p[{i},{j}]: {zp[i, j].X}")

    for i in range(nitems):
        for j in range(n_orient):
            for k in range(n_axes):
                if r[i, j, k].X == 1:
                    print(f"r[{i},{j},{k}]: {r[i, j, k].X}")

    for i in range(nitems):
        if rho[i].X == 1:
            print(f"rho[{i}]: {rho[i].X}")

    for i in range(nitems):
        if g[i].X == 1:
            print(f"g[{i}]: {g[i].X}")

    for i in range(nitems):
        for j in range(nitems):
            if i != j:
                if beta1[i, j].X == 1:
                    print(f"beta1[{i},{j}]: {beta1[i, j].X}")

    for i in range(nitems):
        for j in range(nitems):
            if i != j:
                if beta2[i, j].X == 1:
                    print(f"beta2[{i},{j}]: {beta2[i, j].X}")

    for i in range(nitems):
        if gamma[i].X == 1:
            print(f"gamma[{i}]: {gamma[i].X}")

    for i in range(nitems):
        print(f"x[{i}] = {x[i].X}")
        print(f"z[{i}] = {z[i].X}")
        print(f"xprime[{i}] = {xprime[i].X}")
        print(f"zprime[{i}] = {zprime[i].X}")
>>>>>>> Stashed changes
