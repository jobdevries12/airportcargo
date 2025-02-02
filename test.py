import pickle
from gurobipy import Model, GRB, quicksum
import os

# Specify the full path to your Gurobi license file
gurobi_license_path = "/Users/mariannapiperigou/Documents/gurobi.lic"

# Set the environment variable
os.environ["GRB_LICENSE_FILE"] = gurobi_license_path

#Data Extraction
with open("G1/B.pickle", "rb") as file1:
    bins = pickle.load(file1)

with open("G1/I.pickle", "rb") as file2:
    items = pickle.load(file2)


'''
Parameter Definition
'''
mbins = sum(entry[1][2] for entry in bins.values()) #number of bins
nitems = 10  #number of items
li = [values[0] for values in items.values()] #length of item
hi = [values[1] for values in items.values()] #height of item
ai = [li[i] * hi[i] for i in range(len(li))] #area of iteam
Lj = [300, 300, 300, 300, 192, 192, 192, 192] #Length of bin
L = max(Lj)
Hj = [155, 155, 155, 155, 155, 155, 155, 155] #height of bin
H = max(Hj)
Aj = [Lj[i] * Hj[i] for i in range(len(Lj))] #area of bin
Cj = [200, 200, 200 ,200, 150, 150, 150, 150] #cost of bin
a = [-1, -1, -1, -1, 42, 42, 42, 42] #corner shape of bin
b = [-1, -1, -1, -1, 53, 53, 53, 53] #corner shape of bin

# Orientation parameters
lip = [values[2] for values in items.values()] #item can be rotated by pi/2 or no
lplus = [1 if lip[i] == 1 else 0 for i in range(nitems)]  # 1 if can rotate along length, 0 otherwise
hplus = [1 if lip[i] == 1 else 0 for i in range(nitems)]  # 1 if can rotate along height, 0 otherwise


'''
How to tackle the cut???
'''
# Indices for bins without a cut (0-3) and bins with a cut (4-7)
indices_no_cut = range(4)  # 0, 1, 2, 3
indices_with_cut = range(4, 8)  # 4, 5, 6, 7

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
model.setParam('TimeLimit', 20)
model.setParam('Method', 2)
'''
Variables Definition
'''
p_ij = model.addVars(nitems, mbins, vtype=GRB.BINARY, name="p_ij") #if box i in container j
u_j = model.addVars(mbins, vtype=GRB.BINARY, name='u_j') #if container j is used

xp = model.addVars(nitems, nitems, vtype=GRB.BINARY, name="x_p") #if box i is to the right of box k
zp = model.addVars(nitems, nitems, vtype=GRB.BINARY, name="z_p") #if box i is above box k

# Define variables (coordinates)
x = model.addVars(nitems, vtype=GRB.CONTINUOUS, name="xi")  # Bottom-left x-coordinate
z = model.addVars(nitems, vtype=GRB.CONTINUOUS, name="zi")  # Bottom-left z-coordinate
xprime = model.addVars(nitems, vtype=GRB.CONTINUOUS, name="x_i_prime")  # Top-right x-coordinate
zprime = model.addVars(nitems, vtype=GRB.CONTINUOUS, name="z_i_prime")  # Top-right z-coordinate

r = model.addVars(nitems, 2, 2, vtype=GRB.BINARY, name="r")
rho = model.addVars(nitems, vtype=GRB.BINARY, name='rho')

g = model.addVars(nitems, vtype=GRB.BINARY, name= 'g') #1 if item i lies on the ground of the bin

beta1 = model.addVars(nitems, nitems, vtype=GRB.BINARY, name='beta1') # 1 if vertex 1 of item i is supported by item j
beta2 = model.addVars(nitems, nitems, vtype=GRB.BINARY, name='beta2') # 1 if vertex 2 of item i is supported by item j
gamma = model.addVars(nitems, vtype=GRB.BINARY, name='gamma') # 1 if vertex 1 of item i is supported by the cut of the bin where it is placed


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
    model.addConstr(r[i, 0, 0] <= lplus[i], name=f"OrientationLength_{i}")
    model.addConstr(r[i, 1, 1] <= hplus[i], name=f"OrientationHeight_{i}")

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
#constraint 16 from lecture
for i in range(nitems):
    for j in range(mbins):
        model.addConstr(
            u_j[j] >= p_ij[i, j]
        )

#Constraint for stability
for i in range(nitems):
    model.addConstr(
        gamma[i] + quicksum(beta1[i, j] for j in range(mbins)) + quicksum(beta2[i, j] for j in range(mbins)) + 2*g[i] >= 2
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


model.optimize()

if model.status == GRB.INFEASIBLE:
    print("The model is infeasible. Computing IIS...")
    model.computeIIS()
    model.write("infeasible.ilp")



import matplotlib.pyplot as plt
import matplotlib.patches as patches


def visualize_with_overlap(nitems, mbins, Lj, Hj, xi, zi, x_i_prime, z_i_prime, p_ij):
    fig, axs = plt.subplots(1, mbins, figsize=(15, 5))

    # Ensure the model is optimized before visualization
    if model.status == GRB.OPTIMAL or model.status == GRB.SUBOPTIMAL:
        for j in range(mbins):
            axs[j].set_xlim(0, Lj[j])
            axs[j].set_ylim(0, Hj[j])
            axs[j].set_title(f"Bin {j}")
            axs[j].set_aspect('equal')

            bin_items = []
            for i in range(nitems):
                if p_ij[i, j].X > 0.5:  # Only visualize items that are assigned to bin j
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
                axs[j].text(x + w / 2, z + h / 2, f"Item {item}", ha='center', va='center')

        plt.tight_layout()
        plt.show()
    else:
        print("Model is not yet solved, please run the model optimization first.")


# Call the function
if model.status == GRB.OPTIMAL or model.status == GRB.SUBOPTIMAL:
    visualize_with_overlap(nitems, mbins, Lj, Hj, x, z, xprime, zprime, p_ij)
else:
    print("Model didn't find a solution within the time limit.")
