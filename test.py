import pickle
from gurobipy import Model, GRB, quicksum, gurobi
import os
import numpy as np
<<<<<<< Updated upstream

# Specify the full path to your Gurobi license file
#gurobi_license_path = "/Users/mariannapiperigou/Documents/gurobi.lic"  # marianna
#gurobi_license_path = "C:/Users/Jacob/OneDrive - Delft University of Technology/Documents/gurobi.lic"  # job
os.environ['GRB_LICENSE_FILE'] = '/Users/julia/Documents/OperationsOptimisation/gurobi.lic'
# Set the environment variable
#os.environ["GRB_LICENSE_FILE"] = gurobi_license_path

=======
os.environ['GRB_LICENSE_FILE'] = '/Users/julia/Documents/OperationsOptimisation/gurobi.lic'
>>>>>>> Stashed changes
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

<<<<<<< Updated upstream
mbins = int(sum(entry[1][2] for entry in bins.values())/2)  # number of bins -- should be halved i think
<<<<<<< Updated upstream
nitems = 10                                      # number of items --> why???
=======
nitems = 8                               # number of items --> why???
>>>>>>> Stashed changes
=======
mbins = len(bins)  # number of bins -- should be halved i think
nitems = 25# number of items --> why???
n_axes = 2                                      # number of axes
n_orients = 2                                   # number of different sides/orientations of an item
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
<<<<<<< Updated upstream
model.setParam('TimeLimit', 60*5)
=======
model.setParam('TimeLimit', 3*60)
>>>>>>> Stashed changes
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

<<<<<<< Updated upstream
#constraint 3: area of items not larger than area of bin
for i in range(nitems):
    for j in range(mbins):
        model.addConstr(
            ai[i] * p_ij[i, j] <= Aj[j] * u_j[j],
            name=f"AreaConstraint_{i}_{j}"
        )
=======
""" Geometric Constraints """

# constraint 3: area (instead of mass) of items not larger than area of bin
###for i in range(nitems):
    #for j in range(mbins):
     #   model.addConstr(
      #      ai[i] * p_ij[i, j] <= Aj[j] * u_j[j],
       #     name=f"AreaConstraint_{i}_{j}")

for j in range(mbins):
    model.addConstr(quicksum(ai[i] * p_ij[i, j] for i in range(nitems)) <= Aj[j] * u_j[j], name=f"TotalBinCapacity_{j}")
>>>>>>> Stashed changes

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
<<<<<<< Updated upstream
    model.addConstr(r[i, 0, 0] + r[i, 1, 0] == 1, name=f"LengthAlign_{i}")
=======
    for d in range(n_orients):
        model.addConstr(quicksum(r[i, c, d] for c in range(n_axes)) == 1, name=f"SidesAlign_{d}_{i}")
    # Each axis aligns with exactly one side
    for c in range(n_axes):
        model.addConstr(quicksum(r[i, c, d] for d in range(n_orients)) == 1, name=f"AxesAlign_{c}_{i}")
    """model.addConstr(r[i, 0, 0] + r[i, 1, 0] == 1, name=f"LengthAlign_{i}")
>>>>>>> Stashed changes
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

<<<<<<< Updated upstream
#constraint 16 from lecture (flagging constraint)
for i in range(nitems):
    for j in range(mbins):
=======
beta1 = model.addVars(nitems, nitems, vtype=GRB.BINARY, name='beta1')   # 1 if vertex 1 of item i is supported by item j
beta2 = model.addVars(nitems, nitems, vtype=GRB.BINARY, name='beta2')   # 1 if vertex 2 of item i is supported by item j
gamma = model.addVars(nitems, vtype=GRB.BINARY, name='gamma')           # 1 if vertex 1 of item i is supported by the cut of the bin where it is placed

v = model.addVars(nitems, nitems, vtype=GRB.CONTINUOUS, name='k')   # represents absolute value of z_hi[k] - z_lo[i]
m = model.addVars(nitems, nitems, vtype=GRB.BINARY, name='m')   # 1 if k overlaps or exceeds i in height, i.e. z_hi[k] > z_lo[i]

#'''Vertical Stability Constraints'''

#Constraint 26 adapted for 2D stability
for i in range(nitems):
    model.addConstr(
        gamma[i] + quicksum(beta1[i, j] for j in range(nitems) if i!= j) + quicksum(beta2[i, j] for j in range(nitems)if i!= j) + 2*g[i]
        >= 2
    )

for i in range(nitems):
    # If item i is on the ground, z_lo must be small or equal to 0
    model.addConstr(z_lo[i] <= (1 - g[i]) * H, name=f"GroundConstraint_{i}")
    for k in range(nitems):
        if i != k:
            # 28&29: absolute value of height diff items i and k, v_ik can only be one if z coords at same height
            model.addConstr(z_hi[k] - z_lo[i] <= v[i,k], name=f"AbsZ1{i}_{k}")
            model.addConstr(z_lo[i] - z_hi[k] <= v[i,k], name=f"AbsZ2{i}_{k}")

            # 30: forces m, m must be 0 if z_lo[i] larger than z_hi[k] as this means i above k
            model.addConstr(v[i, k] <= z_hi[k] - z_lo[i] + 2 * H * (1 - m[i, k]), name=f"{k}_Below_{i}")
            # 31: forces m, m must be 1 if z_hi[k] larger than z_lo[i] as this means i not above k
            model.addConstr(v[i, k] <= z_lo[i] - z_hi[k] + 2 * H * m[i, k], name=f"{k}_NotBelow_{i}")

            # 32: forces h, if v_ik is 0 then h must be zero (i and k height compatibility)
            model.addConstr(h[i, k] <= v[i, k], name=f"HeightCompat_{i}_{k}")
            # 33: forces h, if v_ik is larger than 0, then h must be 1 (no height comp.)
            model.addConstr(v[i, k] <= h[i, k] * H, name=f"NoHeightCompat_{i}_{k}")

            ####34 might need revising
            # 34.1: forces o, if horizontal overlap, then o is 0 because neither i or k are to the right of each other
            model.addConstr(o[i, k] <= xp[i, k] + xp[k, i], name=f"Overlap_{i}_{k}")
            # 34.2: forces o, if no horizontal overlap, o must be 1 because either i or k are to right of each other
            model.addConstr(xp[i, k] + xp[k, i] <= 2 * o[i, k], name=f"NoOverlap_{i}_{k}")

            # 35.1: forces s=1, if hor. inters. (o=0) and/or suitable height (h=0), s must be 1 (k supports i)
            model.addConstr((1-s[i, k]) <= h[i, k] + o[i, k], name=f"{j}_Supports_{i}_{k}")
            # 35.2:forces s=0, if hor. inters. (o=1) and suitable height (h=1), s must be 0 (k no support i)
            model.addConstr(h[i, k] + o[i, k] <= 2 * (1 - s[i, k]), name=f"{j}_DoesntSupport_{i}_{k}")

            # 37.1: if beta1 = 1, s automatically also becomes 1
            model.addConstr(beta1[i,k] <= s[i,k], name=f"beta1_{i}_{k}")
            # 37.2: if beta2 = 1, s automatically also becomes 1
            model.addConstr(beta2[i,k] <= s[i,k], name=f"beta2_{i}_{k}")
            # 38: forces s, if beta1 or beta2 is 1, s is flagged to 1
            model.addConstr(beta1[i, k] <= s[i, k], name=f"sBeta1Flag_{i}_{k}")
            model.addConstr(beta2[i, k] <= s[i, k], name=f"sBeta2Flag_{i}_{k}")
            ####39 might need revising, might not even be necessary due to 34
            # 39.1: if i supports j (beta=1), x_l[i]>x_l[j] (eta1 = 0) and x_r[k] > x_l[i]
            model.addConstr(eta1[i, k] <= 1 - beta1[i, k], name=f"Eta1Overlap_{i}_{k}")
            model.addConstr(eta2[i, k] <= 1 - beta2[i, k], name=f"Eta2Overlap_{i}_{k}")
            # 43: forces eta1 to be 1 if x_l[i] is smaller than x_l[k] (k cannot support vertex 1 of i)
            model.addConstr(x_l[k] <= x_l[i] + eta1[i, k] * L, name=f"Eta1Flag_{i}_{k}")
            # 45: forces eta2 to be 1 if x_r[k] is smaller than x_r[i] (k cannot support vertex 2 of i)
            model.addConstr(x_r[i] <= x_r[k] + eta2[i, k] * L, name=f"Eta2Flag_{i}_{k}")
            for j in range(mbins):
                # 36: ensures s holds only for item i and item k in same bin j
                model.addConstr(p_ij[i, j] - p_ij[k, j] <= 1-s[i, k], name=f"{k}_Supports_{i}_InSameBin1_{j}")
                model.addConstr(p_ij[k, j] - p_ij[i, j] <= 1-s[i, k], name=f"{k}_Supports_{i}_InSameBin2_{j}")

            # 49: forces gamma to be 0 if there is no cut (technically gamma could be 1 for no-cut ULD as a=b=-1)
                model.addConstr(p_ij[i,j] + gamma[i] <= a[j] + b[j] + 3, name=f"Gamma0ForItems_{i}_inBins_{j}_NoCut_{k}")
        for j in indices_with_cut:
            # 47: forces gamma if item i is on cut and if i is in bin j
            model.addConstr(z_lo[i] + b[j] / a[j] * x_l[i] - b[j] <= (1 - gamma[i]) * M + (1 - p_ij[i, j])*M,
                            name=f"{i}_OnCutIn{j}_{k}")
            #model.addConstr(z_lo[i] + (b[j] / a[j]) * x_l[i] - b[j] >= -M * (1 - p_ij[i,j])
""" Other Constraints """
#constraint 16 from lecture (flagging constraint)
for i in range(nitems):
    for j in range(mbins):
        model.addConstr(u_j[j] >= p_ij[i, j], name=f"FlaggingConstraint_{i}_{j}")

#some items might be fragile and, as such, no other box can be stacked on top of them
for i in range(nitems):
    for k in range(nitems):
>>>>>>> Stashed changes
        model.addConstr(
            u_j[j] >= p_ij[i, j]
        )

<<<<<<< Updated upstream
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

=======
# Ensure that a ULD cannot contain both perishable and radioactive items


for j in range(mbins):
    model.addConstr(
        quicksum(p_ij[i, j] * perishable[i] for i in range(nitems)) +
        quicksum(p_ij[i, j] * radioactive[i] for i in range(nitems)) <= 1,
        name=f"Perishable_radioactive_{j}"
    )
>>>>>>> Stashed changes
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
<<<<<<< Updated upstream
=======
    for constr in model.getConstrs():
        if constr.IISConstr:
            print(f"Infeasible Constraint: {constr.ConstrName}")
    model.write("infeasible.lp")
    model.computeIIS()
>>>>>>> Stashed changes
    model.write("infeasible.ilp")

else:
    model.write("feasible.lp")

import matplotlib.pyplot as plt
import numpy as np

def visualize_with_overlap(items, nitems, mbins, Lj, Hj, xi, zi, x_i_prime, z_i_prime, p_ij, bins_with_cut, a, b):
    fig, axs = plt.subplots(1, mbins, figsize=(15, 5))
    print("Solution count:", model.SolCount)
    # Ensure the model is optimized before visualization
    #if model.status in [GRB.OPTIMAL, GRB.SUBOPTIMAL, GRB.TIME_LIMIT]:
    if model.SolCount >= 0:
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
<<<<<<< Updated upstream
if model.status == GRB.OPTIMAL or model.status == GRB.SUBOPTIMAL:
    visualize_with_overlap(items, nitems, mbins, Lj, Hj, x, z, xprime, zprime, p_ij, indices_with_cut, a, b)
=======
if model.status == GRB.OPTIMAL or model.status == GRB.SUBOPTIMAL or model.status == GRB.TIME_LIMIT:
    print("Visualize")
    visualize_with_overlap(items, nitems, mbins, Lj, Hj, x_l, z_lo, x_r, z_hi, p_ij, indices_with_cut, a, b)
>>>>>>> Stashed changes
else:
    print("Model didn't find a solution within the time limit, no visualization.")
    print(model.status)

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
        #if u_j[j].X == 1:
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
