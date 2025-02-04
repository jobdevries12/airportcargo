\ Model 2DBPP
\ LP format - for model browsing. Use MPS format to capture full model detail.
\ Signature: 0xff20ee2af7d84391
Minimize
  200 u_j[0] + 200 u_j[1] + 150 u_j[2] + 150 u_j[3] + 0 x_p[0,0]
   + 0 x_p[1,1] + 0 x_p[2,2] + 0 x_p[3,3] + 0 x_p[4,4] + 0 z_p[0,0]
   + 0 z_p[1,1] + 0 z_p[2,2] + 0 z_p[3,3] + 0 z_p[4,4] + 0 rho[0]
   + 0 rho[1] + 0 rho[2] + 0 rho[3] + 0 rho[4] + 0 beta1[0,0]
   + 0 beta1[1,1] + 0 beta1[2,2] + 0 beta1[3,3] + 0 beta1[4,4]
   + 0 beta2[0,0] + 0 beta2[1,1] + 0 beta2[2,2] + 0 beta2[3,3]
   + 0 beta2[4,4]
Subject To
 GroundSupport_UpperBound_0: zi[0] + 100000 g[0] <= 100000
 CutSupport_UpperBound_0_2: 1.261904761904762 xi[0] + zi[0]
   + 100000 gamma[0] <= 100053
 CutSupport_LowerBound_0_2: 1.261904761904762 xi[0] + zi[0]
   - 100000 gamma[0] >= -99947
 CutSupport_UpperBound_0_3: 1.261904761904762 xi[0] + zi[0]
   + 100000 gamma[0] <= 100053
 CutSupport_LowerBound_0_3: 1.261904761904762 xi[0] + zi[0]
   - 100000 gamma[0] >= -99947
 Beta1_0_1: xi[0] - x_i_prime[1] + 100000 beta1[0,1] <= 100000
 Beta1_Z_Upper_0_1: zi[0] - z_i_prime[1] + 100000 beta1[0,1] <= 100000
 Beta1_Z_Lower_0_1: zi[0] - z_i_prime[1] - 100000 beta1[0,1] >= -100000
 Beta2_0_1: - xi[1] + x_i_prime[0] - 100000 beta2[0,1] >= -100000
 Beta2_Z_Upper_0_1: - zi[1] + z_i_prime[0] + 100000 beta2[0,1] <= 100000
 Beta2_Z_Lower_0_1: - zi[1] + z_i_prime[0] - 100000 beta2[0,1] >= -100000
 Beta1_0_2: xi[0] - x_i_prime[2] + 100000 beta1[0,2] <= 100000
 Beta1_Z_Upper_0_2: zi[0] - z_i_prime[2] + 100000 beta1[0,2] <= 100000
 Beta1_Z_Lower_0_2: zi[0] - z_i_prime[2] - 100000 beta1[0,2] >= -100000
 Beta2_0_2: - xi[2] + x_i_prime[0] - 100000 beta2[0,2] >= -100000
 Beta2_Z_Upper_0_2: - zi[2] + z_i_prime[0] + 100000 beta2[0,2] <= 100000
 Beta2_Z_Lower_0_2: - zi[2] + z_i_prime[0] - 100000 beta2[0,2] >= -100000
 Beta1_0_3: xi[0] - x_i_prime[3] + 100000 beta1[0,3] <= 100000
 Beta1_Z_Upper_0_3: zi[0] - z_i_prime[3] + 100000 beta1[0,3] <= 100000
 Beta1_Z_Lower_0_3: zi[0] - z_i_prime[3] - 100000 beta1[0,3] >= -100000
 Beta2_0_3: - xi[3] + x_i_prime[0] - 100000 beta2[0,3] >= -100000
 Beta2_Z_Upper_0_3: - zi[3] + z_i_prime[0] + 100000 beta2[0,3] <= 100000
 Beta2_Z_Lower_0_3: - zi[3] + z_i_prime[0] - 100000 beta2[0,3] >= -100000
 Beta1_0_4: xi[0] - x_i_prime[4] + 100000 beta1[0,4] <= 100000
 Beta1_Z_Upper_0_4: zi[0] - z_i_prime[4] + 100000 beta1[0,4] <= 100000
 Beta1_Z_Lower_0_4: zi[0] - z_i_prime[4] - 100000 beta1[0,4] >= -100000
 Beta2_0_4: - xi[4] + x_i_prime[0] - 100000 beta2[0,4] >= -100000
 Beta2_Z_Upper_0_4: - zi[4] + z_i_prime[0] + 100000 beta2[0,4] <= 100000
 Beta2_Z_Lower_0_4: - zi[4] + z_i_prime[0] - 100000 beta2[0,4] >= -100000
 GroundSupport_UpperBound_1: zi[1] + 100000 g[1] <= 100000
 CutSupport_UpperBound_1_2: 1.261904761904762 xi[1] + zi[1]
   + 100000 gamma[1] <= 100053
 CutSupport_LowerBound_1_2: 1.261904761904762 xi[1] + zi[1]
   - 100000 gamma[1] >= -99947
 CutSupport_UpperBound_1_3: 1.261904761904762 xi[1] + zi[1]
   + 100000 gamma[1] <= 100053
 CutSupport_LowerBound_1_3: 1.261904761904762 xi[1] + zi[1]
   - 100000 gamma[1] >= -99947
 Beta1_1_0: xi[1] - x_i_prime[0] + 100000 beta1[1,0] <= 100000
 Beta1_Z_Upper_1_0: zi[1] - z_i_prime[0] + 100000 beta1[1,0] <= 100000
 Beta1_Z_Lower_1_0: zi[1] - z_i_prime[0] - 100000 beta1[1,0] >= -100000
 Beta2_1_0: - xi[0] + x_i_prime[1] - 100000 beta2[1,0] >= -100000
 Beta2_Z_Upper_1_0: - zi[0] + z_i_prime[1] + 100000 beta2[1,0] <= 100000
 Beta2_Z_Lower_1_0: - zi[0] + z_i_prime[1] - 100000 beta2[1,0] >= -100000
 Beta1_1_2: xi[1] - x_i_prime[2] + 100000 beta1[1,2] <= 100000
 Beta1_Z_Upper_1_2: zi[1] - z_i_prime[2] + 100000 beta1[1,2] <= 100000
 Beta1_Z_Lower_1_2: zi[1] - z_i_prime[2] - 100000 beta1[1,2] >= -100000
 Beta2_1_2: - xi[2] + x_i_prime[1] - 100000 beta2[1,2] >= -100000
 Beta2_Z_Upper_1_2: - zi[2] + z_i_prime[1] + 100000 beta2[1,2] <= 100000
 Beta2_Z_Lower_1_2: - zi[2] + z_i_prime[1] - 100000 beta2[1,2] >= -100000
 Beta1_1_3: xi[1] - x_i_prime[3] + 100000 beta1[1,3] <= 100000
 Beta1_Z_Upper_1_3: zi[1] - z_i_prime[3] + 100000 beta1[1,3] <= 100000
 Beta1_Z_Lower_1_3: zi[1] - z_i_prime[3] - 100000 beta1[1,3] >= -100000
 Beta2_1_3: - xi[3] + x_i_prime[1] - 100000 beta2[1,3] >= -100000
 Beta2_Z_Upper_1_3: - zi[3] + z_i_prime[1] + 100000 beta2[1,3] <= 100000
 Beta2_Z_Lower_1_3: - zi[3] + z_i_prime[1] - 100000 beta2[1,3] >= -100000
 Beta1_1_4: xi[1] - x_i_prime[4] + 100000 beta1[1,4] <= 100000
 Beta1_Z_Upper_1_4: zi[1] - z_i_prime[4] + 100000 beta1[1,4] <= 100000
 Beta1_Z_Lower_1_4: zi[1] - z_i_prime[4] - 100000 beta1[1,4] >= -100000
 Beta2_1_4: - xi[4] + x_i_prime[1] - 100000 beta2[1,4] >= -100000
 Beta2_Z_Upper_1_4: - zi[4] + z_i_prime[1] + 100000 beta2[1,4] <= 100000
 Beta2_Z_Lower_1_4: - zi[4] + z_i_prime[1] - 100000 beta2[1,4] >= -100000
 GroundSupport_UpperBound_2: zi[2] + 100000 g[2] <= 100000
 CutSupport_UpperBound_2_2: 1.261904761904762 xi[2] + zi[2]
   + 100000 gamma[2] <= 100053
 CutSupport_LowerBound_2_2: 1.261904761904762 xi[2] + zi[2]
   - 100000 gamma[2] >= -99947
 CutSupport_UpperBound_2_3: 1.261904761904762 xi[2] + zi[2]
   + 100000 gamma[2] <= 100053
 CutSupport_LowerBound_2_3: 1.261904761904762 xi[2] + zi[2]
   - 100000 gamma[2] >= -99947
 Beta1_2_0: xi[2] - x_i_prime[0] + 100000 beta1[2,0] <= 100000
 Beta1_Z_Upper_2_0: zi[2] - z_i_prime[0] + 100000 beta1[2,0] <= 100000
 Beta1_Z_Lower_2_0: zi[2] - z_i_prime[0] - 100000 beta1[2,0] >= -100000
 Beta2_2_0: - xi[0] + x_i_prime[2] - 100000 beta2[2,0] >= -100000
 Beta2_Z_Upper_2_0: - zi[0] + z_i_prime[2] + 100000 beta2[2,0] <= 100000
 Beta2_Z_Lower_2_0: - zi[0] + z_i_prime[2] - 100000 beta2[2,0] >= -100000
 Beta1_2_1: xi[2] - x_i_prime[1] + 100000 beta1[2,1] <= 100000
 Beta1_Z_Upper_2_1: zi[2] - z_i_prime[1] + 100000 beta1[2,1] <= 100000
 Beta1_Z_Lower_2_1: zi[2] - z_i_prime[1] - 100000 beta1[2,1] >= -100000
 Beta2_2_1: - xi[1] + x_i_prime[2] - 100000 beta2[2,1] >= -100000
 Beta2_Z_Upper_2_1: - zi[1] + z_i_prime[2] + 100000 beta2[2,1] <= 100000
 Beta2_Z_Lower_2_1: - zi[1] + z_i_prime[2] - 100000 beta2[2,1] >= -100000
 Beta1_2_3: xi[2] - x_i_prime[3] + 100000 beta1[2,3] <= 100000
 Beta1_Z_Upper_2_3: zi[2] - z_i_prime[3] + 100000 beta1[2,3] <= 100000
 Beta1_Z_Lower_2_3: zi[2] - z_i_prime[3] - 100000 beta1[2,3] >= -100000
 Beta2_2_3: - xi[3] + x_i_prime[2] - 100000 beta2[2,3] >= -100000
 Beta2_Z_Upper_2_3: - zi[3] + z_i_prime[2] + 100000 beta2[2,3] <= 100000
 Beta2_Z_Lower_2_3: - zi[3] + z_i_prime[2] - 100000 beta2[2,3] >= -100000
 Beta1_2_4: xi[2] - x_i_prime[4] + 100000 beta1[2,4] <= 100000
 Beta1_Z_Upper_2_4: zi[2] - z_i_prime[4] + 100000 beta1[2,4] <= 100000
 Beta1_Z_Lower_2_4: zi[2] - z_i_prime[4] - 100000 beta1[2,4] >= -100000
 Beta2_2_4: - xi[4] + x_i_prime[2] - 100000 beta2[2,4] >= -100000
 Beta2_Z_Upper_2_4: - zi[4] + z_i_prime[2] + 100000 beta2[2,4] <= 100000
 Beta2_Z_Lower_2_4: - zi[4] + z_i_prime[2] - 100000 beta2[2,4] >= -100000
 GroundSupport_UpperBound_3: zi[3] + 100000 g[3] <= 100000
 CutSupport_UpperBound_3_2: 1.261904761904762 xi[3] + zi[3]
   + 100000 gamma[3] <= 100053
 CutSupport_LowerBound_3_2: 1.261904761904762 xi[3] + zi[3]
   - 100000 gamma[3] >= -99947
 CutSupport_UpperBound_3_3: 1.261904761904762 xi[3] + zi[3]
   + 100000 gamma[3] <= 100053
 CutSupport_LowerBound_3_3: 1.261904761904762 xi[3] + zi[3]
   - 100000 gamma[3] >= -99947
 Beta1_3_0: xi[3] - x_i_prime[0] + 100000 beta1[3,0] <= 100000
 Beta1_Z_Upper_3_0: zi[3] - z_i_prime[0] + 100000 beta1[3,0] <= 100000
 Beta1_Z_Lower_3_0: zi[3] - z_i_prime[0] - 100000 beta1[3,0] >= -100000
 Beta2_3_0: - xi[0] + x_i_prime[3] - 100000 beta2[3,0] >= -100000
 Beta2_Z_Upper_3_0: - zi[0] + z_i_prime[3] + 100000 beta2[3,0] <= 100000
 Beta2_Z_Lower_3_0: - zi[0] + z_i_prime[3] - 100000 beta2[3,0] >= -100000
 Beta1_3_1: xi[3] - x_i_prime[1] + 100000 beta1[3,1] <= 100000
 Beta1_Z_Upper_3_1: zi[3] - z_i_prime[1] + 100000 beta1[3,1] <= 100000
 Beta1_Z_Lower_3_1: zi[3] - z_i_prime[1] - 100000 beta1[3,1] >= -100000
 Beta2_3_1: - xi[1] + x_i_prime[3] - 100000 beta2[3,1] >= -100000
 Beta2_Z_Upper_3_1: - zi[1] + z_i_prime[3] + 100000 beta2[3,1] <= 100000
 Beta2_Z_Lower_3_1: - zi[1] + z_i_prime[3] - 100000 beta2[3,1] >= -100000
 Beta1_3_2: xi[3] - x_i_prime[2] + 100000 beta1[3,2] <= 100000
 Beta1_Z_Upper_3_2: zi[3] - z_i_prime[2] + 100000 beta1[3,2] <= 100000
 Beta1_Z_Lower_3_2: zi[3] - z_i_prime[2] - 100000 beta1[3,2] >= -100000
 Beta2_3_2: - xi[2] + x_i_prime[3] - 100000 beta2[3,2] >= -100000
 Beta2_Z_Upper_3_2: - zi[2] + z_i_prime[3] + 100000 beta2[3,2] <= 100000
 Beta2_Z_Lower_3_2: - zi[2] + z_i_prime[3] - 100000 beta2[3,2] >= -100000
 Beta1_3_4: xi[3] - x_i_prime[4] + 100000 beta1[3,4] <= 100000
 Beta1_Z_Upper_3_4: zi[3] - z_i_prime[4] + 100000 beta1[3,4] <= 100000
 Beta1_Z_Lower_3_4: zi[3] - z_i_prime[4] - 100000 beta1[3,4] >= -100000
 Beta2_3_4: - xi[4] + x_i_prime[3] - 100000 beta2[3,4] >= -100000
 Beta2_Z_Upper_3_4: - zi[4] + z_i_prime[3] + 100000 beta2[3,4] <= 100000
 Beta2_Z_Lower_3_4: - zi[4] + z_i_prime[3] - 100000 beta2[3,4] >= -100000
 GroundSupport_UpperBound_4: zi[4] + 100000 g[4] <= 100000
 CutSupport_UpperBound_4_2: 1.261904761904762 xi[4] + zi[4]
   + 100000 gamma[4] <= 100053
 CutSupport_LowerBound_4_2: 1.261904761904762 xi[4] + zi[4]
   - 100000 gamma[4] >= -99947
 CutSupport_UpperBound_4_3: 1.261904761904762 xi[4] + zi[4]
   + 100000 gamma[4] <= 100053
 CutSupport_LowerBound_4_3: 1.261904761904762 xi[4] + zi[4]
   - 100000 gamma[4] >= -99947
 Beta1_4_0: xi[4] - x_i_prime[0] + 100000 beta1[4,0] <= 100000
 Beta1_Z_Upper_4_0: zi[4] - z_i_prime[0] + 100000 beta1[4,0] <= 100000
 Beta1_Z_Lower_4_0: zi[4] - z_i_prime[0] - 100000 beta1[4,0] >= -100000
 Beta2_4_0: - xi[0] + x_i_prime[4] - 100000 beta2[4,0] >= -100000
 Beta2_Z_Upper_4_0: - zi[0] + z_i_prime[4] + 100000 beta2[4,0] <= 100000
 Beta2_Z_Lower_4_0: - zi[0] + z_i_prime[4] - 100000 beta2[4,0] >= -100000
 Beta1_4_1: xi[4] - x_i_prime[1] + 100000 beta1[4,1] <= 100000
 Beta1_Z_Upper_4_1: zi[4] - z_i_prime[1] + 100000 beta1[4,1] <= 100000
 Beta1_Z_Lower_4_1: zi[4] - z_i_prime[1] - 100000 beta1[4,1] >= -100000
 Beta2_4_1: - xi[1] + x_i_prime[4] - 100000 beta2[4,1] >= -100000
 Beta2_Z_Upper_4_1: - zi[1] + z_i_prime[4] + 100000 beta2[4,1] <= 100000
 Beta2_Z_Lower_4_1: - zi[1] + z_i_prime[4] - 100000 beta2[4,1] >= -100000
 Beta1_4_2: xi[4] - x_i_prime[2] + 100000 beta1[4,2] <= 100000
 Beta1_Z_Upper_4_2: zi[4] - z_i_prime[2] + 100000 beta1[4,2] <= 100000
 Beta1_Z_Lower_4_2: zi[4] - z_i_prime[2] - 100000 beta1[4,2] >= -100000
 Beta2_4_2: - xi[2] + x_i_prime[4] - 100000 beta2[4,2] >= -100000
 Beta2_Z_Upper_4_2: - zi[2] + z_i_prime[4] + 100000 beta2[4,2] <= 100000
 Beta2_Z_Lower_4_2: - zi[2] + z_i_prime[4] - 100000 beta2[4,2] >= -100000
 Beta1_4_3: xi[4] - x_i_prime[3] + 100000 beta1[4,3] <= 100000
 Beta1_Z_Upper_4_3: zi[4] - z_i_prime[3] + 100000 beta1[4,3] <= 100000
 Beta1_Z_Lower_4_3: zi[4] - z_i_prime[3] - 100000 beta1[4,3] >= -100000
 Beta2_4_3: - xi[3] + x_i_prime[4] - 100000 beta2[4,3] >= -100000
 Beta2_Z_Upper_4_3: - zi[3] + z_i_prime[4] + 100000 beta2[4,3] <= 100000
 Beta2_Z_Lower_4_3: - zi[3] + z_i_prime[4] - 100000 beta2[4,3] >= -100000
 AreaConstraint_0_0: 1820 p_ij[0,0] - 46500 u_j[0] <= 0
 AreaConstraint_0_1: 1820 p_ij[0,1] - 46500 u_j[1] <= 0
 AreaConstraint_0_2: 1820 p_ij[0,2] - 29760 u_j[2] <= 0
 AreaConstraint_0_3: 1820 p_ij[0,3] - 29760 u_j[3] <= 0
 AreaConstraint_1_0: 2240 p_ij[1,0] - 46500 u_j[0] <= 0
 AreaConstraint_1_1: 2240 p_ij[1,1] - 46500 u_j[1] <= 0
 AreaConstraint_1_2: 2240 p_ij[1,2] - 29760 u_j[2] <= 0
 AreaConstraint_1_3: 2240 p_ij[1,3] - 29760 u_j[3] <= 0
 AreaConstraint_2_0: 1696 p_ij[2,0] - 46500 u_j[0] <= 0
 AreaConstraint_2_1: 1696 p_ij[2,1] - 46500 u_j[1] <= 0
 AreaConstraint_2_2: 1696 p_ij[2,2] - 29760 u_j[2] <= 0
 AreaConstraint_2_3: 1696 p_ij[2,3] - 29760 u_j[3] <= 0
 AreaConstraint_3_0: 3168 p_ij[3,0] - 46500 u_j[0] <= 0
 AreaConstraint_3_1: 3168 p_ij[3,1] - 46500 u_j[1] <= 0
 AreaConstraint_3_2: 3168 p_ij[3,2] - 29760 u_j[2] <= 0
 AreaConstraint_3_3: 3168 p_ij[3,3] - 29760 u_j[3] <= 0
 AreaConstraint_4_0: 2640 p_ij[4,0] - 46500 u_j[0] <= 0
 AreaConstraint_4_1: 2640 p_ij[4,1] - 46500 u_j[1] <= 0
 AreaConstraint_4_2: 2640 p_ij[4,2] - 29760 u_j[2] <= 0
 AreaConstraint_4_3: 2640 p_ij[4,3] - 29760 u_j[3] <= 0
 OneItemOneBin_0: p_ij[0,0] + p_ij[0,1] + p_ij[0,2] + p_ij[0,3] = 1
 OneItemOneBin_1: p_ij[1,0] + p_ij[1,1] + p_ij[1,2] + p_ij[1,3] = 1
 OneItemOneBin_2: p_ij[2,0] + p_ij[2,1] + p_ij[2,2] + p_ij[2,3] = 1
 OneItemOneBin_3: p_ij[3,0] + p_ij[3,1] + p_ij[3,2] + p_ij[3,3] = 1
 OneItemOneBin_4: p_ij[4,0] + p_ij[4,1] + p_ij[4,2] + p_ij[4,3] = 1
 ItemFitToBin_0: - 300 p_ij[0,0] - 300 p_ij[0,1] - 192 p_ij[0,2]
   - 192 p_ij[0,3] + x_i_prime[0] <= 0
 ItemFitToBin_1: - 300 p_ij[1,0] - 300 p_ij[1,1] - 192 p_ij[1,2]
   - 192 p_ij[1,3] + x_i_prime[1] <= 0
 ItemFitToBin_2: - 300 p_ij[2,0] - 300 p_ij[2,1] - 192 p_ij[2,2]
   - 192 p_ij[2,3] + x_i_prime[2] <= 0
 ItemFitToBin_3: - 300 p_ij[3,0] - 300 p_ij[3,1] - 192 p_ij[3,2]
   - 192 p_ij[3,3] + x_i_prime[3] <= 0
 ItemFitToBin_4: - 300 p_ij[4,0] - 300 p_ij[4,1] - 192 p_ij[4,2]
   - 192 p_ij[4,3] + x_i_prime[4] <= 0
 ItemFitToBinZ_0: - 155 p_ij[0,0] - 155 p_ij[0,1] - 155 p_ij[0,2]
   - 155 p_ij[0,3] + z_i_prime[0] <= 0
 ItemFitToBinZ_1: - 155 p_ij[1,0] - 155 p_ij[1,1] - 155 p_ij[1,2]
   - 155 p_ij[1,3] + z_i_prime[1] <= 0
 ItemFitToBinZ_2: - 155 p_ij[2,0] - 155 p_ij[2,1] - 155 p_ij[2,2]
   - 155 p_ij[2,3] + z_i_prime[2] <= 0
 ItemFitToBinZ_3: - 155 p_ij[3,0] - 155 p_ij[3,1] - 155 p_ij[3,2]
   - 155 p_ij[3,3] + z_i_prime[3] <= 0
 ItemFitToBinZ_4: - 155 p_ij[4,0] - 155 p_ij[4,1] - 155 p_ij[4,2]
   - 155 p_ij[4,3] + z_i_prime[4] <= 0
 TransformX_0: - xi[0] + x_i_prime[0] - 65 r[0,0,0] - 28 r[0,0,1] = 0
 TransformZ_0: - zi[0] + z_i_prime[0] - 65 r[0,1,0] - 28 r[0,1,1] = 0
 TransformX_1: - xi[1] + x_i_prime[1] - 64 r[1,0,0] - 35 r[1,0,1] = 0
 TransformZ_1: - zi[1] + z_i_prime[1] - 64 r[1,1,0] - 35 r[1,1,1] = 0
 TransformX_2: - xi[2] + x_i_prime[2] - 53 r[2,0,0] - 32 r[2,0,1] = 0
 TransformZ_2: - zi[2] + z_i_prime[2] - 53 r[2,1,0] - 32 r[2,1,1] = 0
 TransformX_3: - xi[3] + x_i_prime[3] - 88 r[3,0,0] - 36 r[3,0,1] = 0
 TransformZ_3: - zi[3] + z_i_prime[3] - 88 r[3,1,0] - 36 r[3,1,1] = 0
 TransformX_4: - xi[4] + x_i_prime[4] - 88 r[4,0,0] - 30 r[4,0,1] = 0
 TransformZ_4: - zi[4] + z_i_prime[4] - 88 r[4,1,0] - 30 r[4,1,1] = 0
 LengthAlign_0: r[0,0,0] + r[0,1,0] = 1
 HeightAlign_0: r[0,0,1] + r[0,1,1] = 1
 XAxisAlign_0: r[0,0,0] + r[0,0,1] = 1
 ZAxisAlign_0: r[0,1,0] + r[0,1,1] = 1
 LengthAlign_1: r[1,0,0] + r[1,1,0] = 1
 HeightAlign_1: r[1,0,1] + r[1,1,1] = 1
 XAxisAlign_1: r[1,0,0] + r[1,0,1] = 1
 ZAxisAlign_1: r[1,1,0] + r[1,1,1] = 1
 LengthAlign_2: r[2,0,0] + r[2,1,0] = 1
 HeightAlign_2: r[2,0,1] + r[2,1,1] = 1
 XAxisAlign_2: r[2,0,0] + r[2,0,1] = 1
 ZAxisAlign_2: r[2,1,0] + r[2,1,1] = 1
 LengthAlign_3: r[3,0,0] + r[3,1,0] = 1
 HeightAlign_3: r[3,0,1] + r[3,1,1] = 1
 XAxisAlign_3: r[3,0,0] + r[3,0,1] = 1
 ZAxisAlign_3: r[3,1,0] + r[3,1,1] = 1
 LengthAlign_4: r[4,0,0] + r[4,1,0] = 1
 HeightAlign_4: r[4,0,1] + r[4,1,1] = 1
 XAxisAlign_4: r[4,0,0] + r[4,0,1] = 1
 ZAxisAlign_4: r[4,1,0] + r[4,1,1] = 1
 OverlapInBin_0_1_toBin_0: - p_ij[0,0] - p_ij[1,0] + x_p[0,1] + x_p[1,0]
   + z_p[0,1] + z_p[1,0] >= -1
 OverlapInBin_0_1_toBin_1: - p_ij[0,1] - p_ij[1,1] + x_p[0,1] + x_p[1,0]
   + z_p[0,1] + z_p[1,0] >= -1
 OverlapInBin_0_1_toBin_2: - p_ij[0,2] - p_ij[1,2] + x_p[0,1] + x_p[1,0]
   + z_p[0,1] + z_p[1,0] >= -1
 OverlapInBin_0_1_toBin_3: - p_ij[0,3] - p_ij[1,3] + x_p[0,1] + x_p[1,0]
   + z_p[0,1] + z_p[1,0] >= -1
 OverlapInBin_0_2_toBin_0: - p_ij[0,0] - p_ij[2,0] + x_p[0,2] + x_p[2,0]
   + z_p[0,2] + z_p[2,0] >= -1
 OverlapInBin_0_2_toBin_1: - p_ij[0,1] - p_ij[2,1] + x_p[0,2] + x_p[2,0]
   + z_p[0,2] + z_p[2,0] >= -1
 OverlapInBin_0_2_toBin_2: - p_ij[0,2] - p_ij[2,2] + x_p[0,2] + x_p[2,0]
   + z_p[0,2] + z_p[2,0] >= -1
 OverlapInBin_0_2_toBin_3: - p_ij[0,3] - p_ij[2,3] + x_p[0,2] + x_p[2,0]
   + z_p[0,2] + z_p[2,0] >= -1
 OverlapInBin_0_3_toBin_0: - p_ij[0,0] - p_ij[3,0] + x_p[0,3] + x_p[3,0]
   + z_p[0,3] + z_p[3,0] >= -1
 OverlapInBin_0_3_toBin_1: - p_ij[0,1] - p_ij[3,1] + x_p[0,3] + x_p[3,0]
   + z_p[0,3] + z_p[3,0] >= -1
 OverlapInBin_0_3_toBin_2: - p_ij[0,2] - p_ij[3,2] + x_p[0,3] + x_p[3,0]
   + z_p[0,3] + z_p[3,0] >= -1
 OverlapInBin_0_3_toBin_3: - p_ij[0,3] - p_ij[3,3] + x_p[0,3] + x_p[3,0]
   + z_p[0,3] + z_p[3,0] >= -1
 OverlapInBin_0_4_toBin_0: - p_ij[0,0] - p_ij[4,0] + x_p[0,4] + x_p[4,0]
   + z_p[0,4] + z_p[4,0] >= -1
 OverlapInBin_0_4_toBin_1: - p_ij[0,1] - p_ij[4,1] + x_p[0,4] + x_p[4,0]
   + z_p[0,4] + z_p[4,0] >= -1
 OverlapInBin_0_4_toBin_2: - p_ij[0,2] - p_ij[4,2] + x_p[0,4] + x_p[4,0]
   + z_p[0,4] + z_p[4,0] >= -1
 OverlapInBin_0_4_toBin_3: - p_ij[0,3] - p_ij[4,3] + x_p[0,4] + x_p[4,0]
   + z_p[0,4] + z_p[4,0] >= -1
 OverlapInBin_1_0_toBin_0: - p_ij[0,0] - p_ij[1,0] + x_p[0,1] + x_p[1,0]
   + z_p[0,1] + z_p[1,0] >= -1
 OverlapInBin_1_0_toBin_1: - p_ij[0,1] - p_ij[1,1] + x_p[0,1] + x_p[1,0]
   + z_p[0,1] + z_p[1,0] >= -1
 OverlapInBin_1_0_toBin_2: - p_ij[0,2] - p_ij[1,2] + x_p[0,1] + x_p[1,0]
   + z_p[0,1] + z_p[1,0] >= -1
 OverlapInBin_1_0_toBin_3: - p_ij[0,3] - p_ij[1,3] + x_p[0,1] + x_p[1,0]
   + z_p[0,1] + z_p[1,0] >= -1
 OverlapInBin_1_2_toBin_0: - p_ij[1,0] - p_ij[2,0] + x_p[1,2] + x_p[2,1]
   + z_p[1,2] + z_p[2,1] >= -1
 OverlapInBin_1_2_toBin_1: - p_ij[1,1] - p_ij[2,1] + x_p[1,2] + x_p[2,1]
   + z_p[1,2] + z_p[2,1] >= -1
 OverlapInBin_1_2_toBin_2: - p_ij[1,2] - p_ij[2,2] + x_p[1,2] + x_p[2,1]
   + z_p[1,2] + z_p[2,1] >= -1
 OverlapInBin_1_2_toBin_3: - p_ij[1,3] - p_ij[2,3] + x_p[1,2] + x_p[2,1]
   + z_p[1,2] + z_p[2,1] >= -1
 OverlapInBin_1_3_toBin_0: - p_ij[1,0] - p_ij[3,0] + x_p[1,3] + x_p[3,1]
   + z_p[1,3] + z_p[3,1] >= -1
 OverlapInBin_1_3_toBin_1: - p_ij[1,1] - p_ij[3,1] + x_p[1,3] + x_p[3,1]
   + z_p[1,3] + z_p[3,1] >= -1
 OverlapInBin_1_3_toBin_2: - p_ij[1,2] - p_ij[3,2] + x_p[1,3] + x_p[3,1]
   + z_p[1,3] + z_p[3,1] >= -1
 OverlapInBin_1_3_toBin_3: - p_ij[1,3] - p_ij[3,3] + x_p[1,3] + x_p[3,1]
   + z_p[1,3] + z_p[3,1] >= -1
 OverlapInBin_1_4_toBin_0: - p_ij[1,0] - p_ij[4,0] + x_p[1,4] + x_p[4,1]
   + z_p[1,4] + z_p[4,1] >= -1
 OverlapInBin_1_4_toBin_1: - p_ij[1,1] - p_ij[4,1] + x_p[1,4] + x_p[4,1]
   + z_p[1,4] + z_p[4,1] >= -1
 OverlapInBin_1_4_toBin_2: - p_ij[1,2] - p_ij[4,2] + x_p[1,4] + x_p[4,1]
   + z_p[1,4] + z_p[4,1] >= -1
 OverlapInBin_1_4_toBin_3: - p_ij[1,3] - p_ij[4,3] + x_p[1,4] + x_p[4,1]
   + z_p[1,4] + z_p[4,1] >= -1
 OverlapInBin_2_0_toBin_0: - p_ij[0,0] - p_ij[2,0] + x_p[0,2] + x_p[2,0]
   + z_p[0,2] + z_p[2,0] >= -1
 OverlapInBin_2_0_toBin_1: - p_ij[0,1] - p_ij[2,1] + x_p[0,2] + x_p[2,0]
   + z_p[0,2] + z_p[2,0] >= -1
 OverlapInBin_2_0_toBin_2: - p_ij[0,2] - p_ij[2,2] + x_p[0,2] + x_p[2,0]
   + z_p[0,2] + z_p[2,0] >= -1
 OverlapInBin_2_0_toBin_3: - p_ij[0,3] - p_ij[2,3] + x_p[0,2] + x_p[2,0]
   + z_p[0,2] + z_p[2,0] >= -1
 OverlapInBin_2_1_toBin_0: - p_ij[1,0] - p_ij[2,0] + x_p[1,2] + x_p[2,1]
   + z_p[1,2] + z_p[2,1] >= -1
 OverlapInBin_2_1_toBin_1: - p_ij[1,1] - p_ij[2,1] + x_p[1,2] + x_p[2,1]
   + z_p[1,2] + z_p[2,1] >= -1
 OverlapInBin_2_1_toBin_2: - p_ij[1,2] - p_ij[2,2] + x_p[1,2] + x_p[2,1]
   + z_p[1,2] + z_p[2,1] >= -1
 OverlapInBin_2_1_toBin_3: - p_ij[1,3] - p_ij[2,3] + x_p[1,2] + x_p[2,1]
   + z_p[1,2] + z_p[2,1] >= -1
 OverlapInBin_2_3_toBin_0: - p_ij[2,0] - p_ij[3,0] + x_p[2,3] + x_p[3,2]
   + z_p[2,3] + z_p[3,2] >= -1
 OverlapInBin_2_3_toBin_1: - p_ij[2,1] - p_ij[3,1] + x_p[2,3] + x_p[3,2]
   + z_p[2,3] + z_p[3,2] >= -1
 OverlapInBin_2_3_toBin_2: - p_ij[2,2] - p_ij[3,2] + x_p[2,3] + x_p[3,2]
   + z_p[2,3] + z_p[3,2] >= -1
 OverlapInBin_2_3_toBin_3: - p_ij[2,3] - p_ij[3,3] + x_p[2,3] + x_p[3,2]
   + z_p[2,3] + z_p[3,2] >= -1
 OverlapInBin_2_4_toBin_0: - p_ij[2,0] - p_ij[4,0] + x_p[2,4] + x_p[4,2]
   + z_p[2,4] + z_p[4,2] >= -1
 OverlapInBin_2_4_toBin_1: - p_ij[2,1] - p_ij[4,1] + x_p[2,4] + x_p[4,2]
   + z_p[2,4] + z_p[4,2] >= -1
 OverlapInBin_2_4_toBin_2: - p_ij[2,2] - p_ij[4,2] + x_p[2,4] + x_p[4,2]
   + z_p[2,4] + z_p[4,2] >= -1
 OverlapInBin_2_4_toBin_3: - p_ij[2,3] - p_ij[4,3] + x_p[2,4] + x_p[4,2]
   + z_p[2,4] + z_p[4,2] >= -1
 OverlapInBin_3_0_toBin_0: - p_ij[0,0] - p_ij[3,0] + x_p[0,3] + x_p[3,0]
   + z_p[0,3] + z_p[3,0] >= -1
 OverlapInBin_3_0_toBin_1: - p_ij[0,1] - p_ij[3,1] + x_p[0,3] + x_p[3,0]
   + z_p[0,3] + z_p[3,0] >= -1
 OverlapInBin_3_0_toBin_2: - p_ij[0,2] - p_ij[3,2] + x_p[0,3] + x_p[3,0]
   + z_p[0,3] + z_p[3,0] >= -1
 OverlapInBin_3_0_toBin_3: - p_ij[0,3] - p_ij[3,3] + x_p[0,3] + x_p[3,0]
   + z_p[0,3] + z_p[3,0] >= -1
 OverlapInBin_3_1_toBin_0: - p_ij[1,0] - p_ij[3,0] + x_p[1,3] + x_p[3,1]
   + z_p[1,3] + z_p[3,1] >= -1
 OverlapInBin_3_1_toBin_1: - p_ij[1,1] - p_ij[3,1] + x_p[1,3] + x_p[3,1]
   + z_p[1,3] + z_p[3,1] >= -1
 OverlapInBin_3_1_toBin_2: - p_ij[1,2] - p_ij[3,2] + x_p[1,3] + x_p[3,1]
   + z_p[1,3] + z_p[3,1] >= -1
 OverlapInBin_3_1_toBin_3: - p_ij[1,3] - p_ij[3,3] + x_p[1,3] + x_p[3,1]
   + z_p[1,3] + z_p[3,1] >= -1
 OverlapInBin_3_2_toBin_0: - p_ij[2,0] - p_ij[3,0] + x_p[2,3] + x_p[3,2]
   + z_p[2,3] + z_p[3,2] >= -1
 OverlapInBin_3_2_toBin_1: - p_ij[2,1] - p_ij[3,1] + x_p[2,3] + x_p[3,2]
   + z_p[2,3] + z_p[3,2] >= -1
 OverlapInBin_3_2_toBin_2: - p_ij[2,2] - p_ij[3,2] + x_p[2,3] + x_p[3,2]
   + z_p[2,3] + z_p[3,2] >= -1
 OverlapInBin_3_2_toBin_3: - p_ij[2,3] - p_ij[3,3] + x_p[2,3] + x_p[3,2]
   + z_p[2,3] + z_p[3,2] >= -1
 OverlapInBin_3_4_toBin_0: - p_ij[3,0] - p_ij[4,0] + x_p[3,4] + x_p[4,3]
   + z_p[3,4] + z_p[4,3] >= -1
 OverlapInBin_3_4_toBin_1: - p_ij[3,1] - p_ij[4,1] + x_p[3,4] + x_p[4,3]
   + z_p[3,4] + z_p[4,3] >= -1
 OverlapInBin_3_4_toBin_2: - p_ij[3,2] - p_ij[4,2] + x_p[3,4] + x_p[4,3]
   + z_p[3,4] + z_p[4,3] >= -1
 OverlapInBin_3_4_toBin_3: - p_ij[3,3] - p_ij[4,3] + x_p[3,4] + x_p[4,3]
   + z_p[3,4] + z_p[4,3] >= -1
 OverlapInBin_4_0_toBin_0: - p_ij[0,0] - p_ij[4,0] + x_p[0,4] + x_p[4,0]
   + z_p[0,4] + z_p[4,0] >= -1
 OverlapInBin_4_0_toBin_1: - p_ij[0,1] - p_ij[4,1] + x_p[0,4] + x_p[4,0]
   + z_p[0,4] + z_p[4,0] >= -1
 OverlapInBin_4_0_toBin_2: - p_ij[0,2] - p_ij[4,2] + x_p[0,4] + x_p[4,0]
   + z_p[0,4] + z_p[4,0] >= -1
 OverlapInBin_4_0_toBin_3: - p_ij[0,3] - p_ij[4,3] + x_p[0,4] + x_p[4,0]
   + z_p[0,4] + z_p[4,0] >= -1
 OverlapInBin_4_1_toBin_0: - p_ij[1,0] - p_ij[4,0] + x_p[1,4] + x_p[4,1]
   + z_p[1,4] + z_p[4,1] >= -1
 OverlapInBin_4_1_toBin_1: - p_ij[1,1] - p_ij[4,1] + x_p[1,4] + x_p[4,1]
   + z_p[1,4] + z_p[4,1] >= -1
 OverlapInBin_4_1_toBin_2: - p_ij[1,2] - p_ij[4,2] + x_p[1,4] + x_p[4,1]
   + z_p[1,4] + z_p[4,1] >= -1
 OverlapInBin_4_1_toBin_3: - p_ij[1,3] - p_ij[4,3] + x_p[1,4] + x_p[4,1]
   + z_p[1,4] + z_p[4,1] >= -1
 OverlapInBin_4_2_toBin_0: - p_ij[2,0] - p_ij[4,0] + x_p[2,4] + x_p[4,2]
   + z_p[2,4] + z_p[4,2] >= -1
 OverlapInBin_4_2_toBin_1: - p_ij[2,1] - p_ij[4,1] + x_p[2,4] + x_p[4,2]
   + z_p[2,4] + z_p[4,2] >= -1
 OverlapInBin_4_2_toBin_2: - p_ij[2,2] - p_ij[4,2] + x_p[2,4] + x_p[4,2]
   + z_p[2,4] + z_p[4,2] >= -1
 OverlapInBin_4_2_toBin_3: - p_ij[2,3] - p_ij[4,3] + x_p[2,4] + x_p[4,2]
   + z_p[2,4] + z_p[4,2] >= -1
 OverlapInBin_4_3_toBin_0: - p_ij[3,0] - p_ij[4,0] + x_p[3,4] + x_p[4,3]
   + z_p[3,4] + z_p[4,3] >= -1
 OverlapInBin_4_3_toBin_1: - p_ij[3,1] - p_ij[4,1] + x_p[3,4] + x_p[4,3]
   + z_p[3,4] + z_p[4,3] >= -1
 OverlapInBin_4_3_toBin_2: - p_ij[3,2] - p_ij[4,2] + x_p[3,4] + x_p[4,3]
   + z_p[3,4] + z_p[4,3] >= -1
 OverlapInBin_4_3_toBin_3: - p_ij[3,3] - p_ij[4,3] + x_p[3,4] + x_p[4,3]
   + z_p[3,4] + z_p[4,3] >= -1
 Overlap_0_1: x_p[0,1] + x_p[1,0] + z_p[0,1] + z_p[1,0] >= 1
 RightLeft_0_1: x_p[0,1] + x_p[1,0] <= 1
 LeftRight_0_1: x_p[0,1] + x_p[1,0] <= 1
 AboveBelow_0_1: z_p[0,1] + z_p[1,0] <= 1
 BelowAbove_0_1: z_p[0,1] + z_p[1,0] <= 1
 Overlap_0_2: x_p[0,2] + x_p[2,0] + z_p[0,2] + z_p[2,0] >= 1
 RightLeft_0_2: x_p[0,2] + x_p[2,0] <= 1
 LeftRight_0_2: x_p[0,2] + x_p[2,0] <= 1
 AboveBelow_0_2: z_p[0,2] + z_p[2,0] <= 1
 BelowAbove_0_2: z_p[0,2] + z_p[2,0] <= 1
 Overlap_0_3: x_p[0,3] + x_p[3,0] + z_p[0,3] + z_p[3,0] >= 1
 RightLeft_0_3: x_p[0,3] + x_p[3,0] <= 1
 LeftRight_0_3: x_p[0,3] + x_p[3,0] <= 1
 AboveBelow_0_3: z_p[0,3] + z_p[3,0] <= 1
 BelowAbove_0_3: z_p[0,3] + z_p[3,0] <= 1
 Overlap_0_4: x_p[0,4] + x_p[4,0] + z_p[0,4] + z_p[4,0] >= 1
 RightLeft_0_4: x_p[0,4] + x_p[4,0] <= 1
 LeftRight_0_4: x_p[0,4] + x_p[4,0] <= 1
 AboveBelow_0_4: z_p[0,4] + z_p[4,0] <= 1
 BelowAbove_0_4: z_p[0,4] + z_p[4,0] <= 1
 Overlap_1_0: x_p[0,1] + x_p[1,0] + z_p[0,1] + z_p[1,0] >= 1
 RightLeft_1_0: x_p[0,1] + x_p[1,0] <= 1
 LeftRight_1_0: x_p[0,1] + x_p[1,0] <= 1
 AboveBelow_1_0: z_p[0,1] + z_p[1,0] <= 1
 BelowAbove_1_0: z_p[0,1] + z_p[1,0] <= 1
 Overlap_1_2: x_p[1,2] + x_p[2,1] + z_p[1,2] + z_p[2,1] >= 1
 RightLeft_1_2: x_p[1,2] + x_p[2,1] <= 1
 LeftRight_1_2: x_p[1,2] + x_p[2,1] <= 1
 AboveBelow_1_2: z_p[1,2] + z_p[2,1] <= 1
 BelowAbove_1_2: z_p[1,2] + z_p[2,1] <= 1
 Overlap_1_3: x_p[1,3] + x_p[3,1] + z_p[1,3] + z_p[3,1] >= 1
 RightLeft_1_3: x_p[1,3] + x_p[3,1] <= 1
 LeftRight_1_3: x_p[1,3] + x_p[3,1] <= 1
 AboveBelow_1_3: z_p[1,3] + z_p[3,1] <= 1
 BelowAbove_1_3: z_p[1,3] + z_p[3,1] <= 1
 Overlap_1_4: x_p[1,4] + x_p[4,1] + z_p[1,4] + z_p[4,1] >= 1
 RightLeft_1_4: x_p[1,4] + x_p[4,1] <= 1
 LeftRight_1_4: x_p[1,4] + x_p[4,1] <= 1
 AboveBelow_1_4: z_p[1,4] + z_p[4,1] <= 1
 BelowAbove_1_4: z_p[1,4] + z_p[4,1] <= 1
 Overlap_2_0: x_p[0,2] + x_p[2,0] + z_p[0,2] + z_p[2,0] >= 1
 RightLeft_2_0: x_p[0,2] + x_p[2,0] <= 1
 LeftRight_2_0: x_p[0,2] + x_p[2,0] <= 1
 AboveBelow_2_0: z_p[0,2] + z_p[2,0] <= 1
 BelowAbove_2_0: z_p[0,2] + z_p[2,0] <= 1
 Overlap_2_1: x_p[1,2] + x_p[2,1] + z_p[1,2] + z_p[2,1] >= 1
 RightLeft_2_1: x_p[1,2] + x_p[2,1] <= 1
 LeftRight_2_1: x_p[1,2] + x_p[2,1] <= 1
 AboveBelow_2_1: z_p[1,2] + z_p[2,1] <= 1
 BelowAbove_2_1: z_p[1,2] + z_p[2,1] <= 1
 Overlap_2_3: x_p[2,3] + x_p[3,2] + z_p[2,3] + z_p[3,2] >= 1
 RightLeft_2_3: x_p[2,3] + x_p[3,2] <= 1
 LeftRight_2_3: x_p[2,3] + x_p[3,2] <= 1
 AboveBelow_2_3: z_p[2,3] + z_p[3,2] <= 1
 BelowAbove_2_3: z_p[2,3] + z_p[3,2] <= 1
 Overlap_2_4: x_p[2,4] + x_p[4,2] + z_p[2,4] + z_p[4,2] >= 1
 RightLeft_2_4: x_p[2,4] + x_p[4,2] <= 1
 LeftRight_2_4: x_p[2,4] + x_p[4,2] <= 1
 AboveBelow_2_4: z_p[2,4] + z_p[4,2] <= 1
 BelowAbove_2_4: z_p[2,4] + z_p[4,2] <= 1
 Overlap_3_0: x_p[0,3] + x_p[3,0] + z_p[0,3] + z_p[3,0] >= 1
 RightLeft_3_0: x_p[0,3] + x_p[3,0] <= 1
 LeftRight_3_0: x_p[0,3] + x_p[3,0] <= 1
 AboveBelow_3_0: z_p[0,3] + z_p[3,0] <= 1
 BelowAbove_3_0: z_p[0,3] + z_p[3,0] <= 1
 Overlap_3_1: x_p[1,3] + x_p[3,1] + z_p[1,3] + z_p[3,1] >= 1
 RightLeft_3_1: x_p[1,3] + x_p[3,1] <= 1
 LeftRight_3_1: x_p[1,3] + x_p[3,1] <= 1
 AboveBelow_3_1: z_p[1,3] + z_p[3,1] <= 1
 BelowAbove_3_1: z_p[1,3] + z_p[3,1] <= 1
 Overlap_3_2: x_p[2,3] + x_p[3,2] + z_p[2,3] + z_p[3,2] >= 1
 RightLeft_3_2: x_p[2,3] + x_p[3,2] <= 1
 LeftRight_3_2: x_p[2,3] + x_p[3,2] <= 1
 AboveBelow_3_2: z_p[2,3] + z_p[3,2] <= 1
 BelowAbove_3_2: z_p[2,3] + z_p[3,2] <= 1
 Overlap_3_4: x_p[3,4] + x_p[4,3] + z_p[3,4] + z_p[4,3] >= 1
 RightLeft_3_4: x_p[3,4] + x_p[4,3] <= 1
 LeftRight_3_4: x_p[3,4] + x_p[4,3] <= 1
 AboveBelow_3_4: z_p[3,4] + z_p[4,3] <= 1
 BelowAbove_3_4: z_p[3,4] + z_p[4,3] <= 1
 Overlap_4_0: x_p[0,4] + x_p[4,0] + z_p[0,4] + z_p[4,0] >= 1
 RightLeft_4_0: x_p[0,4] + x_p[4,0] <= 1
 LeftRight_4_0: x_p[0,4] + x_p[4,0] <= 1
 AboveBelow_4_0: z_p[0,4] + z_p[4,0] <= 1
 BelowAbove_4_0: z_p[0,4] + z_p[4,0] <= 1
 Overlap_4_1: x_p[1,4] + x_p[4,1] + z_p[1,4] + z_p[4,1] >= 1
 RightLeft_4_1: x_p[1,4] + x_p[4,1] <= 1
 LeftRight_4_1: x_p[1,4] + x_p[4,1] <= 1
 AboveBelow_4_1: z_p[1,4] + z_p[4,1] <= 1
 BelowAbove_4_1: z_p[1,4] + z_p[4,1] <= 1
 Overlap_4_2: x_p[2,4] + x_p[4,2] + z_p[2,4] + z_p[4,2] >= 1
 RightLeft_4_2: x_p[2,4] + x_p[4,2] <= 1
 LeftRight_4_2: x_p[2,4] + x_p[4,2] <= 1
 AboveBelow_4_2: z_p[2,4] + z_p[4,2] <= 1
 BelowAbove_4_2: z_p[2,4] + z_p[4,2] <= 1
 Overlap_4_3: x_p[3,4] + x_p[4,3] + z_p[3,4] + z_p[4,3] >= 1
 RightLeft_4_3: x_p[3,4] + x_p[4,3] <= 1
 LeftRight_4_3: x_p[3,4] + x_p[4,3] <= 1
 AboveBelow_4_3: z_p[3,4] + z_p[4,3] <= 1
 BelowAbove_4_3: z_p[3,4] + z_p[4,3] <= 1
 Constraint14_0_1: 300 x_p[0,1] - xi[0] + x_i_prime[1] <= 300
 Constraint14_0_2: 300 x_p[0,2] - xi[0] + x_i_prime[2] <= 300
 Constraint14_0_3: 300 x_p[0,3] - xi[0] + x_i_prime[3] <= 300
 Constraint14_0_4: 300 x_p[0,4] - xi[0] + x_i_prime[4] <= 300
 Constraint14_1_0: 300 x_p[1,0] - xi[1] + x_i_prime[0] <= 300
 Constraint14_1_2: 300 x_p[1,2] - xi[1] + x_i_prime[2] <= 300
 Constraint14_1_3: 300 x_p[1,3] - xi[1] + x_i_prime[3] <= 300
 Constraint14_1_4: 300 x_p[1,4] - xi[1] + x_i_prime[4] <= 300
 Constraint14_2_0: 300 x_p[2,0] - xi[2] + x_i_prime[0] <= 300
 Constraint14_2_1: 300 x_p[2,1] - xi[2] + x_i_prime[1] <= 300
 Constraint14_2_3: 300 x_p[2,3] - xi[2] + x_i_prime[3] <= 300
 Constraint14_2_4: 300 x_p[2,4] - xi[2] + x_i_prime[4] <= 300
 Constraint14_3_0: 300 x_p[3,0] - xi[3] + x_i_prime[0] <= 300
 Constraint14_3_1: 300 x_p[3,1] - xi[3] + x_i_prime[1] <= 300
 Constraint14_3_2: 300 x_p[3,2] - xi[3] + x_i_prime[2] <= 300
 Constraint14_3_4: 300 x_p[3,4] - xi[3] + x_i_prime[4] <= 300
 Constraint14_4_0: 300 x_p[4,0] - xi[4] + x_i_prime[0] <= 300
 Constraint14_4_1: 300 x_p[4,1] - xi[4] + x_i_prime[1] <= 300
 Constraint14_4_2: 300 x_p[4,2] - xi[4] + x_i_prime[2] <= 300
 Constraint14_4_3: 300 x_p[4,3] - xi[4] + x_i_prime[3] <= 300
 Constraint15_0_1: - 300 x_p[0,1] + xi[0] - x_i_prime[1] <= -1
 Constraint15_0_2: - 300 x_p[0,2] + xi[0] - x_i_prime[2] <= -1
 Constraint15_0_3: - 300 x_p[0,3] + xi[0] - x_i_prime[3] <= -1
 Constraint15_0_4: - 300 x_p[0,4] + xi[0] - x_i_prime[4] <= -1
 Constraint15_1_0: - 300 x_p[1,0] + xi[1] - x_i_prime[0] <= -1
 Constraint15_1_2: - 300 x_p[1,2] + xi[1] - x_i_prime[2] <= -1
 Constraint15_1_3: - 300 x_p[1,3] + xi[1] - x_i_prime[3] <= -1
 Constraint15_1_4: - 300 x_p[1,4] + xi[1] - x_i_prime[4] <= -1
 Constraint15_2_0: - 300 x_p[2,0] + xi[2] - x_i_prime[0] <= -1
 Constraint15_2_1: - 300 x_p[2,1] + xi[2] - x_i_prime[1] <= -1
 Constraint15_2_3: - 300 x_p[2,3] + xi[2] - x_i_prime[3] <= -1
 Constraint15_2_4: - 300 x_p[2,4] + xi[2] - x_i_prime[4] <= -1
 Constraint15_3_0: - 300 x_p[3,0] + xi[3] - x_i_prime[0] <= -1
 Constraint15_3_1: - 300 x_p[3,1] + xi[3] - x_i_prime[1] <= -1
 Constraint15_3_2: - 300 x_p[3,2] + xi[3] - x_i_prime[2] <= -1
 Constraint15_3_4: - 300 x_p[3,4] + xi[3] - x_i_prime[4] <= -1
 Constraint15_4_0: - 300 x_p[4,0] + xi[4] - x_i_prime[0] <= -1
 Constraint15_4_1: - 300 x_p[4,1] + xi[4] - x_i_prime[1] <= -1
 Constraint15_4_2: - 300 x_p[4,2] + xi[4] - x_i_prime[2] <= -1
 Constraint15_4_3: - 300 x_p[4,3] + xi[4] - x_i_prime[3] <= -1
 Constraint16_0_1: 155 z_p[0,1] - zi[0] + z_i_prime[1] <= 155
 Constraint16_0_2: 155 z_p[0,2] - zi[0] + z_i_prime[2] <= 155
 Constraint16_0_3: 155 z_p[0,3] - zi[0] + z_i_prime[3] <= 155
 Constraint16_0_4: 155 z_p[0,4] - zi[0] + z_i_prime[4] <= 155
 Constraint16_1_0: 155 z_p[1,0] - zi[1] + z_i_prime[0] <= 155
 Constraint16_1_2: 155 z_p[1,2] - zi[1] + z_i_prime[2] <= 155
 Constraint16_1_3: 155 z_p[1,3] - zi[1] + z_i_prime[3] <= 155
 Constraint16_1_4: 155 z_p[1,4] - zi[1] + z_i_prime[4] <= 155
 Constraint16_2_0: 155 z_p[2,0] - zi[2] + z_i_prime[0] <= 155
 Constraint16_2_1: 155 z_p[2,1] - zi[2] + z_i_prime[1] <= 155
 Constraint16_2_3: 155 z_p[2,3] - zi[2] + z_i_prime[3] <= 155
 Constraint16_2_4: 155 z_p[2,4] - zi[2] + z_i_prime[4] <= 155
 Constraint16_3_0: 155 z_p[3,0] - zi[3] + z_i_prime[0] <= 155
 Constraint16_3_1: 155 z_p[3,1] - zi[3] + z_i_prime[1] <= 155
 Constraint16_3_2: 155 z_p[3,2] - zi[3] + z_i_prime[2] <= 155
 Constraint16_3_4: 155 z_p[3,4] - zi[3] + z_i_prime[4] <= 155
 Constraint16_4_0: 155 z_p[4,0] - zi[4] + z_i_prime[0] <= 155
 Constraint16_4_1: 155 z_p[4,1] - zi[4] + z_i_prime[1] <= 155
 Constraint16_4_2: 155 z_p[4,2] - zi[4] + z_i_prime[2] <= 155
 Constraint16_4_3: 155 z_p[4,3] - zi[4] + z_i_prime[3] <= 155
 OrientationLength_0: r[0,0,0] <= 1
 OrientationHeight_0: r[0,1,1] <= 1
 OrientationLength_1: r[1,0,0] <= 0
 OrientationHeight_1: r[1,1,1] <= 0
 OrientationLength_2: r[2,0,0] <= 1
 OrientationHeight_2: r[2,1,1] <= 1
 OrientationLength_3: r[3,0,0] <= 1
 OrientationHeight_3: r[3,1,1] <= 1
 OrientationLength_4: r[4,0,0] <= 1
 OrientationHeight_4: r[4,1,1] <= 1
 Flagging_0_0: - p_ij[0,0] + u_j[0] >= 0
 Flagging_0_1: - p_ij[0,1] + u_j[1] >= 0
 Flagging_0_2: - p_ij[0,2] + u_j[2] >= 0
 Flagging_0_3: - p_ij[0,3] + u_j[3] >= 0
 Flagging_1_0: - p_ij[1,0] + u_j[0] >= 0
 Flagging_1_1: - p_ij[1,1] + u_j[1] >= 0
 Flagging_1_2: - p_ij[1,2] + u_j[2] >= 0
 Flagging_1_3: - p_ij[1,3] + u_j[3] >= 0
 Flagging_2_0: - p_ij[2,0] + u_j[0] >= 0
 Flagging_2_1: - p_ij[2,1] + u_j[1] >= 0
 Flagging_2_2: - p_ij[2,2] + u_j[2] >= 0
 Flagging_2_3: - p_ij[2,3] + u_j[3] >= 0
 Flagging_3_0: - p_ij[3,0] + u_j[0] >= 0
 Flagging_3_1: - p_ij[3,1] + u_j[1] >= 0
 Flagging_3_2: - p_ij[3,2] + u_j[2] >= 0
 Flagging_3_3: - p_ij[3,3] + u_j[3] >= 0
 Flagging_4_0: - p_ij[4,0] + u_j[0] >= 0
 Flagging_4_1: - p_ij[4,1] + u_j[1] >= 0
 Flagging_4_2: - p_ij[4,2] + u_j[2] >= 0
 Flagging_4_3: - p_ij[4,3] + u_j[3] >= 0
 Stability_0_3: 2 g[0] + beta1[0,1] + beta1[0,2] + beta1[0,3] + beta1[0,4]
   + beta2[0,1] + beta2[0,2] + beta2[0,3] + beta2[0,4] + gamma[0] >= 2
 Stability_1_3: 2 g[1] + beta1[1,0] + beta1[1,2] + beta1[1,3] + beta1[1,4]
   + beta2[1,0] + beta2[1,2] + beta2[1,3] + beta2[1,4] + gamma[1] >= 2
 Stability_2_3: 2 g[2] + beta1[2,0] + beta1[2,1] + beta1[2,3] + beta1[2,4]
   + beta2[2,0] + beta2[2,1] + beta2[2,3] + beta2[2,4] + gamma[2] >= 2
 Stability_3_3: 2 g[3] + beta1[3,0] + beta1[3,1] + beta1[3,2] + beta1[3,4]
   + beta2[3,0] + beta2[3,1] + beta2[3,2] + beta2[3,4] + gamma[3] >= 2
 Stability_4_3: 2 g[4] + beta1[4,0] + beta1[4,1] + beta1[4,2] + beta1[4,3]
   + beta2[4,0] + beta2[4,1] + beta2[4,2] + beta2[4,3] + gamma[4] >= 2
Bounds
Binaries
 p_ij[0,0] p_ij[0,1] p_ij[0,2] p_ij[0,3] p_ij[1,0] p_ij[1,1] p_ij[1,2]
 p_ij[1,3] p_ij[2,0] p_ij[2,1] p_ij[2,2] p_ij[2,3] p_ij[3,0] p_ij[3,1]
 p_ij[3,2] p_ij[3,3] p_ij[4,0] p_ij[4,1] p_ij[4,2] p_ij[4,3] u_j[0] u_j[1]
 u_j[2] u_j[3] x_p[0,0] x_p[0,1] x_p[0,2] x_p[0,3] x_p[0,4] x_p[1,0]
 x_p[1,1] x_p[1,2] x_p[1,3] x_p[1,4] x_p[2,0] x_p[2,1] x_p[2,2] x_p[2,3]
 x_p[2,4] x_p[3,0] x_p[3,1] x_p[3,2] x_p[3,3] x_p[3,4] x_p[4,0] x_p[4,1]
 x_p[4,2] x_p[4,3] x_p[4,4] z_p[0,0] z_p[0,1] z_p[0,2] z_p[0,3] z_p[0,4]
 z_p[1,0] z_p[1,1] z_p[1,2] z_p[1,3] z_p[1,4] z_p[2,0] z_p[2,1] z_p[2,2]
 z_p[2,3] z_p[2,4] z_p[3,0] z_p[3,1] z_p[3,2] z_p[3,3] z_p[3,4] z_p[4,0]
 z_p[4,1] z_p[4,2] z_p[4,3] z_p[4,4] r[0,0,0] r[0,0,1] r[0,1,0] r[0,1,1]
 r[1,0,0] r[1,0,1] r[1,1,0] r[1,1,1] r[2,0,0] r[2,0,1] r[2,1,0] r[2,1,1]
 r[3,0,0] r[3,0,1] r[3,1,0] r[3,1,1] r[4,0,0] r[4,0,1] r[4,1,0] r[4,1,1]
 rho[0] rho[1] rho[2] rho[3] rho[4] g[0] g[1] g[2] g[3] g[4] beta1[0,0]
 beta1[0,1] beta1[0,2] beta1[0,3] beta1[0,4] beta1[1,0] beta1[1,1]
 beta1[1,2] beta1[1,3] beta1[1,4] beta1[2,0] beta1[2,1] beta1[2,2]
 beta1[2,3] beta1[2,4] beta1[3,0] beta1[3,1] beta1[3,2] beta1[3,3]
 beta1[3,4] beta1[4,0] beta1[4,1] beta1[4,2] beta1[4,3] beta1[4,4]
 beta2[0,0] beta2[0,1] beta2[0,2] beta2[0,3] beta2[0,4] beta2[1,0]
 beta2[1,1] beta2[1,2] beta2[1,3] beta2[1,4] beta2[2,0] beta2[2,1]
 beta2[2,2] beta2[2,3] beta2[2,4] beta2[3,0] beta2[3,1] beta2[3,2]
 beta2[3,3] beta2[3,4] beta2[4,0] beta2[4,1] beta2[4,2] beta2[4,3]
 beta2[4,4] gamma[0] gamma[1] gamma[2] gamma[3] gamma[4]
End
