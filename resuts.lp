\ Model 2DBPP
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  200 u_j[0] + 200 u_j[1] + 150 u_j[2] + 150 u_j[3] + 0 x_p[0,0]
   + 0 x_p[1,1] + 0 x_p[2,2] + 0 x_p[3,3] + 0 x_p[4,4] + 0 x_p[5,5]
   + 0 x_p[6,6] + 0 x_p[7,7] + 0 z_p[0,0] + 0 z_p[1,1] + 0 z_p[2,2]
   + 0 z_p[3,3] + 0 z_p[4,4] + 0 z_p[5,5] + 0 z_p[6,6] + 0 z_p[7,7]
   + 0 rho[0] + 0 rho[1] + 0 rho[2] + 0 rho[3] + 0 rho[4] + 0 rho[5]
   + 0 rho[6] + 0 rho[7] + 0 beta1[0,0] + 0 beta1[1,1] + 0 beta1[2,2]
   + 0 beta1[3,3] + 0 beta1[4,4] + 0 beta1[5,5] + 0 beta1[6,6]
   + 0 beta1[7,7] + 0 beta2[0,0] + 0 beta2[1,1] + 0 beta2[2,2]
   + 0 beta2[3,3] + 0 beta2[4,4] + 0 beta2[5,5] + 0 beta2[6,6]
   + 0 beta2[7,7] + 0 v_i_k[0,0] + 0 v_i_k[0,1] + 0 v_i_k[0,2]
   + 0 v_i_k[0,3] + 0 v_i_k[0,4] + 0 v_i_k[0,5] + 0 v_i_k[0,6]
   + 0 v_i_k[0,7] + 0 v_i_k[1,0] + 0 v_i_k[1,1] + 0 v_i_k[1,2]
   + 0 v_i_k[1,3] + 0 v_i_k[1,4] + 0 v_i_k[1,5] + 0 v_i_k[1,6]
   + 0 v_i_k[1,7] + 0 v_i_k[2,0] + 0 v_i_k[2,1] + 0 v_i_k[2,2]
   + 0 v_i_k[2,3] + 0 v_i_k[2,4] + 0 v_i_k[2,5] + 0 v_i_k[2,6]
   + 0 v_i_k[2,7] + 0 v_i_k[3,0] + 0 v_i_k[3,1] + 0 v_i_k[3,2]
   + 0 v_i_k[3,3] + 0 v_i_k[3,4] + 0 v_i_k[3,5] + 0 v_i_k[3,6]
   + 0 v_i_k[3,7] + 0 v_i_k[4,0] + 0 v_i_k[4,1] + 0 v_i_k[4,2]
   + 0 v_i_k[4,3] + 0 v_i_k[4,4] + 0 v_i_k[4,5] + 0 v_i_k[4,6]
   + 0 v_i_k[4,7] + 0 v_i_k[5,0] + 0 v_i_k[5,1] + 0 v_i_k[5,2]
   + 0 v_i_k[5,3] + 0 v_i_k[5,4] + 0 v_i_k[5,5] + 0 v_i_k[5,6]
   + 0 v_i_k[5,7] + 0 v_i_k[6,0] + 0 v_i_k[6,1] + 0 v_i_k[6,2]
   + 0 v_i_k[6,3] + 0 v_i_k[6,4] + 0 v_i_k[6,5] + 0 v_i_k[6,6]
   + 0 v_i_k[6,7] + 0 v_i_k[7,0] + 0 v_i_k[7,1] + 0 v_i_k[7,2]
   + 0 v_i_k[7,3] + 0 v_i_k[7,4] + 0 v_i_k[7,5] + 0 v_i_k[7,6]
   + 0 v_i_k[7,7] + 0 h_i_k[0,0] + 0 h_i_k[0,1] + 0 h_i_k[0,2]
   + 0 h_i_k[0,3] + 0 h_i_k[0,4] + 0 h_i_k[0,5] + 0 h_i_k[0,6]
   + 0 h_i_k[0,7] + 0 h_i_k[1,0] + 0 h_i_k[1,1] + 0 h_i_k[1,2]
   + 0 h_i_k[1,3] + 0 h_i_k[1,4] + 0 h_i_k[1,5] + 0 h_i_k[1,6]
   + 0 h_i_k[1,7] + 0 h_i_k[2,0] + 0 h_i_k[2,1] + 0 h_i_k[2,2]
   + 0 h_i_k[2,3] + 0 h_i_k[2,4] + 0 h_i_k[2,5] + 0 h_i_k[2,6]
   + 0 h_i_k[2,7] + 0 h_i_k[3,0] + 0 h_i_k[3,1] + 0 h_i_k[3,2]
   + 0 h_i_k[3,3] + 0 h_i_k[3,4] + 0 h_i_k[3,5] + 0 h_i_k[3,6]
   + 0 h_i_k[3,7] + 0 h_i_k[4,0] + 0 h_i_k[4,1] + 0 h_i_k[4,2]
   + 0 h_i_k[4,3] + 0 h_i_k[4,4] + 0 h_i_k[4,5] + 0 h_i_k[4,6]
   + 0 h_i_k[4,7] + 0 h_i_k[5,0] + 0 h_i_k[5,1] + 0 h_i_k[5,2]
   + 0 h_i_k[5,3] + 0 h_i_k[5,4] + 0 h_i_k[5,5] + 0 h_i_k[5,6]
   + 0 h_i_k[5,7] + 0 h_i_k[6,0] + 0 h_i_k[6,1] + 0 h_i_k[6,2]
   + 0 h_i_k[6,3] + 0 h_i_k[6,4] + 0 h_i_k[6,5] + 0 h_i_k[6,6]
   + 0 h_i_k[6,7] + 0 h_i_k[7,0] + 0 h_i_k[7,1] + 0 h_i_k[7,2]
   + 0 h_i_k[7,3] + 0 h_i_k[7,4] + 0 h_i_k[7,5] + 0 h_i_k[7,6]
   + 0 h_i_k[7,7] + 0 n1_i_k[0,0] + 0 n1_i_k[0,1] + 0 n1_i_k[0,2]
   + 0 n1_i_k[0,3] + 0 n1_i_k[0,4] + 0 n1_i_k[0,5] + 0 n1_i_k[0,6]
   + 0 n1_i_k[0,7] + 0 n1_i_k[1,0] + 0 n1_i_k[1,1] + 0 n1_i_k[1,2]
   + 0 n1_i_k[1,3] + 0 n1_i_k[1,4] + 0 n1_i_k[1,5] + 0 n1_i_k[1,6]
   + 0 n1_i_k[1,7] + 0 n1_i_k[2,0] + 0 n1_i_k[2,1] + 0 n1_i_k[2,2]
   + 0 n1_i_k[2,3] + 0 n1_i_k[2,4] + 0 n1_i_k[2,5] + 0 n1_i_k[2,6]
   + 0 n1_i_k[2,7] + 0 n1_i_k[3,0] + 0 n1_i_k[3,1] + 0 n1_i_k[3,2]
   + 0 n1_i_k[3,3] + 0 n1_i_k[3,4] + 0 n1_i_k[3,5] + 0 n1_i_k[3,6]
   + 0 n1_i_k[3,7] + 0 n1_i_k[4,0] + 0 n1_i_k[4,1] + 0 n1_i_k[4,2]
   + 0 n1_i_k[4,3] + 0 n1_i_k[4,4] + 0 n1_i_k[4,5] + 0 n1_i_k[4,6]
   + 0 n1_i_k[4,7] + 0 n1_i_k[5,0] + 0 n1_i_k[5,1] + 0 n1_i_k[5,2]
   + 0 n1_i_k[5,3] + 0 n1_i_k[5,4] + 0 n1_i_k[5,5] + 0 n1_i_k[5,6]
   + 0 n1_i_k[5,7] + 0 n1_i_k[6,0] + 0 n1_i_k[6,1] + 0 n1_i_k[6,2]
   + 0 n1_i_k[6,3] + 0 n1_i_k[6,4] + 0 n1_i_k[6,5] + 0 n1_i_k[6,6]
   + 0 n1_i_k[6,7] + 0 n1_i_k[7,0] + 0 n1_i_k[7,1] + 0 n1_i_k[7,2]
   + 0 n1_i_k[7,3] + 0 n1_i_k[7,4] + 0 n1_i_k[7,5] + 0 n1_i_k[7,6]
   + 0 n1_i_k[7,7] + 0 n2_i_k[0,0] + 0 n2_i_k[0,1] + 0 n2_i_k[0,2]
   + 0 n2_i_k[0,3] + 0 n2_i_k[0,4] + 0 n2_i_k[0,5] + 0 n2_i_k[0,6]
   + 0 n2_i_k[0,7] + 0 n2_i_k[1,0] + 0 n2_i_k[1,1] + 0 n2_i_k[1,2]
   + 0 n2_i_k[1,3] + 0 n2_i_k[1,4] + 0 n2_i_k[1,5] + 0 n2_i_k[1,6]
   + 0 n2_i_k[1,7] + 0 n2_i_k[2,0] + 0 n2_i_k[2,1] + 0 n2_i_k[2,2]
   + 0 n2_i_k[2,3] + 0 n2_i_k[2,4] + 0 n2_i_k[2,5] + 0 n2_i_k[2,6]
   + 0 n2_i_k[2,7] + 0 n2_i_k[3,0] + 0 n2_i_k[3,1] + 0 n2_i_k[3,2]
   + 0 n2_i_k[3,3] + 0 n2_i_k[3,4] + 0 n2_i_k[3,5] + 0 n2_i_k[3,6]
   + 0 n2_i_k[3,7] + 0 n2_i_k[4,0] + 0 n2_i_k[4,1] + 0 n2_i_k[4,2]
   + 0 n2_i_k[4,3] + 0 n2_i_k[4,4] + 0 n2_i_k[4,5] + 0 n2_i_k[4,6]
   + 0 n2_i_k[4,7] + 0 n2_i_k[5,0] + 0 n2_i_k[5,1] + 0 n2_i_k[5,2]
   + 0 n2_i_k[5,3] + 0 n2_i_k[5,4] + 0 n2_i_k[5,5] + 0 n2_i_k[5,6]
   + 0 n2_i_k[5,7] + 0 n2_i_k[6,0] + 0 n2_i_k[6,1] + 0 n2_i_k[6,2]
   + 0 n2_i_k[6,3] + 0 n2_i_k[6,4] + 0 n2_i_k[6,5] + 0 n2_i_k[6,6]
   + 0 n2_i_k[6,7] + 0 n2_i_k[7,0] + 0 n2_i_k[7,1] + 0 n2_i_k[7,2]
   + 0 n2_i_k[7,3] + 0 n2_i_k[7,4] + 0 n2_i_k[7,5] + 0 n2_i_k[7,6]
   + 0 n2_i_k[7,7]
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
 Beta1_0_5: xi[0] - x_i_prime[5] + 100000 beta1[0,5] <= 100000
 Beta1_Z_Upper_0_5: zi[0] - z_i_prime[5] + 100000 beta1[0,5] <= 100000
 Beta1_Z_Lower_0_5: zi[0] - z_i_prime[5] - 100000 beta1[0,5] >= -100000
 Beta2_0_5: - xi[5] + x_i_prime[0] - 100000 beta2[0,5] >= -100000
 Beta2_Z_Upper_0_5: - zi[5] + z_i_prime[0] + 100000 beta2[0,5] <= 100000
 Beta2_Z_Lower_0_5: - zi[5] + z_i_prime[0] - 100000 beta2[0,5] >= -100000
 Beta1_0_6: xi[0] - x_i_prime[6] + 100000 beta1[0,6] <= 100000
 Beta1_Z_Upper_0_6: zi[0] - z_i_prime[6] + 100000 beta1[0,6] <= 100000
 Beta1_Z_Lower_0_6: zi[0] - z_i_prime[6] - 100000 beta1[0,6] >= -100000
 Beta2_0_6: - xi[6] + x_i_prime[0] - 100000 beta2[0,6] >= -100000
 Beta2_Z_Upper_0_6: - zi[6] + z_i_prime[0] + 100000 beta2[0,6] <= 100000
 Beta2_Z_Lower_0_6: - zi[6] + z_i_prime[0] - 100000 beta2[0,6] >= -100000
 Beta1_0_7: xi[0] - x_i_prime[7] + 100000 beta1[0,7] <= 100000
 Beta1_Z_Upper_0_7: zi[0] - z_i_prime[7] + 100000 beta1[0,7] <= 100000
 Beta1_Z_Lower_0_7: zi[0] - z_i_prime[7] - 100000 beta1[0,7] >= -100000
 Beta2_0_7: - xi[7] + x_i_prime[0] - 100000 beta2[0,7] >= -100000
 Beta2_Z_Upper_0_7: - zi[7] + z_i_prime[0] + 100000 beta2[0,7] <= 100000
 Beta2_Z_Lower_0_7: - zi[7] + z_i_prime[0] - 100000 beta2[0,7] >= -100000
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
 Beta1_1_5: xi[1] - x_i_prime[5] + 100000 beta1[1,5] <= 100000
 Beta1_Z_Upper_1_5: zi[1] - z_i_prime[5] + 100000 beta1[1,5] <= 100000
 Beta1_Z_Lower_1_5: zi[1] - z_i_prime[5] - 100000 beta1[1,5] >= -100000
 Beta2_1_5: - xi[5] + x_i_prime[1] - 100000 beta2[1,5] >= -100000
 Beta2_Z_Upper_1_5: - zi[5] + z_i_prime[1] + 100000 beta2[1,5] <= 100000
 Beta2_Z_Lower_1_5: - zi[5] + z_i_prime[1] - 100000 beta2[1,5] >= -100000
 Beta1_1_6: xi[1] - x_i_prime[6] + 100000 beta1[1,6] <= 100000
 Beta1_Z_Upper_1_6: zi[1] - z_i_prime[6] + 100000 beta1[1,6] <= 100000
 Beta1_Z_Lower_1_6: zi[1] - z_i_prime[6] - 100000 beta1[1,6] >= -100000
 Beta2_1_6: - xi[6] + x_i_prime[1] - 100000 beta2[1,6] >= -100000
 Beta2_Z_Upper_1_6: - zi[6] + z_i_prime[1] + 100000 beta2[1,6] <= 100000
 Beta2_Z_Lower_1_6: - zi[6] + z_i_prime[1] - 100000 beta2[1,6] >= -100000
 Beta1_1_7: xi[1] - x_i_prime[7] + 100000 beta1[1,7] <= 100000
 Beta1_Z_Upper_1_7: zi[1] - z_i_prime[7] + 100000 beta1[1,7] <= 100000
 Beta1_Z_Lower_1_7: zi[1] - z_i_prime[7] - 100000 beta1[1,7] >= -100000
 Beta2_1_7: - xi[7] + x_i_prime[1] - 100000 beta2[1,7] >= -100000
 Beta2_Z_Upper_1_7: - zi[7] + z_i_prime[1] + 100000 beta2[1,7] <= 100000
 Beta2_Z_Lower_1_7: - zi[7] + z_i_prime[1] - 100000 beta2[1,7] >= -100000
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
 Beta1_2_5: xi[2] - x_i_prime[5] + 100000 beta1[2,5] <= 100000
 Beta1_Z_Upper_2_5: zi[2] - z_i_prime[5] + 100000 beta1[2,5] <= 100000
 Beta1_Z_Lower_2_5: zi[2] - z_i_prime[5] - 100000 beta1[2,5] >= -100000
 Beta2_2_5: - xi[5] + x_i_prime[2] - 100000 beta2[2,5] >= -100000
 Beta2_Z_Upper_2_5: - zi[5] + z_i_prime[2] + 100000 beta2[2,5] <= 100000
 Beta2_Z_Lower_2_5: - zi[5] + z_i_prime[2] - 100000 beta2[2,5] >= -100000
 Beta1_2_6: xi[2] - x_i_prime[6] + 100000 beta1[2,6] <= 100000
 Beta1_Z_Upper_2_6: zi[2] - z_i_prime[6] + 100000 beta1[2,6] <= 100000
 Beta1_Z_Lower_2_6: zi[2] - z_i_prime[6] - 100000 beta1[2,6] >= -100000
 Beta2_2_6: - xi[6] + x_i_prime[2] - 100000 beta2[2,6] >= -100000
 Beta2_Z_Upper_2_6: - zi[6] + z_i_prime[2] + 100000 beta2[2,6] <= 100000
 Beta2_Z_Lower_2_6: - zi[6] + z_i_prime[2] - 100000 beta2[2,6] >= -100000
 Beta1_2_7: xi[2] - x_i_prime[7] + 100000 beta1[2,7] <= 100000
 Beta1_Z_Upper_2_7: zi[2] - z_i_prime[7] + 100000 beta1[2,7] <= 100000
 Beta1_Z_Lower_2_7: zi[2] - z_i_prime[7] - 100000 beta1[2,7] >= -100000
 Beta2_2_7: - xi[7] + x_i_prime[2] - 100000 beta2[2,7] >= -100000
 Beta2_Z_Upper_2_7: - zi[7] + z_i_prime[2] + 100000 beta2[2,7] <= 100000
 Beta2_Z_Lower_2_7: - zi[7] + z_i_prime[2] - 100000 beta2[2,7] >= -100000
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
 Beta1_3_5: xi[3] - x_i_prime[5] + 100000 beta1[3,5] <= 100000
 Beta1_Z_Upper_3_5: zi[3] - z_i_prime[5] + 100000 beta1[3,5] <= 100000
 Beta1_Z_Lower_3_5: zi[3] - z_i_prime[5] - 100000 beta1[3,5] >= -100000
 Beta2_3_5: - xi[5] + x_i_prime[3] - 100000 beta2[3,5] >= -100000
 Beta2_Z_Upper_3_5: - zi[5] + z_i_prime[3] + 100000 beta2[3,5] <= 100000
 Beta2_Z_Lower_3_5: - zi[5] + z_i_prime[3] - 100000 beta2[3,5] >= -100000
 Beta1_3_6: xi[3] - x_i_prime[6] + 100000 beta1[3,6] <= 100000
 Beta1_Z_Upper_3_6: zi[3] - z_i_prime[6] + 100000 beta1[3,6] <= 100000
 Beta1_Z_Lower_3_6: zi[3] - z_i_prime[6] - 100000 beta1[3,6] >= -100000
 Beta2_3_6: - xi[6] + x_i_prime[3] - 100000 beta2[3,6] >= -100000
 Beta2_Z_Upper_3_6: - zi[6] + z_i_prime[3] + 100000 beta2[3,6] <= 100000
 Beta2_Z_Lower_3_6: - zi[6] + z_i_prime[3] - 100000 beta2[3,6] >= -100000
 Beta1_3_7: xi[3] - x_i_prime[7] + 100000 beta1[3,7] <= 100000
 Beta1_Z_Upper_3_7: zi[3] - z_i_prime[7] + 100000 beta1[3,7] <= 100000
 Beta1_Z_Lower_3_7: zi[3] - z_i_prime[7] - 100000 beta1[3,7] >= -100000
 Beta2_3_7: - xi[7] + x_i_prime[3] - 100000 beta2[3,7] >= -100000
 Beta2_Z_Upper_3_7: - zi[7] + z_i_prime[3] + 100000 beta2[3,7] <= 100000
 Beta2_Z_Lower_3_7: - zi[7] + z_i_prime[3] - 100000 beta2[3,7] >= -100000
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
 Beta1_4_5: xi[4] - x_i_prime[5] + 100000 beta1[4,5] <= 100000
 Beta1_Z_Upper_4_5: zi[4] - z_i_prime[5] + 100000 beta1[4,5] <= 100000
 Beta1_Z_Lower_4_5: zi[4] - z_i_prime[5] - 100000 beta1[4,5] >= -100000
 Beta2_4_5: - xi[5] + x_i_prime[4] - 100000 beta2[4,5] >= -100000
 Beta2_Z_Upper_4_5: - zi[5] + z_i_prime[4] + 100000 beta2[4,5] <= 100000
 Beta2_Z_Lower_4_5: - zi[5] + z_i_prime[4] - 100000 beta2[4,5] >= -100000
 Beta1_4_6: xi[4] - x_i_prime[6] + 100000 beta1[4,6] <= 100000
 Beta1_Z_Upper_4_6: zi[4] - z_i_prime[6] + 100000 beta1[4,6] <= 100000
 Beta1_Z_Lower_4_6: zi[4] - z_i_prime[6] - 100000 beta1[4,6] >= -100000
 Beta2_4_6: - xi[6] + x_i_prime[4] - 100000 beta2[4,6] >= -100000
 Beta2_Z_Upper_4_6: - zi[6] + z_i_prime[4] + 100000 beta2[4,6] <= 100000
 Beta2_Z_Lower_4_6: - zi[6] + z_i_prime[4] - 100000 beta2[4,6] >= -100000
 Beta1_4_7: xi[4] - x_i_prime[7] + 100000 beta1[4,7] <= 100000
 Beta1_Z_Upper_4_7: zi[4] - z_i_prime[7] + 100000 beta1[4,7] <= 100000
 Beta1_Z_Lower_4_7: zi[4] - z_i_prime[7] - 100000 beta1[4,7] >= -100000
 Beta2_4_7: - xi[7] + x_i_prime[4] - 100000 beta2[4,7] >= -100000
 Beta2_Z_Upper_4_7: - zi[7] + z_i_prime[4] + 100000 beta2[4,7] <= 100000
 Beta2_Z_Lower_4_7: - zi[7] + z_i_prime[4] - 100000 beta2[4,7] >= -100000
 GroundSupport_UpperBound_5: zi[5] + 100000 g[5] <= 100000
 CutSupport_UpperBound_5_2: 1.261904761904762 xi[5] + zi[5]
   + 100000 gamma[5] <= 100053
 CutSupport_LowerBound_5_2: 1.261904761904762 xi[5] + zi[5]
   - 100000 gamma[5] >= -99947
 CutSupport_UpperBound_5_3: 1.261904761904762 xi[5] + zi[5]
   + 100000 gamma[5] <= 100053
 CutSupport_LowerBound_5_3: 1.261904761904762 xi[5] + zi[5]
   - 100000 gamma[5] >= -99947
 Beta1_5_0: xi[5] - x_i_prime[0] + 100000 beta1[5,0] <= 100000
 Beta1_Z_Upper_5_0: zi[5] - z_i_prime[0] + 100000 beta1[5,0] <= 100000
 Beta1_Z_Lower_5_0: zi[5] - z_i_prime[0] - 100000 beta1[5,0] >= -100000
 Beta2_5_0: - xi[0] + x_i_prime[5] - 100000 beta2[5,0] >= -100000
 Beta2_Z_Upper_5_0: - zi[0] + z_i_prime[5] + 100000 beta2[5,0] <= 100000
 Beta2_Z_Lower_5_0: - zi[0] + z_i_prime[5] - 100000 beta2[5,0] >= -100000
 Beta1_5_1: xi[5] - x_i_prime[1] + 100000 beta1[5,1] <= 100000
 Beta1_Z_Upper_5_1: zi[5] - z_i_prime[1] + 100000 beta1[5,1] <= 100000
 Beta1_Z_Lower_5_1: zi[5] - z_i_prime[1] - 100000 beta1[5,1] >= -100000
 Beta2_5_1: - xi[1] + x_i_prime[5] - 100000 beta2[5,1] >= -100000
 Beta2_Z_Upper_5_1: - zi[1] + z_i_prime[5] + 100000 beta2[5,1] <= 100000
 Beta2_Z_Lower_5_1: - zi[1] + z_i_prime[5] - 100000 beta2[5,1] >= -100000
 Beta1_5_2: xi[5] - x_i_prime[2] + 100000 beta1[5,2] <= 100000
 Beta1_Z_Upper_5_2: zi[5] - z_i_prime[2] + 100000 beta1[5,2] <= 100000
 Beta1_Z_Lower_5_2: zi[5] - z_i_prime[2] - 100000 beta1[5,2] >= -100000
 Beta2_5_2: - xi[2] + x_i_prime[5] - 100000 beta2[5,2] >= -100000
 Beta2_Z_Upper_5_2: - zi[2] + z_i_prime[5] + 100000 beta2[5,2] <= 100000
 Beta2_Z_Lower_5_2: - zi[2] + z_i_prime[5] - 100000 beta2[5,2] >= -100000
 Beta1_5_3: xi[5] - x_i_prime[3] + 100000 beta1[5,3] <= 100000
 Beta1_Z_Upper_5_3: zi[5] - z_i_prime[3] + 100000 beta1[5,3] <= 100000
 Beta1_Z_Lower_5_3: zi[5] - z_i_prime[3] - 100000 beta1[5,3] >= -100000
 Beta2_5_3: - xi[3] + x_i_prime[5] - 100000 beta2[5,3] >= -100000
 Beta2_Z_Upper_5_3: - zi[3] + z_i_prime[5] + 100000 beta2[5,3] <= 100000
 Beta2_Z_Lower_5_3: - zi[3] + z_i_prime[5] - 100000 beta2[5,3] >= -100000
 Beta1_5_4: xi[5] - x_i_prime[4] + 100000 beta1[5,4] <= 100000
 Beta1_Z_Upper_5_4: zi[5] - z_i_prime[4] + 100000 beta1[5,4] <= 100000
 Beta1_Z_Lower_5_4: zi[5] - z_i_prime[4] - 100000 beta1[5,4] >= -100000
 Beta2_5_4: - xi[4] + x_i_prime[5] - 100000 beta2[5,4] >= -100000
 Beta2_Z_Upper_5_4: - zi[4] + z_i_prime[5] + 100000 beta2[5,4] <= 100000
 Beta2_Z_Lower_5_4: - zi[4] + z_i_prime[5] - 100000 beta2[5,4] >= -100000
 Beta1_5_6: xi[5] - x_i_prime[6] + 100000 beta1[5,6] <= 100000
 Beta1_Z_Upper_5_6: zi[5] - z_i_prime[6] + 100000 beta1[5,6] <= 100000
 Beta1_Z_Lower_5_6: zi[5] - z_i_prime[6] - 100000 beta1[5,6] >= -100000
 Beta2_5_6: - xi[6] + x_i_prime[5] - 100000 beta2[5,6] >= -100000
 Beta2_Z_Upper_5_6: - zi[6] + z_i_prime[5] + 100000 beta2[5,6] <= 100000
 Beta2_Z_Lower_5_6: - zi[6] + z_i_prime[5] - 100000 beta2[5,6] >= -100000
 Beta1_5_7: xi[5] - x_i_prime[7] + 100000 beta1[5,7] <= 100000
 Beta1_Z_Upper_5_7: zi[5] - z_i_prime[7] + 100000 beta1[5,7] <= 100000
 Beta1_Z_Lower_5_7: zi[5] - z_i_prime[7] - 100000 beta1[5,7] >= -100000
 Beta2_5_7: - xi[7] + x_i_prime[5] - 100000 beta2[5,7] >= -100000
 Beta2_Z_Upper_5_7: - zi[7] + z_i_prime[5] + 100000 beta2[5,7] <= 100000
 Beta2_Z_Lower_5_7: - zi[7] + z_i_prime[5] - 100000 beta2[5,7] >= -100000
 GroundSupport_UpperBound_6: zi[6] + 100000 g[6] <= 100000
 CutSupport_UpperBound_6_2: 1.261904761904762 xi[6] + zi[6]
   + 100000 gamma[6] <= 100053
 CutSupport_LowerBound_6_2: 1.261904761904762 xi[6] + zi[6]
   - 100000 gamma[6] >= -99947
 CutSupport_UpperBound_6_3: 1.261904761904762 xi[6] + zi[6]
   + 100000 gamma[6] <= 100053
 CutSupport_LowerBound_6_3: 1.261904761904762 xi[6] + zi[6]
   - 100000 gamma[6] >= -99947
 Beta1_6_0: xi[6] - x_i_prime[0] + 100000 beta1[6,0] <= 100000
 Beta1_Z_Upper_6_0: zi[6] - z_i_prime[0] + 100000 beta1[6,0] <= 100000
 Beta1_Z_Lower_6_0: zi[6] - z_i_prime[0] - 100000 beta1[6,0] >= -100000
 Beta2_6_0: - xi[0] + x_i_prime[6] - 100000 beta2[6,0] >= -100000
 Beta2_Z_Upper_6_0: - zi[0] + z_i_prime[6] + 100000 beta2[6,0] <= 100000
 Beta2_Z_Lower_6_0: - zi[0] + z_i_prime[6] - 100000 beta2[6,0] >= -100000
 Beta1_6_1: xi[6] - x_i_prime[1] + 100000 beta1[6,1] <= 100000
 Beta1_Z_Upper_6_1: zi[6] - z_i_prime[1] + 100000 beta1[6,1] <= 100000
 Beta1_Z_Lower_6_1: zi[6] - z_i_prime[1] - 100000 beta1[6,1] >= -100000
 Beta2_6_1: - xi[1] + x_i_prime[6] - 100000 beta2[6,1] >= -100000
 Beta2_Z_Upper_6_1: - zi[1] + z_i_prime[6] + 100000 beta2[6,1] <= 100000
 Beta2_Z_Lower_6_1: - zi[1] + z_i_prime[6] - 100000 beta2[6,1] >= -100000
 Beta1_6_2: xi[6] - x_i_prime[2] + 100000 beta1[6,2] <= 100000
 Beta1_Z_Upper_6_2: zi[6] - z_i_prime[2] + 100000 beta1[6,2] <= 100000
 Beta1_Z_Lower_6_2: zi[6] - z_i_prime[2] - 100000 beta1[6,2] >= -100000
 Beta2_6_2: - xi[2] + x_i_prime[6] - 100000 beta2[6,2] >= -100000
 Beta2_Z_Upper_6_2: - zi[2] + z_i_prime[6] + 100000 beta2[6,2] <= 100000
 Beta2_Z_Lower_6_2: - zi[2] + z_i_prime[6] - 100000 beta2[6,2] >= -100000
 Beta1_6_3: xi[6] - x_i_prime[3] + 100000 beta1[6,3] <= 100000
 Beta1_Z_Upper_6_3: zi[6] - z_i_prime[3] + 100000 beta1[6,3] <= 100000
 Beta1_Z_Lower_6_3: zi[6] - z_i_prime[3] - 100000 beta1[6,3] >= -100000
 Beta2_6_3: - xi[3] + x_i_prime[6] - 100000 beta2[6,3] >= -100000
 Beta2_Z_Upper_6_3: - zi[3] + z_i_prime[6] + 100000 beta2[6,3] <= 100000
 Beta2_Z_Lower_6_3: - zi[3] + z_i_prime[6] - 100000 beta2[6,3] >= -100000
 Beta1_6_4: xi[6] - x_i_prime[4] + 100000 beta1[6,4] <= 100000
 Beta1_Z_Upper_6_4: zi[6] - z_i_prime[4] + 100000 beta1[6,4] <= 100000
 Beta1_Z_Lower_6_4: zi[6] - z_i_prime[4] - 100000 beta1[6,4] >= -100000
 Beta2_6_4: - xi[4] + x_i_prime[6] - 100000 beta2[6,4] >= -100000
 Beta2_Z_Upper_6_4: - zi[4] + z_i_prime[6] + 100000 beta2[6,4] <= 100000
 Beta2_Z_Lower_6_4: - zi[4] + z_i_prime[6] - 100000 beta2[6,4] >= -100000
 Beta1_6_5: xi[6] - x_i_prime[5] + 100000 beta1[6,5] <= 100000
 Beta1_Z_Upper_6_5: zi[6] - z_i_prime[5] + 100000 beta1[6,5] <= 100000
 Beta1_Z_Lower_6_5: zi[6] - z_i_prime[5] - 100000 beta1[6,5] >= -100000
 Beta2_6_5: - xi[5] + x_i_prime[6] - 100000 beta2[6,5] >= -100000
 Beta2_Z_Upper_6_5: - zi[5] + z_i_prime[6] + 100000 beta2[6,5] <= 100000
 Beta2_Z_Lower_6_5: - zi[5] + z_i_prime[6] - 100000 beta2[6,5] >= -100000
 Beta1_6_7: xi[6] - x_i_prime[7] + 100000 beta1[6,7] <= 100000
 Beta1_Z_Upper_6_7: zi[6] - z_i_prime[7] + 100000 beta1[6,7] <= 100000
 Beta1_Z_Lower_6_7: zi[6] - z_i_prime[7] - 100000 beta1[6,7] >= -100000
 Beta2_6_7: - xi[7] + x_i_prime[6] - 100000 beta2[6,7] >= -100000
 Beta2_Z_Upper_6_7: - zi[7] + z_i_prime[6] + 100000 beta2[6,7] <= 100000
 Beta2_Z_Lower_6_7: - zi[7] + z_i_prime[6] - 100000 beta2[6,7] >= -100000
 GroundSupport_UpperBound_7: zi[7] + 100000 g[7] <= 100000
 CutSupport_UpperBound_7_2: 1.261904761904762 xi[7] + zi[7]
   + 100000 gamma[7] <= 100053
 CutSupport_LowerBound_7_2: 1.261904761904762 xi[7] + zi[7]
   - 100000 gamma[7] >= -99947
 CutSupport_UpperBound_7_3: 1.261904761904762 xi[7] + zi[7]
   + 100000 gamma[7] <= 100053
 CutSupport_LowerBound_7_3: 1.261904761904762 xi[7] + zi[7]
   - 100000 gamma[7] >= -99947
 Beta1_7_0: xi[7] - x_i_prime[0] + 100000 beta1[7,0] <= 100000
 Beta1_Z_Upper_7_0: zi[7] - z_i_prime[0] + 100000 beta1[7,0] <= 100000
 Beta1_Z_Lower_7_0: zi[7] - z_i_prime[0] - 100000 beta1[7,0] >= -100000
 Beta2_7_0: - xi[0] + x_i_prime[7] - 100000 beta2[7,0] >= -100000
 Beta2_Z_Upper_7_0: - zi[0] + z_i_prime[7] + 100000 beta2[7,0] <= 100000
 Beta2_Z_Lower_7_0: - zi[0] + z_i_prime[7] - 100000 beta2[7,0] >= -100000
 Beta1_7_1: xi[7] - x_i_prime[1] + 100000 beta1[7,1] <= 100000
 Beta1_Z_Upper_7_1: zi[7] - z_i_prime[1] + 100000 beta1[7,1] <= 100000
 Beta1_Z_Lower_7_1: zi[7] - z_i_prime[1] - 100000 beta1[7,1] >= -100000
 Beta2_7_1: - xi[1] + x_i_prime[7] - 100000 beta2[7,1] >= -100000
 Beta2_Z_Upper_7_1: - zi[1] + z_i_prime[7] + 100000 beta2[7,1] <= 100000
 Beta2_Z_Lower_7_1: - zi[1] + z_i_prime[7] - 100000 beta2[7,1] >= -100000
 Beta1_7_2: xi[7] - x_i_prime[2] + 100000 beta1[7,2] <= 100000
 Beta1_Z_Upper_7_2: zi[7] - z_i_prime[2] + 100000 beta1[7,2] <= 100000
 Beta1_Z_Lower_7_2: zi[7] - z_i_prime[2] - 100000 beta1[7,2] >= -100000
 Beta2_7_2: - xi[2] + x_i_prime[7] - 100000 beta2[7,2] >= -100000
 Beta2_Z_Upper_7_2: - zi[2] + z_i_prime[7] + 100000 beta2[7,2] <= 100000
 Beta2_Z_Lower_7_2: - zi[2] + z_i_prime[7] - 100000 beta2[7,2] >= -100000
 Beta1_7_3: xi[7] - x_i_prime[3] + 100000 beta1[7,3] <= 100000
 Beta1_Z_Upper_7_3: zi[7] - z_i_prime[3] + 100000 beta1[7,3] <= 100000
 Beta1_Z_Lower_7_3: zi[7] - z_i_prime[3] - 100000 beta1[7,3] >= -100000
 Beta2_7_3: - xi[3] + x_i_prime[7] - 100000 beta2[7,3] >= -100000
 Beta2_Z_Upper_7_3: - zi[3] + z_i_prime[7] + 100000 beta2[7,3] <= 100000
 Beta2_Z_Lower_7_3: - zi[3] + z_i_prime[7] - 100000 beta2[7,3] >= -100000
 Beta1_7_4: xi[7] - x_i_prime[4] + 100000 beta1[7,4] <= 100000
 Beta1_Z_Upper_7_4: zi[7] - z_i_prime[4] + 100000 beta1[7,4] <= 100000
 Beta1_Z_Lower_7_4: zi[7] - z_i_prime[4] - 100000 beta1[7,4] >= -100000
 Beta2_7_4: - xi[4] + x_i_prime[7] - 100000 beta2[7,4] >= -100000
 Beta2_Z_Upper_7_4: - zi[4] + z_i_prime[7] + 100000 beta2[7,4] <= 100000
 Beta2_Z_Lower_7_4: - zi[4] + z_i_prime[7] - 100000 beta2[7,4] >= -100000
 Beta1_7_5: xi[7] - x_i_prime[5] + 100000 beta1[7,5] <= 100000
 Beta1_Z_Upper_7_5: zi[7] - z_i_prime[5] + 100000 beta1[7,5] <= 100000
 Beta1_Z_Lower_7_5: zi[7] - z_i_prime[5] - 100000 beta1[7,5] >= -100000
 Beta2_7_5: - xi[5] + x_i_prime[7] - 100000 beta2[7,5] >= -100000
 Beta2_Z_Upper_7_5: - zi[5] + z_i_prime[7] + 100000 beta2[7,5] <= 100000
 Beta2_Z_Lower_7_5: - zi[5] + z_i_prime[7] - 100000 beta2[7,5] >= -100000
 Beta1_7_6: xi[7] - x_i_prime[6] + 100000 beta1[7,6] <= 100000
 Beta1_Z_Upper_7_6: zi[7] - z_i_prime[6] + 100000 beta1[7,6] <= 100000
 Beta1_Z_Lower_7_6: zi[7] - z_i_prime[6] - 100000 beta1[7,6] >= -100000
 Beta2_7_6: - xi[6] + x_i_prime[7] - 100000 beta2[7,6] >= -100000
 Beta2_Z_Upper_7_6: - zi[6] + z_i_prime[7] + 100000 beta2[7,6] <= 100000
 Beta2_Z_Lower_7_6: - zi[6] + z_i_prime[7] - 100000 beta2[7,6] >= -100000
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
 AreaConstraint_5_0: 2494 p_ij[5,0] - 46500 u_j[0] <= 0
 AreaConstraint_5_1: 2494 p_ij[5,1] - 46500 u_j[1] <= 0
 AreaConstraint_5_2: 2494 p_ij[5,2] - 29760 u_j[2] <= 0
 AreaConstraint_5_3: 2494 p_ij[5,3] - 29760 u_j[3] <= 0
 AreaConstraint_6_0: 2262 p_ij[6,0] - 46500 u_j[0] <= 0
 AreaConstraint_6_1: 2262 p_ij[6,1] - 46500 u_j[1] <= 0
 AreaConstraint_6_2: 2262 p_ij[6,2] - 29760 u_j[2] <= 0
 AreaConstraint_6_3: 2262 p_ij[6,3] - 29760 u_j[3] <= 0
 AreaConstraint_7_0: 2838 p_ij[7,0] - 46500 u_j[0] <= 0
 AreaConstraint_7_1: 2838 p_ij[7,1] - 46500 u_j[1] <= 0
 AreaConstraint_7_2: 2838 p_ij[7,2] - 29760 u_j[2] <= 0
 AreaConstraint_7_3: 2838 p_ij[7,3] - 29760 u_j[3] <= 0
 OneItemOneBin_0: p_ij[0,0] + p_ij[0,1] + p_ij[0,2] + p_ij[0,3] = 1
 OneItemOneBin_1: p_ij[1,0] + p_ij[1,1] + p_ij[1,2] + p_ij[1,3] = 1
 OneItemOneBin_2: p_ij[2,0] + p_ij[2,1] + p_ij[2,2] + p_ij[2,3] = 1
 OneItemOneBin_3: p_ij[3,0] + p_ij[3,1] + p_ij[3,2] + p_ij[3,3] = 1
 OneItemOneBin_4: p_ij[4,0] + p_ij[4,1] + p_ij[4,2] + p_ij[4,3] = 1
 OneItemOneBin_5: p_ij[5,0] + p_ij[5,1] + p_ij[5,2] + p_ij[5,3] = 1
 OneItemOneBin_6: p_ij[6,0] + p_ij[6,1] + p_ij[6,2] + p_ij[6,3] = 1
 OneItemOneBin_7: p_ij[7,0] + p_ij[7,1] + p_ij[7,2] + p_ij[7,3] = 1
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
 ItemFitToBin_5: - 300 p_ij[5,0] - 300 p_ij[5,1] - 192 p_ij[5,2]
   - 192 p_ij[5,3] + x_i_prime[5] <= 0
 ItemFitToBin_6: - 300 p_ij[6,0] - 300 p_ij[6,1] - 192 p_ij[6,2]
   - 192 p_ij[6,3] + x_i_prime[6] <= 0
 ItemFitToBin_7: - 300 p_ij[7,0] - 300 p_ij[7,1] - 192 p_ij[7,2]
   - 192 p_ij[7,3] + x_i_prime[7] <= 0
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
 ItemFitToBinZ_5: - 155 p_ij[5,0] - 155 p_ij[5,1] - 155 p_ij[5,2]
   - 155 p_ij[5,3] + z_i_prime[5] <= 0
 ItemFitToBinZ_6: - 155 p_ij[6,0] - 155 p_ij[6,1] - 155 p_ij[6,2]
   - 155 p_ij[6,3] + z_i_prime[6] <= 0
 ItemFitToBinZ_7: - 155 p_ij[7,0] - 155 p_ij[7,1] - 155 p_ij[7,2]
   - 155 p_ij[7,3] + z_i_prime[7] <= 0
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
 TransformX_5: - xi[5] + x_i_prime[5] - 86 r[5,0,0] - 29 r[5,0,1] = 0
 TransformZ_5: - zi[5] + z_i_prime[5] - 86 r[5,1,0] - 29 r[5,1,1] = 0
 TransformX_6: - xi[6] + x_i_prime[6] - 78 r[6,0,0] - 29 r[6,0,1] = 0
 TransformZ_6: - zi[6] + z_i_prime[6] - 78 r[6,1,0] - 29 r[6,1,1] = 0
 TransformX_7: - xi[7] + x_i_prime[7] - 66 r[7,0,0] - 43 r[7,0,1] = 0
 TransformZ_7: - zi[7] + z_i_prime[7] - 66 r[7,1,0] - 43 r[7,1,1] = 0
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
 LengthAlign_5: r[5,0,0] + r[5,1,0] = 1
 HeightAlign_5: r[5,0,1] + r[5,1,1] = 1
 XAxisAlign_5: r[5,0,0] + r[5,0,1] = 1
 ZAxisAlign_5: r[5,1,0] + r[5,1,1] = 1
 LengthAlign_6: r[6,0,0] + r[6,1,0] = 1
 HeightAlign_6: r[6,0,1] + r[6,1,1] = 1
 XAxisAlign_6: r[6,0,0] + r[6,0,1] = 1
 ZAxisAlign_6: r[6,1,0] + r[6,1,1] = 1
 LengthAlign_7: r[7,0,0] + r[7,1,0] = 1
 HeightAlign_7: r[7,0,1] + r[7,1,1] = 1
 XAxisAlign_7: r[7,0,0] + r[7,0,1] = 1
 ZAxisAlign_7: r[7,1,0] + r[7,1,1] = 1
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
 OverlapInBin_0_5_toBin_0: - p_ij[0,0] - p_ij[5,0] + x_p[0,5] + x_p[5,0]
   + z_p[0,5] + z_p[5,0] >= -1
 OverlapInBin_0_5_toBin_1: - p_ij[0,1] - p_ij[5,1] + x_p[0,5] + x_p[5,0]
   + z_p[0,5] + z_p[5,0] >= -1
 OverlapInBin_0_5_toBin_2: - p_ij[0,2] - p_ij[5,2] + x_p[0,5] + x_p[5,0]
   + z_p[0,5] + z_p[5,0] >= -1
 OverlapInBin_0_5_toBin_3: - p_ij[0,3] - p_ij[5,3] + x_p[0,5] + x_p[5,0]
   + z_p[0,5] + z_p[5,0] >= -1
 OverlapInBin_0_6_toBin_0: - p_ij[0,0] - p_ij[6,0] + x_p[0,6] + x_p[6,0]
   + z_p[0,6] + z_p[6,0] >= -1
 OverlapInBin_0_6_toBin_1: - p_ij[0,1] - p_ij[6,1] + x_p[0,6] + x_p[6,0]
   + z_p[0,6] + z_p[6,0] >= -1
 OverlapInBin_0_6_toBin_2: - p_ij[0,2] - p_ij[6,2] + x_p[0,6] + x_p[6,0]
   + z_p[0,6] + z_p[6,0] >= -1
 OverlapInBin_0_6_toBin_3: - p_ij[0,3] - p_ij[6,3] + x_p[0,6] + x_p[6,0]
   + z_p[0,6] + z_p[6,0] >= -1
 OverlapInBin_0_7_toBin_0: - p_ij[0,0] - p_ij[7,0] + x_p[0,7] + x_p[7,0]
   + z_p[0,7] + z_p[7,0] >= -1
 OverlapInBin_0_7_toBin_1: - p_ij[0,1] - p_ij[7,1] + x_p[0,7] + x_p[7,0]
   + z_p[0,7] + z_p[7,0] >= -1
 OverlapInBin_0_7_toBin_2: - p_ij[0,2] - p_ij[7,2] + x_p[0,7] + x_p[7,0]
   + z_p[0,7] + z_p[7,0] >= -1
 OverlapInBin_0_7_toBin_3: - p_ij[0,3] - p_ij[7,3] + x_p[0,7] + x_p[7,0]
   + z_p[0,7] + z_p[7,0] >= -1
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
 OverlapInBin_1_5_toBin_0: - p_ij[1,0] - p_ij[5,0] + x_p[1,5] + x_p[5,1]
   + z_p[1,5] + z_p[5,1] >= -1
 OverlapInBin_1_5_toBin_1: - p_ij[1,1] - p_ij[5,1] + x_p[1,5] + x_p[5,1]
   + z_p[1,5] + z_p[5,1] >= -1
 OverlapInBin_1_5_toBin_2: - p_ij[1,2] - p_ij[5,2] + x_p[1,5] + x_p[5,1]
   + z_p[1,5] + z_p[5,1] >= -1
 OverlapInBin_1_5_toBin_3: - p_ij[1,3] - p_ij[5,3] + x_p[1,5] + x_p[5,1]
   + z_p[1,5] + z_p[5,1] >= -1
 OverlapInBin_1_6_toBin_0: - p_ij[1,0] - p_ij[6,0] + x_p[1,6] + x_p[6,1]
   + z_p[1,6] + z_p[6,1] >= -1
 OverlapInBin_1_6_toBin_1: - p_ij[1,1] - p_ij[6,1] + x_p[1,6] + x_p[6,1]
   + z_p[1,6] + z_p[6,1] >= -1
 OverlapInBin_1_6_toBin_2: - p_ij[1,2] - p_ij[6,2] + x_p[1,6] + x_p[6,1]
   + z_p[1,6] + z_p[6,1] >= -1
 OverlapInBin_1_6_toBin_3: - p_ij[1,3] - p_ij[6,3] + x_p[1,6] + x_p[6,1]
   + z_p[1,6] + z_p[6,1] >= -1
 OverlapInBin_1_7_toBin_0: - p_ij[1,0] - p_ij[7,0] + x_p[1,7] + x_p[7,1]
   + z_p[1,7] + z_p[7,1] >= -1
 OverlapInBin_1_7_toBin_1: - p_ij[1,1] - p_ij[7,1] + x_p[1,7] + x_p[7,1]
   + z_p[1,7] + z_p[7,1] >= -1
 OverlapInBin_1_7_toBin_2: - p_ij[1,2] - p_ij[7,2] + x_p[1,7] + x_p[7,1]
   + z_p[1,7] + z_p[7,1] >= -1
 OverlapInBin_1_7_toBin_3: - p_ij[1,3] - p_ij[7,3] + x_p[1,7] + x_p[7,1]
   + z_p[1,7] + z_p[7,1] >= -1
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
 OverlapInBin_2_5_toBin_0: - p_ij[2,0] - p_ij[5,0] + x_p[2,5] + x_p[5,2]
   + z_p[2,5] + z_p[5,2] >= -1
 OverlapInBin_2_5_toBin_1: - p_ij[2,1] - p_ij[5,1] + x_p[2,5] + x_p[5,2]
   + z_p[2,5] + z_p[5,2] >= -1
 OverlapInBin_2_5_toBin_2: - p_ij[2,2] - p_ij[5,2] + x_p[2,5] + x_p[5,2]
   + z_p[2,5] + z_p[5,2] >= -1
 OverlapInBin_2_5_toBin_3: - p_ij[2,3] - p_ij[5,3] + x_p[2,5] + x_p[5,2]
   + z_p[2,5] + z_p[5,2] >= -1
 OverlapInBin_2_6_toBin_0: - p_ij[2,0] - p_ij[6,0] + x_p[2,6] + x_p[6,2]
   + z_p[2,6] + z_p[6,2] >= -1
 OverlapInBin_2_6_toBin_1: - p_ij[2,1] - p_ij[6,1] + x_p[2,6] + x_p[6,2]
   + z_p[2,6] + z_p[6,2] >= -1
 OverlapInBin_2_6_toBin_2: - p_ij[2,2] - p_ij[6,2] + x_p[2,6] + x_p[6,2]
   + z_p[2,6] + z_p[6,2] >= -1
 OverlapInBin_2_6_toBin_3: - p_ij[2,3] - p_ij[6,3] + x_p[2,6] + x_p[6,2]
   + z_p[2,6] + z_p[6,2] >= -1
 OverlapInBin_2_7_toBin_0: - p_ij[2,0] - p_ij[7,0] + x_p[2,7] + x_p[7,2]
   + z_p[2,7] + z_p[7,2] >= -1
 OverlapInBin_2_7_toBin_1: - p_ij[2,1] - p_ij[7,1] + x_p[2,7] + x_p[7,2]
   + z_p[2,7] + z_p[7,2] >= -1
 OverlapInBin_2_7_toBin_2: - p_ij[2,2] - p_ij[7,2] + x_p[2,7] + x_p[7,2]
   + z_p[2,7] + z_p[7,2] >= -1
 OverlapInBin_2_7_toBin_3: - p_ij[2,3] - p_ij[7,3] + x_p[2,7] + x_p[7,2]
   + z_p[2,7] + z_p[7,2] >= -1
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
 OverlapInBin_3_5_toBin_0: - p_ij[3,0] - p_ij[5,0] + x_p[3,5] + x_p[5,3]
   + z_p[3,5] + z_p[5,3] >= -1
 OverlapInBin_3_5_toBin_1: - p_ij[3,1] - p_ij[5,1] + x_p[3,5] + x_p[5,3]
   + z_p[3,5] + z_p[5,3] >= -1
 OverlapInBin_3_5_toBin_2: - p_ij[3,2] - p_ij[5,2] + x_p[3,5] + x_p[5,3]
   + z_p[3,5] + z_p[5,3] >= -1
 OverlapInBin_3_5_toBin_3: - p_ij[3,3] - p_ij[5,3] + x_p[3,5] + x_p[5,3]
   + z_p[3,5] + z_p[5,3] >= -1
 OverlapInBin_3_6_toBin_0: - p_ij[3,0] - p_ij[6,0] + x_p[3,6] + x_p[6,3]
   + z_p[3,6] + z_p[6,3] >= -1
 OverlapInBin_3_6_toBin_1: - p_ij[3,1] - p_ij[6,1] + x_p[3,6] + x_p[6,3]
   + z_p[3,6] + z_p[6,3] >= -1
 OverlapInBin_3_6_toBin_2: - p_ij[3,2] - p_ij[6,2] + x_p[3,6] + x_p[6,3]
   + z_p[3,6] + z_p[6,3] >= -1
 OverlapInBin_3_6_toBin_3: - p_ij[3,3] - p_ij[6,3] + x_p[3,6] + x_p[6,3]
   + z_p[3,6] + z_p[6,3] >= -1
 OverlapInBin_3_7_toBin_0: - p_ij[3,0] - p_ij[7,0] + x_p[3,7] + x_p[7,3]
   + z_p[3,7] + z_p[7,3] >= -1
 OverlapInBin_3_7_toBin_1: - p_ij[3,1] - p_ij[7,1] + x_p[3,7] + x_p[7,3]
   + z_p[3,7] + z_p[7,3] >= -1
 OverlapInBin_3_7_toBin_2: - p_ij[3,2] - p_ij[7,2] + x_p[3,7] + x_p[7,3]
   + z_p[3,7] + z_p[7,3] >= -1
 OverlapInBin_3_7_toBin_3: - p_ij[3,3] - p_ij[7,3] + x_p[3,7] + x_p[7,3]
   + z_p[3,7] + z_p[7,3] >= -1
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
 OverlapInBin_4_5_toBin_0: - p_ij[4,0] - p_ij[5,0] + x_p[4,5] + x_p[5,4]
   + z_p[4,5] + z_p[5,4] >= -1
 OverlapInBin_4_5_toBin_1: - p_ij[4,1] - p_ij[5,1] + x_p[4,5] + x_p[5,4]
   + z_p[4,5] + z_p[5,4] >= -1
 OverlapInBin_4_5_toBin_2: - p_ij[4,2] - p_ij[5,2] + x_p[4,5] + x_p[5,4]
   + z_p[4,5] + z_p[5,4] >= -1
 OverlapInBin_4_5_toBin_3: - p_ij[4,3] - p_ij[5,3] + x_p[4,5] + x_p[5,4]
   + z_p[4,5] + z_p[5,4] >= -1
 OverlapInBin_4_6_toBin_0: - p_ij[4,0] - p_ij[6,0] + x_p[4,6] + x_p[6,4]
   + z_p[4,6] + z_p[6,4] >= -1
 OverlapInBin_4_6_toBin_1: - p_ij[4,1] - p_ij[6,1] + x_p[4,6] + x_p[6,4]
   + z_p[4,6] + z_p[6,4] >= -1
 OverlapInBin_4_6_toBin_2: - p_ij[4,2] - p_ij[6,2] + x_p[4,6] + x_p[6,4]
   + z_p[4,6] + z_p[6,4] >= -1
 OverlapInBin_4_6_toBin_3: - p_ij[4,3] - p_ij[6,3] + x_p[4,6] + x_p[6,4]
   + z_p[4,6] + z_p[6,4] >= -1
 OverlapInBin_4_7_toBin_0: - p_ij[4,0] - p_ij[7,0] + x_p[4,7] + x_p[7,4]
   + z_p[4,7] + z_p[7,4] >= -1
 OverlapInBin_4_7_toBin_1: - p_ij[4,1] - p_ij[7,1] + x_p[4,7] + x_p[7,4]
   + z_p[4,7] + z_p[7,4] >= -1
 OverlapInBin_4_7_toBin_2: - p_ij[4,2] - p_ij[7,2] + x_p[4,7] + x_p[7,4]
   + z_p[4,7] + z_p[7,4] >= -1
 OverlapInBin_4_7_toBin_3: - p_ij[4,3] - p_ij[7,3] + x_p[4,7] + x_p[7,4]
   + z_p[4,7] + z_p[7,4] >= -1
 OverlapInBin_5_0_toBin_0: - p_ij[0,0] - p_ij[5,0] + x_p[0,5] + x_p[5,0]
   + z_p[0,5] + z_p[5,0] >= -1
 OverlapInBin_5_0_toBin_1: - p_ij[0,1] - p_ij[5,1] + x_p[0,5] + x_p[5,0]
   + z_p[0,5] + z_p[5,0] >= -1
 OverlapInBin_5_0_toBin_2: - p_ij[0,2] - p_ij[5,2] + x_p[0,5] + x_p[5,0]
   + z_p[0,5] + z_p[5,0] >= -1
 OverlapInBin_5_0_toBin_3: - p_ij[0,3] - p_ij[5,3] + x_p[0,5] + x_p[5,0]
   + z_p[0,5] + z_p[5,0] >= -1
 OverlapInBin_5_1_toBin_0: - p_ij[1,0] - p_ij[5,0] + x_p[1,5] + x_p[5,1]
   + z_p[1,5] + z_p[5,1] >= -1
 OverlapInBin_5_1_toBin_1: - p_ij[1,1] - p_ij[5,1] + x_p[1,5] + x_p[5,1]
   + z_p[1,5] + z_p[5,1] >= -1
 OverlapInBin_5_1_toBin_2: - p_ij[1,2] - p_ij[5,2] + x_p[1,5] + x_p[5,1]
   + z_p[1,5] + z_p[5,1] >= -1
 OverlapInBin_5_1_toBin_3: - p_ij[1,3] - p_ij[5,3] + x_p[1,5] + x_p[5,1]
   + z_p[1,5] + z_p[5,1] >= -1
 OverlapInBin_5_2_toBin_0: - p_ij[2,0] - p_ij[5,0] + x_p[2,5] + x_p[5,2]
   + z_p[2,5] + z_p[5,2] >= -1
 OverlapInBin_5_2_toBin_1: - p_ij[2,1] - p_ij[5,1] + x_p[2,5] + x_p[5,2]
   + z_p[2,5] + z_p[5,2] >= -1
 OverlapInBin_5_2_toBin_2: - p_ij[2,2] - p_ij[5,2] + x_p[2,5] + x_p[5,2]
   + z_p[2,5] + z_p[5,2] >= -1
 OverlapInBin_5_2_toBin_3: - p_ij[2,3] - p_ij[5,3] + x_p[2,5] + x_p[5,2]
   + z_p[2,5] + z_p[5,2] >= -1
 OverlapInBin_5_3_toBin_0: - p_ij[3,0] - p_ij[5,0] + x_p[3,5] + x_p[5,3]
   + z_p[3,5] + z_p[5,3] >= -1
 OverlapInBin_5_3_toBin_1: - p_ij[3,1] - p_ij[5,1] + x_p[3,5] + x_p[5,3]
   + z_p[3,5] + z_p[5,3] >= -1
 OverlapInBin_5_3_toBin_2: - p_ij[3,2] - p_ij[5,2] + x_p[3,5] + x_p[5,3]
   + z_p[3,5] + z_p[5,3] >= -1
 OverlapInBin_5_3_toBin_3: - p_ij[3,3] - p_ij[5,3] + x_p[3,5] + x_p[5,3]
   + z_p[3,5] + z_p[5,3] >= -1
 OverlapInBin_5_4_toBin_0: - p_ij[4,0] - p_ij[5,0] + x_p[4,5] + x_p[5,4]
   + z_p[4,5] + z_p[5,4] >= -1
 OverlapInBin_5_4_toBin_1: - p_ij[4,1] - p_ij[5,1] + x_p[4,5] + x_p[5,4]
   + z_p[4,5] + z_p[5,4] >= -1
 OverlapInBin_5_4_toBin_2: - p_ij[4,2] - p_ij[5,2] + x_p[4,5] + x_p[5,4]
   + z_p[4,5] + z_p[5,4] >= -1
 OverlapInBin_5_4_toBin_3: - p_ij[4,3] - p_ij[5,3] + x_p[4,5] + x_p[5,4]
   + z_p[4,5] + z_p[5,4] >= -1
 OverlapInBin_5_6_toBin_0: - p_ij[5,0] - p_ij[6,0] + x_p[5,6] + x_p[6,5]
   + z_p[5,6] + z_p[6,5] >= -1
 OverlapInBin_5_6_toBin_1: - p_ij[5,1] - p_ij[6,1] + x_p[5,6] + x_p[6,5]
   + z_p[5,6] + z_p[6,5] >= -1
 OverlapInBin_5_6_toBin_2: - p_ij[5,2] - p_ij[6,2] + x_p[5,6] + x_p[6,5]
   + z_p[5,6] + z_p[6,5] >= -1
 OverlapInBin_5_6_toBin_3: - p_ij[5,3] - p_ij[6,3] + x_p[5,6] + x_p[6,5]
   + z_p[5,6] + z_p[6,5] >= -1
 OverlapInBin_5_7_toBin_0: - p_ij[5,0] - p_ij[7,0] + x_p[5,7] + x_p[7,5]
   + z_p[5,7] + z_p[7,5] >= -1
 OverlapInBin_5_7_toBin_1: - p_ij[5,1] - p_ij[7,1] + x_p[5,7] + x_p[7,5]
   + z_p[5,7] + z_p[7,5] >= -1
 OverlapInBin_5_7_toBin_2: - p_ij[5,2] - p_ij[7,2] + x_p[5,7] + x_p[7,5]
   + z_p[5,7] + z_p[7,5] >= -1
 OverlapInBin_5_7_toBin_3: - p_ij[5,3] - p_ij[7,3] + x_p[5,7] + x_p[7,5]
   + z_p[5,7] + z_p[7,5] >= -1
 OverlapInBin_6_0_toBin_0: - p_ij[0,0] - p_ij[6,0] + x_p[0,6] + x_p[6,0]
   + z_p[0,6] + z_p[6,0] >= -1
 OverlapInBin_6_0_toBin_1: - p_ij[0,1] - p_ij[6,1] + x_p[0,6] + x_p[6,0]
   + z_p[0,6] + z_p[6,0] >= -1
 OverlapInBin_6_0_toBin_2: - p_ij[0,2] - p_ij[6,2] + x_p[0,6] + x_p[6,0]
   + z_p[0,6] + z_p[6,0] >= -1
 OverlapInBin_6_0_toBin_3: - p_ij[0,3] - p_ij[6,3] + x_p[0,6] + x_p[6,0]
   + z_p[0,6] + z_p[6,0] >= -1
 OverlapInBin_6_1_toBin_0: - p_ij[1,0] - p_ij[6,0] + x_p[1,6] + x_p[6,1]
   + z_p[1,6] + z_p[6,1] >= -1
 OverlapInBin_6_1_toBin_1: - p_ij[1,1] - p_ij[6,1] + x_p[1,6] + x_p[6,1]
   + z_p[1,6] + z_p[6,1] >= -1
 OverlapInBin_6_1_toBin_2: - p_ij[1,2] - p_ij[6,2] + x_p[1,6] + x_p[6,1]
   + z_p[1,6] + z_p[6,1] >= -1
 OverlapInBin_6_1_toBin_3: - p_ij[1,3] - p_ij[6,3] + x_p[1,6] + x_p[6,1]
   + z_p[1,6] + z_p[6,1] >= -1
 OverlapInBin_6_2_toBin_0: - p_ij[2,0] - p_ij[6,0] + x_p[2,6] + x_p[6,2]
   + z_p[2,6] + z_p[6,2] >= -1
 OverlapInBin_6_2_toBin_1: - p_ij[2,1] - p_ij[6,1] + x_p[2,6] + x_p[6,2]
   + z_p[2,6] + z_p[6,2] >= -1
 OverlapInBin_6_2_toBin_2: - p_ij[2,2] - p_ij[6,2] + x_p[2,6] + x_p[6,2]
   + z_p[2,6] + z_p[6,2] >= -1
 OverlapInBin_6_2_toBin_3: - p_ij[2,3] - p_ij[6,3] + x_p[2,6] + x_p[6,2]
   + z_p[2,6] + z_p[6,2] >= -1
 OverlapInBin_6_3_toBin_0: - p_ij[3,0] - p_ij[6,0] + x_p[3,6] + x_p[6,3]
   + z_p[3,6] + z_p[6,3] >= -1
 OverlapInBin_6_3_toBin_1: - p_ij[3,1] - p_ij[6,1] + x_p[3,6] + x_p[6,3]
   + z_p[3,6] + z_p[6,3] >= -1
 OverlapInBin_6_3_toBin_2: - p_ij[3,2] - p_ij[6,2] + x_p[3,6] + x_p[6,3]
   + z_p[3,6] + z_p[6,3] >= -1
 OverlapInBin_6_3_toBin_3: - p_ij[3,3] - p_ij[6,3] + x_p[3,6] + x_p[6,3]
   + z_p[3,6] + z_p[6,3] >= -1
 OverlapInBin_6_4_toBin_0: - p_ij[4,0] - p_ij[6,0] + x_p[4,6] + x_p[6,4]
   + z_p[4,6] + z_p[6,4] >= -1
 OverlapInBin_6_4_toBin_1: - p_ij[4,1] - p_ij[6,1] + x_p[4,6] + x_p[6,4]
   + z_p[4,6] + z_p[6,4] >= -1
 OverlapInBin_6_4_toBin_2: - p_ij[4,2] - p_ij[6,2] + x_p[4,6] + x_p[6,4]
   + z_p[4,6] + z_p[6,4] >= -1
 OverlapInBin_6_4_toBin_3: - p_ij[4,3] - p_ij[6,3] + x_p[4,6] + x_p[6,4]
   + z_p[4,6] + z_p[6,4] >= -1
 OverlapInBin_6_5_toBin_0: - p_ij[5,0] - p_ij[6,0] + x_p[5,6] + x_p[6,5]
   + z_p[5,6] + z_p[6,5] >= -1
 OverlapInBin_6_5_toBin_1: - p_ij[5,1] - p_ij[6,1] + x_p[5,6] + x_p[6,5]
   + z_p[5,6] + z_p[6,5] >= -1
 OverlapInBin_6_5_toBin_2: - p_ij[5,2] - p_ij[6,2] + x_p[5,6] + x_p[6,5]
   + z_p[5,6] + z_p[6,5] >= -1
 OverlapInBin_6_5_toBin_3: - p_ij[5,3] - p_ij[6,3] + x_p[5,6] + x_p[6,5]
   + z_p[5,6] + z_p[6,5] >= -1
 OverlapInBin_6_7_toBin_0: - p_ij[6,0] - p_ij[7,0] + x_p[6,7] + x_p[7,6]
   + z_p[6,7] + z_p[7,6] >= -1
 OverlapInBin_6_7_toBin_1: - p_ij[6,1] - p_ij[7,1] + x_p[6,7] + x_p[7,6]
   + z_p[6,7] + z_p[7,6] >= -1
 OverlapInBin_6_7_toBin_2: - p_ij[6,2] - p_ij[7,2] + x_p[6,7] + x_p[7,6]
   + z_p[6,7] + z_p[7,6] >= -1
 OverlapInBin_6_7_toBin_3: - p_ij[6,3] - p_ij[7,3] + x_p[6,7] + x_p[7,6]
   + z_p[6,7] + z_p[7,6] >= -1
 OverlapInBin_7_0_toBin_0: - p_ij[0,0] - p_ij[7,0] + x_p[0,7] + x_p[7,0]
   + z_p[0,7] + z_p[7,0] >= -1
 OverlapInBin_7_0_toBin_1: - p_ij[0,1] - p_ij[7,1] + x_p[0,7] + x_p[7,0]
   + z_p[0,7] + z_p[7,0] >= -1
 OverlapInBin_7_0_toBin_2: - p_ij[0,2] - p_ij[7,2] + x_p[0,7] + x_p[7,0]
   + z_p[0,7] + z_p[7,0] >= -1
 OverlapInBin_7_0_toBin_3: - p_ij[0,3] - p_ij[7,3] + x_p[0,7] + x_p[7,0]
   + z_p[0,7] + z_p[7,0] >= -1
 OverlapInBin_7_1_toBin_0: - p_ij[1,0] - p_ij[7,0] + x_p[1,7] + x_p[7,1]
   + z_p[1,7] + z_p[7,1] >= -1
 OverlapInBin_7_1_toBin_1: - p_ij[1,1] - p_ij[7,1] + x_p[1,7] + x_p[7,1]
   + z_p[1,7] + z_p[7,1] >= -1
 OverlapInBin_7_1_toBin_2: - p_ij[1,2] - p_ij[7,2] + x_p[1,7] + x_p[7,1]
   + z_p[1,7] + z_p[7,1] >= -1
 OverlapInBin_7_1_toBin_3: - p_ij[1,3] - p_ij[7,3] + x_p[1,7] + x_p[7,1]
   + z_p[1,7] + z_p[7,1] >= -1
 OverlapInBin_7_2_toBin_0: - p_ij[2,0] - p_ij[7,0] + x_p[2,7] + x_p[7,2]
   + z_p[2,7] + z_p[7,2] >= -1
 OverlapInBin_7_2_toBin_1: - p_ij[2,1] - p_ij[7,1] + x_p[2,7] + x_p[7,2]
   + z_p[2,7] + z_p[7,2] >= -1
 OverlapInBin_7_2_toBin_2: - p_ij[2,2] - p_ij[7,2] + x_p[2,7] + x_p[7,2]
   + z_p[2,7] + z_p[7,2] >= -1
 OverlapInBin_7_2_toBin_3: - p_ij[2,3] - p_ij[7,3] + x_p[2,7] + x_p[7,2]
   + z_p[2,7] + z_p[7,2] >= -1
 OverlapInBin_7_3_toBin_0: - p_ij[3,0] - p_ij[7,0] + x_p[3,7] + x_p[7,3]
   + z_p[3,7] + z_p[7,3] >= -1
 OverlapInBin_7_3_toBin_1: - p_ij[3,1] - p_ij[7,1] + x_p[3,7] + x_p[7,3]
   + z_p[3,7] + z_p[7,3] >= -1
 OverlapInBin_7_3_toBin_2: - p_ij[3,2] - p_ij[7,2] + x_p[3,7] + x_p[7,3]
   + z_p[3,7] + z_p[7,3] >= -1
 OverlapInBin_7_3_toBin_3: - p_ij[3,3] - p_ij[7,3] + x_p[3,7] + x_p[7,3]
   + z_p[3,7] + z_p[7,3] >= -1
 OverlapInBin_7_4_toBin_0: - p_ij[4,0] - p_ij[7,0] + x_p[4,7] + x_p[7,4]
   + z_p[4,7] + z_p[7,4] >= -1
 OverlapInBin_7_4_toBin_1: - p_ij[4,1] - p_ij[7,1] + x_p[4,7] + x_p[7,4]
   + z_p[4,7] + z_p[7,4] >= -1
 OverlapInBin_7_4_toBin_2: - p_ij[4,2] - p_ij[7,2] + x_p[4,7] + x_p[7,4]
   + z_p[4,7] + z_p[7,4] >= -1
 OverlapInBin_7_4_toBin_3: - p_ij[4,3] - p_ij[7,3] + x_p[4,7] + x_p[7,4]
   + z_p[4,7] + z_p[7,4] >= -1
 OverlapInBin_7_5_toBin_0: - p_ij[5,0] - p_ij[7,0] + x_p[5,7] + x_p[7,5]
   + z_p[5,7] + z_p[7,5] >= -1
 OverlapInBin_7_5_toBin_1: - p_ij[5,1] - p_ij[7,1] + x_p[5,7] + x_p[7,5]
   + z_p[5,7] + z_p[7,5] >= -1
 OverlapInBin_7_5_toBin_2: - p_ij[5,2] - p_ij[7,2] + x_p[5,7] + x_p[7,5]
   + z_p[5,7] + z_p[7,5] >= -1
 OverlapInBin_7_5_toBin_3: - p_ij[5,3] - p_ij[7,3] + x_p[5,7] + x_p[7,5]
   + z_p[5,7] + z_p[7,5] >= -1
 OverlapInBin_7_6_toBin_0: - p_ij[6,0] - p_ij[7,0] + x_p[6,7] + x_p[7,6]
   + z_p[6,7] + z_p[7,6] >= -1
 OverlapInBin_7_6_toBin_1: - p_ij[6,1] - p_ij[7,1] + x_p[6,7] + x_p[7,6]
   + z_p[6,7] + z_p[7,6] >= -1
 OverlapInBin_7_6_toBin_2: - p_ij[6,2] - p_ij[7,2] + x_p[6,7] + x_p[7,6]
   + z_p[6,7] + z_p[7,6] >= -1
 OverlapInBin_7_6_toBin_3: - p_ij[6,3] - p_ij[7,3] + x_p[6,7] + x_p[7,6]
   + z_p[6,7] + z_p[7,6] >= -1
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
 Overlap_0_5: x_p[0,5] + x_p[5,0] + z_p[0,5] + z_p[5,0] >= 1
 RightLeft_0_5: x_p[0,5] + x_p[5,0] <= 1
 LeftRight_0_5: x_p[0,5] + x_p[5,0] <= 1
 AboveBelow_0_5: z_p[0,5] + z_p[5,0] <= 1
 BelowAbove_0_5: z_p[0,5] + z_p[5,0] <= 1
 Overlap_0_6: x_p[0,6] + x_p[6,0] + z_p[0,6] + z_p[6,0] >= 1
 RightLeft_0_6: x_p[0,6] + x_p[6,0] <= 1
 LeftRight_0_6: x_p[0,6] + x_p[6,0] <= 1
 AboveBelow_0_6: z_p[0,6] + z_p[6,0] <= 1
 BelowAbove_0_6: z_p[0,6] + z_p[6,0] <= 1
 Overlap_0_7: x_p[0,7] + x_p[7,0] + z_p[0,7] + z_p[7,0] >= 1
 RightLeft_0_7: x_p[0,7] + x_p[7,0] <= 1
 LeftRight_0_7: x_p[0,7] + x_p[7,0] <= 1
 AboveBelow_0_7: z_p[0,7] + z_p[7,0] <= 1
 BelowAbove_0_7: z_p[0,7] + z_p[7,0] <= 1
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
 Overlap_1_5: x_p[1,5] + x_p[5,1] + z_p[1,5] + z_p[5,1] >= 1
 RightLeft_1_5: x_p[1,5] + x_p[5,1] <= 1
 LeftRight_1_5: x_p[1,5] + x_p[5,1] <= 1
 AboveBelow_1_5: z_p[1,5] + z_p[5,1] <= 1
 BelowAbove_1_5: z_p[1,5] + z_p[5,1] <= 1
 Overlap_1_6: x_p[1,6] + x_p[6,1] + z_p[1,6] + z_p[6,1] >= 1
 RightLeft_1_6: x_p[1,6] + x_p[6,1] <= 1
 LeftRight_1_6: x_p[1,6] + x_p[6,1] <= 1
 AboveBelow_1_6: z_p[1,6] + z_p[6,1] <= 1
 BelowAbove_1_6: z_p[1,6] + z_p[6,1] <= 1
 Overlap_1_7: x_p[1,7] + x_p[7,1] + z_p[1,7] + z_p[7,1] >= 1
 RightLeft_1_7: x_p[1,7] + x_p[7,1] <= 1
 LeftRight_1_7: x_p[1,7] + x_p[7,1] <= 1
 AboveBelow_1_7: z_p[1,7] + z_p[7,1] <= 1
 BelowAbove_1_7: z_p[1,7] + z_p[7,1] <= 1
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
 Overlap_2_5: x_p[2,5] + x_p[5,2] + z_p[2,5] + z_p[5,2] >= 1
 RightLeft_2_5: x_p[2,5] + x_p[5,2] <= 1
 LeftRight_2_5: x_p[2,5] + x_p[5,2] <= 1
 AboveBelow_2_5: z_p[2,5] + z_p[5,2] <= 1
 BelowAbove_2_5: z_p[2,5] + z_p[5,2] <= 1
 Overlap_2_6: x_p[2,6] + x_p[6,2] + z_p[2,6] + z_p[6,2] >= 1
 RightLeft_2_6: x_p[2,6] + x_p[6,2] <= 1
 LeftRight_2_6: x_p[2,6] + x_p[6,2] <= 1
 AboveBelow_2_6: z_p[2,6] + z_p[6,2] <= 1
 BelowAbove_2_6: z_p[2,6] + z_p[6,2] <= 1
 Overlap_2_7: x_p[2,7] + x_p[7,2] + z_p[2,7] + z_p[7,2] >= 1
 RightLeft_2_7: x_p[2,7] + x_p[7,2] <= 1
 LeftRight_2_7: x_p[2,7] + x_p[7,2] <= 1
 AboveBelow_2_7: z_p[2,7] + z_p[7,2] <= 1
 BelowAbove_2_7: z_p[2,7] + z_p[7,2] <= 1
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
 Overlap_3_5: x_p[3,5] + x_p[5,3] + z_p[3,5] + z_p[5,3] >= 1
 RightLeft_3_5: x_p[3,5] + x_p[5,3] <= 1
 LeftRight_3_5: x_p[3,5] + x_p[5,3] <= 1
 AboveBelow_3_5: z_p[3,5] + z_p[5,3] <= 1
 BelowAbove_3_5: z_p[3,5] + z_p[5,3] <= 1
 Overlap_3_6: x_p[3,6] + x_p[6,3] + z_p[3,6] + z_p[6,3] >= 1
 RightLeft_3_6: x_p[3,6] + x_p[6,3] <= 1
 LeftRight_3_6: x_p[3,6] + x_p[6,3] <= 1
 AboveBelow_3_6: z_p[3,6] + z_p[6,3] <= 1
 BelowAbove_3_6: z_p[3,6] + z_p[6,3] <= 1
 Overlap_3_7: x_p[3,7] + x_p[7,3] + z_p[3,7] + z_p[7,3] >= 1
 RightLeft_3_7: x_p[3,7] + x_p[7,3] <= 1
 LeftRight_3_7: x_p[3,7] + x_p[7,3] <= 1
 AboveBelow_3_7: z_p[3,7] + z_p[7,3] <= 1
 BelowAbove_3_7: z_p[3,7] + z_p[7,3] <= 1
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
 Overlap_4_5: x_p[4,5] + x_p[5,4] + z_p[4,5] + z_p[5,4] >= 1
 RightLeft_4_5: x_p[4,5] + x_p[5,4] <= 1
 LeftRight_4_5: x_p[4,5] + x_p[5,4] <= 1
 AboveBelow_4_5: z_p[4,5] + z_p[5,4] <= 1
 BelowAbove_4_5: z_p[4,5] + z_p[5,4] <= 1
 Overlap_4_6: x_p[4,6] + x_p[6,4] + z_p[4,6] + z_p[6,4] >= 1
 RightLeft_4_6: x_p[4,6] + x_p[6,4] <= 1
 LeftRight_4_6: x_p[4,6] + x_p[6,4] <= 1
 AboveBelow_4_6: z_p[4,6] + z_p[6,4] <= 1
 BelowAbove_4_6: z_p[4,6] + z_p[6,4] <= 1
 Overlap_4_7: x_p[4,7] + x_p[7,4] + z_p[4,7] + z_p[7,4] >= 1
 RightLeft_4_7: x_p[4,7] + x_p[7,4] <= 1
 LeftRight_4_7: x_p[4,7] + x_p[7,4] <= 1
 AboveBelow_4_7: z_p[4,7] + z_p[7,4] <= 1
 BelowAbove_4_7: z_p[4,7] + z_p[7,4] <= 1
 Overlap_5_0: x_p[0,5] + x_p[5,0] + z_p[0,5] + z_p[5,0] >= 1
 RightLeft_5_0: x_p[0,5] + x_p[5,0] <= 1
 LeftRight_5_0: x_p[0,5] + x_p[5,0] <= 1
 AboveBelow_5_0: z_p[0,5] + z_p[5,0] <= 1
 BelowAbove_5_0: z_p[0,5] + z_p[5,0] <= 1
 Overlap_5_1: x_p[1,5] + x_p[5,1] + z_p[1,5] + z_p[5,1] >= 1
 RightLeft_5_1: x_p[1,5] + x_p[5,1] <= 1
 LeftRight_5_1: x_p[1,5] + x_p[5,1] <= 1
 AboveBelow_5_1: z_p[1,5] + z_p[5,1] <= 1
 BelowAbove_5_1: z_p[1,5] + z_p[5,1] <= 1
 Overlap_5_2: x_p[2,5] + x_p[5,2] + z_p[2,5] + z_p[5,2] >= 1
 RightLeft_5_2: x_p[2,5] + x_p[5,2] <= 1
 LeftRight_5_2: x_p[2,5] + x_p[5,2] <= 1
 AboveBelow_5_2: z_p[2,5] + z_p[5,2] <= 1
 BelowAbove_5_2: z_p[2,5] + z_p[5,2] <= 1
 Overlap_5_3: x_p[3,5] + x_p[5,3] + z_p[3,5] + z_p[5,3] >= 1
 RightLeft_5_3: x_p[3,5] + x_p[5,3] <= 1
 LeftRight_5_3: x_p[3,5] + x_p[5,3] <= 1
 AboveBelow_5_3: z_p[3,5] + z_p[5,3] <= 1
 BelowAbove_5_3: z_p[3,5] + z_p[5,3] <= 1
 Overlap_5_4: x_p[4,5] + x_p[5,4] + z_p[4,5] + z_p[5,4] >= 1
 RightLeft_5_4: x_p[4,5] + x_p[5,4] <= 1
 LeftRight_5_4: x_p[4,5] + x_p[5,4] <= 1
 AboveBelow_5_4: z_p[4,5] + z_p[5,4] <= 1
 BelowAbove_5_4: z_p[4,5] + z_p[5,4] <= 1
 Overlap_5_6: x_p[5,6] + x_p[6,5] + z_p[5,6] + z_p[6,5] >= 1
 RightLeft_5_6: x_p[5,6] + x_p[6,5] <= 1
 LeftRight_5_6: x_p[5,6] + x_p[6,5] <= 1
 AboveBelow_5_6: z_p[5,6] + z_p[6,5] <= 1
 BelowAbove_5_6: z_p[5,6] + z_p[6,5] <= 1
 Overlap_5_7: x_p[5,7] + x_p[7,5] + z_p[5,7] + z_p[7,5] >= 1
 RightLeft_5_7: x_p[5,7] + x_p[7,5] <= 1
 LeftRight_5_7: x_p[5,7] + x_p[7,5] <= 1
 AboveBelow_5_7: z_p[5,7] + z_p[7,5] <= 1
 BelowAbove_5_7: z_p[5,7] + z_p[7,5] <= 1
 Overlap_6_0: x_p[0,6] + x_p[6,0] + z_p[0,6] + z_p[6,0] >= 1
 RightLeft_6_0: x_p[0,6] + x_p[6,0] <= 1
 LeftRight_6_0: x_p[0,6] + x_p[6,0] <= 1
 AboveBelow_6_0: z_p[0,6] + z_p[6,0] <= 1
 BelowAbove_6_0: z_p[0,6] + z_p[6,0] <= 1
 Overlap_6_1: x_p[1,6] + x_p[6,1] + z_p[1,6] + z_p[6,1] >= 1
 RightLeft_6_1: x_p[1,6] + x_p[6,1] <= 1
 LeftRight_6_1: x_p[1,6] + x_p[6,1] <= 1
 AboveBelow_6_1: z_p[1,6] + z_p[6,1] <= 1
 BelowAbove_6_1: z_p[1,6] + z_p[6,1] <= 1
 Overlap_6_2: x_p[2,6] + x_p[6,2] + z_p[2,6] + z_p[6,2] >= 1
 RightLeft_6_2: x_p[2,6] + x_p[6,2] <= 1
 LeftRight_6_2: x_p[2,6] + x_p[6,2] <= 1
 AboveBelow_6_2: z_p[2,6] + z_p[6,2] <= 1
 BelowAbove_6_2: z_p[2,6] + z_p[6,2] <= 1
 Overlap_6_3: x_p[3,6] + x_p[6,3] + z_p[3,6] + z_p[6,3] >= 1
 RightLeft_6_3: x_p[3,6] + x_p[6,3] <= 1
 LeftRight_6_3: x_p[3,6] + x_p[6,3] <= 1
 AboveBelow_6_3: z_p[3,6] + z_p[6,3] <= 1
 BelowAbove_6_3: z_p[3,6] + z_p[6,3] <= 1
 Overlap_6_4: x_p[4,6] + x_p[6,4] + z_p[4,6] + z_p[6,4] >= 1
 RightLeft_6_4: x_p[4,6] + x_p[6,4] <= 1
 LeftRight_6_4: x_p[4,6] + x_p[6,4] <= 1
 AboveBelow_6_4: z_p[4,6] + z_p[6,4] <= 1
 BelowAbove_6_4: z_p[4,6] + z_p[6,4] <= 1
 Overlap_6_5: x_p[5,6] + x_p[6,5] + z_p[5,6] + z_p[6,5] >= 1
 RightLeft_6_5: x_p[5,6] + x_p[6,5] <= 1
 LeftRight_6_5: x_p[5,6] + x_p[6,5] <= 1
 AboveBelow_6_5: z_p[5,6] + z_p[6,5] <= 1
 BelowAbove_6_5: z_p[5,6] + z_p[6,5] <= 1
 Overlap_6_7: x_p[6,7] + x_p[7,6] + z_p[6,7] + z_p[7,6] >= 1
 RightLeft_6_7: x_p[6,7] + x_p[7,6] <= 1
 LeftRight_6_7: x_p[6,7] + x_p[7,6] <= 1
 AboveBelow_6_7: z_p[6,7] + z_p[7,6] <= 1
 BelowAbove_6_7: z_p[6,7] + z_p[7,6] <= 1
 Overlap_7_0: x_p[0,7] + x_p[7,0] + z_p[0,7] + z_p[7,0] >= 1
 RightLeft_7_0: x_p[0,7] + x_p[7,0] <= 1
 LeftRight_7_0: x_p[0,7] + x_p[7,0] <= 1
 AboveBelow_7_0: z_p[0,7] + z_p[7,0] <= 1
 BelowAbove_7_0: z_p[0,7] + z_p[7,0] <= 1
 Overlap_7_1: x_p[1,7] + x_p[7,1] + z_p[1,7] + z_p[7,1] >= 1
 RightLeft_7_1: x_p[1,7] + x_p[7,1] <= 1
 LeftRight_7_1: x_p[1,7] + x_p[7,1] <= 1
 AboveBelow_7_1: z_p[1,7] + z_p[7,1] <= 1
 BelowAbove_7_1: z_p[1,7] + z_p[7,1] <= 1
 Overlap_7_2: x_p[2,7] + x_p[7,2] + z_p[2,7] + z_p[7,2] >= 1
 RightLeft_7_2: x_p[2,7] + x_p[7,2] <= 1
 LeftRight_7_2: x_p[2,7] + x_p[7,2] <= 1
 AboveBelow_7_2: z_p[2,7] + z_p[7,2] <= 1
 BelowAbove_7_2: z_p[2,7] + z_p[7,2] <= 1
 Overlap_7_3: x_p[3,7] + x_p[7,3] + z_p[3,7] + z_p[7,3] >= 1
 RightLeft_7_3: x_p[3,7] + x_p[7,3] <= 1
 LeftRight_7_3: x_p[3,7] + x_p[7,3] <= 1
 AboveBelow_7_3: z_p[3,7] + z_p[7,3] <= 1
 BelowAbove_7_3: z_p[3,7] + z_p[7,3] <= 1
 Overlap_7_4: x_p[4,7] + x_p[7,4] + z_p[4,7] + z_p[7,4] >= 1
 RightLeft_7_4: x_p[4,7] + x_p[7,4] <= 1
 LeftRight_7_4: x_p[4,7] + x_p[7,4] <= 1
 AboveBelow_7_4: z_p[4,7] + z_p[7,4] <= 1
 BelowAbove_7_4: z_p[4,7] + z_p[7,4] <= 1
 Overlap_7_5: x_p[5,7] + x_p[7,5] + z_p[5,7] + z_p[7,5] >= 1
 RightLeft_7_5: x_p[5,7] + x_p[7,5] <= 1
 LeftRight_7_5: x_p[5,7] + x_p[7,5] <= 1
 AboveBelow_7_5: z_p[5,7] + z_p[7,5] <= 1
 BelowAbove_7_5: z_p[5,7] + z_p[7,5] <= 1
 Overlap_7_6: x_p[6,7] + x_p[7,6] + z_p[6,7] + z_p[7,6] >= 1
 RightLeft_7_6: x_p[6,7] + x_p[7,6] <= 1
 LeftRight_7_6: x_p[6,7] + x_p[7,6] <= 1
 AboveBelow_7_6: z_p[6,7] + z_p[7,6] <= 1
 BelowAbove_7_6: z_p[6,7] + z_p[7,6] <= 1
 Constraint14_0_1: 300 x_p[0,1] - xi[0] + x_i_prime[1] <= 300
 Constraint14_0_2: 300 x_p[0,2] - xi[0] + x_i_prime[2] <= 300
 Constraint14_0_3: 300 x_p[0,3] - xi[0] + x_i_prime[3] <= 300
 Constraint14_0_4: 300 x_p[0,4] - xi[0] + x_i_prime[4] <= 300
 Constraint14_0_5: 300 x_p[0,5] - xi[0] + x_i_prime[5] <= 300
 Constraint14_0_6: 300 x_p[0,6] - xi[0] + x_i_prime[6] <= 300
 Constraint14_0_7: 300 x_p[0,7] - xi[0] + x_i_prime[7] <= 300
 Constraint14_1_0: 300 x_p[1,0] - xi[1] + x_i_prime[0] <= 300
 Constraint14_1_2: 300 x_p[1,2] - xi[1] + x_i_prime[2] <= 300
 Constraint14_1_3: 300 x_p[1,3] - xi[1] + x_i_prime[3] <= 300
 Constraint14_1_4: 300 x_p[1,4] - xi[1] + x_i_prime[4] <= 300
 Constraint14_1_5: 300 x_p[1,5] - xi[1] + x_i_prime[5] <= 300
 Constraint14_1_6: 300 x_p[1,6] - xi[1] + x_i_prime[6] <= 300
 Constraint14_1_7: 300 x_p[1,7] - xi[1] + x_i_prime[7] <= 300
 Constraint14_2_0: 300 x_p[2,0] - xi[2] + x_i_prime[0] <= 300
 Constraint14_2_1: 300 x_p[2,1] - xi[2] + x_i_prime[1] <= 300
 Constraint14_2_3: 300 x_p[2,3] - xi[2] + x_i_prime[3] <= 300
 Constraint14_2_4: 300 x_p[2,4] - xi[2] + x_i_prime[4] <= 300
 Constraint14_2_5: 300 x_p[2,5] - xi[2] + x_i_prime[5] <= 300
 Constraint14_2_6: 300 x_p[2,6] - xi[2] + x_i_prime[6] <= 300
 Constraint14_2_7: 300 x_p[2,7] - xi[2] + x_i_prime[7] <= 300
 Constraint14_3_0: 300 x_p[3,0] - xi[3] + x_i_prime[0] <= 300
 Constraint14_3_1: 300 x_p[3,1] - xi[3] + x_i_prime[1] <= 300
 Constraint14_3_2: 300 x_p[3,2] - xi[3] + x_i_prime[2] <= 300
 Constraint14_3_4: 300 x_p[3,4] - xi[3] + x_i_prime[4] <= 300
 Constraint14_3_5: 300 x_p[3,5] - xi[3] + x_i_prime[5] <= 300
 Constraint14_3_6: 300 x_p[3,6] - xi[3] + x_i_prime[6] <= 300
 Constraint14_3_7: 300 x_p[3,7] - xi[3] + x_i_prime[7] <= 300
 Constraint14_4_0: 300 x_p[4,0] - xi[4] + x_i_prime[0] <= 300
 Constraint14_4_1: 300 x_p[4,1] - xi[4] + x_i_prime[1] <= 300
 Constraint14_4_2: 300 x_p[4,2] - xi[4] + x_i_prime[2] <= 300
 Constraint14_4_3: 300 x_p[4,3] - xi[4] + x_i_prime[3] <= 300
 Constraint14_4_5: 300 x_p[4,5] - xi[4] + x_i_prime[5] <= 300
 Constraint14_4_6: 300 x_p[4,6] - xi[4] + x_i_prime[6] <= 300
 Constraint14_4_7: 300 x_p[4,7] - xi[4] + x_i_prime[7] <= 300
 Constraint14_5_0: 300 x_p[5,0] - xi[5] + x_i_prime[0] <= 300
 Constraint14_5_1: 300 x_p[5,1] - xi[5] + x_i_prime[1] <= 300
 Constraint14_5_2: 300 x_p[5,2] - xi[5] + x_i_prime[2] <= 300
 Constraint14_5_3: 300 x_p[5,3] - xi[5] + x_i_prime[3] <= 300
 Constraint14_5_4: 300 x_p[5,4] - xi[5] + x_i_prime[4] <= 300
 Constraint14_5_6: 300 x_p[5,6] - xi[5] + x_i_prime[6] <= 300
 Constraint14_5_7: 300 x_p[5,7] - xi[5] + x_i_prime[7] <= 300
 Constraint14_6_0: 300 x_p[6,0] - xi[6] + x_i_prime[0] <= 300
 Constraint14_6_1: 300 x_p[6,1] - xi[6] + x_i_prime[1] <= 300
 Constraint14_6_2: 300 x_p[6,2] - xi[6] + x_i_prime[2] <= 300
 Constraint14_6_3: 300 x_p[6,3] - xi[6] + x_i_prime[3] <= 300
 Constraint14_6_4: 300 x_p[6,4] - xi[6] + x_i_prime[4] <= 300
 Constraint14_6_5: 300 x_p[6,5] - xi[6] + x_i_prime[5] <= 300
 Constraint14_6_7: 300 x_p[6,7] - xi[6] + x_i_prime[7] <= 300
 Constraint14_7_0: 300 x_p[7,0] - xi[7] + x_i_prime[0] <= 300
 Constraint14_7_1: 300 x_p[7,1] - xi[7] + x_i_prime[1] <= 300
 Constraint14_7_2: 300 x_p[7,2] - xi[7] + x_i_prime[2] <= 300
 Constraint14_7_3: 300 x_p[7,3] - xi[7] + x_i_prime[3] <= 300
 Constraint14_7_4: 300 x_p[7,4] - xi[7] + x_i_prime[4] <= 300
 Constraint14_7_5: 300 x_p[7,5] - xi[7] + x_i_prime[5] <= 300
 Constraint14_7_6: 300 x_p[7,6] - xi[7] + x_i_prime[6] <= 300
 Constraint15_0_1: - 300 x_p[0,1] + xi[0] - x_i_prime[1] <= -1
 Constraint15_0_2: - 300 x_p[0,2] + xi[0] - x_i_prime[2] <= -1
 Constraint15_0_3: - 300 x_p[0,3] + xi[0] - x_i_prime[3] <= -1
 Constraint15_0_4: - 300 x_p[0,4] + xi[0] - x_i_prime[4] <= -1
 Constraint15_0_5: - 300 x_p[0,5] + xi[0] - x_i_prime[5] <= -1
 Constraint15_0_6: - 300 x_p[0,6] + xi[0] - x_i_prime[6] <= -1
 Constraint15_0_7: - 300 x_p[0,7] + xi[0] - x_i_prime[7] <= -1
 Constraint15_1_0: - 300 x_p[1,0] + xi[1] - x_i_prime[0] <= -1
 Constraint15_1_2: - 300 x_p[1,2] + xi[1] - x_i_prime[2] <= -1
 Constraint15_1_3: - 300 x_p[1,3] + xi[1] - x_i_prime[3] <= -1
 Constraint15_1_4: - 300 x_p[1,4] + xi[1] - x_i_prime[4] <= -1
 Constraint15_1_5: - 300 x_p[1,5] + xi[1] - x_i_prime[5] <= -1
 Constraint15_1_6: - 300 x_p[1,6] + xi[1] - x_i_prime[6] <= -1
 Constraint15_1_7: - 300 x_p[1,7] + xi[1] - x_i_prime[7] <= -1
 Constraint15_2_0: - 300 x_p[2,0] + xi[2] - x_i_prime[0] <= -1
 Constraint15_2_1: - 300 x_p[2,1] + xi[2] - x_i_prime[1] <= -1
 Constraint15_2_3: - 300 x_p[2,3] + xi[2] - x_i_prime[3] <= -1
 Constraint15_2_4: - 300 x_p[2,4] + xi[2] - x_i_prime[4] <= -1
 Constraint15_2_5: - 300 x_p[2,5] + xi[2] - x_i_prime[5] <= -1
 Constraint15_2_6: - 300 x_p[2,6] + xi[2] - x_i_prime[6] <= -1
 Constraint15_2_7: - 300 x_p[2,7] + xi[2] - x_i_prime[7] <= -1
 Constraint15_3_0: - 300 x_p[3,0] + xi[3] - x_i_prime[0] <= -1
 Constraint15_3_1: - 300 x_p[3,1] + xi[3] - x_i_prime[1] <= -1
 Constraint15_3_2: - 300 x_p[3,2] + xi[3] - x_i_prime[2] <= -1
 Constraint15_3_4: - 300 x_p[3,4] + xi[3] - x_i_prime[4] <= -1
 Constraint15_3_5: - 300 x_p[3,5] + xi[3] - x_i_prime[5] <= -1
 Constraint15_3_6: - 300 x_p[3,6] + xi[3] - x_i_prime[6] <= -1
 Constraint15_3_7: - 300 x_p[3,7] + xi[3] - x_i_prime[7] <= -1
 Constraint15_4_0: - 300 x_p[4,0] + xi[4] - x_i_prime[0] <= -1
 Constraint15_4_1: - 300 x_p[4,1] + xi[4] - x_i_prime[1] <= -1
 Constraint15_4_2: - 300 x_p[4,2] + xi[4] - x_i_prime[2] <= -1
 Constraint15_4_3: - 300 x_p[4,3] + xi[4] - x_i_prime[3] <= -1
 Constraint15_4_5: - 300 x_p[4,5] + xi[4] - x_i_prime[5] <= -1
 Constraint15_4_6: - 300 x_p[4,6] + xi[4] - x_i_prime[6] <= -1
 Constraint15_4_7: - 300 x_p[4,7] + xi[4] - x_i_prime[7] <= -1
 Constraint15_5_0: - 300 x_p[5,0] + xi[5] - x_i_prime[0] <= -1
 Constraint15_5_1: - 300 x_p[5,1] + xi[5] - x_i_prime[1] <= -1
 Constraint15_5_2: - 300 x_p[5,2] + xi[5] - x_i_prime[2] <= -1
 Constraint15_5_3: - 300 x_p[5,3] + xi[5] - x_i_prime[3] <= -1
 Constraint15_5_4: - 300 x_p[5,4] + xi[5] - x_i_prime[4] <= -1
 Constraint15_5_6: - 300 x_p[5,6] + xi[5] - x_i_prime[6] <= -1
 Constraint15_5_7: - 300 x_p[5,7] + xi[5] - x_i_prime[7] <= -1
 Constraint15_6_0: - 300 x_p[6,0] + xi[6] - x_i_prime[0] <= -1
 Constraint15_6_1: - 300 x_p[6,1] + xi[6] - x_i_prime[1] <= -1
 Constraint15_6_2: - 300 x_p[6,2] + xi[6] - x_i_prime[2] <= -1
 Constraint15_6_3: - 300 x_p[6,3] + xi[6] - x_i_prime[3] <= -1
 Constraint15_6_4: - 300 x_p[6,4] + xi[6] - x_i_prime[4] <= -1
 Constraint15_6_5: - 300 x_p[6,5] + xi[6] - x_i_prime[5] <= -1
 Constraint15_6_7: - 300 x_p[6,7] + xi[6] - x_i_prime[7] <= -1
 Constraint15_7_0: - 300 x_p[7,0] + xi[7] - x_i_prime[0] <= -1
 Constraint15_7_1: - 300 x_p[7,1] + xi[7] - x_i_prime[1] <= -1
 Constraint15_7_2: - 300 x_p[7,2] + xi[7] - x_i_prime[2] <= -1
 Constraint15_7_3: - 300 x_p[7,3] + xi[7] - x_i_prime[3] <= -1
 Constraint15_7_4: - 300 x_p[7,4] + xi[7] - x_i_prime[4] <= -1
 Constraint15_7_5: - 300 x_p[7,5] + xi[7] - x_i_prime[5] <= -1
 Constraint15_7_6: - 300 x_p[7,6] + xi[7] - x_i_prime[6] <= -1
 Constraint16_0_1: 155 z_p[0,1] - zi[0] + z_i_prime[1] <= 155
 Constraint16_0_2: 155 z_p[0,2] - zi[0] + z_i_prime[2] <= 155
 Constraint16_0_3: 155 z_p[0,3] - zi[0] + z_i_prime[3] <= 155
 Constraint16_0_4: 155 z_p[0,4] - zi[0] + z_i_prime[4] <= 155
 Constraint16_0_5: 155 z_p[0,5] - zi[0] + z_i_prime[5] <= 155
 Constraint16_0_6: 155 z_p[0,6] - zi[0] + z_i_prime[6] <= 155
 Constraint16_0_7: 155 z_p[0,7] - zi[0] + z_i_prime[7] <= 155
 Constraint16_1_0: 155 z_p[1,0] - zi[1] + z_i_prime[0] <= 155
 Constraint16_1_2: 155 z_p[1,2] - zi[1] + z_i_prime[2] <= 155
 Constraint16_1_3: 155 z_p[1,3] - zi[1] + z_i_prime[3] <= 155
 Constraint16_1_4: 155 z_p[1,4] - zi[1] + z_i_prime[4] <= 155
 Constraint16_1_5: 155 z_p[1,5] - zi[1] + z_i_prime[5] <= 155
 Constraint16_1_6: 155 z_p[1,6] - zi[1] + z_i_prime[6] <= 155
 Constraint16_1_7: 155 z_p[1,7] - zi[1] + z_i_prime[7] <= 155
 Constraint16_2_0: 155 z_p[2,0] - zi[2] + z_i_prime[0] <= 155
 Constraint16_2_1: 155 z_p[2,1] - zi[2] + z_i_prime[1] <= 155
 Constraint16_2_3: 155 z_p[2,3] - zi[2] + z_i_prime[3] <= 155
 Constraint16_2_4: 155 z_p[2,4] - zi[2] + z_i_prime[4] <= 155
 Constraint16_2_5: 155 z_p[2,5] - zi[2] + z_i_prime[5] <= 155
 Constraint16_2_6: 155 z_p[2,6] - zi[2] + z_i_prime[6] <= 155
 Constraint16_2_7: 155 z_p[2,7] - zi[2] + z_i_prime[7] <= 155
 Constraint16_3_0: 155 z_p[3,0] - zi[3] + z_i_prime[0] <= 155
 Constraint16_3_1: 155 z_p[3,1] - zi[3] + z_i_prime[1] <= 155
 Constraint16_3_2: 155 z_p[3,2] - zi[3] + z_i_prime[2] <= 155
 Constraint16_3_4: 155 z_p[3,4] - zi[3] + z_i_prime[4] <= 155
 Constraint16_3_5: 155 z_p[3,5] - zi[3] + z_i_prime[5] <= 155
 Constraint16_3_6: 155 z_p[3,6] - zi[3] + z_i_prime[6] <= 155
 Constraint16_3_7: 155 z_p[3,7] - zi[3] + z_i_prime[7] <= 155
 Constraint16_4_0: 155 z_p[4,0] - zi[4] + z_i_prime[0] <= 155
 Constraint16_4_1: 155 z_p[4,1] - zi[4] + z_i_prime[1] <= 155
 Constraint16_4_2: 155 z_p[4,2] - zi[4] + z_i_prime[2] <= 155
 Constraint16_4_3: 155 z_p[4,3] - zi[4] + z_i_prime[3] <= 155
 Constraint16_4_5: 155 z_p[4,5] - zi[4] + z_i_prime[5] <= 155
 Constraint16_4_6: 155 z_p[4,6] - zi[4] + z_i_prime[6] <= 155
 Constraint16_4_7: 155 z_p[4,7] - zi[4] + z_i_prime[7] <= 155
 Constraint16_5_0: 155 z_p[5,0] - zi[5] + z_i_prime[0] <= 155
 Constraint16_5_1: 155 z_p[5,1] - zi[5] + z_i_prime[1] <= 155
 Constraint16_5_2: 155 z_p[5,2] - zi[5] + z_i_prime[2] <= 155
 Constraint16_5_3: 155 z_p[5,3] - zi[5] + z_i_prime[3] <= 155
 Constraint16_5_4: 155 z_p[5,4] - zi[5] + z_i_prime[4] <= 155
 Constraint16_5_6: 155 z_p[5,6] - zi[5] + z_i_prime[6] <= 155
 Constraint16_5_7: 155 z_p[5,7] - zi[5] + z_i_prime[7] <= 155
 Constraint16_6_0: 155 z_p[6,0] - zi[6] + z_i_prime[0] <= 155
 Constraint16_6_1: 155 z_p[6,1] - zi[6] + z_i_prime[1] <= 155
 Constraint16_6_2: 155 z_p[6,2] - zi[6] + z_i_prime[2] <= 155
 Constraint16_6_3: 155 z_p[6,3] - zi[6] + z_i_prime[3] <= 155
 Constraint16_6_4: 155 z_p[6,4] - zi[6] + z_i_prime[4] <= 155
 Constraint16_6_5: 155 z_p[6,5] - zi[6] + z_i_prime[5] <= 155
 Constraint16_6_7: 155 z_p[6,7] - zi[6] + z_i_prime[7] <= 155
 Constraint16_7_0: 155 z_p[7,0] - zi[7] + z_i_prime[0] <= 155
 Constraint16_7_1: 155 z_p[7,1] - zi[7] + z_i_prime[1] <= 155
 Constraint16_7_2: 155 z_p[7,2] - zi[7] + z_i_prime[2] <= 155
 Constraint16_7_3: 155 z_p[7,3] - zi[7] + z_i_prime[3] <= 155
 Constraint16_7_4: 155 z_p[7,4] - zi[7] + z_i_prime[4] <= 155
 Constraint16_7_5: 155 z_p[7,5] - zi[7] + z_i_prime[5] <= 155
 Constraint16_7_6: 155 z_p[7,6] - zi[7] + z_i_prime[6] <= 155
 OrientationLength_0: r[0,0,1] <= 1
 OrientationHeight_0: r[0,1,0] <= 1
 OrientationLength_1: r[1,0,1] <= 0
 OrientationHeight_1: r[1,1,0] <= 0
 OrientationLength_2: r[2,0,1] <= 1
 OrientationHeight_2: r[2,1,0] <= 1
 OrientationLength_3: r[3,0,1] <= 1
 OrientationHeight_3: r[3,1,0] <= 1
 OrientationLength_4: r[4,0,1] <= 1
 OrientationHeight_4: r[4,1,0] <= 1
 OrientationLength_5: r[5,0,1] <= 1
 OrientationHeight_5: r[5,1,0] <= 1
 OrientationLength_6: r[6,0,1] <= 1
 OrientationHeight_6: r[6,1,0] <= 1
 OrientationLength_7: r[7,0,1] <= 0
 OrientationHeight_7: r[7,1,0] <= 0
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
 Flagging_5_0: - p_ij[5,0] + u_j[0] >= 0
 Flagging_5_1: - p_ij[5,1] + u_j[1] >= 0
 Flagging_5_2: - p_ij[5,2] + u_j[2] >= 0
 Flagging_5_3: - p_ij[5,3] + u_j[3] >= 0
 Flagging_6_0: - p_ij[6,0] + u_j[0] >= 0
 Flagging_6_1: - p_ij[6,1] + u_j[1] >= 0
 Flagging_6_2: - p_ij[6,2] + u_j[2] >= 0
 Flagging_6_3: - p_ij[6,3] + u_j[3] >= 0
 Flagging_7_0: - p_ij[7,0] + u_j[0] >= 0
 Flagging_7_1: - p_ij[7,1] + u_j[1] >= 0
 Flagging_7_2: - p_ij[7,2] + u_j[2] >= 0
 Flagging_7_3: - p_ij[7,3] + u_j[3] >= 0
 Stability_0_3: 2 g[0] + beta1[0,1] + beta1[0,2] + beta1[0,3] + beta1[0,4]
   + beta1[0,5] + beta1[0,6] + beta1[0,7] + beta2[0,1] + beta2[0,2]
   + beta2[0,3] + beta2[0,4] + beta2[0,5] + beta2[0,6] + beta2[0,7]
   + gamma[0] >= 2
 Stability_1_3: 2 g[1] + beta1[1,0] + beta1[1,2] + beta1[1,3] + beta1[1,4]
   + beta1[1,5] + beta1[1,6] + beta1[1,7] + beta2[1,0] + beta2[1,2]
   + beta2[1,3] + beta2[1,4] + beta2[1,5] + beta2[1,6] + beta2[1,7]
   + gamma[1] >= 2
 Stability_2_3: 2 g[2] + beta1[2,0] + beta1[2,1] + beta1[2,3] + beta1[2,4]
   + beta1[2,5] + beta1[2,6] + beta1[2,7] + beta2[2,0] + beta2[2,1]
   + beta2[2,3] + beta2[2,4] + beta2[2,5] + beta2[2,6] + beta2[2,7]
   + gamma[2] >= 2
 Stability_3_3: 2 g[3] + beta1[3,0] + beta1[3,1] + beta1[3,2] + beta1[3,4]
   + beta1[3,5] + beta1[3,6] + beta1[3,7] + beta2[3,0] + beta2[3,1]
   + beta2[3,2] + beta2[3,4] + beta2[3,5] + beta2[3,6] + beta2[3,7]
   + gamma[3] >= 2
 Stability_4_3: 2 g[4] + beta1[4,0] + beta1[4,1] + beta1[4,2] + beta1[4,3]
   + beta1[4,5] + beta1[4,6] + beta1[4,7] + beta2[4,0] + beta2[4,1]
   + beta2[4,2] + beta2[4,3] + beta2[4,5] + beta2[4,6] + beta2[4,7]
   + gamma[4] >= 2
 Stability_5_3: 2 g[5] + beta1[5,0] + beta1[5,1] + beta1[5,2] + beta1[5,3]
   + beta1[5,4] + beta1[5,6] + beta1[5,7] + beta2[5,0] + beta2[5,1]
   + beta2[5,2] + beta2[5,3] + beta2[5,4] + beta2[5,6] + beta2[5,7]
   + gamma[5] >= 2
 Stability_6_3: 2 g[6] + beta1[6,0] + beta1[6,1] + beta1[6,2] + beta1[6,3]
   + beta1[6,4] + beta1[6,5] + beta1[6,7] + beta2[6,0] + beta2[6,1]
   + beta2[6,2] + beta2[6,3] + beta2[6,4] + beta2[6,5] + beta2[6,7]
   + gamma[6] >= 2
 Stability_7_3: 2 g[7] + beta1[7,0] + beta1[7,1] + beta1[7,2] + beta1[7,3]
   + beta1[7,4] + beta1[7,5] + beta1[7,6] + beta2[7,0] + beta2[7,1]
   + beta2[7,2] + beta2[7,3] + beta2[7,4] + beta2[7,5] + beta2[7,6]
   + gamma[7] >= 2
Bounds
Binaries
 p_ij[0,0] p_ij[0,1] p_ij[0,2] p_ij[0,3] p_ij[1,0] p_ij[1,1] p_ij[1,2]
 p_ij[1,3] p_ij[2,0] p_ij[2,1] p_ij[2,2] p_ij[2,3] p_ij[3,0] p_ij[3,1]
 p_ij[3,2] p_ij[3,3] p_ij[4,0] p_ij[4,1] p_ij[4,2] p_ij[4,3] p_ij[5,0]
 p_ij[5,1] p_ij[5,2] p_ij[5,3] p_ij[6,0] p_ij[6,1] p_ij[6,2] p_ij[6,3]
 p_ij[7,0] p_ij[7,1] p_ij[7,2] p_ij[7,3] u_j[0] u_j[1] u_j[2] u_j[3]
 x_p[0,0] x_p[0,1] x_p[0,2] x_p[0,3] x_p[0,4] x_p[0,5] x_p[0,6] x_p[0,7]
 x_p[1,0] x_p[1,1] x_p[1,2] x_p[1,3] x_p[1,4] x_p[1,5] x_p[1,6] x_p[1,7]
 x_p[2,0] x_p[2,1] x_p[2,2] x_p[2,3] x_p[2,4] x_p[2,5] x_p[2,6] x_p[2,7]
 x_p[3,0] x_p[3,1] x_p[3,2] x_p[3,3] x_p[3,4] x_p[3,5] x_p[3,6] x_p[3,7]
 x_p[4,0] x_p[4,1] x_p[4,2] x_p[4,3] x_p[4,4] x_p[4,5] x_p[4,6] x_p[4,7]
 x_p[5,0] x_p[5,1] x_p[5,2] x_p[5,3] x_p[5,4] x_p[5,5] x_p[5,6] x_p[5,7]
 x_p[6,0] x_p[6,1] x_p[6,2] x_p[6,3] x_p[6,4] x_p[6,5] x_p[6,6] x_p[6,7]
 x_p[7,0] x_p[7,1] x_p[7,2] x_p[7,3] x_p[7,4] x_p[7,5] x_p[7,6] x_p[7,7]
 z_p[0,0] z_p[0,1] z_p[0,2] z_p[0,3] z_p[0,4] z_p[0,5] z_p[0,6] z_p[0,7]
 z_p[1,0] z_p[1,1] z_p[1,2] z_p[1,3] z_p[1,4] z_p[1,5] z_p[1,6] z_p[1,7]
 z_p[2,0] z_p[2,1] z_p[2,2] z_p[2,3] z_p[2,4] z_p[2,5] z_p[2,6] z_p[2,7]
 z_p[3,0] z_p[3,1] z_p[3,2] z_p[3,3] z_p[3,4] z_p[3,5] z_p[3,6] z_p[3,7]
 z_p[4,0] z_p[4,1] z_p[4,2] z_p[4,3] z_p[4,4] z_p[4,5] z_p[4,6] z_p[4,7]
 z_p[5,0] z_p[5,1] z_p[5,2] z_p[5,3] z_p[5,4] z_p[5,5] z_p[5,6] z_p[5,7]
 z_p[6,0] z_p[6,1] z_p[6,2] z_p[6,3] z_p[6,4] z_p[6,5] z_p[6,6] z_p[6,7]
 z_p[7,0] z_p[7,1] z_p[7,2] z_p[7,3] z_p[7,4] z_p[7,5] z_p[7,6] z_p[7,7]
 r[0,0,0] r[0,0,1] r[0,1,0] r[0,1,1] r[1,0,0] r[1,0,1] r[1,1,0] r[1,1,1]
 r[2,0,0] r[2,0,1] r[2,1,0] r[2,1,1] r[3,0,0] r[3,0,1] r[3,1,0] r[3,1,1]
 r[4,0,0] r[4,0,1] r[4,1,0] r[4,1,1] r[5,0,0] r[5,0,1] r[5,1,0] r[5,1,1]
 r[6,0,0] r[6,0,1] r[6,1,0] r[6,1,1] r[7,0,0] r[7,0,1] r[7,1,0] r[7,1,1]
 rho[0] rho[1] rho[2] rho[3] rho[4] rho[5] rho[6] rho[7] g[0] g[1] g[2]
 g[3] g[4] g[5] g[6] g[7] beta1[0,0] beta1[0,1] beta1[0,2] beta1[0,3]
 beta1[0,4] beta1[0,5] beta1[0,6] beta1[0,7] beta1[1,0] beta1[1,1]
 beta1[1,2] beta1[1,3] beta1[1,4] beta1[1,5] beta1[1,6] beta1[1,7]
 beta1[2,0] beta1[2,1] beta1[2,2] beta1[2,3] beta1[2,4] beta1[2,5]
 beta1[2,6] beta1[2,7] beta1[3,0] beta1[3,1] beta1[3,2] beta1[3,3]
 beta1[3,4] beta1[3,5] beta1[3,6] beta1[3,7] beta1[4,0] beta1[4,1]
 beta1[4,2] beta1[4,3] beta1[4,4] beta1[4,5] beta1[4,6] beta1[4,7]
 beta1[5,0] beta1[5,1] beta1[5,2] beta1[5,3] beta1[5,4] beta1[5,5]
 beta1[5,6] beta1[5,7] beta1[6,0] beta1[6,1] beta1[6,2] beta1[6,3]
 beta1[6,4] beta1[6,5] beta1[6,6] beta1[6,7] beta1[7,0] beta1[7,1]
 beta1[7,2] beta1[7,3] beta1[7,4] beta1[7,5] beta1[7,6] beta1[7,7]
 beta2[0,0] beta2[0,1] beta2[0,2] beta2[0,3] beta2[0,4] beta2[0,5]
 beta2[0,6] beta2[0,7] beta2[1,0] beta2[1,1] beta2[1,2] beta2[1,3]
 beta2[1,4] beta2[1,5] beta2[1,6] beta2[1,7] beta2[2,0] beta2[2,1]
 beta2[2,2] beta2[2,3] beta2[2,4] beta2[2,5] beta2[2,6] beta2[2,7]
 beta2[3,0] beta2[3,1] beta2[3,2] beta2[3,3] beta2[3,4] beta2[3,5]
 beta2[3,6] beta2[3,7] beta2[4,0] beta2[4,1] beta2[4,2] beta2[4,3]
 beta2[4,4] beta2[4,5] beta2[4,6] beta2[4,7] beta2[5,0] beta2[5,1]
 beta2[5,2] beta2[5,3] beta2[5,4] beta2[5,5] beta2[5,6] beta2[5,7]
 beta2[6,0] beta2[6,1] beta2[6,2] beta2[6,3] beta2[6,4] beta2[6,5]
 beta2[6,6] beta2[6,7] beta2[7,0] beta2[7,1] beta2[7,2] beta2[7,3]
 beta2[7,4] beta2[7,5] beta2[7,6] beta2[7,7] gamma[0] gamma[1] gamma[2]
 gamma[3] gamma[4] gamma[5] gamma[6] gamma[7] h_i_k[0,0] h_i_k[0,1]
 h_i_k[0,2] h_i_k[0,3] h_i_k[0,4] h_i_k[0,5] h_i_k[0,6] h_i_k[0,7]
 h_i_k[1,0] h_i_k[1,1] h_i_k[1,2] h_i_k[1,3] h_i_k[1,4] h_i_k[1,5]
 h_i_k[1,6] h_i_k[1,7] h_i_k[2,0] h_i_k[2,1] h_i_k[2,2] h_i_k[2,3]
 h_i_k[2,4] h_i_k[2,5] h_i_k[2,6] h_i_k[2,7] h_i_k[3,0] h_i_k[3,1]
 h_i_k[3,2] h_i_k[3,3] h_i_k[3,4] h_i_k[3,5] h_i_k[3,6] h_i_k[3,7]
 h_i_k[4,0] h_i_k[4,1] h_i_k[4,2] h_i_k[4,3] h_i_k[4,4] h_i_k[4,5]
 h_i_k[4,6] h_i_k[4,7] h_i_k[5,0] h_i_k[5,1] h_i_k[5,2] h_i_k[5,3]
 h_i_k[5,4] h_i_k[5,5] h_i_k[5,6] h_i_k[5,7] h_i_k[6,0] h_i_k[6,1]
 h_i_k[6,2] h_i_k[6,3] h_i_k[6,4] h_i_k[6,5] h_i_k[6,6] h_i_k[6,7]
 h_i_k[7,0] h_i_k[7,1] h_i_k[7,2] h_i_k[7,3] h_i_k[7,4] h_i_k[7,5]
 h_i_k[7,6] h_i_k[7,7] n1_i_k[0,0] n1_i_k[0,1] n1_i_k[0,2] n1_i_k[0,3]
 n1_i_k[0,4] n1_i_k[0,5] n1_i_k[0,6] n1_i_k[0,7] n1_i_k[1,0] n1_i_k[1,1]
 n1_i_k[1,2] n1_i_k[1,3] n1_i_k[1,4] n1_i_k[1,5] n1_i_k[1,6] n1_i_k[1,7]
 n1_i_k[2,0] n1_i_k[2,1] n1_i_k[2,2] n1_i_k[2,3] n1_i_k[2,4] n1_i_k[2,5]
 n1_i_k[2,6] n1_i_k[2,7] n1_i_k[3,0] n1_i_k[3,1] n1_i_k[3,2] n1_i_k[3,3]
 n1_i_k[3,4] n1_i_k[3,5] n1_i_k[3,6] n1_i_k[3,7] n1_i_k[4,0] n1_i_k[4,1]
 n1_i_k[4,2] n1_i_k[4,3] n1_i_k[4,4] n1_i_k[4,5] n1_i_k[4,6] n1_i_k[4,7]
 n1_i_k[5,0] n1_i_k[5,1] n1_i_k[5,2] n1_i_k[5,3] n1_i_k[5,4] n1_i_k[5,5]
 n1_i_k[5,6] n1_i_k[5,7] n1_i_k[6,0] n1_i_k[6,1] n1_i_k[6,2] n1_i_k[6,3]
 n1_i_k[6,4] n1_i_k[6,5] n1_i_k[6,6] n1_i_k[6,7] n1_i_k[7,0] n1_i_k[7,1]
 n1_i_k[7,2] n1_i_k[7,3] n1_i_k[7,4] n1_i_k[7,5] n1_i_k[7,6] n1_i_k[7,7]
 n2_i_k[0,0] n2_i_k[0,1] n2_i_k[0,2] n2_i_k[0,3] n2_i_k[0,4] n2_i_k[0,5]
 n2_i_k[0,6] n2_i_k[0,7] n2_i_k[1,0] n2_i_k[1,1] n2_i_k[1,2] n2_i_k[1,3]
 n2_i_k[1,4] n2_i_k[1,5] n2_i_k[1,6] n2_i_k[1,7] n2_i_k[2,0] n2_i_k[2,1]
 n2_i_k[2,2] n2_i_k[2,3] n2_i_k[2,4] n2_i_k[2,5] n2_i_k[2,6] n2_i_k[2,7]
 n2_i_k[3,0] n2_i_k[3,1] n2_i_k[3,2] n2_i_k[3,3] n2_i_k[3,4] n2_i_k[3,5]
 n2_i_k[3,6] n2_i_k[3,7] n2_i_k[4,0] n2_i_k[4,1] n2_i_k[4,2] n2_i_k[4,3]
 n2_i_k[4,4] n2_i_k[4,5] n2_i_k[4,6] n2_i_k[4,7] n2_i_k[5,0] n2_i_k[5,1]
 n2_i_k[5,2] n2_i_k[5,3] n2_i_k[5,4] n2_i_k[5,5] n2_i_k[5,6] n2_i_k[5,7]
 n2_i_k[6,0] n2_i_k[6,1] n2_i_k[6,2] n2_i_k[6,3] n2_i_k[6,4] n2_i_k[6,5]
 n2_i_k[6,6] n2_i_k[6,7] n2_i_k[7,0] n2_i_k[7,1] n2_i_k[7,2] n2_i_k[7,3]
 n2_i_k[7,4] n2_i_k[7,5] n2_i_k[7,6] n2_i_k[7,7]
End
