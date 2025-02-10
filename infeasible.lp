\ Model 2DBPP
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  200 u_j[0] + 200 u_j[1] + 150 u_j[2] + 150 u_j[3] + 0 x_p[0,0]
   + 0 x_p[1,1] + 0 x_p[2,2] + 0 x_p[3,3] + 0 x_p[4,4] + 0 x_p[5,5]
   + 0 x_p[6,6] + 0 x_p[7,7] + 0 x_p[8,8] + 0 x_p[9,9] + 0 x_p[10,10]
   + 0 x_p[11,11] + 0 x_p[12,12] + 0 z_p[0,0] + 0 z_p[1,1] + 0 z_p[2,2]
   + 0 z_p[3,3] + 0 z_p[4,4] + 0 z_p[5,5] + 0 z_p[6,6] + 0 z_p[7,7]
   + 0 z_p[8,8] + 0 z_p[9,9] + 0 z_p[10,10] + 0 z_p[11,11] + 0 z_p[12,12]
   + 0 h[0,0] + 0 h[1,1] + 0 h[2,2] + 0 h[3,3] + 0 h[4,4] + 0 h[5,5]
   + 0 h[6,6] + 0 h[7,7] + 0 h[8,8] + 0 h[9,9] + 0 h[10,10] + 0 h[11,11]
   + 0 h[12,12] + 0 o[0,0] + 0 o[1,1] + 0 o[2,2] + 0 o[3,3] + 0 o[4,4]
   + 0 o[5,5] + 0 o[6,6] + 0 o[7,7] + 0 o[8,8] + 0 o[9,9] + 0 o[10,10]
   + 0 o[11,11] + 0 o[12,12] + 0 eta1[0,0] + 0 eta1[1,1] + 0 eta1[2,2]
   + 0 eta1[3,3] + 0 eta1[4,4] + 0 eta1[5,5] + 0 eta1[6,6] + 0 eta1[7,7]
   + 0 eta1[8,8] + 0 eta1[9,9] + 0 eta1[10,10] + 0 eta1[11,11]
   + 0 eta1[12,12] + 0 eta2[0,0] + 0 eta2[1,1] + 0 eta2[2,2] + 0 eta2[3,3]
   + 0 eta2[4,4] + 0 eta2[5,5] + 0 eta2[6,6] + 0 eta2[7,7] + 0 eta2[8,8]
   + 0 eta2[9,9] + 0 eta2[10,10] + 0 eta2[11,11] + 0 eta2[12,12]
   + 0 beta1[0,0] + 0 beta1[1,1] + 0 beta1[2,2] + 0 beta1[3,3]
   + 0 beta1[4,4] + 0 beta1[5,5] + 0 beta1[6,6] + 0 beta1[7,7]
   + 0 beta1[8,8] + 0 beta1[9,9] + 0 beta1[10,10] + 0 beta1[11,11]
   + 0 beta1[12,12] + 0 beta2[0,0] + 0 beta2[1,1] + 0 beta2[2,2]
   + 0 beta2[3,3] + 0 beta2[4,4] + 0 beta2[5,5] + 0 beta2[6,6]
   + 0 beta2[7,7] + 0 beta2[8,8] + 0 beta2[9,9] + 0 beta2[10,10]
   + 0 beta2[11,11] + 0 beta2[12,12] + 0 k[0,0] + 0 k[1,1] + 0 k[2,2]
   + 0 k[3,3] + 0 k[4,4] + 0 k[5,5] + 0 k[6,6] + 0 k[7,7] + 0 k[8,8]
   + 0 k[9,9] + 0 k[10,10] + 0 k[11,11] + 0 k[12,12] + 0 m[0,0] + 0 m[1,1]
   + 0 m[2,2] + 0 m[3,3] + 0 m[4,4] + 0 m[5,5] + 0 m[6,6] + 0 m[7,7]
   + 0 m[8,8] + 0 m[9,9] + 0 m[10,10] + 0 m[11,11] + 0 m[12,12]
Subject To
 TotalBinCapacity_0: 1820 p_ij[0,0] + 2240 p_ij[1,0] + 1696 p_ij[2,0]
   + 3168 p_ij[3,0] + 2640 p_ij[4,0] + 2494 p_ij[5,0] + 2262 p_ij[6,0]
   + 2838 p_ij[7,0] + 2418 p_ij[8,0] + 1692 p_ij[9,0] + 2772 p_ij[10,0]
   + 3256 p_ij[11,0] + 2009 p_ij[12,0] - 46500 u_j[0] <= 0
 TotalBinCapacity_1: 1820 p_ij[0,1] + 2240 p_ij[1,1] + 1696 p_ij[2,1]
   + 3168 p_ij[3,1] + 2640 p_ij[4,1] + 2494 p_ij[5,1] + 2262 p_ij[6,1]
   + 2838 p_ij[7,1] + 2418 p_ij[8,1] + 1692 p_ij[9,1] + 2772 p_ij[10,1]
   + 3256 p_ij[11,1] + 2009 p_ij[12,1] - 46500 u_j[1] <= 0
 TotalBinCapacity_2: 1820 p_ij[0,2] + 2240 p_ij[1,2] + 1696 p_ij[2,2]
   + 3168 p_ij[3,2] + 2640 p_ij[4,2] + 2494 p_ij[5,2] + 2262 p_ij[6,2]
   + 2838 p_ij[7,2] + 2418 p_ij[8,2] + 1692 p_ij[9,2] + 2772 p_ij[10,2]
   + 3256 p_ij[11,2] + 2009 p_ij[12,2] - 29760 u_j[2] <= 0
 TotalBinCapacity_3: 1820 p_ij[0,3] + 2240 p_ij[1,3] + 1696 p_ij[2,3]
   + 3168 p_ij[3,3] + 2640 p_ij[4,3] + 2494 p_ij[5,3] + 2262 p_ij[6,3]
   + 2838 p_ij[7,3] + 2418 p_ij[8,3] + 1692 p_ij[9,3] + 2772 p_ij[10,3]
   + 3256 p_ij[11,3] + 2009 p_ij[12,3] - 29760 u_j[3] <= 0
 OneItemOneBin_0: p_ij[0,0] + p_ij[0,1] + p_ij[0,2] + p_ij[0,3] = 1
 OneItemOneBin_1: p_ij[1,0] + p_ij[1,1] + p_ij[1,2] + p_ij[1,3] = 1
 OneItemOneBin_2: p_ij[2,0] + p_ij[2,1] + p_ij[2,2] + p_ij[2,3] = 1
 OneItemOneBin_3: p_ij[3,0] + p_ij[3,1] + p_ij[3,2] + p_ij[3,3] = 1
 OneItemOneBin_4: p_ij[4,0] + p_ij[4,1] + p_ij[4,2] + p_ij[4,3] = 1
 OneItemOneBin_5: p_ij[5,0] + p_ij[5,1] + p_ij[5,2] + p_ij[5,3] = 1
 OneItemOneBin_6: p_ij[6,0] + p_ij[6,1] + p_ij[6,2] + p_ij[6,3] = 1
 OneItemOneBin_7: p_ij[7,0] + p_ij[7,1] + p_ij[7,2] + p_ij[7,3] = 1
 OneItemOneBin_8: p_ij[8,0] + p_ij[8,1] + p_ij[8,2] + p_ij[8,3] = 1
 OneItemOneBin_9: p_ij[9,0] + p_ij[9,1] + p_ij[9,2] + p_ij[9,3] = 1
 OneItemOneBin_10: p_ij[10,0] + p_ij[10,1] + p_ij[10,2] + p_ij[10,3] = 1
 OneItemOneBin_11: p_ij[11,0] + p_ij[11,1] + p_ij[11,2] + p_ij[11,3] = 1
 OneItemOneBin_12: p_ij[12,0] + p_ij[12,1] + p_ij[12,2] + p_ij[12,3] = 1
 ItemFitToBin_0: - 300 p_ij[0,0] - 300 p_ij[0,1] - 192 p_ij[0,2]
   - 192 p_ij[0,3] + xi_r[0] <= 0
 ItemFitToBin_1: - 300 p_ij[1,0] - 300 p_ij[1,1] - 192 p_ij[1,2]
   - 192 p_ij[1,3] + xi_r[1] <= 0
 ItemFitToBin_2: - 300 p_ij[2,0] - 300 p_ij[2,1] - 192 p_ij[2,2]
   - 192 p_ij[2,3] + xi_r[2] <= 0
 ItemFitToBin_3: - 300 p_ij[3,0] - 300 p_ij[3,1] - 192 p_ij[3,2]
   - 192 p_ij[3,3] + xi_r[3] <= 0
 ItemFitToBin_4: - 300 p_ij[4,0] - 300 p_ij[4,1] - 192 p_ij[4,2]
   - 192 p_ij[4,3] + xi_r[4] <= 0
 ItemFitToBin_5: - 300 p_ij[5,0] - 300 p_ij[5,1] - 192 p_ij[5,2]
   - 192 p_ij[5,3] + xi_r[5] <= 0
 ItemFitToBin_6: - 300 p_ij[6,0] - 300 p_ij[6,1] - 192 p_ij[6,2]
   - 192 p_ij[6,3] + xi_r[6] <= 0
 ItemFitToBin_7: - 300 p_ij[7,0] - 300 p_ij[7,1] - 192 p_ij[7,2]
   - 192 p_ij[7,3] + xi_r[7] <= 0
 ItemFitToBin_8: - 300 p_ij[8,0] - 300 p_ij[8,1] - 192 p_ij[8,2]
   - 192 p_ij[8,3] + xi_r[8] <= 0
 ItemFitToBin_9: - 300 p_ij[9,0] - 300 p_ij[9,1] - 192 p_ij[9,2]
   - 192 p_ij[9,3] + xi_r[9] <= 0
 ItemFitToBin_10: - 300 p_ij[10,0] - 300 p_ij[10,1] - 192 p_ij[10,2]
   - 192 p_ij[10,3] + xi_r[10] <= 0
 ItemFitToBin_11: - 300 p_ij[11,0] - 300 p_ij[11,1] - 192 p_ij[11,2]
   - 192 p_ij[11,3] + xi_r[11] <= 0
 ItemFitToBin_12: - 300 p_ij[12,0] - 300 p_ij[12,1] - 192 p_ij[12,2]
   - 192 p_ij[12,3] + xi_r[12] <= 0
 ItemFitToBinZ_0: - 155 p_ij[0,0] - 155 p_ij[0,1] - 155 p_ij[0,2]
   - 155 p_ij[0,3] + zi_hi[0] <= 0
 ItemFitToBinZ_1: - 155 p_ij[1,0] - 155 p_ij[1,1] - 155 p_ij[1,2]
   - 155 p_ij[1,3] + zi_hi[1] <= 0
 ItemFitToBinZ_2: - 155 p_ij[2,0] - 155 p_ij[2,1] - 155 p_ij[2,2]
   - 155 p_ij[2,3] + zi_hi[2] <= 0
 ItemFitToBinZ_3: - 155 p_ij[3,0] - 155 p_ij[3,1] - 155 p_ij[3,2]
   - 155 p_ij[3,3] + zi_hi[3] <= 0
 ItemFitToBinZ_4: - 155 p_ij[4,0] - 155 p_ij[4,1] - 155 p_ij[4,2]
   - 155 p_ij[4,3] + zi_hi[4] <= 0
 ItemFitToBinZ_5: - 155 p_ij[5,0] - 155 p_ij[5,1] - 155 p_ij[5,2]
   - 155 p_ij[5,3] + zi_hi[5] <= 0
 ItemFitToBinZ_6: - 155 p_ij[6,0] - 155 p_ij[6,1] - 155 p_ij[6,2]
   - 155 p_ij[6,3] + zi_hi[6] <= 0
 ItemFitToBinZ_7: - 155 p_ij[7,0] - 155 p_ij[7,1] - 155 p_ij[7,2]
   - 155 p_ij[7,3] + zi_hi[7] <= 0
 ItemFitToBinZ_8: - 155 p_ij[8,0] - 155 p_ij[8,1] - 155 p_ij[8,2]
   - 155 p_ij[8,3] + zi_hi[8] <= 0
 ItemFitToBinZ_9: - 155 p_ij[9,0] - 155 p_ij[9,1] - 155 p_ij[9,2]
   - 155 p_ij[9,3] + zi_hi[9] <= 0
 ItemFitToBinZ_10: - 155 p_ij[10,0] - 155 p_ij[10,1] - 155 p_ij[10,2]
   - 155 p_ij[10,3] + zi_hi[10] <= 0
 ItemFitToBinZ_11: - 155 p_ij[11,0] - 155 p_ij[11,1] - 155 p_ij[11,2]
   - 155 p_ij[11,3] + zi_hi[11] <= 0
 ItemFitToBinZ_12: - 155 p_ij[12,0] - 155 p_ij[12,1] - 155 p_ij[12,2]
   - 155 p_ij[12,3] + zi_hi[12] <= 0
 TransformX_0: - xi_l[0] + xi_r[0] - 65 r[0,0,0] - 28 r[0,0,1] = 0
 TransformZ_0: - zi[0] + zi_hi[0] - 65 r[0,1,0] - 28 r[0,1,1] = 0
 TransformX_1: - xi_l[1] + xi_r[1] - 64 r[1,0,0] - 35 r[1,0,1] = 0
 TransformZ_1: - zi[1] + zi_hi[1] - 64 r[1,1,0] - 35 r[1,1,1] = 0
 TransformX_2: - xi_l[2] + xi_r[2] - 53 r[2,0,0] - 32 r[2,0,1] = 0
 TransformZ_2: - zi[2] + zi_hi[2] - 53 r[2,1,0] - 32 r[2,1,1] = 0
 TransformX_3: - xi_l[3] + xi_r[3] - 88 r[3,0,0] - 36 r[3,0,1] = 0
 TransformZ_3: - zi[3] + zi_hi[3] - 88 r[3,1,0] - 36 r[3,1,1] = 0
 TransformX_4: - xi_l[4] + xi_r[4] - 88 r[4,0,0] - 30 r[4,0,1] = 0
 TransformZ_4: - zi[4] + zi_hi[4] - 88 r[4,1,0] - 30 r[4,1,1] = 0
 TransformX_5: - xi_l[5] + xi_r[5] - 86 r[5,0,0] - 29 r[5,0,1] = 0
 TransformZ_5: - zi[5] + zi_hi[5] - 86 r[5,1,0] - 29 r[5,1,1] = 0
 TransformX_6: - xi_l[6] + xi_r[6] - 78 r[6,0,0] - 29 r[6,0,1] = 0
 TransformZ_6: - zi[6] + zi_hi[6] - 78 r[6,1,0] - 29 r[6,1,1] = 0
 TransformX_7: - xi_l[7] + xi_r[7] - 66 r[7,0,0] - 43 r[7,0,1] = 0
 TransformZ_7: - zi[7] + zi_hi[7] - 66 r[7,1,0] - 43 r[7,1,1] = 0
 TransformX_8: - xi_l[8] + xi_r[8] - 78 r[8,0,0] - 31 r[8,0,1] = 0
 TransformZ_8: - zi[8] + zi_hi[8] - 78 r[8,1,0] - 31 r[8,1,1] = 0
 TransformX_9: - xi_l[9] + xi_r[9] - 47 r[9,0,0] - 36 r[9,0,1] = 0
 TransformZ_9: - zi[9] + zi_hi[9] - 47 r[9,1,0] - 36 r[9,1,1] = 0
 TransformX_10: - xi_l[10] + xi_r[10] - 77 r[10,0,0] - 36 r[10,0,1] = 0
 TransformZ_10: - zi[10] + zi_hi[10] - 77 r[10,1,0] - 36 r[10,1,1] = 0
 TransformX_11: - xi_l[11] + xi_r[11] - 74 r[11,0,0] - 44 r[11,0,1] = 0
 TransformZ_11: - zi[11] + zi_hi[11] - 74 r[11,1,0] - 44 r[11,1,1] = 0
 TransformX_12: - xi_l[12] + xi_r[12] - 49 r[12,0,0] - 41 r[12,0,1] = 0
 TransformZ_12: - zi[12] + zi_hi[12] - 49 r[12,1,0] - 41 r[12,1,1] = 0
 SidesAlign_0_0: r[0,0,0] + r[0,1,0] = 1
 SidesAlign_1_0: r[0,0,1] + r[0,1,1] = 1
 AxesAlign_0_0: r[0,0,0] + r[0,0,1] = 1
 AxesAlign_1_0: r[0,1,0] + r[0,1,1] = 1
 SidesAlign_0_1: r[1,0,0] + r[1,1,0] = 1
 SidesAlign_1_1: r[1,0,1] + r[1,1,1] = 1
 AxesAlign_0_1: r[1,0,0] + r[1,0,1] = 1
 AxesAlign_1_1: r[1,1,0] + r[1,1,1] = 1
 SidesAlign_0_2: r[2,0,0] + r[2,1,0] = 1
 SidesAlign_1_2: r[2,0,1] + r[2,1,1] = 1
 AxesAlign_0_2: r[2,0,0] + r[2,0,1] = 1
 AxesAlign_1_2: r[2,1,0] + r[2,1,1] = 1
 SidesAlign_0_3: r[3,0,0] + r[3,1,0] = 1
 SidesAlign_1_3: r[3,0,1] + r[3,1,1] = 1
 AxesAlign_0_3: r[3,0,0] + r[3,0,1] = 1
 AxesAlign_1_3: r[3,1,0] + r[3,1,1] = 1
 SidesAlign_0_4: r[4,0,0] + r[4,1,0] = 1
 SidesAlign_1_4: r[4,0,1] + r[4,1,1] = 1
 AxesAlign_0_4: r[4,0,0] + r[4,0,1] = 1
 AxesAlign_1_4: r[4,1,0] + r[4,1,1] = 1
 SidesAlign_0_5: r[5,0,0] + r[5,1,0] = 1
 SidesAlign_1_5: r[5,0,1] + r[5,1,1] = 1
 AxesAlign_0_5: r[5,0,0] + r[5,0,1] = 1
 AxesAlign_1_5: r[5,1,0] + r[5,1,1] = 1
 SidesAlign_0_6: r[6,0,0] + r[6,1,0] = 1
 SidesAlign_1_6: r[6,0,1] + r[6,1,1] = 1
 AxesAlign_0_6: r[6,0,0] + r[6,0,1] = 1
 AxesAlign_1_6: r[6,1,0] + r[6,1,1] = 1
 SidesAlign_0_7: r[7,0,0] + r[7,1,0] = 1
 SidesAlign_1_7: r[7,0,1] + r[7,1,1] = 1
 AxesAlign_0_7: r[7,0,0] + r[7,0,1] = 1
 AxesAlign_1_7: r[7,1,0] + r[7,1,1] = 1
 SidesAlign_0_8: r[8,0,0] + r[8,1,0] = 1
 SidesAlign_1_8: r[8,0,1] + r[8,1,1] = 1
 AxesAlign_0_8: r[8,0,0] + r[8,0,1] = 1
 AxesAlign_1_8: r[8,1,0] + r[8,1,1] = 1
 SidesAlign_0_9: r[9,0,0] + r[9,1,0] = 1
 SidesAlign_1_9: r[9,0,1] + r[9,1,1] = 1
 AxesAlign_0_9: r[9,0,0] + r[9,0,1] = 1
 AxesAlign_1_9: r[9,1,0] + r[9,1,1] = 1
 SidesAlign_0_10: r[10,0,0] + r[10,1,0] = 1
 SidesAlign_1_10: r[10,0,1] + r[10,1,1] = 1
 AxesAlign_0_10: r[10,0,0] + r[10,0,1] = 1
 AxesAlign_1_10: r[10,1,0] + r[10,1,1] = 1
 SidesAlign_0_11: r[11,0,0] + r[11,1,0] = 1
 SidesAlign_1_11: r[11,0,1] + r[11,1,1] = 1
 AxesAlign_0_11: r[11,0,0] + r[11,0,1] = 1
 AxesAlign_1_11: r[11,1,0] + r[11,1,1] = 1
 SidesAlign_0_12: r[12,0,0] + r[12,1,0] = 1
 SidesAlign_1_12: r[12,0,1] + r[12,1,1] = 1
 AxesAlign_0_12: r[12,0,0] + r[12,0,1] = 1
 AxesAlign_1_12: r[12,1,0] + r[12,1,1] = 1
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
 OverlapInBin_0_8_toBin_0: - p_ij[0,0] - p_ij[8,0] + x_p[0,8] + x_p[8,0]
   + z_p[0,8] + z_p[8,0] >= -1
 OverlapInBin_0_8_toBin_1: - p_ij[0,1] - p_ij[8,1] + x_p[0,8] + x_p[8,0]
   + z_p[0,8] + z_p[8,0] >= -1
 OverlapInBin_0_8_toBin_2: - p_ij[0,2] - p_ij[8,2] + x_p[0,8] + x_p[8,0]
   + z_p[0,8] + z_p[8,0] >= -1
 OverlapInBin_0_8_toBin_3: - p_ij[0,3] - p_ij[8,3] + x_p[0,8] + x_p[8,0]
   + z_p[0,8] + z_p[8,0] >= -1
 OverlapInBin_0_9_toBin_0: - p_ij[0,0] - p_ij[9,0] + x_p[0,9] + x_p[9,0]
   + z_p[0,9] + z_p[9,0] >= -1
 OverlapInBin_0_9_toBin_1: - p_ij[0,1] - p_ij[9,1] + x_p[0,9] + x_p[9,0]
   + z_p[0,9] + z_p[9,0] >= -1
 OverlapInBin_0_9_toBin_2: - p_ij[0,2] - p_ij[9,2] + x_p[0,9] + x_p[9,0]
   + z_p[0,9] + z_p[9,0] >= -1
 OverlapInBin_0_9_toBin_3: - p_ij[0,3] - p_ij[9,3] + x_p[0,9] + x_p[9,0]
   + z_p[0,9] + z_p[9,0] >= -1
 OverlapInBin_0_10_toBin_0: - p_ij[0,0] - p_ij[10,0] + x_p[0,10]
   + x_p[10,0] + z_p[0,10] + z_p[10,0] >= -1
 OverlapInBin_0_10_toBin_1: - p_ij[0,1] - p_ij[10,1] + x_p[0,10]
   + x_p[10,0] + z_p[0,10] + z_p[10,0] >= -1
 OverlapInBin_0_10_toBin_2: - p_ij[0,2] - p_ij[10,2] + x_p[0,10]
   + x_p[10,0] + z_p[0,10] + z_p[10,0] >= -1
 OverlapInBin_0_10_toBin_3: - p_ij[0,3] - p_ij[10,3] + x_p[0,10]
   + x_p[10,0] + z_p[0,10] + z_p[10,0] >= -1
 OverlapInBin_0_11_toBin_0: - p_ij[0,0] - p_ij[11,0] + x_p[0,11]
   + x_p[11,0] + z_p[0,11] + z_p[11,0] >= -1
 OverlapInBin_0_11_toBin_1: - p_ij[0,1] - p_ij[11,1] + x_p[0,11]
   + x_p[11,0] + z_p[0,11] + z_p[11,0] >= -1
 OverlapInBin_0_11_toBin_2: - p_ij[0,2] - p_ij[11,2] + x_p[0,11]
   + x_p[11,0] + z_p[0,11] + z_p[11,0] >= -1
 OverlapInBin_0_11_toBin_3: - p_ij[0,3] - p_ij[11,3] + x_p[0,11]
   + x_p[11,0] + z_p[0,11] + z_p[11,0] >= -1
 OverlapInBin_0_12_toBin_0: - p_ij[0,0] - p_ij[12,0] + x_p[0,12]
   + x_p[12,0] + z_p[0,12] + z_p[12,0] >= -1
 OverlapInBin_0_12_toBin_1: - p_ij[0,1] - p_ij[12,1] + x_p[0,12]
   + x_p[12,0] + z_p[0,12] + z_p[12,0] >= -1
 OverlapInBin_0_12_toBin_2: - p_ij[0,2] - p_ij[12,2] + x_p[0,12]
   + x_p[12,0] + z_p[0,12] + z_p[12,0] >= -1
 OverlapInBin_0_12_toBin_3: - p_ij[0,3] - p_ij[12,3] + x_p[0,12]
   + x_p[12,0] + z_p[0,12] + z_p[12,0] >= -1
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
 OverlapInBin_1_8_toBin_0: - p_ij[1,0] - p_ij[8,0] + x_p[1,8] + x_p[8,1]
   + z_p[1,8] + z_p[8,1] >= -1
 OverlapInBin_1_8_toBin_1: - p_ij[1,1] - p_ij[8,1] + x_p[1,8] + x_p[8,1]
   + z_p[1,8] + z_p[8,1] >= -1
 OverlapInBin_1_8_toBin_2: - p_ij[1,2] - p_ij[8,2] + x_p[1,8] + x_p[8,1]
   + z_p[1,8] + z_p[8,1] >= -1
 OverlapInBin_1_8_toBin_3: - p_ij[1,3] - p_ij[8,3] + x_p[1,8] + x_p[8,1]
   + z_p[1,8] + z_p[8,1] >= -1
 OverlapInBin_1_9_toBin_0: - p_ij[1,0] - p_ij[9,0] + x_p[1,9] + x_p[9,1]
   + z_p[1,9] + z_p[9,1] >= -1
 OverlapInBin_1_9_toBin_1: - p_ij[1,1] - p_ij[9,1] + x_p[1,9] + x_p[9,1]
   + z_p[1,9] + z_p[9,1] >= -1
 OverlapInBin_1_9_toBin_2: - p_ij[1,2] - p_ij[9,2] + x_p[1,9] + x_p[9,1]
   + z_p[1,9] + z_p[9,1] >= -1
 OverlapInBin_1_9_toBin_3: - p_ij[1,3] - p_ij[9,3] + x_p[1,9] + x_p[9,1]
   + z_p[1,9] + z_p[9,1] >= -1
 OverlapInBin_1_10_toBin_0: - p_ij[1,0] - p_ij[10,0] + x_p[1,10]
   + x_p[10,1] + z_p[1,10] + z_p[10,1] >= -1
 OverlapInBin_1_10_toBin_1: - p_ij[1,1] - p_ij[10,1] + x_p[1,10]
   + x_p[10,1] + z_p[1,10] + z_p[10,1] >= -1
 OverlapInBin_1_10_toBin_2: - p_ij[1,2] - p_ij[10,2] + x_p[1,10]
   + x_p[10,1] + z_p[1,10] + z_p[10,1] >= -1
 OverlapInBin_1_10_toBin_3: - p_ij[1,3] - p_ij[10,3] + x_p[1,10]
   + x_p[10,1] + z_p[1,10] + z_p[10,1] >= -1
 OverlapInBin_1_11_toBin_0: - p_ij[1,0] - p_ij[11,0] + x_p[1,11]
   + x_p[11,1] + z_p[1,11] + z_p[11,1] >= -1
 OverlapInBin_1_11_toBin_1: - p_ij[1,1] - p_ij[11,1] + x_p[1,11]
   + x_p[11,1] + z_p[1,11] + z_p[11,1] >= -1
 OverlapInBin_1_11_toBin_2: - p_ij[1,2] - p_ij[11,2] + x_p[1,11]
   + x_p[11,1] + z_p[1,11] + z_p[11,1] >= -1
 OverlapInBin_1_11_toBin_3: - p_ij[1,3] - p_ij[11,3] + x_p[1,11]
   + x_p[11,1] + z_p[1,11] + z_p[11,1] >= -1
 OverlapInBin_1_12_toBin_0: - p_ij[1,0] - p_ij[12,0] + x_p[1,12]
   + x_p[12,1] + z_p[1,12] + z_p[12,1] >= -1
 OverlapInBin_1_12_toBin_1: - p_ij[1,1] - p_ij[12,1] + x_p[1,12]
   + x_p[12,1] + z_p[1,12] + z_p[12,1] >= -1
 OverlapInBin_1_12_toBin_2: - p_ij[1,2] - p_ij[12,2] + x_p[1,12]
   + x_p[12,1] + z_p[1,12] + z_p[12,1] >= -1
 OverlapInBin_1_12_toBin_3: - p_ij[1,3] - p_ij[12,3] + x_p[1,12]
   + x_p[12,1] + z_p[1,12] + z_p[12,1] >= -1
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
 OverlapInBin_2_8_toBin_0: - p_ij[2,0] - p_ij[8,0] + x_p[2,8] + x_p[8,2]
   + z_p[2,8] + z_p[8,2] >= -1
 OverlapInBin_2_8_toBin_1: - p_ij[2,1] - p_ij[8,1] + x_p[2,8] + x_p[8,2]
   + z_p[2,8] + z_p[8,2] >= -1
 OverlapInBin_2_8_toBin_2: - p_ij[2,2] - p_ij[8,2] + x_p[2,8] + x_p[8,2]
   + z_p[2,8] + z_p[8,2] >= -1
 OverlapInBin_2_8_toBin_3: - p_ij[2,3] - p_ij[8,3] + x_p[2,8] + x_p[8,2]
   + z_p[2,8] + z_p[8,2] >= -1
 OverlapInBin_2_9_toBin_0: - p_ij[2,0] - p_ij[9,0] + x_p[2,9] + x_p[9,2]
   + z_p[2,9] + z_p[9,2] >= -1
 OverlapInBin_2_9_toBin_1: - p_ij[2,1] - p_ij[9,1] + x_p[2,9] + x_p[9,2]
   + z_p[2,9] + z_p[9,2] >= -1
 OverlapInBin_2_9_toBin_2: - p_ij[2,2] - p_ij[9,2] + x_p[2,9] + x_p[9,2]
   + z_p[2,9] + z_p[9,2] >= -1
 OverlapInBin_2_9_toBin_3: - p_ij[2,3] - p_ij[9,3] + x_p[2,9] + x_p[9,2]
   + z_p[2,9] + z_p[9,2] >= -1
 OverlapInBin_2_10_toBin_0: - p_ij[2,0] - p_ij[10,0] + x_p[2,10]
   + x_p[10,2] + z_p[2,10] + z_p[10,2] >= -1
 OverlapInBin_2_10_toBin_1: - p_ij[2,1] - p_ij[10,1] + x_p[2,10]
   + x_p[10,2] + z_p[2,10] + z_p[10,2] >= -1
 OverlapInBin_2_10_toBin_2: - p_ij[2,2] - p_ij[10,2] + x_p[2,10]
   + x_p[10,2] + z_p[2,10] + z_p[10,2] >= -1
 OverlapInBin_2_10_toBin_3: - p_ij[2,3] - p_ij[10,3] + x_p[2,10]
   + x_p[10,2] + z_p[2,10] + z_p[10,2] >= -1
 OverlapInBin_2_11_toBin_0: - p_ij[2,0] - p_ij[11,0] + x_p[2,11]
   + x_p[11,2] + z_p[2,11] + z_p[11,2] >= -1
 OverlapInBin_2_11_toBin_1: - p_ij[2,1] - p_ij[11,1] + x_p[2,11]
   + x_p[11,2] + z_p[2,11] + z_p[11,2] >= -1
 OverlapInBin_2_11_toBin_2: - p_ij[2,2] - p_ij[11,2] + x_p[2,11]
   + x_p[11,2] + z_p[2,11] + z_p[11,2] >= -1
 OverlapInBin_2_11_toBin_3: - p_ij[2,3] - p_ij[11,3] + x_p[2,11]
   + x_p[11,2] + z_p[2,11] + z_p[11,2] >= -1
 OverlapInBin_2_12_toBin_0: - p_ij[2,0] - p_ij[12,0] + x_p[2,12]
   + x_p[12,2] + z_p[2,12] + z_p[12,2] >= -1
 OverlapInBin_2_12_toBin_1: - p_ij[2,1] - p_ij[12,1] + x_p[2,12]
   + x_p[12,2] + z_p[2,12] + z_p[12,2] >= -1
 OverlapInBin_2_12_toBin_2: - p_ij[2,2] - p_ij[12,2] + x_p[2,12]
   + x_p[12,2] + z_p[2,12] + z_p[12,2] >= -1
 OverlapInBin_2_12_toBin_3: - p_ij[2,3] - p_ij[12,3] + x_p[2,12]
   + x_p[12,2] + z_p[2,12] + z_p[12,2] >= -1
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
 OverlapInBin_3_8_toBin_0: - p_ij[3,0] - p_ij[8,0] + x_p[3,8] + x_p[8,3]
   + z_p[3,8] + z_p[8,3] >= -1
 OverlapInBin_3_8_toBin_1: - p_ij[3,1] - p_ij[8,1] + x_p[3,8] + x_p[8,3]
   + z_p[3,8] + z_p[8,3] >= -1
 OverlapInBin_3_8_toBin_2: - p_ij[3,2] - p_ij[8,2] + x_p[3,8] + x_p[8,3]
   + z_p[3,8] + z_p[8,3] >= -1
 OverlapInBin_3_8_toBin_3: - p_ij[3,3] - p_ij[8,3] + x_p[3,8] + x_p[8,3]
   + z_p[3,8] + z_p[8,3] >= -1
 OverlapInBin_3_9_toBin_0: - p_ij[3,0] - p_ij[9,0] + x_p[3,9] + x_p[9,3]
   + z_p[3,9] + z_p[9,3] >= -1
 OverlapInBin_3_9_toBin_1: - p_ij[3,1] - p_ij[9,1] + x_p[3,9] + x_p[9,3]
   + z_p[3,9] + z_p[9,3] >= -1
 OverlapInBin_3_9_toBin_2: - p_ij[3,2] - p_ij[9,2] + x_p[3,9] + x_p[9,3]
   + z_p[3,9] + z_p[9,3] >= -1
 OverlapInBin_3_9_toBin_3: - p_ij[3,3] - p_ij[9,3] + x_p[3,9] + x_p[9,3]
   + z_p[3,9] + z_p[9,3] >= -1
 OverlapInBin_3_10_toBin_0: - p_ij[3,0] - p_ij[10,0] + x_p[3,10]
   + x_p[10,3] + z_p[3,10] + z_p[10,3] >= -1
 OverlapInBin_3_10_toBin_1: - p_ij[3,1] - p_ij[10,1] + x_p[3,10]
   + x_p[10,3] + z_p[3,10] + z_p[10,3] >= -1
 OverlapInBin_3_10_toBin_2: - p_ij[3,2] - p_ij[10,2] + x_p[3,10]
   + x_p[10,3] + z_p[3,10] + z_p[10,3] >= -1
 OverlapInBin_3_10_toBin_3: - p_ij[3,3] - p_ij[10,3] + x_p[3,10]
   + x_p[10,3] + z_p[3,10] + z_p[10,3] >= -1
 OverlapInBin_3_11_toBin_0: - p_ij[3,0] - p_ij[11,0] + x_p[3,11]
   + x_p[11,3] + z_p[3,11] + z_p[11,3] >= -1
 OverlapInBin_3_11_toBin_1: - p_ij[3,1] - p_ij[11,1] + x_p[3,11]
   + x_p[11,3] + z_p[3,11] + z_p[11,3] >= -1
 OverlapInBin_3_11_toBin_2: - p_ij[3,2] - p_ij[11,2] + x_p[3,11]
   + x_p[11,3] + z_p[3,11] + z_p[11,3] >= -1
 OverlapInBin_3_11_toBin_3: - p_ij[3,3] - p_ij[11,3] + x_p[3,11]
   + x_p[11,3] + z_p[3,11] + z_p[11,3] >= -1
 OverlapInBin_3_12_toBin_0: - p_ij[3,0] - p_ij[12,0] + x_p[3,12]
   + x_p[12,3] + z_p[3,12] + z_p[12,3] >= -1
 OverlapInBin_3_12_toBin_1: - p_ij[3,1] - p_ij[12,1] + x_p[3,12]
   + x_p[12,3] + z_p[3,12] + z_p[12,3] >= -1
 OverlapInBin_3_12_toBin_2: - p_ij[3,2] - p_ij[12,2] + x_p[3,12]
   + x_p[12,3] + z_p[3,12] + z_p[12,3] >= -1
 OverlapInBin_3_12_toBin_3: - p_ij[3,3] - p_ij[12,3] + x_p[3,12]
   + x_p[12,3] + z_p[3,12] + z_p[12,3] >= -1
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
 OverlapInBin_4_8_toBin_0: - p_ij[4,0] - p_ij[8,0] + x_p[4,8] + x_p[8,4]
   + z_p[4,8] + z_p[8,4] >= -1
 OverlapInBin_4_8_toBin_1: - p_ij[4,1] - p_ij[8,1] + x_p[4,8] + x_p[8,4]
   + z_p[4,8] + z_p[8,4] >= -1
 OverlapInBin_4_8_toBin_2: - p_ij[4,2] - p_ij[8,2] + x_p[4,8] + x_p[8,4]
   + z_p[4,8] + z_p[8,4] >= -1
 OverlapInBin_4_8_toBin_3: - p_ij[4,3] - p_ij[8,3] + x_p[4,8] + x_p[8,4]
   + z_p[4,8] + z_p[8,4] >= -1
 OverlapInBin_4_9_toBin_0: - p_ij[4,0] - p_ij[9,0] + x_p[4,9] + x_p[9,4]
   + z_p[4,9] + z_p[9,4] >= -1
 OverlapInBin_4_9_toBin_1: - p_ij[4,1] - p_ij[9,1] + x_p[4,9] + x_p[9,4]
   + z_p[4,9] + z_p[9,4] >= -1
 OverlapInBin_4_9_toBin_2: - p_ij[4,2] - p_ij[9,2] + x_p[4,9] + x_p[9,4]
   + z_p[4,9] + z_p[9,4] >= -1
 OverlapInBin_4_9_toBin_3: - p_ij[4,3] - p_ij[9,3] + x_p[4,9] + x_p[9,4]
   + z_p[4,9] + z_p[9,4] >= -1
 OverlapInBin_4_10_toBin_0: - p_ij[4,0] - p_ij[10,0] + x_p[4,10]
   + x_p[10,4] + z_p[4,10] + z_p[10,4] >= -1
 OverlapInBin_4_10_toBin_1: - p_ij[4,1] - p_ij[10,1] + x_p[4,10]
   + x_p[10,4] + z_p[4,10] + z_p[10,4] >= -1
 OverlapInBin_4_10_toBin_2: - p_ij[4,2] - p_ij[10,2] + x_p[4,10]
   + x_p[10,4] + z_p[4,10] + z_p[10,4] >= -1
 OverlapInBin_4_10_toBin_3: - p_ij[4,3] - p_ij[10,3] + x_p[4,10]
   + x_p[10,4] + z_p[4,10] + z_p[10,4] >= -1
 OverlapInBin_4_11_toBin_0: - p_ij[4,0] - p_ij[11,0] + x_p[4,11]
   + x_p[11,4] + z_p[4,11] + z_p[11,4] >= -1
 OverlapInBin_4_11_toBin_1: - p_ij[4,1] - p_ij[11,1] + x_p[4,11]
   + x_p[11,4] + z_p[4,11] + z_p[11,4] >= -1
 OverlapInBin_4_11_toBin_2: - p_ij[4,2] - p_ij[11,2] + x_p[4,11]
   + x_p[11,4] + z_p[4,11] + z_p[11,4] >= -1
 OverlapInBin_4_11_toBin_3: - p_ij[4,3] - p_ij[11,3] + x_p[4,11]
   + x_p[11,4] + z_p[4,11] + z_p[11,4] >= -1
 OverlapInBin_4_12_toBin_0: - p_ij[4,0] - p_ij[12,0] + x_p[4,12]
   + x_p[12,4] + z_p[4,12] + z_p[12,4] >= -1
 OverlapInBin_4_12_toBin_1: - p_ij[4,1] - p_ij[12,1] + x_p[4,12]
   + x_p[12,4] + z_p[4,12] + z_p[12,4] >= -1
 OverlapInBin_4_12_toBin_2: - p_ij[4,2] - p_ij[12,2] + x_p[4,12]
   + x_p[12,4] + z_p[4,12] + z_p[12,4] >= -1
 OverlapInBin_4_12_toBin_3: - p_ij[4,3] - p_ij[12,3] + x_p[4,12]
   + x_p[12,4] + z_p[4,12] + z_p[12,4] >= -1
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
 OverlapInBin_5_8_toBin_0: - p_ij[5,0] - p_ij[8,0] + x_p[5,8] + x_p[8,5]
   + z_p[5,8] + z_p[8,5] >= -1
 OverlapInBin_5_8_toBin_1: - p_ij[5,1] - p_ij[8,1] + x_p[5,8] + x_p[8,5]
   + z_p[5,8] + z_p[8,5] >= -1
 OverlapInBin_5_8_toBin_2: - p_ij[5,2] - p_ij[8,2] + x_p[5,8] + x_p[8,5]
   + z_p[5,8] + z_p[8,5] >= -1
 OverlapInBin_5_8_toBin_3: - p_ij[5,3] - p_ij[8,3] + x_p[5,8] + x_p[8,5]
   + z_p[5,8] + z_p[8,5] >= -1
 OverlapInBin_5_9_toBin_0: - p_ij[5,0] - p_ij[9,0] + x_p[5,9] + x_p[9,5]
   + z_p[5,9] + z_p[9,5] >= -1
 OverlapInBin_5_9_toBin_1: - p_ij[5,1] - p_ij[9,1] + x_p[5,9] + x_p[9,5]
   + z_p[5,9] + z_p[9,5] >= -1
 OverlapInBin_5_9_toBin_2: - p_ij[5,2] - p_ij[9,2] + x_p[5,9] + x_p[9,5]
   + z_p[5,9] + z_p[9,5] >= -1
 OverlapInBin_5_9_toBin_3: - p_ij[5,3] - p_ij[9,3] + x_p[5,9] + x_p[9,5]
   + z_p[5,9] + z_p[9,5] >= -1
 OverlapInBin_5_10_toBin_0: - p_ij[5,0] - p_ij[10,0] + x_p[5,10]
   + x_p[10,5] + z_p[5,10] + z_p[10,5] >= -1
 OverlapInBin_5_10_toBin_1: - p_ij[5,1] - p_ij[10,1] + x_p[5,10]
   + x_p[10,5] + z_p[5,10] + z_p[10,5] >= -1
 OverlapInBin_5_10_toBin_2: - p_ij[5,2] - p_ij[10,2] + x_p[5,10]
   + x_p[10,5] + z_p[5,10] + z_p[10,5] >= -1
 OverlapInBin_5_10_toBin_3: - p_ij[5,3] - p_ij[10,3] + x_p[5,10]
   + x_p[10,5] + z_p[5,10] + z_p[10,5] >= -1
 OverlapInBin_5_11_toBin_0: - p_ij[5,0] - p_ij[11,0] + x_p[5,11]
   + x_p[11,5] + z_p[5,11] + z_p[11,5] >= -1
 OverlapInBin_5_11_toBin_1: - p_ij[5,1] - p_ij[11,1] + x_p[5,11]
   + x_p[11,5] + z_p[5,11] + z_p[11,5] >= -1
 OverlapInBin_5_11_toBin_2: - p_ij[5,2] - p_ij[11,2] + x_p[5,11]
   + x_p[11,5] + z_p[5,11] + z_p[11,5] >= -1
 OverlapInBin_5_11_toBin_3: - p_ij[5,3] - p_ij[11,3] + x_p[5,11]
   + x_p[11,5] + z_p[5,11] + z_p[11,5] >= -1
 OverlapInBin_5_12_toBin_0: - p_ij[5,0] - p_ij[12,0] + x_p[5,12]
   + x_p[12,5] + z_p[5,12] + z_p[12,5] >= -1
 OverlapInBin_5_12_toBin_1: - p_ij[5,1] - p_ij[12,1] + x_p[5,12]
   + x_p[12,5] + z_p[5,12] + z_p[12,5] >= -1
 OverlapInBin_5_12_toBin_2: - p_ij[5,2] - p_ij[12,2] + x_p[5,12]
   + x_p[12,5] + z_p[5,12] + z_p[12,5] >= -1
 OverlapInBin_5_12_toBin_3: - p_ij[5,3] - p_ij[12,3] + x_p[5,12]
   + x_p[12,5] + z_p[5,12] + z_p[12,5] >= -1
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
 OverlapInBin_6_8_toBin_0: - p_ij[6,0] - p_ij[8,0] + x_p[6,8] + x_p[8,6]
   + z_p[6,8] + z_p[8,6] >= -1
 OverlapInBin_6_8_toBin_1: - p_ij[6,1] - p_ij[8,1] + x_p[6,8] + x_p[8,6]
   + z_p[6,8] + z_p[8,6] >= -1
 OverlapInBin_6_8_toBin_2: - p_ij[6,2] - p_ij[8,2] + x_p[6,8] + x_p[8,6]
   + z_p[6,8] + z_p[8,6] >= -1
 OverlapInBin_6_8_toBin_3: - p_ij[6,3] - p_ij[8,3] + x_p[6,8] + x_p[8,6]
   + z_p[6,8] + z_p[8,6] >= -1
 OverlapInBin_6_9_toBin_0: - p_ij[6,0] - p_ij[9,0] + x_p[6,9] + x_p[9,6]
   + z_p[6,9] + z_p[9,6] >= -1
 OverlapInBin_6_9_toBin_1: - p_ij[6,1] - p_ij[9,1] + x_p[6,9] + x_p[9,6]
   + z_p[6,9] + z_p[9,6] >= -1
 OverlapInBin_6_9_toBin_2: - p_ij[6,2] - p_ij[9,2] + x_p[6,9] + x_p[9,6]
   + z_p[6,9] + z_p[9,6] >= -1
 OverlapInBin_6_9_toBin_3: - p_ij[6,3] - p_ij[9,3] + x_p[6,9] + x_p[9,6]
   + z_p[6,9] + z_p[9,6] >= -1
 OverlapInBin_6_10_toBin_0: - p_ij[6,0] - p_ij[10,0] + x_p[6,10]
   + x_p[10,6] + z_p[6,10] + z_p[10,6] >= -1
 OverlapInBin_6_10_toBin_1: - p_ij[6,1] - p_ij[10,1] + x_p[6,10]
   + x_p[10,6] + z_p[6,10] + z_p[10,6] >= -1
 OverlapInBin_6_10_toBin_2: - p_ij[6,2] - p_ij[10,2] + x_p[6,10]
   + x_p[10,6] + z_p[6,10] + z_p[10,6] >= -1
 OverlapInBin_6_10_toBin_3: - p_ij[6,3] - p_ij[10,3] + x_p[6,10]
   + x_p[10,6] + z_p[6,10] + z_p[10,6] >= -1
 OverlapInBin_6_11_toBin_0: - p_ij[6,0] - p_ij[11,0] + x_p[6,11]
   + x_p[11,6] + z_p[6,11] + z_p[11,6] >= -1
 OverlapInBin_6_11_toBin_1: - p_ij[6,1] - p_ij[11,1] + x_p[6,11]
   + x_p[11,6] + z_p[6,11] + z_p[11,6] >= -1
 OverlapInBin_6_11_toBin_2: - p_ij[6,2] - p_ij[11,2] + x_p[6,11]
   + x_p[11,6] + z_p[6,11] + z_p[11,6] >= -1
 OverlapInBin_6_11_toBin_3: - p_ij[6,3] - p_ij[11,3] + x_p[6,11]
   + x_p[11,6] + z_p[6,11] + z_p[11,6] >= -1
 OverlapInBin_6_12_toBin_0: - p_ij[6,0] - p_ij[12,0] + x_p[6,12]
   + x_p[12,6] + z_p[6,12] + z_p[12,6] >= -1
 OverlapInBin_6_12_toBin_1: - p_ij[6,1] - p_ij[12,1] + x_p[6,12]
   + x_p[12,6] + z_p[6,12] + z_p[12,6] >= -1
 OverlapInBin_6_12_toBin_2: - p_ij[6,2] - p_ij[12,2] + x_p[6,12]
   + x_p[12,6] + z_p[6,12] + z_p[12,6] >= -1
 OverlapInBin_6_12_toBin_3: - p_ij[6,3] - p_ij[12,3] + x_p[6,12]
   + x_p[12,6] + z_p[6,12] + z_p[12,6] >= -1
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
 OverlapInBin_7_8_toBin_0: - p_ij[7,0] - p_ij[8,0] + x_p[7,8] + x_p[8,7]
   + z_p[7,8] + z_p[8,7] >= -1
 OverlapInBin_7_8_toBin_1: - p_ij[7,1] - p_ij[8,1] + x_p[7,8] + x_p[8,7]
   + z_p[7,8] + z_p[8,7] >= -1
 OverlapInBin_7_8_toBin_2: - p_ij[7,2] - p_ij[8,2] + x_p[7,8] + x_p[8,7]
   + z_p[7,8] + z_p[8,7] >= -1
 OverlapInBin_7_8_toBin_3: - p_ij[7,3] - p_ij[8,3] + x_p[7,8] + x_p[8,7]
   + z_p[7,8] + z_p[8,7] >= -1
 OverlapInBin_7_9_toBin_0: - p_ij[7,0] - p_ij[9,0] + x_p[7,9] + x_p[9,7]
   + z_p[7,9] + z_p[9,7] >= -1
 OverlapInBin_7_9_toBin_1: - p_ij[7,1] - p_ij[9,1] + x_p[7,9] + x_p[9,7]
   + z_p[7,9] + z_p[9,7] >= -1
 OverlapInBin_7_9_toBin_2: - p_ij[7,2] - p_ij[9,2] + x_p[7,9] + x_p[9,7]
   + z_p[7,9] + z_p[9,7] >= -1
 OverlapInBin_7_9_toBin_3: - p_ij[7,3] - p_ij[9,3] + x_p[7,9] + x_p[9,7]
   + z_p[7,9] + z_p[9,7] >= -1
 OverlapInBin_7_10_toBin_0: - p_ij[7,0] - p_ij[10,0] + x_p[7,10]
   + x_p[10,7] + z_p[7,10] + z_p[10,7] >= -1
 OverlapInBin_7_10_toBin_1: - p_ij[7,1] - p_ij[10,1] + x_p[7,10]
   + x_p[10,7] + z_p[7,10] + z_p[10,7] >= -1
 OverlapInBin_7_10_toBin_2: - p_ij[7,2] - p_ij[10,2] + x_p[7,10]
   + x_p[10,7] + z_p[7,10] + z_p[10,7] >= -1
 OverlapInBin_7_10_toBin_3: - p_ij[7,3] - p_ij[10,3] + x_p[7,10]
   + x_p[10,7] + z_p[7,10] + z_p[10,7] >= -1
 OverlapInBin_7_11_toBin_0: - p_ij[7,0] - p_ij[11,0] + x_p[7,11]
   + x_p[11,7] + z_p[7,11] + z_p[11,7] >= -1
 OverlapInBin_7_11_toBin_1: - p_ij[7,1] - p_ij[11,1] + x_p[7,11]
   + x_p[11,7] + z_p[7,11] + z_p[11,7] >= -1
 OverlapInBin_7_11_toBin_2: - p_ij[7,2] - p_ij[11,2] + x_p[7,11]
   + x_p[11,7] + z_p[7,11] + z_p[11,7] >= -1
 OverlapInBin_7_11_toBin_3: - p_ij[7,3] - p_ij[11,3] + x_p[7,11]
   + x_p[11,7] + z_p[7,11] + z_p[11,7] >= -1
 OverlapInBin_7_12_toBin_0: - p_ij[7,0] - p_ij[12,0] + x_p[7,12]
   + x_p[12,7] + z_p[7,12] + z_p[12,7] >= -1
 OverlapInBin_7_12_toBin_1: - p_ij[7,1] - p_ij[12,1] + x_p[7,12]
   + x_p[12,7] + z_p[7,12] + z_p[12,7] >= -1
 OverlapInBin_7_12_toBin_2: - p_ij[7,2] - p_ij[12,2] + x_p[7,12]
   + x_p[12,7] + z_p[7,12] + z_p[12,7] >= -1
 OverlapInBin_7_12_toBin_3: - p_ij[7,3] - p_ij[12,3] + x_p[7,12]
   + x_p[12,7] + z_p[7,12] + z_p[12,7] >= -1
 OverlapInBin_8_0_toBin_0: - p_ij[0,0] - p_ij[8,0] + x_p[0,8] + x_p[8,0]
   + z_p[0,8] + z_p[8,0] >= -1
 OverlapInBin_8_0_toBin_1: - p_ij[0,1] - p_ij[8,1] + x_p[0,8] + x_p[8,0]
   + z_p[0,8] + z_p[8,0] >= -1
 OverlapInBin_8_0_toBin_2: - p_ij[0,2] - p_ij[8,2] + x_p[0,8] + x_p[8,0]
   + z_p[0,8] + z_p[8,0] >= -1
 OverlapInBin_8_0_toBin_3: - p_ij[0,3] - p_ij[8,3] + x_p[0,8] + x_p[8,0]
   + z_p[0,8] + z_p[8,0] >= -1
 OverlapInBin_8_1_toBin_0: - p_ij[1,0] - p_ij[8,0] + x_p[1,8] + x_p[8,1]
   + z_p[1,8] + z_p[8,1] >= -1
 OverlapInBin_8_1_toBin_1: - p_ij[1,1] - p_ij[8,1] + x_p[1,8] + x_p[8,1]
   + z_p[1,8] + z_p[8,1] >= -1
 OverlapInBin_8_1_toBin_2: - p_ij[1,2] - p_ij[8,2] + x_p[1,8] + x_p[8,1]
   + z_p[1,8] + z_p[8,1] >= -1
 OverlapInBin_8_1_toBin_3: - p_ij[1,3] - p_ij[8,3] + x_p[1,8] + x_p[8,1]
   + z_p[1,8] + z_p[8,1] >= -1
 OverlapInBin_8_2_toBin_0: - p_ij[2,0] - p_ij[8,0] + x_p[2,8] + x_p[8,2]
   + z_p[2,8] + z_p[8,2] >= -1
 OverlapInBin_8_2_toBin_1: - p_ij[2,1] - p_ij[8,1] + x_p[2,8] + x_p[8,2]
   + z_p[2,8] + z_p[8,2] >= -1
 OverlapInBin_8_2_toBin_2: - p_ij[2,2] - p_ij[8,2] + x_p[2,8] + x_p[8,2]
   + z_p[2,8] + z_p[8,2] >= -1
 OverlapInBin_8_2_toBin_3: - p_ij[2,3] - p_ij[8,3] + x_p[2,8] + x_p[8,2]
   + z_p[2,8] + z_p[8,2] >= -1
 OverlapInBin_8_3_toBin_0: - p_ij[3,0] - p_ij[8,0] + x_p[3,8] + x_p[8,3]
   + z_p[3,8] + z_p[8,3] >= -1
 OverlapInBin_8_3_toBin_1: - p_ij[3,1] - p_ij[8,1] + x_p[3,8] + x_p[8,3]
   + z_p[3,8] + z_p[8,3] >= -1
 OverlapInBin_8_3_toBin_2: - p_ij[3,2] - p_ij[8,2] + x_p[3,8] + x_p[8,3]
   + z_p[3,8] + z_p[8,3] >= -1
 OverlapInBin_8_3_toBin_3: - p_ij[3,3] - p_ij[8,3] + x_p[3,8] + x_p[8,3]
   + z_p[3,8] + z_p[8,3] >= -1
 OverlapInBin_8_4_toBin_0: - p_ij[4,0] - p_ij[8,0] + x_p[4,8] + x_p[8,4]
   + z_p[4,8] + z_p[8,4] >= -1
 OverlapInBin_8_4_toBin_1: - p_ij[4,1] - p_ij[8,1] + x_p[4,8] + x_p[8,4]
   + z_p[4,8] + z_p[8,4] >= -1
 OverlapInBin_8_4_toBin_2: - p_ij[4,2] - p_ij[8,2] + x_p[4,8] + x_p[8,4]
   + z_p[4,8] + z_p[8,4] >= -1
 OverlapInBin_8_4_toBin_3: - p_ij[4,3] - p_ij[8,3] + x_p[4,8] + x_p[8,4]
   + z_p[4,8] + z_p[8,4] >= -1
 OverlapInBin_8_5_toBin_0: - p_ij[5,0] - p_ij[8,0] + x_p[5,8] + x_p[8,5]
   + z_p[5,8] + z_p[8,5] >= -1
 OverlapInBin_8_5_toBin_1: - p_ij[5,1] - p_ij[8,1] + x_p[5,8] + x_p[8,5]
   + z_p[5,8] + z_p[8,5] >= -1
 OverlapInBin_8_5_toBin_2: - p_ij[5,2] - p_ij[8,2] + x_p[5,8] + x_p[8,5]
   + z_p[5,8] + z_p[8,5] >= -1
 OverlapInBin_8_5_toBin_3: - p_ij[5,3] - p_ij[8,3] + x_p[5,8] + x_p[8,5]
   + z_p[5,8] + z_p[8,5] >= -1
 OverlapInBin_8_6_toBin_0: - p_ij[6,0] - p_ij[8,0] + x_p[6,8] + x_p[8,6]
   + z_p[6,8] + z_p[8,6] >= -1
 OverlapInBin_8_6_toBin_1: - p_ij[6,1] - p_ij[8,1] + x_p[6,8] + x_p[8,6]
   + z_p[6,8] + z_p[8,6] >= -1
 OverlapInBin_8_6_toBin_2: - p_ij[6,2] - p_ij[8,2] + x_p[6,8] + x_p[8,6]
   + z_p[6,8] + z_p[8,6] >= -1
 OverlapInBin_8_6_toBin_3: - p_ij[6,3] - p_ij[8,3] + x_p[6,8] + x_p[8,6]
   + z_p[6,8] + z_p[8,6] >= -1
 OverlapInBin_8_7_toBin_0: - p_ij[7,0] - p_ij[8,0] + x_p[7,8] + x_p[8,7]
   + z_p[7,8] + z_p[8,7] >= -1
 OverlapInBin_8_7_toBin_1: - p_ij[7,1] - p_ij[8,1] + x_p[7,8] + x_p[8,7]
   + z_p[7,8] + z_p[8,7] >= -1
 OverlapInBin_8_7_toBin_2: - p_ij[7,2] - p_ij[8,2] + x_p[7,8] + x_p[8,7]
   + z_p[7,8] + z_p[8,7] >= -1
 OverlapInBin_8_7_toBin_3: - p_ij[7,3] - p_ij[8,3] + x_p[7,8] + x_p[8,7]
   + z_p[7,8] + z_p[8,7] >= -1
 OverlapInBin_8_9_toBin_0: - p_ij[8,0] - p_ij[9,0] + x_p[8,9] + x_p[9,8]
   + z_p[8,9] + z_p[9,8] >= -1
 OverlapInBin_8_9_toBin_1: - p_ij[8,1] - p_ij[9,1] + x_p[8,9] + x_p[9,8]
   + z_p[8,9] + z_p[9,8] >= -1
 OverlapInBin_8_9_toBin_2: - p_ij[8,2] - p_ij[9,2] + x_p[8,9] + x_p[9,8]
   + z_p[8,9] + z_p[9,8] >= -1
 OverlapInBin_8_9_toBin_3: - p_ij[8,3] - p_ij[9,3] + x_p[8,9] + x_p[9,8]
   + z_p[8,9] + z_p[9,8] >= -1
 OverlapInBin_8_10_toBin_0: - p_ij[8,0] - p_ij[10,0] + x_p[8,10]
   + x_p[10,8] + z_p[8,10] + z_p[10,8] >= -1
 OverlapInBin_8_10_toBin_1: - p_ij[8,1] - p_ij[10,1] + x_p[8,10]
   + x_p[10,8] + z_p[8,10] + z_p[10,8] >= -1
 OverlapInBin_8_10_toBin_2: - p_ij[8,2] - p_ij[10,2] + x_p[8,10]
   + x_p[10,8] + z_p[8,10] + z_p[10,8] >= -1
 OverlapInBin_8_10_toBin_3: - p_ij[8,3] - p_ij[10,3] + x_p[8,10]
   + x_p[10,8] + z_p[8,10] + z_p[10,8] >= -1
 OverlapInBin_8_11_toBin_0: - p_ij[8,0] - p_ij[11,0] + x_p[8,11]
   + x_p[11,8] + z_p[8,11] + z_p[11,8] >= -1
 OverlapInBin_8_11_toBin_1: - p_ij[8,1] - p_ij[11,1] + x_p[8,11]
   + x_p[11,8] + z_p[8,11] + z_p[11,8] >= -1
 OverlapInBin_8_11_toBin_2: - p_ij[8,2] - p_ij[11,2] + x_p[8,11]
   + x_p[11,8] + z_p[8,11] + z_p[11,8] >= -1
 OverlapInBin_8_11_toBin_3: - p_ij[8,3] - p_ij[11,3] + x_p[8,11]
   + x_p[11,8] + z_p[8,11] + z_p[11,8] >= -1
 OverlapInBin_8_12_toBin_0: - p_ij[8,0] - p_ij[12,0] + x_p[8,12]
   + x_p[12,8] + z_p[8,12] + z_p[12,8] >= -1
 OverlapInBin_8_12_toBin_1: - p_ij[8,1] - p_ij[12,1] + x_p[8,12]
   + x_p[12,8] + z_p[8,12] + z_p[12,8] >= -1
 OverlapInBin_8_12_toBin_2: - p_ij[8,2] - p_ij[12,2] + x_p[8,12]
   + x_p[12,8] + z_p[8,12] + z_p[12,8] >= -1
 OverlapInBin_8_12_toBin_3: - p_ij[8,3] - p_ij[12,3] + x_p[8,12]
   + x_p[12,8] + z_p[8,12] + z_p[12,8] >= -1
 OverlapInBin_9_0_toBin_0: - p_ij[0,0] - p_ij[9,0] + x_p[0,9] + x_p[9,0]
   + z_p[0,9] + z_p[9,0] >= -1
 OverlapInBin_9_0_toBin_1: - p_ij[0,1] - p_ij[9,1] + x_p[0,9] + x_p[9,0]
   + z_p[0,9] + z_p[9,0] >= -1
 OverlapInBin_9_0_toBin_2: - p_ij[0,2] - p_ij[9,2] + x_p[0,9] + x_p[9,0]
   + z_p[0,9] + z_p[9,0] >= -1
 OverlapInBin_9_0_toBin_3: - p_ij[0,3] - p_ij[9,3] + x_p[0,9] + x_p[9,0]
   + z_p[0,9] + z_p[9,0] >= -1
 OverlapInBin_9_1_toBin_0: - p_ij[1,0] - p_ij[9,0] + x_p[1,9] + x_p[9,1]
   + z_p[1,9] + z_p[9,1] >= -1
 OverlapInBin_9_1_toBin_1: - p_ij[1,1] - p_ij[9,1] + x_p[1,9] + x_p[9,1]
   + z_p[1,9] + z_p[9,1] >= -1
 OverlapInBin_9_1_toBin_2: - p_ij[1,2] - p_ij[9,2] + x_p[1,9] + x_p[9,1]
   + z_p[1,9] + z_p[9,1] >= -1
 OverlapInBin_9_1_toBin_3: - p_ij[1,3] - p_ij[9,3] + x_p[1,9] + x_p[9,1]
   + z_p[1,9] + z_p[9,1] >= -1
 OverlapInBin_9_2_toBin_0: - p_ij[2,0] - p_ij[9,0] + x_p[2,9] + x_p[9,2]
   + z_p[2,9] + z_p[9,2] >= -1
 OverlapInBin_9_2_toBin_1: - p_ij[2,1] - p_ij[9,1] + x_p[2,9] + x_p[9,2]
   + z_p[2,9] + z_p[9,2] >= -1
 OverlapInBin_9_2_toBin_2: - p_ij[2,2] - p_ij[9,2] + x_p[2,9] + x_p[9,2]
   + z_p[2,9] + z_p[9,2] >= -1
 OverlapInBin_9_2_toBin_3: - p_ij[2,3] - p_ij[9,3] + x_p[2,9] + x_p[9,2]
   + z_p[2,9] + z_p[9,2] >= -1
 OverlapInBin_9_3_toBin_0: - p_ij[3,0] - p_ij[9,0] + x_p[3,9] + x_p[9,3]
   + z_p[3,9] + z_p[9,3] >= -1
 OverlapInBin_9_3_toBin_1: - p_ij[3,1] - p_ij[9,1] + x_p[3,9] + x_p[9,3]
   + z_p[3,9] + z_p[9,3] >= -1
 OverlapInBin_9_3_toBin_2: - p_ij[3,2] - p_ij[9,2] + x_p[3,9] + x_p[9,3]
   + z_p[3,9] + z_p[9,3] >= -1
 OverlapInBin_9_3_toBin_3: - p_ij[3,3] - p_ij[9,3] + x_p[3,9] + x_p[9,3]
   + z_p[3,9] + z_p[9,3] >= -1
 OverlapInBin_9_4_toBin_0: - p_ij[4,0] - p_ij[9,0] + x_p[4,9] + x_p[9,4]
   + z_p[4,9] + z_p[9,4] >= -1
 OverlapInBin_9_4_toBin_1: - p_ij[4,1] - p_ij[9,1] + x_p[4,9] + x_p[9,4]
   + z_p[4,9] + z_p[9,4] >= -1
 OverlapInBin_9_4_toBin_2: - p_ij[4,2] - p_ij[9,2] + x_p[4,9] + x_p[9,4]
   + z_p[4,9] + z_p[9,4] >= -1
 OverlapInBin_9_4_toBin_3: - p_ij[4,3] - p_ij[9,3] + x_p[4,9] + x_p[9,4]
   + z_p[4,9] + z_p[9,4] >= -1
 OverlapInBin_9_5_toBin_0: - p_ij[5,0] - p_ij[9,0] + x_p[5,9] + x_p[9,5]
   + z_p[5,9] + z_p[9,5] >= -1
 OverlapInBin_9_5_toBin_1: - p_ij[5,1] - p_ij[9,1] + x_p[5,9] + x_p[9,5]
   + z_p[5,9] + z_p[9,5] >= -1
 OverlapInBin_9_5_toBin_2: - p_ij[5,2] - p_ij[9,2] + x_p[5,9] + x_p[9,5]
   + z_p[5,9] + z_p[9,5] >= -1
 OverlapInBin_9_5_toBin_3: - p_ij[5,3] - p_ij[9,3] + x_p[5,9] + x_p[9,5]
   + z_p[5,9] + z_p[9,5] >= -1
 OverlapInBin_9_6_toBin_0: - p_ij[6,0] - p_ij[9,0] + x_p[6,9] + x_p[9,6]
   + z_p[6,9] + z_p[9,6] >= -1
 OverlapInBin_9_6_toBin_1: - p_ij[6,1] - p_ij[9,1] + x_p[6,9] + x_p[9,6]
   + z_p[6,9] + z_p[9,6] >= -1
 OverlapInBin_9_6_toBin_2: - p_ij[6,2] - p_ij[9,2] + x_p[6,9] + x_p[9,6]
   + z_p[6,9] + z_p[9,6] >= -1
 OverlapInBin_9_6_toBin_3: - p_ij[6,3] - p_ij[9,3] + x_p[6,9] + x_p[9,6]
   + z_p[6,9] + z_p[9,6] >= -1
 OverlapInBin_9_7_toBin_0: - p_ij[7,0] - p_ij[9,0] + x_p[7,9] + x_p[9,7]
   + z_p[7,9] + z_p[9,7] >= -1
 OverlapInBin_9_7_toBin_1: - p_ij[7,1] - p_ij[9,1] + x_p[7,9] + x_p[9,7]
   + z_p[7,9] + z_p[9,7] >= -1
 OverlapInBin_9_7_toBin_2: - p_ij[7,2] - p_ij[9,2] + x_p[7,9] + x_p[9,7]
   + z_p[7,9] + z_p[9,7] >= -1
 OverlapInBin_9_7_toBin_3: - p_ij[7,3] - p_ij[9,3] + x_p[7,9] + x_p[9,7]
   + z_p[7,9] + z_p[9,7] >= -1
 OverlapInBin_9_8_toBin_0: - p_ij[8,0] - p_ij[9,0] + x_p[8,9] + x_p[9,8]
   + z_p[8,9] + z_p[9,8] >= -1
 OverlapInBin_9_8_toBin_1: - p_ij[8,1] - p_ij[9,1] + x_p[8,9] + x_p[9,8]
   + z_p[8,9] + z_p[9,8] >= -1
 OverlapInBin_9_8_toBin_2: - p_ij[8,2] - p_ij[9,2] + x_p[8,9] + x_p[9,8]
   + z_p[8,9] + z_p[9,8] >= -1
 OverlapInBin_9_8_toBin_3: - p_ij[8,3] - p_ij[9,3] + x_p[8,9] + x_p[9,8]
   + z_p[8,9] + z_p[9,8] >= -1
 OverlapInBin_9_10_toBin_0: - p_ij[9,0] - p_ij[10,0] + x_p[9,10]
   + x_p[10,9] + z_p[9,10] + z_p[10,9] >= -1
 OverlapInBin_9_10_toBin_1: - p_ij[9,1] - p_ij[10,1] + x_p[9,10]
   + x_p[10,9] + z_p[9,10] + z_p[10,9] >= -1
 OverlapInBin_9_10_toBin_2: - p_ij[9,2] - p_ij[10,2] + x_p[9,10]
   + x_p[10,9] + z_p[9,10] + z_p[10,9] >= -1
 OverlapInBin_9_10_toBin_3: - p_ij[9,3] - p_ij[10,3] + x_p[9,10]
   + x_p[10,9] + z_p[9,10] + z_p[10,9] >= -1
 OverlapInBin_9_11_toBin_0: - p_ij[9,0] - p_ij[11,0] + x_p[9,11]
   + x_p[11,9] + z_p[9,11] + z_p[11,9] >= -1
 OverlapInBin_9_11_toBin_1: - p_ij[9,1] - p_ij[11,1] + x_p[9,11]
   + x_p[11,9] + z_p[9,11] + z_p[11,9] >= -1
 OverlapInBin_9_11_toBin_2: - p_ij[9,2] - p_ij[11,2] + x_p[9,11]
   + x_p[11,9] + z_p[9,11] + z_p[11,9] >= -1
 OverlapInBin_9_11_toBin_3: - p_ij[9,3] - p_ij[11,3] + x_p[9,11]
   + x_p[11,9] + z_p[9,11] + z_p[11,9] >= -1
 OverlapInBin_9_12_toBin_0: - p_ij[9,0] - p_ij[12,0] + x_p[9,12]
   + x_p[12,9] + z_p[9,12] + z_p[12,9] >= -1
 OverlapInBin_9_12_toBin_1: - p_ij[9,1] - p_ij[12,1] + x_p[9,12]
   + x_p[12,9] + z_p[9,12] + z_p[12,9] >= -1
 OverlapInBin_9_12_toBin_2: - p_ij[9,2] - p_ij[12,2] + x_p[9,12]
   + x_p[12,9] + z_p[9,12] + z_p[12,9] >= -1
 OverlapInBin_9_12_toBin_3: - p_ij[9,3] - p_ij[12,3] + x_p[9,12]
   + x_p[12,9] + z_p[9,12] + z_p[12,9] >= -1
 OverlapInBin_10_0_toBin_0: - p_ij[0,0] - p_ij[10,0] + x_p[0,10]
   + x_p[10,0] + z_p[0,10] + z_p[10,0] >= -1
 OverlapInBin_10_0_toBin_1: - p_ij[0,1] - p_ij[10,1] + x_p[0,10]
   + x_p[10,0] + z_p[0,10] + z_p[10,0] >= -1
 OverlapInBin_10_0_toBin_2: - p_ij[0,2] - p_ij[10,2] + x_p[0,10]
   + x_p[10,0] + z_p[0,10] + z_p[10,0] >= -1
 OverlapInBin_10_0_toBin_3: - p_ij[0,3] - p_ij[10,3] + x_p[0,10]
   + x_p[10,0] + z_p[0,10] + z_p[10,0] >= -1
 OverlapInBin_10_1_toBin_0: - p_ij[1,0] - p_ij[10,0] + x_p[1,10]
   + x_p[10,1] + z_p[1,10] + z_p[10,1] >= -1
 OverlapInBin_10_1_toBin_1: - p_ij[1,1] - p_ij[10,1] + x_p[1,10]
   + x_p[10,1] + z_p[1,10] + z_p[10,1] >= -1
 OverlapInBin_10_1_toBin_2: - p_ij[1,2] - p_ij[10,2] + x_p[1,10]
   + x_p[10,1] + z_p[1,10] + z_p[10,1] >= -1
 OverlapInBin_10_1_toBin_3: - p_ij[1,3] - p_ij[10,3] + x_p[1,10]
   + x_p[10,1] + z_p[1,10] + z_p[10,1] >= -1
 OverlapInBin_10_2_toBin_0: - p_ij[2,0] - p_ij[10,0] + x_p[2,10]
   + x_p[10,2] + z_p[2,10] + z_p[10,2] >= -1
 OverlapInBin_10_2_toBin_1: - p_ij[2,1] - p_ij[10,1] + x_p[2,10]
   + x_p[10,2] + z_p[2,10] + z_p[10,2] >= -1
 OverlapInBin_10_2_toBin_2: - p_ij[2,2] - p_ij[10,2] + x_p[2,10]
   + x_p[10,2] + z_p[2,10] + z_p[10,2] >= -1
 OverlapInBin_10_2_toBin_3: - p_ij[2,3] - p_ij[10,3] + x_p[2,10]
   + x_p[10,2] + z_p[2,10] + z_p[10,2] >= -1
 OverlapInBin_10_3_toBin_0: - p_ij[3,0] - p_ij[10,0] + x_p[3,10]
   + x_p[10,3] + z_p[3,10] + z_p[10,3] >= -1
 OverlapInBin_10_3_toBin_1: - p_ij[3,1] - p_ij[10,1] + x_p[3,10]
   + x_p[10,3] + z_p[3,10] + z_p[10,3] >= -1
 OverlapInBin_10_3_toBin_2: - p_ij[3,2] - p_ij[10,2] + x_p[3,10]
   + x_p[10,3] + z_p[3,10] + z_p[10,3] >= -1
 OverlapInBin_10_3_toBin_3: - p_ij[3,3] - p_ij[10,3] + x_p[3,10]
   + x_p[10,3] + z_p[3,10] + z_p[10,3] >= -1
 OverlapInBin_10_4_toBin_0: - p_ij[4,0] - p_ij[10,0] + x_p[4,10]
   + x_p[10,4] + z_p[4,10] + z_p[10,4] >= -1
 OverlapInBin_10_4_toBin_1: - p_ij[4,1] - p_ij[10,1] + x_p[4,10]
   + x_p[10,4] + z_p[4,10] + z_p[10,4] >= -1
 OverlapInBin_10_4_toBin_2: - p_ij[4,2] - p_ij[10,2] + x_p[4,10]
   + x_p[10,4] + z_p[4,10] + z_p[10,4] >= -1
 OverlapInBin_10_4_toBin_3: - p_ij[4,3] - p_ij[10,3] + x_p[4,10]
   + x_p[10,4] + z_p[4,10] + z_p[10,4] >= -1
 OverlapInBin_10_5_toBin_0: - p_ij[5,0] - p_ij[10,0] + x_p[5,10]
   + x_p[10,5] + z_p[5,10] + z_p[10,5] >= -1
 OverlapInBin_10_5_toBin_1: - p_ij[5,1] - p_ij[10,1] + x_p[5,10]
   + x_p[10,5] + z_p[5,10] + z_p[10,5] >= -1
 OverlapInBin_10_5_toBin_2: - p_ij[5,2] - p_ij[10,2] + x_p[5,10]
   + x_p[10,5] + z_p[5,10] + z_p[10,5] >= -1
 OverlapInBin_10_5_toBin_3: - p_ij[5,3] - p_ij[10,3] + x_p[5,10]
   + x_p[10,5] + z_p[5,10] + z_p[10,5] >= -1
 OverlapInBin_10_6_toBin_0: - p_ij[6,0] - p_ij[10,0] + x_p[6,10]
   + x_p[10,6] + z_p[6,10] + z_p[10,6] >= -1
 OverlapInBin_10_6_toBin_1: - p_ij[6,1] - p_ij[10,1] + x_p[6,10]
   + x_p[10,6] + z_p[6,10] + z_p[10,6] >= -1
 OverlapInBin_10_6_toBin_2: - p_ij[6,2] - p_ij[10,2] + x_p[6,10]
   + x_p[10,6] + z_p[6,10] + z_p[10,6] >= -1
 OverlapInBin_10_6_toBin_3: - p_ij[6,3] - p_ij[10,3] + x_p[6,10]
   + x_p[10,6] + z_p[6,10] + z_p[10,6] >= -1
 OverlapInBin_10_7_toBin_0: - p_ij[7,0] - p_ij[10,0] + x_p[7,10]
   + x_p[10,7] + z_p[7,10] + z_p[10,7] >= -1
 OverlapInBin_10_7_toBin_1: - p_ij[7,1] - p_ij[10,1] + x_p[7,10]
   + x_p[10,7] + z_p[7,10] + z_p[10,7] >= -1
 OverlapInBin_10_7_toBin_2: - p_ij[7,2] - p_ij[10,2] + x_p[7,10]
   + x_p[10,7] + z_p[7,10] + z_p[10,7] >= -1
 OverlapInBin_10_7_toBin_3: - p_ij[7,3] - p_ij[10,3] + x_p[7,10]
   + x_p[10,7] + z_p[7,10] + z_p[10,7] >= -1
 OverlapInBin_10_8_toBin_0: - p_ij[8,0] - p_ij[10,0] + x_p[8,10]
   + x_p[10,8] + z_p[8,10] + z_p[10,8] >= -1
 OverlapInBin_10_8_toBin_1: - p_ij[8,1] - p_ij[10,1] + x_p[8,10]
   + x_p[10,8] + z_p[8,10] + z_p[10,8] >= -1
 OverlapInBin_10_8_toBin_2: - p_ij[8,2] - p_ij[10,2] + x_p[8,10]
   + x_p[10,8] + z_p[8,10] + z_p[10,8] >= -1
 OverlapInBin_10_8_toBin_3: - p_ij[8,3] - p_ij[10,3] + x_p[8,10]
   + x_p[10,8] + z_p[8,10] + z_p[10,8] >= -1
 OverlapInBin_10_9_toBin_0: - p_ij[9,0] - p_ij[10,0] + x_p[9,10]
   + x_p[10,9] + z_p[9,10] + z_p[10,9] >= -1
 OverlapInBin_10_9_toBin_1: - p_ij[9,1] - p_ij[10,1] + x_p[9,10]
   + x_p[10,9] + z_p[9,10] + z_p[10,9] >= -1
 OverlapInBin_10_9_toBin_2: - p_ij[9,2] - p_ij[10,2] + x_p[9,10]
   + x_p[10,9] + z_p[9,10] + z_p[10,9] >= -1
 OverlapInBin_10_9_toBin_3: - p_ij[9,3] - p_ij[10,3] + x_p[9,10]
   + x_p[10,9] + z_p[9,10] + z_p[10,9] >= -1
 OverlapInBin_10_11_toBin_0: - p_ij[10,0] - p_ij[11,0] + x_p[10,11]
   + x_p[11,10] + z_p[10,11] + z_p[11,10] >= -1
 OverlapInBin_10_11_toBin_1: - p_ij[10,1] - p_ij[11,1] + x_p[10,11]
   + x_p[11,10] + z_p[10,11] + z_p[11,10] >= -1
 OverlapInBin_10_11_toBin_2: - p_ij[10,2] - p_ij[11,2] + x_p[10,11]
   + x_p[11,10] + z_p[10,11] + z_p[11,10] >= -1
 OverlapInBin_10_11_toBin_3: - p_ij[10,3] - p_ij[11,3] + x_p[10,11]
   + x_p[11,10] + z_p[10,11] + z_p[11,10] >= -1
 OverlapInBin_10_12_toBin_0: - p_ij[10,0] - p_ij[12,0] + x_p[10,12]
   + x_p[12,10] + z_p[10,12] + z_p[12,10] >= -1
 OverlapInBin_10_12_toBin_1: - p_ij[10,1] - p_ij[12,1] + x_p[10,12]
   + x_p[12,10] + z_p[10,12] + z_p[12,10] >= -1
 OverlapInBin_10_12_toBin_2: - p_ij[10,2] - p_ij[12,2] + x_p[10,12]
   + x_p[12,10] + z_p[10,12] + z_p[12,10] >= -1
 OverlapInBin_10_12_toBin_3: - p_ij[10,3] - p_ij[12,3] + x_p[10,12]
   + x_p[12,10] + z_p[10,12] + z_p[12,10] >= -1
 OverlapInBin_11_0_toBin_0: - p_ij[0,0] - p_ij[11,0] + x_p[0,11]
   + x_p[11,0] + z_p[0,11] + z_p[11,0] >= -1
 OverlapInBin_11_0_toBin_1: - p_ij[0,1] - p_ij[11,1] + x_p[0,11]
   + x_p[11,0] + z_p[0,11] + z_p[11,0] >= -1
 OverlapInBin_11_0_toBin_2: - p_ij[0,2] - p_ij[11,2] + x_p[0,11]
   + x_p[11,0] + z_p[0,11] + z_p[11,0] >= -1
 OverlapInBin_11_0_toBin_3: - p_ij[0,3] - p_ij[11,3] + x_p[0,11]
   + x_p[11,0] + z_p[0,11] + z_p[11,0] >= -1
 OverlapInBin_11_1_toBin_0: - p_ij[1,0] - p_ij[11,0] + x_p[1,11]
   + x_p[11,1] + z_p[1,11] + z_p[11,1] >= -1
 OverlapInBin_11_1_toBin_1: - p_ij[1,1] - p_ij[11,1] + x_p[1,11]
   + x_p[11,1] + z_p[1,11] + z_p[11,1] >= -1
 OverlapInBin_11_1_toBin_2: - p_ij[1,2] - p_ij[11,2] + x_p[1,11]
   + x_p[11,1] + z_p[1,11] + z_p[11,1] >= -1
 OverlapInBin_11_1_toBin_3: - p_ij[1,3] - p_ij[11,3] + x_p[1,11]
   + x_p[11,1] + z_p[1,11] + z_p[11,1] >= -1
 OverlapInBin_11_2_toBin_0: - p_ij[2,0] - p_ij[11,0] + x_p[2,11]
   + x_p[11,2] + z_p[2,11] + z_p[11,2] >= -1
 OverlapInBin_11_2_toBin_1: - p_ij[2,1] - p_ij[11,1] + x_p[2,11]
   + x_p[11,2] + z_p[2,11] + z_p[11,2] >= -1
 OverlapInBin_11_2_toBin_2: - p_ij[2,2] - p_ij[11,2] + x_p[2,11]
   + x_p[11,2] + z_p[2,11] + z_p[11,2] >= -1
 OverlapInBin_11_2_toBin_3: - p_ij[2,3] - p_ij[11,3] + x_p[2,11]
   + x_p[11,2] + z_p[2,11] + z_p[11,2] >= -1
 OverlapInBin_11_3_toBin_0: - p_ij[3,0] - p_ij[11,0] + x_p[3,11]
   + x_p[11,3] + z_p[3,11] + z_p[11,3] >= -1
 OverlapInBin_11_3_toBin_1: - p_ij[3,1] - p_ij[11,1] + x_p[3,11]
   + x_p[11,3] + z_p[3,11] + z_p[11,3] >= -1
 OverlapInBin_11_3_toBin_2: - p_ij[3,2] - p_ij[11,2] + x_p[3,11]
   + x_p[11,3] + z_p[3,11] + z_p[11,3] >= -1
 OverlapInBin_11_3_toBin_3: - p_ij[3,3] - p_ij[11,3] + x_p[3,11]
   + x_p[11,3] + z_p[3,11] + z_p[11,3] >= -1
 OverlapInBin_11_4_toBin_0: - p_ij[4,0] - p_ij[11,0] + x_p[4,11]
   + x_p[11,4] + z_p[4,11] + z_p[11,4] >= -1
 OverlapInBin_11_4_toBin_1: - p_ij[4,1] - p_ij[11,1] + x_p[4,11]
   + x_p[11,4] + z_p[4,11] + z_p[11,4] >= -1
 OverlapInBin_11_4_toBin_2: - p_ij[4,2] - p_ij[11,2] + x_p[4,11]
   + x_p[11,4] + z_p[4,11] + z_p[11,4] >= -1
 OverlapInBin_11_4_toBin_3: - p_ij[4,3] - p_ij[11,3] + x_p[4,11]
   + x_p[11,4] + z_p[4,11] + z_p[11,4] >= -1
 OverlapInBin_11_5_toBin_0: - p_ij[5,0] - p_ij[11,0] + x_p[5,11]
   + x_p[11,5] + z_p[5,11] + z_p[11,5] >= -1
 OverlapInBin_11_5_toBin_1: - p_ij[5,1] - p_ij[11,1] + x_p[5,11]
   + x_p[11,5] + z_p[5,11] + z_p[11,5] >= -1
 OverlapInBin_11_5_toBin_2: - p_ij[5,2] - p_ij[11,2] + x_p[5,11]
   + x_p[11,5] + z_p[5,11] + z_p[11,5] >= -1
 OverlapInBin_11_5_toBin_3: - p_ij[5,3] - p_ij[11,3] + x_p[5,11]
   + x_p[11,5] + z_p[5,11] + z_p[11,5] >= -1
 OverlapInBin_11_6_toBin_0: - p_ij[6,0] - p_ij[11,0] + x_p[6,11]
   + x_p[11,6] + z_p[6,11] + z_p[11,6] >= -1
 OverlapInBin_11_6_toBin_1: - p_ij[6,1] - p_ij[11,1] + x_p[6,11]
   + x_p[11,6] + z_p[6,11] + z_p[11,6] >= -1
 OverlapInBin_11_6_toBin_2: - p_ij[6,2] - p_ij[11,2] + x_p[6,11]
   + x_p[11,6] + z_p[6,11] + z_p[11,6] >= -1
 OverlapInBin_11_6_toBin_3: - p_ij[6,3] - p_ij[11,3] + x_p[6,11]
   + x_p[11,6] + z_p[6,11] + z_p[11,6] >= -1
 OverlapInBin_11_7_toBin_0: - p_ij[7,0] - p_ij[11,0] + x_p[7,11]
   + x_p[11,7] + z_p[7,11] + z_p[11,7] >= -1
 OverlapInBin_11_7_toBin_1: - p_ij[7,1] - p_ij[11,1] + x_p[7,11]
   + x_p[11,7] + z_p[7,11] + z_p[11,7] >= -1
 OverlapInBin_11_7_toBin_2: - p_ij[7,2] - p_ij[11,2] + x_p[7,11]
   + x_p[11,7] + z_p[7,11] + z_p[11,7] >= -1
 OverlapInBin_11_7_toBin_3: - p_ij[7,3] - p_ij[11,3] + x_p[7,11]
   + x_p[11,7] + z_p[7,11] + z_p[11,7] >= -1
 OverlapInBin_11_8_toBin_0: - p_ij[8,0] - p_ij[11,0] + x_p[8,11]
   + x_p[11,8] + z_p[8,11] + z_p[11,8] >= -1
 OverlapInBin_11_8_toBin_1: - p_ij[8,1] - p_ij[11,1] + x_p[8,11]
   + x_p[11,8] + z_p[8,11] + z_p[11,8] >= -1
 OverlapInBin_11_8_toBin_2: - p_ij[8,2] - p_ij[11,2] + x_p[8,11]
   + x_p[11,8] + z_p[8,11] + z_p[11,8] >= -1
 OverlapInBin_11_8_toBin_3: - p_ij[8,3] - p_ij[11,3] + x_p[8,11]
   + x_p[11,8] + z_p[8,11] + z_p[11,8] >= -1
 OverlapInBin_11_9_toBin_0: - p_ij[9,0] - p_ij[11,0] + x_p[9,11]
   + x_p[11,9] + z_p[9,11] + z_p[11,9] >= -1
 OverlapInBin_11_9_toBin_1: - p_ij[9,1] - p_ij[11,1] + x_p[9,11]
   + x_p[11,9] + z_p[9,11] + z_p[11,9] >= -1
 OverlapInBin_11_9_toBin_2: - p_ij[9,2] - p_ij[11,2] + x_p[9,11]
   + x_p[11,9] + z_p[9,11] + z_p[11,9] >= -1
 OverlapInBin_11_9_toBin_3: - p_ij[9,3] - p_ij[11,3] + x_p[9,11]
   + x_p[11,9] + z_p[9,11] + z_p[11,9] >= -1
 OverlapInBin_11_10_toBin_0: - p_ij[10,0] - p_ij[11,0] + x_p[10,11]
   + x_p[11,10] + z_p[10,11] + z_p[11,10] >= -1
 OverlapInBin_11_10_toBin_1: - p_ij[10,1] - p_ij[11,1] + x_p[10,11]
   + x_p[11,10] + z_p[10,11] + z_p[11,10] >= -1
 OverlapInBin_11_10_toBin_2: - p_ij[10,2] - p_ij[11,2] + x_p[10,11]
   + x_p[11,10] + z_p[10,11] + z_p[11,10] >= -1
 OverlapInBin_11_10_toBin_3: - p_ij[10,3] - p_ij[11,3] + x_p[10,11]
   + x_p[11,10] + z_p[10,11] + z_p[11,10] >= -1
 OverlapInBin_11_12_toBin_0: - p_ij[11,0] - p_ij[12,0] + x_p[11,12]
   + x_p[12,11] + z_p[11,12] + z_p[12,11] >= -1
 OverlapInBin_11_12_toBin_1: - p_ij[11,1] - p_ij[12,1] + x_p[11,12]
   + x_p[12,11] + z_p[11,12] + z_p[12,11] >= -1
 OverlapInBin_11_12_toBin_2: - p_ij[11,2] - p_ij[12,2] + x_p[11,12]
   + x_p[12,11] + z_p[11,12] + z_p[12,11] >= -1
 OverlapInBin_11_12_toBin_3: - p_ij[11,3] - p_ij[12,3] + x_p[11,12]
   + x_p[12,11] + z_p[11,12] + z_p[12,11] >= -1
 OverlapInBin_12_0_toBin_0: - p_ij[0,0] - p_ij[12,0] + x_p[0,12]
   + x_p[12,0] + z_p[0,12] + z_p[12,0] >= -1
 OverlapInBin_12_0_toBin_1: - p_ij[0,1] - p_ij[12,1] + x_p[0,12]
   + x_p[12,0] + z_p[0,12] + z_p[12,0] >= -1
 OverlapInBin_12_0_toBin_2: - p_ij[0,2] - p_ij[12,2] + x_p[0,12]
   + x_p[12,0] + z_p[0,12] + z_p[12,0] >= -1
 OverlapInBin_12_0_toBin_3: - p_ij[0,3] - p_ij[12,3] + x_p[0,12]
   + x_p[12,0] + z_p[0,12] + z_p[12,0] >= -1
 OverlapInBin_12_1_toBin_0: - p_ij[1,0] - p_ij[12,0] + x_p[1,12]
   + x_p[12,1] + z_p[1,12] + z_p[12,1] >= -1
 OverlapInBin_12_1_toBin_1: - p_ij[1,1] - p_ij[12,1] + x_p[1,12]
   + x_p[12,1] + z_p[1,12] + z_p[12,1] >= -1
 OverlapInBin_12_1_toBin_2: - p_ij[1,2] - p_ij[12,2] + x_p[1,12]
   + x_p[12,1] + z_p[1,12] + z_p[12,1] >= -1
 OverlapInBin_12_1_toBin_3: - p_ij[1,3] - p_ij[12,3] + x_p[1,12]
   + x_p[12,1] + z_p[1,12] + z_p[12,1] >= -1
 OverlapInBin_12_2_toBin_0: - p_ij[2,0] - p_ij[12,0] + x_p[2,12]
   + x_p[12,2] + z_p[2,12] + z_p[12,2] >= -1
 OverlapInBin_12_2_toBin_1: - p_ij[2,1] - p_ij[12,1] + x_p[2,12]
   + x_p[12,2] + z_p[2,12] + z_p[12,2] >= -1
 OverlapInBin_12_2_toBin_2: - p_ij[2,2] - p_ij[12,2] + x_p[2,12]
   + x_p[12,2] + z_p[2,12] + z_p[12,2] >= -1
 OverlapInBin_12_2_toBin_3: - p_ij[2,3] - p_ij[12,3] + x_p[2,12]
   + x_p[12,2] + z_p[2,12] + z_p[12,2] >= -1
 OverlapInBin_12_3_toBin_0: - p_ij[3,0] - p_ij[12,0] + x_p[3,12]
   + x_p[12,3] + z_p[3,12] + z_p[12,3] >= -1
 OverlapInBin_12_3_toBin_1: - p_ij[3,1] - p_ij[12,1] + x_p[3,12]
   + x_p[12,3] + z_p[3,12] + z_p[12,3] >= -1
 OverlapInBin_12_3_toBin_2: - p_ij[3,2] - p_ij[12,2] + x_p[3,12]
   + x_p[12,3] + z_p[3,12] + z_p[12,3] >= -1
 OverlapInBin_12_3_toBin_3: - p_ij[3,3] - p_ij[12,3] + x_p[3,12]
   + x_p[12,3] + z_p[3,12] + z_p[12,3] >= -1
 OverlapInBin_12_4_toBin_0: - p_ij[4,0] - p_ij[12,0] + x_p[4,12]
   + x_p[12,4] + z_p[4,12] + z_p[12,4] >= -1
 OverlapInBin_12_4_toBin_1: - p_ij[4,1] - p_ij[12,1] + x_p[4,12]
   + x_p[12,4] + z_p[4,12] + z_p[12,4] >= -1
 OverlapInBin_12_4_toBin_2: - p_ij[4,2] - p_ij[12,2] + x_p[4,12]
   + x_p[12,4] + z_p[4,12] + z_p[12,4] >= -1
 OverlapInBin_12_4_toBin_3: - p_ij[4,3] - p_ij[12,3] + x_p[4,12]
   + x_p[12,4] + z_p[4,12] + z_p[12,4] >= -1
 OverlapInBin_12_5_toBin_0: - p_ij[5,0] - p_ij[12,0] + x_p[5,12]
   + x_p[12,5] + z_p[5,12] + z_p[12,5] >= -1
 OverlapInBin_12_5_toBin_1: - p_ij[5,1] - p_ij[12,1] + x_p[5,12]
   + x_p[12,5] + z_p[5,12] + z_p[12,5] >= -1
 OverlapInBin_12_5_toBin_2: - p_ij[5,2] - p_ij[12,2] + x_p[5,12]
   + x_p[12,5] + z_p[5,12] + z_p[12,5] >= -1
 OverlapInBin_12_5_toBin_3: - p_ij[5,3] - p_ij[12,3] + x_p[5,12]
   + x_p[12,5] + z_p[5,12] + z_p[12,5] >= -1
 OverlapInBin_12_6_toBin_0: - p_ij[6,0] - p_ij[12,0] + x_p[6,12]
   + x_p[12,6] + z_p[6,12] + z_p[12,6] >= -1
 OverlapInBin_12_6_toBin_1: - p_ij[6,1] - p_ij[12,1] + x_p[6,12]
   + x_p[12,6] + z_p[6,12] + z_p[12,6] >= -1
 OverlapInBin_12_6_toBin_2: - p_ij[6,2] - p_ij[12,2] + x_p[6,12]
   + x_p[12,6] + z_p[6,12] + z_p[12,6] >= -1
 OverlapInBin_12_6_toBin_3: - p_ij[6,3] - p_ij[12,3] + x_p[6,12]
   + x_p[12,6] + z_p[6,12] + z_p[12,6] >= -1
 OverlapInBin_12_7_toBin_0: - p_ij[7,0] - p_ij[12,0] + x_p[7,12]
   + x_p[12,7] + z_p[7,12] + z_p[12,7] >= -1
 OverlapInBin_12_7_toBin_1: - p_ij[7,1] - p_ij[12,1] + x_p[7,12]
   + x_p[12,7] + z_p[7,12] + z_p[12,7] >= -1
 OverlapInBin_12_7_toBin_2: - p_ij[7,2] - p_ij[12,2] + x_p[7,12]
   + x_p[12,7] + z_p[7,12] + z_p[12,7] >= -1
 OverlapInBin_12_7_toBin_3: - p_ij[7,3] - p_ij[12,3] + x_p[7,12]
   + x_p[12,7] + z_p[7,12] + z_p[12,7] >= -1
 OverlapInBin_12_8_toBin_0: - p_ij[8,0] - p_ij[12,0] + x_p[8,12]
   + x_p[12,8] + z_p[8,12] + z_p[12,8] >= -1
 OverlapInBin_12_8_toBin_1: - p_ij[8,1] - p_ij[12,1] + x_p[8,12]
   + x_p[12,8] + z_p[8,12] + z_p[12,8] >= -1
 OverlapInBin_12_8_toBin_2: - p_ij[8,2] - p_ij[12,2] + x_p[8,12]
   + x_p[12,8] + z_p[8,12] + z_p[12,8] >= -1
 OverlapInBin_12_8_toBin_3: - p_ij[8,3] - p_ij[12,3] + x_p[8,12]
   + x_p[12,8] + z_p[8,12] + z_p[12,8] >= -1
 OverlapInBin_12_9_toBin_0: - p_ij[9,0] - p_ij[12,0] + x_p[9,12]
   + x_p[12,9] + z_p[9,12] + z_p[12,9] >= -1
 OverlapInBin_12_9_toBin_1: - p_ij[9,1] - p_ij[12,1] + x_p[9,12]
   + x_p[12,9] + z_p[9,12] + z_p[12,9] >= -1
 OverlapInBin_12_9_toBin_2: - p_ij[9,2] - p_ij[12,2] + x_p[9,12]
   + x_p[12,9] + z_p[9,12] + z_p[12,9] >= -1
 OverlapInBin_12_9_toBin_3: - p_ij[9,3] - p_ij[12,3] + x_p[9,12]
   + x_p[12,9] + z_p[9,12] + z_p[12,9] >= -1
 OverlapInBin_12_10_toBin_0: - p_ij[10,0] - p_ij[12,0] + x_p[10,12]
   + x_p[12,10] + z_p[10,12] + z_p[12,10] >= -1
 OverlapInBin_12_10_toBin_1: - p_ij[10,1] - p_ij[12,1] + x_p[10,12]
   + x_p[12,10] + z_p[10,12] + z_p[12,10] >= -1
 OverlapInBin_12_10_toBin_2: - p_ij[10,2] - p_ij[12,2] + x_p[10,12]
   + x_p[12,10] + z_p[10,12] + z_p[12,10] >= -1
 OverlapInBin_12_10_toBin_3: - p_ij[10,3] - p_ij[12,3] + x_p[10,12]
   + x_p[12,10] + z_p[10,12] + z_p[12,10] >= -1
 OverlapInBin_12_11_toBin_0: - p_ij[11,0] - p_ij[12,0] + x_p[11,12]
   + x_p[12,11] + z_p[11,12] + z_p[12,11] >= -1
 OverlapInBin_12_11_toBin_1: - p_ij[11,1] - p_ij[12,1] + x_p[11,12]
   + x_p[12,11] + z_p[11,12] + z_p[12,11] >= -1
 OverlapInBin_12_11_toBin_2: - p_ij[11,2] - p_ij[12,2] + x_p[11,12]
   + x_p[12,11] + z_p[11,12] + z_p[12,11] >= -1
 OverlapInBin_12_11_toBin_3: - p_ij[11,3] - p_ij[12,3] + x_p[11,12]
   + x_p[12,11] + z_p[11,12] + z_p[12,11] >= -1
 Overlap_0_ToRightOf_1: 300 x_p[0,1] - xi_l[0] + xi_r[1] <= 300
 Overlap_0NotToRightOf_1: - 300 x_p[0,1] + xi_l[0] - xi_r[1] <= -1
 Overlap_0_Above_1: 155 z_p[0,1] - zi[0] + zi_hi[1] <= 155
 Overlap_0_ToRightOf_2: 300 x_p[0,2] - xi_l[0] + xi_r[2] <= 300
 Overlap_0NotToRightOf_2: - 300 x_p[0,2] + xi_l[0] - xi_r[2] <= -1
 Overlap_0_Above_2: 155 z_p[0,2] - zi[0] + zi_hi[2] <= 155
 Overlap_0_ToRightOf_3: 300 x_p[0,3] - xi_l[0] + xi_r[3] <= 300
 Overlap_0NotToRightOf_3: - 300 x_p[0,3] + xi_l[0] - xi_r[3] <= -1
 Overlap_0_Above_3: 155 z_p[0,3] - zi[0] + zi_hi[3] <= 155
 Overlap_0_ToRightOf_4: 300 x_p[0,4] - xi_l[0] + xi_r[4] <= 300
 Overlap_0NotToRightOf_4: - 300 x_p[0,4] + xi_l[0] - xi_r[4] <= -1
 Overlap_0_Above_4: 155 z_p[0,4] - zi[0] + zi_hi[4] <= 155
 Overlap_0_ToRightOf_5: 300 x_p[0,5] - xi_l[0] + xi_r[5] <= 300
 Overlap_0NotToRightOf_5: - 300 x_p[0,5] + xi_l[0] - xi_r[5] <= -1
 Overlap_0_Above_5: 155 z_p[0,5] - zi[0] + zi_hi[5] <= 155
 Overlap_0_ToRightOf_6: 300 x_p[0,6] - xi_l[0] + xi_r[6] <= 300
 Overlap_0NotToRightOf_6: - 300 x_p[0,6] + xi_l[0] - xi_r[6] <= -1
 Overlap_0_Above_6: 155 z_p[0,6] - zi[0] + zi_hi[6] <= 155
 Overlap_0_ToRightOf_7: 300 x_p[0,7] - xi_l[0] + xi_r[7] <= 300
 Overlap_0NotToRightOf_7: - 300 x_p[0,7] + xi_l[0] - xi_r[7] <= -1
 Overlap_0_Above_7: 155 z_p[0,7] - zi[0] + zi_hi[7] <= 155
 Overlap_0_ToRightOf_8: 300 x_p[0,8] - xi_l[0] + xi_r[8] <= 300
 Overlap_0NotToRightOf_8: - 300 x_p[0,8] + xi_l[0] - xi_r[8] <= -1
 Overlap_0_Above_8: 155 z_p[0,8] - zi[0] + zi_hi[8] <= 155
 Overlap_0_ToRightOf_9: 300 x_p[0,9] - xi_l[0] + xi_r[9] <= 300
 Overlap_0NotToRightOf_9: - 300 x_p[0,9] + xi_l[0] - xi_r[9] <= -1
 Overlap_0_Above_9: 155 z_p[0,9] - zi[0] + zi_hi[9] <= 155
 Overlap_0_ToRightOf_10: 300 x_p[0,10] - xi_l[0] + xi_r[10] <= 300
 Overlap_0NotToRightOf_10: - 300 x_p[0,10] + xi_l[0] - xi_r[10] <= -1
 Overlap_0_Above_10: 155 z_p[0,10] - zi[0] + zi_hi[10] <= 155
 Overlap_0_ToRightOf_11: 300 x_p[0,11] - xi_l[0] + xi_r[11] <= 300
 Overlap_0NotToRightOf_11: - 300 x_p[0,11] + xi_l[0] - xi_r[11] <= -1
 Overlap_0_Above_11: 155 z_p[0,11] - zi[0] + zi_hi[11] <= 155
 Overlap_0_ToRightOf_12: 300 x_p[0,12] - xi_l[0] + xi_r[12] <= 300
 Overlap_0NotToRightOf_12: - 300 x_p[0,12] + xi_l[0] - xi_r[12] <= -1
 Overlap_0_Above_12: 155 z_p[0,12] - zi[0] + zi_hi[12] <= 155
 Overlap_1_ToRightOf_0: 300 x_p[1,0] - xi_l[1] + xi_r[0] <= 300
 Overlap_1NotToRightOf_0: - 300 x_p[1,0] + xi_l[1] - xi_r[0] <= -1
 Overlap_1_Above_0: 155 z_p[1,0] - zi[1] + zi_hi[0] <= 155
 Overlap_1_ToRightOf_2: 300 x_p[1,2] - xi_l[1] + xi_r[2] <= 300
 Overlap_1NotToRightOf_2: - 300 x_p[1,2] + xi_l[1] - xi_r[2] <= -1
 Overlap_1_Above_2: 155 z_p[1,2] - zi[1] + zi_hi[2] <= 155
 Overlap_1_ToRightOf_3: 300 x_p[1,3] - xi_l[1] + xi_r[3] <= 300
 Overlap_1NotToRightOf_3: - 300 x_p[1,3] + xi_l[1] - xi_r[3] <= -1
 Overlap_1_Above_3: 155 z_p[1,3] - zi[1] + zi_hi[3] <= 155
 Overlap_1_ToRightOf_4: 300 x_p[1,4] - xi_l[1] + xi_r[4] <= 300
 Overlap_1NotToRightOf_4: - 300 x_p[1,4] + xi_l[1] - xi_r[4] <= -1
 Overlap_1_Above_4: 155 z_p[1,4] - zi[1] + zi_hi[4] <= 155
 Overlap_1_ToRightOf_5: 300 x_p[1,5] - xi_l[1] + xi_r[5] <= 300
 Overlap_1NotToRightOf_5: - 300 x_p[1,5] + xi_l[1] - xi_r[5] <= -1
 Overlap_1_Above_5: 155 z_p[1,5] - zi[1] + zi_hi[5] <= 155
 Overlap_1_ToRightOf_6: 300 x_p[1,6] - xi_l[1] + xi_r[6] <= 300
 Overlap_1NotToRightOf_6: - 300 x_p[1,6] + xi_l[1] - xi_r[6] <= -1
 Overlap_1_Above_6: 155 z_p[1,6] - zi[1] + zi_hi[6] <= 155
 Overlap_1_ToRightOf_7: 300 x_p[1,7] - xi_l[1] + xi_r[7] <= 300
 Overlap_1NotToRightOf_7: - 300 x_p[1,7] + xi_l[1] - xi_r[7] <= -1
 Overlap_1_Above_7: 155 z_p[1,7] - zi[1] + zi_hi[7] <= 155
 Overlap_1_ToRightOf_8: 300 x_p[1,8] - xi_l[1] + xi_r[8] <= 300
 Overlap_1NotToRightOf_8: - 300 x_p[1,8] + xi_l[1] - xi_r[8] <= -1
 Overlap_1_Above_8: 155 z_p[1,8] - zi[1] + zi_hi[8] <= 155
 Overlap_1_ToRightOf_9: 300 x_p[1,9] - xi_l[1] + xi_r[9] <= 300
 Overlap_1NotToRightOf_9: - 300 x_p[1,9] + xi_l[1] - xi_r[9] <= -1
 Overlap_1_Above_9: 155 z_p[1,9] - zi[1] + zi_hi[9] <= 155
 Overlap_1_ToRightOf_10: 300 x_p[1,10] - xi_l[1] + xi_r[10] <= 300
 Overlap_1NotToRightOf_10: - 300 x_p[1,10] + xi_l[1] - xi_r[10] <= -1
 Overlap_1_Above_10: 155 z_p[1,10] - zi[1] + zi_hi[10] <= 155
 Overlap_1_ToRightOf_11: 300 x_p[1,11] - xi_l[1] + xi_r[11] <= 300
 Overlap_1NotToRightOf_11: - 300 x_p[1,11] + xi_l[1] - xi_r[11] <= -1
 Overlap_1_Above_11: 155 z_p[1,11] - zi[1] + zi_hi[11] <= 155
 Overlap_1_ToRightOf_12: 300 x_p[1,12] - xi_l[1] + xi_r[12] <= 300
 Overlap_1NotToRightOf_12: - 300 x_p[1,12] + xi_l[1] - xi_r[12] <= -1
 Overlap_1_Above_12: 155 z_p[1,12] - zi[1] + zi_hi[12] <= 155
 Overlap_2_ToRightOf_0: 300 x_p[2,0] - xi_l[2] + xi_r[0] <= 300
 Overlap_2NotToRightOf_0: - 300 x_p[2,0] + xi_l[2] - xi_r[0] <= -1
 Overlap_2_Above_0: 155 z_p[2,0] - zi[2] + zi_hi[0] <= 155
 Overlap_2_ToRightOf_1: 300 x_p[2,1] - xi_l[2] + xi_r[1] <= 300
 Overlap_2NotToRightOf_1: - 300 x_p[2,1] + xi_l[2] - xi_r[1] <= -1
 Overlap_2_Above_1: 155 z_p[2,1] - zi[2] + zi_hi[1] <= 155
 Overlap_2_ToRightOf_3: 300 x_p[2,3] - xi_l[2] + xi_r[3] <= 300
 Overlap_2NotToRightOf_3: - 300 x_p[2,3] + xi_l[2] - xi_r[3] <= -1
 Overlap_2_Above_3: 155 z_p[2,3] - zi[2] + zi_hi[3] <= 155
 Overlap_2_ToRightOf_4: 300 x_p[2,4] - xi_l[2] + xi_r[4] <= 300
 Overlap_2NotToRightOf_4: - 300 x_p[2,4] + xi_l[2] - xi_r[4] <= -1
 Overlap_2_Above_4: 155 z_p[2,4] - zi[2] + zi_hi[4] <= 155
 Overlap_2_ToRightOf_5: 300 x_p[2,5] - xi_l[2] + xi_r[5] <= 300
 Overlap_2NotToRightOf_5: - 300 x_p[2,5] + xi_l[2] - xi_r[5] <= -1
 Overlap_2_Above_5: 155 z_p[2,5] - zi[2] + zi_hi[5] <= 155
 Overlap_2_ToRightOf_6: 300 x_p[2,6] - xi_l[2] + xi_r[6] <= 300
 Overlap_2NotToRightOf_6: - 300 x_p[2,6] + xi_l[2] - xi_r[6] <= -1
 Overlap_2_Above_6: 155 z_p[2,6] - zi[2] + zi_hi[6] <= 155
 Overlap_2_ToRightOf_7: 300 x_p[2,7] - xi_l[2] + xi_r[7] <= 300
 Overlap_2NotToRightOf_7: - 300 x_p[2,7] + xi_l[2] - xi_r[7] <= -1
 Overlap_2_Above_7: 155 z_p[2,7] - zi[2] + zi_hi[7] <= 155
 Overlap_2_ToRightOf_8: 300 x_p[2,8] - xi_l[2] + xi_r[8] <= 300
 Overlap_2NotToRightOf_8: - 300 x_p[2,8] + xi_l[2] - xi_r[8] <= -1
 Overlap_2_Above_8: 155 z_p[2,8] - zi[2] + zi_hi[8] <= 155
 Overlap_2_ToRightOf_9: 300 x_p[2,9] - xi_l[2] + xi_r[9] <= 300
 Overlap_2NotToRightOf_9: - 300 x_p[2,9] + xi_l[2] - xi_r[9] <= -1
 Overlap_2_Above_9: 155 z_p[2,9] - zi[2] + zi_hi[9] <= 155
 Overlap_2_ToRightOf_10: 300 x_p[2,10] - xi_l[2] + xi_r[10] <= 300
 Overlap_2NotToRightOf_10: - 300 x_p[2,10] + xi_l[2] - xi_r[10] <= -1
 Overlap_2_Above_10: 155 z_p[2,10] - zi[2] + zi_hi[10] <= 155
 Overlap_2_ToRightOf_11: 300 x_p[2,11] - xi_l[2] + xi_r[11] <= 300
 Overlap_2NotToRightOf_11: - 300 x_p[2,11] + xi_l[2] - xi_r[11] <= -1
 Overlap_2_Above_11: 155 z_p[2,11] - zi[2] + zi_hi[11] <= 155
 Overlap_2_ToRightOf_12: 300 x_p[2,12] - xi_l[2] + xi_r[12] <= 300
 Overlap_2NotToRightOf_12: - 300 x_p[2,12] + xi_l[2] - xi_r[12] <= -1
 Overlap_2_Above_12: 155 z_p[2,12] - zi[2] + zi_hi[12] <= 155
 Overlap_3_ToRightOf_0: 300 x_p[3,0] - xi_l[3] + xi_r[0] <= 300
 Overlap_3NotToRightOf_0: - 300 x_p[3,0] + xi_l[3] - xi_r[0] <= -1
 Overlap_3_Above_0: 155 z_p[3,0] - zi[3] + zi_hi[0] <= 155
 Overlap_3_ToRightOf_1: 300 x_p[3,1] - xi_l[3] + xi_r[1] <= 300
 Overlap_3NotToRightOf_1: - 300 x_p[3,1] + xi_l[3] - xi_r[1] <= -1
 Overlap_3_Above_1: 155 z_p[3,1] - zi[3] + zi_hi[1] <= 155
 Overlap_3_ToRightOf_2: 300 x_p[3,2] - xi_l[3] + xi_r[2] <= 300
 Overlap_3NotToRightOf_2: - 300 x_p[3,2] + xi_l[3] - xi_r[2] <= -1
 Overlap_3_Above_2: 155 z_p[3,2] - zi[3] + zi_hi[2] <= 155
 Overlap_3_ToRightOf_4: 300 x_p[3,4] - xi_l[3] + xi_r[4] <= 300
 Overlap_3NotToRightOf_4: - 300 x_p[3,4] + xi_l[3] - xi_r[4] <= -1
 Overlap_3_Above_4: 155 z_p[3,4] - zi[3] + zi_hi[4] <= 155
 Overlap_3_ToRightOf_5: 300 x_p[3,5] - xi_l[3] + xi_r[5] <= 300
 Overlap_3NotToRightOf_5: - 300 x_p[3,5] + xi_l[3] - xi_r[5] <= -1
 Overlap_3_Above_5: 155 z_p[3,5] - zi[3] + zi_hi[5] <= 155
 Overlap_3_ToRightOf_6: 300 x_p[3,6] - xi_l[3] + xi_r[6] <= 300
 Overlap_3NotToRightOf_6: - 300 x_p[3,6] + xi_l[3] - xi_r[6] <= -1
 Overlap_3_Above_6: 155 z_p[3,6] - zi[3] + zi_hi[6] <= 155
 Overlap_3_ToRightOf_7: 300 x_p[3,7] - xi_l[3] + xi_r[7] <= 300
 Overlap_3NotToRightOf_7: - 300 x_p[3,7] + xi_l[3] - xi_r[7] <= -1
 Overlap_3_Above_7: 155 z_p[3,7] - zi[3] + zi_hi[7] <= 155
 Overlap_3_ToRightOf_8: 300 x_p[3,8] - xi_l[3] + xi_r[8] <= 300
 Overlap_3NotToRightOf_8: - 300 x_p[3,8] + xi_l[3] - xi_r[8] <= -1
 Overlap_3_Above_8: 155 z_p[3,8] - zi[3] + zi_hi[8] <= 155
 Overlap_3_ToRightOf_9: 300 x_p[3,9] - xi_l[3] + xi_r[9] <= 300
 Overlap_3NotToRightOf_9: - 300 x_p[3,9] + xi_l[3] - xi_r[9] <= -1
 Overlap_3_Above_9: 155 z_p[3,9] - zi[3] + zi_hi[9] <= 155
 Overlap_3_ToRightOf_10: 300 x_p[3,10] - xi_l[3] + xi_r[10] <= 300
 Overlap_3NotToRightOf_10: - 300 x_p[3,10] + xi_l[3] - xi_r[10] <= -1
 Overlap_3_Above_10: 155 z_p[3,10] - zi[3] + zi_hi[10] <= 155
 Overlap_3_ToRightOf_11: 300 x_p[3,11] - xi_l[3] + xi_r[11] <= 300
 Overlap_3NotToRightOf_11: - 300 x_p[3,11] + xi_l[3] - xi_r[11] <= -1
 Overlap_3_Above_11: 155 z_p[3,11] - zi[3] + zi_hi[11] <= 155
 Overlap_3_ToRightOf_12: 300 x_p[3,12] - xi_l[3] + xi_r[12] <= 300
 Overlap_3NotToRightOf_12: - 300 x_p[3,12] + xi_l[3] - xi_r[12] <= -1
 Overlap_3_Above_12: 155 z_p[3,12] - zi[3] + zi_hi[12] <= 155
 Overlap_4_ToRightOf_0: 300 x_p[4,0] - xi_l[4] + xi_r[0] <= 300
 Overlap_4NotToRightOf_0: - 300 x_p[4,0] + xi_l[4] - xi_r[0] <= -1
 Overlap_4_Above_0: 155 z_p[4,0] - zi[4] + zi_hi[0] <= 155
 Overlap_4_ToRightOf_1: 300 x_p[4,1] - xi_l[4] + xi_r[1] <= 300
 Overlap_4NotToRightOf_1: - 300 x_p[4,1] + xi_l[4] - xi_r[1] <= -1
 Overlap_4_Above_1: 155 z_p[4,1] - zi[4] + zi_hi[1] <= 155
 Overlap_4_ToRightOf_2: 300 x_p[4,2] - xi_l[4] + xi_r[2] <= 300
 Overlap_4NotToRightOf_2: - 300 x_p[4,2] + xi_l[4] - xi_r[2] <= -1
 Overlap_4_Above_2: 155 z_p[4,2] - zi[4] + zi_hi[2] <= 155
 Overlap_4_ToRightOf_3: 300 x_p[4,3] - xi_l[4] + xi_r[3] <= 300
 Overlap_4NotToRightOf_3: - 300 x_p[4,3] + xi_l[4] - xi_r[3] <= -1
 Overlap_4_Above_3: 155 z_p[4,3] - zi[4] + zi_hi[3] <= 155
 Overlap_4_ToRightOf_5: 300 x_p[4,5] - xi_l[4] + xi_r[5] <= 300
 Overlap_4NotToRightOf_5: - 300 x_p[4,5] + xi_l[4] - xi_r[5] <= -1
 Overlap_4_Above_5: 155 z_p[4,5] - zi[4] + zi_hi[5] <= 155
 Overlap_4_ToRightOf_6: 300 x_p[4,6] - xi_l[4] + xi_r[6] <= 300
 Overlap_4NotToRightOf_6: - 300 x_p[4,6] + xi_l[4] - xi_r[6] <= -1
 Overlap_4_Above_6: 155 z_p[4,6] - zi[4] + zi_hi[6] <= 155
 Overlap_4_ToRightOf_7: 300 x_p[4,7] - xi_l[4] + xi_r[7] <= 300
 Overlap_4NotToRightOf_7: - 300 x_p[4,7] + xi_l[4] - xi_r[7] <= -1
 Overlap_4_Above_7: 155 z_p[4,7] - zi[4] + zi_hi[7] <= 155
 Overlap_4_ToRightOf_8: 300 x_p[4,8] - xi_l[4] + xi_r[8] <= 300
 Overlap_4NotToRightOf_8: - 300 x_p[4,8] + xi_l[4] - xi_r[8] <= -1
 Overlap_4_Above_8: 155 z_p[4,8] - zi[4] + zi_hi[8] <= 155
 Overlap_4_ToRightOf_9: 300 x_p[4,9] - xi_l[4] + xi_r[9] <= 300
 Overlap_4NotToRightOf_9: - 300 x_p[4,9] + xi_l[4] - xi_r[9] <= -1
 Overlap_4_Above_9: 155 z_p[4,9] - zi[4] + zi_hi[9] <= 155
 Overlap_4_ToRightOf_10: 300 x_p[4,10] - xi_l[4] + xi_r[10] <= 300
 Overlap_4NotToRightOf_10: - 300 x_p[4,10] + xi_l[4] - xi_r[10] <= -1
 Overlap_4_Above_10: 155 z_p[4,10] - zi[4] + zi_hi[10] <= 155
 Overlap_4_ToRightOf_11: 300 x_p[4,11] - xi_l[4] + xi_r[11] <= 300
 Overlap_4NotToRightOf_11: - 300 x_p[4,11] + xi_l[4] - xi_r[11] <= -1
 Overlap_4_Above_11: 155 z_p[4,11] - zi[4] + zi_hi[11] <= 155
 Overlap_4_ToRightOf_12: 300 x_p[4,12] - xi_l[4] + xi_r[12] <= 300
 Overlap_4NotToRightOf_12: - 300 x_p[4,12] + xi_l[4] - xi_r[12] <= -1
 Overlap_4_Above_12: 155 z_p[4,12] - zi[4] + zi_hi[12] <= 155
 Overlap_5_ToRightOf_0: 300 x_p[5,0] - xi_l[5] + xi_r[0] <= 300
 Overlap_5NotToRightOf_0: - 300 x_p[5,0] + xi_l[5] - xi_r[0] <= -1
 Overlap_5_Above_0: 155 z_p[5,0] - zi[5] + zi_hi[0] <= 155
 Overlap_5_ToRightOf_1: 300 x_p[5,1] - xi_l[5] + xi_r[1] <= 300
 Overlap_5NotToRightOf_1: - 300 x_p[5,1] + xi_l[5] - xi_r[1] <= -1
 Overlap_5_Above_1: 155 z_p[5,1] - zi[5] + zi_hi[1] <= 155
 Overlap_5_ToRightOf_2: 300 x_p[5,2] - xi_l[5] + xi_r[2] <= 300
 Overlap_5NotToRightOf_2: - 300 x_p[5,2] + xi_l[5] - xi_r[2] <= -1
 Overlap_5_Above_2: 155 z_p[5,2] - zi[5] + zi_hi[2] <= 155
 Overlap_5_ToRightOf_3: 300 x_p[5,3] - xi_l[5] + xi_r[3] <= 300
 Overlap_5NotToRightOf_3: - 300 x_p[5,3] + xi_l[5] - xi_r[3] <= -1
 Overlap_5_Above_3: 155 z_p[5,3] - zi[5] + zi_hi[3] <= 155
 Overlap_5_ToRightOf_4: 300 x_p[5,4] - xi_l[5] + xi_r[4] <= 300
 Overlap_5NotToRightOf_4: - 300 x_p[5,4] + xi_l[5] - xi_r[4] <= -1
 Overlap_5_Above_4: 155 z_p[5,4] - zi[5] + zi_hi[4] <= 155
 Overlap_5_ToRightOf_6: 300 x_p[5,6] - xi_l[5] + xi_r[6] <= 300
 Overlap_5NotToRightOf_6: - 300 x_p[5,6] + xi_l[5] - xi_r[6] <= -1
 Overlap_5_Above_6: 155 z_p[5,6] - zi[5] + zi_hi[6] <= 155
 Overlap_5_ToRightOf_7: 300 x_p[5,7] - xi_l[5] + xi_r[7] <= 300
 Overlap_5NotToRightOf_7: - 300 x_p[5,7] + xi_l[5] - xi_r[7] <= -1
 Overlap_5_Above_7: 155 z_p[5,7] - zi[5] + zi_hi[7] <= 155
 Overlap_5_ToRightOf_8: 300 x_p[5,8] - xi_l[5] + xi_r[8] <= 300
 Overlap_5NotToRightOf_8: - 300 x_p[5,8] + xi_l[5] - xi_r[8] <= -1
 Overlap_5_Above_8: 155 z_p[5,8] - zi[5] + zi_hi[8] <= 155
 Overlap_5_ToRightOf_9: 300 x_p[5,9] - xi_l[5] + xi_r[9] <= 300
 Overlap_5NotToRightOf_9: - 300 x_p[5,9] + xi_l[5] - xi_r[9] <= -1
 Overlap_5_Above_9: 155 z_p[5,9] - zi[5] + zi_hi[9] <= 155
 Overlap_5_ToRightOf_10: 300 x_p[5,10] - xi_l[5] + xi_r[10] <= 300
 Overlap_5NotToRightOf_10: - 300 x_p[5,10] + xi_l[5] - xi_r[10] <= -1
 Overlap_5_Above_10: 155 z_p[5,10] - zi[5] + zi_hi[10] <= 155
 Overlap_5_ToRightOf_11: 300 x_p[5,11] - xi_l[5] + xi_r[11] <= 300
 Overlap_5NotToRightOf_11: - 300 x_p[5,11] + xi_l[5] - xi_r[11] <= -1
 Overlap_5_Above_11: 155 z_p[5,11] - zi[5] + zi_hi[11] <= 155
 Overlap_5_ToRightOf_12: 300 x_p[5,12] - xi_l[5] + xi_r[12] <= 300
 Overlap_5NotToRightOf_12: - 300 x_p[5,12] + xi_l[5] - xi_r[12] <= -1
 Overlap_5_Above_12: 155 z_p[5,12] - zi[5] + zi_hi[12] <= 155
 Overlap_6_ToRightOf_0: 300 x_p[6,0] - xi_l[6] + xi_r[0] <= 300
 Overlap_6NotToRightOf_0: - 300 x_p[6,0] + xi_l[6] - xi_r[0] <= -1
 Overlap_6_Above_0: 155 z_p[6,0] - zi[6] + zi_hi[0] <= 155
 Overlap_6_ToRightOf_1: 300 x_p[6,1] - xi_l[6] + xi_r[1] <= 300
 Overlap_6NotToRightOf_1: - 300 x_p[6,1] + xi_l[6] - xi_r[1] <= -1
 Overlap_6_Above_1: 155 z_p[6,1] - zi[6] + zi_hi[1] <= 155
 Overlap_6_ToRightOf_2: 300 x_p[6,2] - xi_l[6] + xi_r[2] <= 300
 Overlap_6NotToRightOf_2: - 300 x_p[6,2] + xi_l[6] - xi_r[2] <= -1
 Overlap_6_Above_2: 155 z_p[6,2] - zi[6] + zi_hi[2] <= 155
 Overlap_6_ToRightOf_3: 300 x_p[6,3] - xi_l[6] + xi_r[3] <= 300
 Overlap_6NotToRightOf_3: - 300 x_p[6,3] + xi_l[6] - xi_r[3] <= -1
 Overlap_6_Above_3: 155 z_p[6,3] - zi[6] + zi_hi[3] <= 155
 Overlap_6_ToRightOf_4: 300 x_p[6,4] - xi_l[6] + xi_r[4] <= 300
 Overlap_6NotToRightOf_4: - 300 x_p[6,4] + xi_l[6] - xi_r[4] <= -1
 Overlap_6_Above_4: 155 z_p[6,4] - zi[6] + zi_hi[4] <= 155
 Overlap_6_ToRightOf_5: 300 x_p[6,5] - xi_l[6] + xi_r[5] <= 300
 Overlap_6NotToRightOf_5: - 300 x_p[6,5] + xi_l[6] - xi_r[5] <= -1
 Overlap_6_Above_5: 155 z_p[6,5] - zi[6] + zi_hi[5] <= 155
 Overlap_6_ToRightOf_7: 300 x_p[6,7] - xi_l[6] + xi_r[7] <= 300
 Overlap_6NotToRightOf_7: - 300 x_p[6,7] + xi_l[6] - xi_r[7] <= -1
 Overlap_6_Above_7: 155 z_p[6,7] - zi[6] + zi_hi[7] <= 155
 Overlap_6_ToRightOf_8: 300 x_p[6,8] - xi_l[6] + xi_r[8] <= 300
 Overlap_6NotToRightOf_8: - 300 x_p[6,8] + xi_l[6] - xi_r[8] <= -1
 Overlap_6_Above_8: 155 z_p[6,8] - zi[6] + zi_hi[8] <= 155
 Overlap_6_ToRightOf_9: 300 x_p[6,9] - xi_l[6] + xi_r[9] <= 300
 Overlap_6NotToRightOf_9: - 300 x_p[6,9] + xi_l[6] - xi_r[9] <= -1
 Overlap_6_Above_9: 155 z_p[6,9] - zi[6] + zi_hi[9] <= 155
 Overlap_6_ToRightOf_10: 300 x_p[6,10] - xi_l[6] + xi_r[10] <= 300
 Overlap_6NotToRightOf_10: - 300 x_p[6,10] + xi_l[6] - xi_r[10] <= -1
 Overlap_6_Above_10: 155 z_p[6,10] - zi[6] + zi_hi[10] <= 155
 Overlap_6_ToRightOf_11: 300 x_p[6,11] - xi_l[6] + xi_r[11] <= 300
 Overlap_6NotToRightOf_11: - 300 x_p[6,11] + xi_l[6] - xi_r[11] <= -1
 Overlap_6_Above_11: 155 z_p[6,11] - zi[6] + zi_hi[11] <= 155
 Overlap_6_ToRightOf_12: 300 x_p[6,12] - xi_l[6] + xi_r[12] <= 300
 Overlap_6NotToRightOf_12: - 300 x_p[6,12] + xi_l[6] - xi_r[12] <= -1
 Overlap_6_Above_12: 155 z_p[6,12] - zi[6] + zi_hi[12] <= 155
 Overlap_7_ToRightOf_0: 300 x_p[7,0] - xi_l[7] + xi_r[0] <= 300
 Overlap_7NotToRightOf_0: - 300 x_p[7,0] + xi_l[7] - xi_r[0] <= -1
 Overlap_7_Above_0: 155 z_p[7,0] - zi[7] + zi_hi[0] <= 155
 Overlap_7_ToRightOf_1: 300 x_p[7,1] - xi_l[7] + xi_r[1] <= 300
 Overlap_7NotToRightOf_1: - 300 x_p[7,1] + xi_l[7] - xi_r[1] <= -1
 Overlap_7_Above_1: 155 z_p[7,1] - zi[7] + zi_hi[1] <= 155
 Overlap_7_ToRightOf_2: 300 x_p[7,2] - xi_l[7] + xi_r[2] <= 300
 Overlap_7NotToRightOf_2: - 300 x_p[7,2] + xi_l[7] - xi_r[2] <= -1
 Overlap_7_Above_2: 155 z_p[7,2] - zi[7] + zi_hi[2] <= 155
 Overlap_7_ToRightOf_3: 300 x_p[7,3] - xi_l[7] + xi_r[3] <= 300
 Overlap_7NotToRightOf_3: - 300 x_p[7,3] + xi_l[7] - xi_r[3] <= -1
 Overlap_7_Above_3: 155 z_p[7,3] - zi[7] + zi_hi[3] <= 155
 Overlap_7_ToRightOf_4: 300 x_p[7,4] - xi_l[7] + xi_r[4] <= 300
 Overlap_7NotToRightOf_4: - 300 x_p[7,4] + xi_l[7] - xi_r[4] <= -1
 Overlap_7_Above_4: 155 z_p[7,4] - zi[7] + zi_hi[4] <= 155
 Overlap_7_ToRightOf_5: 300 x_p[7,5] - xi_l[7] + xi_r[5] <= 300
 Overlap_7NotToRightOf_5: - 300 x_p[7,5] + xi_l[7] - xi_r[5] <= -1
 Overlap_7_Above_5: 155 z_p[7,5] - zi[7] + zi_hi[5] <= 155
 Overlap_7_ToRightOf_6: 300 x_p[7,6] - xi_l[7] + xi_r[6] <= 300
 Overlap_7NotToRightOf_6: - 300 x_p[7,6] + xi_l[7] - xi_r[6] <= -1
 Overlap_7_Above_6: 155 z_p[7,6] - zi[7] + zi_hi[6] <= 155
 Overlap_7_ToRightOf_8: 300 x_p[7,8] - xi_l[7] + xi_r[8] <= 300
 Overlap_7NotToRightOf_8: - 300 x_p[7,8] + xi_l[7] - xi_r[8] <= -1
 Overlap_7_Above_8: 155 z_p[7,8] - zi[7] + zi_hi[8] <= 155
 Overlap_7_ToRightOf_9: 300 x_p[7,9] - xi_l[7] + xi_r[9] <= 300
 Overlap_7NotToRightOf_9: - 300 x_p[7,9] + xi_l[7] - xi_r[9] <= -1
 Overlap_7_Above_9: 155 z_p[7,9] - zi[7] + zi_hi[9] <= 155
 Overlap_7_ToRightOf_10: 300 x_p[7,10] - xi_l[7] + xi_r[10] <= 300
 Overlap_7NotToRightOf_10: - 300 x_p[7,10] + xi_l[7] - xi_r[10] <= -1
 Overlap_7_Above_10: 155 z_p[7,10] - zi[7] + zi_hi[10] <= 155
 Overlap_7_ToRightOf_11: 300 x_p[7,11] - xi_l[7] + xi_r[11] <= 300
 Overlap_7NotToRightOf_11: - 300 x_p[7,11] + xi_l[7] - xi_r[11] <= -1
 Overlap_7_Above_11: 155 z_p[7,11] - zi[7] + zi_hi[11] <= 155
 Overlap_7_ToRightOf_12: 300 x_p[7,12] - xi_l[7] + xi_r[12] <= 300
 Overlap_7NotToRightOf_12: - 300 x_p[7,12] + xi_l[7] - xi_r[12] <= -1
 Overlap_7_Above_12: 155 z_p[7,12] - zi[7] + zi_hi[12] <= 155
 Overlap_8_ToRightOf_0: 300 x_p[8,0] - xi_l[8] + xi_r[0] <= 300
 Overlap_8NotToRightOf_0: - 300 x_p[8,0] + xi_l[8] - xi_r[0] <= -1
 Overlap_8_Above_0: 155 z_p[8,0] - zi[8] + zi_hi[0] <= 155
 Overlap_8_ToRightOf_1: 300 x_p[8,1] - xi_l[8] + xi_r[1] <= 300
 Overlap_8NotToRightOf_1: - 300 x_p[8,1] + xi_l[8] - xi_r[1] <= -1
 Overlap_8_Above_1: 155 z_p[8,1] - zi[8] + zi_hi[1] <= 155
 Overlap_8_ToRightOf_2: 300 x_p[8,2] - xi_l[8] + xi_r[2] <= 300
 Overlap_8NotToRightOf_2: - 300 x_p[8,2] + xi_l[8] - xi_r[2] <= -1
 Overlap_8_Above_2: 155 z_p[8,2] - zi[8] + zi_hi[2] <= 155
 Overlap_8_ToRightOf_3: 300 x_p[8,3] - xi_l[8] + xi_r[3] <= 300
 Overlap_8NotToRightOf_3: - 300 x_p[8,3] + xi_l[8] - xi_r[3] <= -1
 Overlap_8_Above_3: 155 z_p[8,3] - zi[8] + zi_hi[3] <= 155
 Overlap_8_ToRightOf_4: 300 x_p[8,4] - xi_l[8] + xi_r[4] <= 300
 Overlap_8NotToRightOf_4: - 300 x_p[8,4] + xi_l[8] - xi_r[4] <= -1
 Overlap_8_Above_4: 155 z_p[8,4] - zi[8] + zi_hi[4] <= 155
 Overlap_8_ToRightOf_5: 300 x_p[8,5] - xi_l[8] + xi_r[5] <= 300
 Overlap_8NotToRightOf_5: - 300 x_p[8,5] + xi_l[8] - xi_r[5] <= -1
 Overlap_8_Above_5: 155 z_p[8,5] - zi[8] + zi_hi[5] <= 155
 Overlap_8_ToRightOf_6: 300 x_p[8,6] - xi_l[8] + xi_r[6] <= 300
 Overlap_8NotToRightOf_6: - 300 x_p[8,6] + xi_l[8] - xi_r[6] <= -1
 Overlap_8_Above_6: 155 z_p[8,6] - zi[8] + zi_hi[6] <= 155
 Overlap_8_ToRightOf_7: 300 x_p[8,7] - xi_l[8] + xi_r[7] <= 300
 Overlap_8NotToRightOf_7: - 300 x_p[8,7] + xi_l[8] - xi_r[7] <= -1
 Overlap_8_Above_7: 155 z_p[8,7] - zi[8] + zi_hi[7] <= 155
 Overlap_8_ToRightOf_9: 300 x_p[8,9] - xi_l[8] + xi_r[9] <= 300
 Overlap_8NotToRightOf_9: - 300 x_p[8,9] + xi_l[8] - xi_r[9] <= -1
 Overlap_8_Above_9: 155 z_p[8,9] - zi[8] + zi_hi[9] <= 155
 Overlap_8_ToRightOf_10: 300 x_p[8,10] - xi_l[8] + xi_r[10] <= 300
 Overlap_8NotToRightOf_10: - 300 x_p[8,10] + xi_l[8] - xi_r[10] <= -1
 Overlap_8_Above_10: 155 z_p[8,10] - zi[8] + zi_hi[10] <= 155
 Overlap_8_ToRightOf_11: 300 x_p[8,11] - xi_l[8] + xi_r[11] <= 300
 Overlap_8NotToRightOf_11: - 300 x_p[8,11] + xi_l[8] - xi_r[11] <= -1
 Overlap_8_Above_11: 155 z_p[8,11] - zi[8] + zi_hi[11] <= 155
 Overlap_8_ToRightOf_12: 300 x_p[8,12] - xi_l[8] + xi_r[12] <= 300
 Overlap_8NotToRightOf_12: - 300 x_p[8,12] + xi_l[8] - xi_r[12] <= -1
 Overlap_8_Above_12: 155 z_p[8,12] - zi[8] + zi_hi[12] <= 155
 Overlap_9_ToRightOf_0: 300 x_p[9,0] - xi_l[9] + xi_r[0] <= 300
 Overlap_9NotToRightOf_0: - 300 x_p[9,0] + xi_l[9] - xi_r[0] <= -1
 Overlap_9_Above_0: 155 z_p[9,0] - zi[9] + zi_hi[0] <= 155
 Overlap_9_ToRightOf_1: 300 x_p[9,1] - xi_l[9] + xi_r[1] <= 300
 Overlap_9NotToRightOf_1: - 300 x_p[9,1] + xi_l[9] - xi_r[1] <= -1
 Overlap_9_Above_1: 155 z_p[9,1] - zi[9] + zi_hi[1] <= 155
 Overlap_9_ToRightOf_2: 300 x_p[9,2] - xi_l[9] + xi_r[2] <= 300
 Overlap_9NotToRightOf_2: - 300 x_p[9,2] + xi_l[9] - xi_r[2] <= -1
 Overlap_9_Above_2: 155 z_p[9,2] - zi[9] + zi_hi[2] <= 155
 Overlap_9_ToRightOf_3: 300 x_p[9,3] - xi_l[9] + xi_r[3] <= 300
 Overlap_9NotToRightOf_3: - 300 x_p[9,3] + xi_l[9] - xi_r[3] <= -1
 Overlap_9_Above_3: 155 z_p[9,3] - zi[9] + zi_hi[3] <= 155
 Overlap_9_ToRightOf_4: 300 x_p[9,4] - xi_l[9] + xi_r[4] <= 300
 Overlap_9NotToRightOf_4: - 300 x_p[9,4] + xi_l[9] - xi_r[4] <= -1
 Overlap_9_Above_4: 155 z_p[9,4] - zi[9] + zi_hi[4] <= 155
 Overlap_9_ToRightOf_5: 300 x_p[9,5] - xi_l[9] + xi_r[5] <= 300
 Overlap_9NotToRightOf_5: - 300 x_p[9,5] + xi_l[9] - xi_r[5] <= -1
 Overlap_9_Above_5: 155 z_p[9,5] - zi[9] + zi_hi[5] <= 155
 Overlap_9_ToRightOf_6: 300 x_p[9,6] - xi_l[9] + xi_r[6] <= 300
 Overlap_9NotToRightOf_6: - 300 x_p[9,6] + xi_l[9] - xi_r[6] <= -1
 Overlap_9_Above_6: 155 z_p[9,6] - zi[9] + zi_hi[6] <= 155
 Overlap_9_ToRightOf_7: 300 x_p[9,7] - xi_l[9] + xi_r[7] <= 300
 Overlap_9NotToRightOf_7: - 300 x_p[9,7] + xi_l[9] - xi_r[7] <= -1
 Overlap_9_Above_7: 155 z_p[9,7] - zi[9] + zi_hi[7] <= 155
 Overlap_9_ToRightOf_8: 300 x_p[9,8] - xi_l[9] + xi_r[8] <= 300
 Overlap_9NotToRightOf_8: - 300 x_p[9,8] + xi_l[9] - xi_r[8] <= -1
 Overlap_9_Above_8: 155 z_p[9,8] - zi[9] + zi_hi[8] <= 155
 Overlap_9_ToRightOf_10: 300 x_p[9,10] - xi_l[9] + xi_r[10] <= 300
 Overlap_9NotToRightOf_10: - 300 x_p[9,10] + xi_l[9] - xi_r[10] <= -1
 Overlap_9_Above_10: 155 z_p[9,10] - zi[9] + zi_hi[10] <= 155
 Overlap_9_ToRightOf_11: 300 x_p[9,11] - xi_l[9] + xi_r[11] <= 300
 Overlap_9NotToRightOf_11: - 300 x_p[9,11] + xi_l[9] - xi_r[11] <= -1
 Overlap_9_Above_11: 155 z_p[9,11] - zi[9] + zi_hi[11] <= 155
 Overlap_9_ToRightOf_12: 300 x_p[9,12] - xi_l[9] + xi_r[12] <= 300
 Overlap_9NotToRightOf_12: - 300 x_p[9,12] + xi_l[9] - xi_r[12] <= -1
 Overlap_9_Above_12: 155 z_p[9,12] - zi[9] + zi_hi[12] <= 155
 Overlap_10_ToRightOf_0: 300 x_p[10,0] - xi_l[10] + xi_r[0] <= 300
 Overlap_10NotToRightOf_0: - 300 x_p[10,0] + xi_l[10] - xi_r[0] <= -1
 Overlap_10_Above_0: 155 z_p[10,0] - zi[10] + zi_hi[0] <= 155
 Overlap_10_ToRightOf_1: 300 x_p[10,1] - xi_l[10] + xi_r[1] <= 300
 Overlap_10NotToRightOf_1: - 300 x_p[10,1] + xi_l[10] - xi_r[1] <= -1
 Overlap_10_Above_1: 155 z_p[10,1] - zi[10] + zi_hi[1] <= 155
 Overlap_10_ToRightOf_2: 300 x_p[10,2] - xi_l[10] + xi_r[2] <= 300
 Overlap_10NotToRightOf_2: - 300 x_p[10,2] + xi_l[10] - xi_r[2] <= -1
 Overlap_10_Above_2: 155 z_p[10,2] - zi[10] + zi_hi[2] <= 155
 Overlap_10_ToRightOf_3: 300 x_p[10,3] - xi_l[10] + xi_r[3] <= 300
 Overlap_10NotToRightOf_3: - 300 x_p[10,3] + xi_l[10] - xi_r[3] <= -1
 Overlap_10_Above_3: 155 z_p[10,3] - zi[10] + zi_hi[3] <= 155
 Overlap_10_ToRightOf_4: 300 x_p[10,4] - xi_l[10] + xi_r[4] <= 300
 Overlap_10NotToRightOf_4: - 300 x_p[10,4] + xi_l[10] - xi_r[4] <= -1
 Overlap_10_Above_4: 155 z_p[10,4] - zi[10] + zi_hi[4] <= 155
 Overlap_10_ToRightOf_5: 300 x_p[10,5] - xi_l[10] + xi_r[5] <= 300
 Overlap_10NotToRightOf_5: - 300 x_p[10,5] + xi_l[10] - xi_r[5] <= -1
 Overlap_10_Above_5: 155 z_p[10,5] - zi[10] + zi_hi[5] <= 155
 Overlap_10_ToRightOf_6: 300 x_p[10,6] - xi_l[10] + xi_r[6] <= 300
 Overlap_10NotToRightOf_6: - 300 x_p[10,6] + xi_l[10] - xi_r[6] <= -1
 Overlap_10_Above_6: 155 z_p[10,6] - zi[10] + zi_hi[6] <= 155
 Overlap_10_ToRightOf_7: 300 x_p[10,7] - xi_l[10] + xi_r[7] <= 300
 Overlap_10NotToRightOf_7: - 300 x_p[10,7] + xi_l[10] - xi_r[7] <= -1
 Overlap_10_Above_7: 155 z_p[10,7] - zi[10] + zi_hi[7] <= 155
 Overlap_10_ToRightOf_8: 300 x_p[10,8] - xi_l[10] + xi_r[8] <= 300
 Overlap_10NotToRightOf_8: - 300 x_p[10,8] + xi_l[10] - xi_r[8] <= -1
 Overlap_10_Above_8: 155 z_p[10,8] - zi[10] + zi_hi[8] <= 155
 Overlap_10_ToRightOf_9: 300 x_p[10,9] - xi_l[10] + xi_r[9] <= 300
 Overlap_10NotToRightOf_9: - 300 x_p[10,9] + xi_l[10] - xi_r[9] <= -1
 Overlap_10_Above_9: 155 z_p[10,9] - zi[10] + zi_hi[9] <= 155
 Overlap_10_ToRightOf_11: 300 x_p[10,11] - xi_l[10] + xi_r[11] <= 300
 Overlap_10NotToRightOf_11: - 300 x_p[10,11] + xi_l[10] - xi_r[11] <= -1
 Overlap_10_Above_11: 155 z_p[10,11] - zi[10] + zi_hi[11] <= 155
 Overlap_10_ToRightOf_12: 300 x_p[10,12] - xi_l[10] + xi_r[12] <= 300
 Overlap_10NotToRightOf_12: - 300 x_p[10,12] + xi_l[10] - xi_r[12] <= -1
 Overlap_10_Above_12: 155 z_p[10,12] - zi[10] + zi_hi[12] <= 155
 Overlap_11_ToRightOf_0: 300 x_p[11,0] - xi_l[11] + xi_r[0] <= 300
 Overlap_11NotToRightOf_0: - 300 x_p[11,0] + xi_l[11] - xi_r[0] <= -1
 Overlap_11_Above_0: 155 z_p[11,0] - zi[11] + zi_hi[0] <= 155
 Overlap_11_ToRightOf_1: 300 x_p[11,1] - xi_l[11] + xi_r[1] <= 300
 Overlap_11NotToRightOf_1: - 300 x_p[11,1] + xi_l[11] - xi_r[1] <= -1
 Overlap_11_Above_1: 155 z_p[11,1] - zi[11] + zi_hi[1] <= 155
 Overlap_11_ToRightOf_2: 300 x_p[11,2] - xi_l[11] + xi_r[2] <= 300
 Overlap_11NotToRightOf_2: - 300 x_p[11,2] + xi_l[11] - xi_r[2] <= -1
 Overlap_11_Above_2: 155 z_p[11,2] - zi[11] + zi_hi[2] <= 155
 Overlap_11_ToRightOf_3: 300 x_p[11,3] - xi_l[11] + xi_r[3] <= 300
 Overlap_11NotToRightOf_3: - 300 x_p[11,3] + xi_l[11] - xi_r[3] <= -1
 Overlap_11_Above_3: 155 z_p[11,3] - zi[11] + zi_hi[3] <= 155
 Overlap_11_ToRightOf_4: 300 x_p[11,4] - xi_l[11] + xi_r[4] <= 300
 Overlap_11NotToRightOf_4: - 300 x_p[11,4] + xi_l[11] - xi_r[4] <= -1
 Overlap_11_Above_4: 155 z_p[11,4] - zi[11] + zi_hi[4] <= 155
 Overlap_11_ToRightOf_5: 300 x_p[11,5] - xi_l[11] + xi_r[5] <= 300
 Overlap_11NotToRightOf_5: - 300 x_p[11,5] + xi_l[11] - xi_r[5] <= -1
 Overlap_11_Above_5: 155 z_p[11,5] - zi[11] + zi_hi[5] <= 155
 Overlap_11_ToRightOf_6: 300 x_p[11,6] - xi_l[11] + xi_r[6] <= 300
 Overlap_11NotToRightOf_6: - 300 x_p[11,6] + xi_l[11] - xi_r[6] <= -1
 Overlap_11_Above_6: 155 z_p[11,6] - zi[11] + zi_hi[6] <= 155
 Overlap_11_ToRightOf_7: 300 x_p[11,7] - xi_l[11] + xi_r[7] <= 300
 Overlap_11NotToRightOf_7: - 300 x_p[11,7] + xi_l[11] - xi_r[7] <= -1
 Overlap_11_Above_7: 155 z_p[11,7] - zi[11] + zi_hi[7] <= 155
 Overlap_11_ToRightOf_8: 300 x_p[11,8] - xi_l[11] + xi_r[8] <= 300
 Overlap_11NotToRightOf_8: - 300 x_p[11,8] + xi_l[11] - xi_r[8] <= -1
 Overlap_11_Above_8: 155 z_p[11,8] - zi[11] + zi_hi[8] <= 155
 Overlap_11_ToRightOf_9: 300 x_p[11,9] - xi_l[11] + xi_r[9] <= 300
 Overlap_11NotToRightOf_9: - 300 x_p[11,9] + xi_l[11] - xi_r[9] <= -1
 Overlap_11_Above_9: 155 z_p[11,9] - zi[11] + zi_hi[9] <= 155
 Overlap_11_ToRightOf_10: 300 x_p[11,10] - xi_l[11] + xi_r[10] <= 300
 Overlap_11NotToRightOf_10: - 300 x_p[11,10] + xi_l[11] - xi_r[10] <= -1
 Overlap_11_Above_10: 155 z_p[11,10] - zi[11] + zi_hi[10] <= 155
 Overlap_11_ToRightOf_12: 300 x_p[11,12] - xi_l[11] + xi_r[12] <= 300
 Overlap_11NotToRightOf_12: - 300 x_p[11,12] + xi_l[11] - xi_r[12] <= -1
 Overlap_11_Above_12: 155 z_p[11,12] - zi[11] + zi_hi[12] <= 155
 Overlap_12_ToRightOf_0: 300 x_p[12,0] - xi_l[12] + xi_r[0] <= 300
 Overlap_12NotToRightOf_0: - 300 x_p[12,0] + xi_l[12] - xi_r[0] <= -1
 Overlap_12_Above_0: 155 z_p[12,0] - zi[12] + zi_hi[0] <= 155
 Overlap_12_ToRightOf_1: 300 x_p[12,1] - xi_l[12] + xi_r[1] <= 300
 Overlap_12NotToRightOf_1: - 300 x_p[12,1] + xi_l[12] - xi_r[1] <= -1
 Overlap_12_Above_1: 155 z_p[12,1] - zi[12] + zi_hi[1] <= 155
 Overlap_12_ToRightOf_2: 300 x_p[12,2] - xi_l[12] + xi_r[2] <= 300
 Overlap_12NotToRightOf_2: - 300 x_p[12,2] + xi_l[12] - xi_r[2] <= -1
 Overlap_12_Above_2: 155 z_p[12,2] - zi[12] + zi_hi[2] <= 155
 Overlap_12_ToRightOf_3: 300 x_p[12,3] - xi_l[12] + xi_r[3] <= 300
 Overlap_12NotToRightOf_3: - 300 x_p[12,3] + xi_l[12] - xi_r[3] <= -1
 Overlap_12_Above_3: 155 z_p[12,3] - zi[12] + zi_hi[3] <= 155
 Overlap_12_ToRightOf_4: 300 x_p[12,4] - xi_l[12] + xi_r[4] <= 300
 Overlap_12NotToRightOf_4: - 300 x_p[12,4] + xi_l[12] - xi_r[4] <= -1
 Overlap_12_Above_4: 155 z_p[12,4] - zi[12] + zi_hi[4] <= 155
 Overlap_12_ToRightOf_5: 300 x_p[12,5] - xi_l[12] + xi_r[5] <= 300
 Overlap_12NotToRightOf_5: - 300 x_p[12,5] + xi_l[12] - xi_r[5] <= -1
 Overlap_12_Above_5: 155 z_p[12,5] - zi[12] + zi_hi[5] <= 155
 Overlap_12_ToRightOf_6: 300 x_p[12,6] - xi_l[12] + xi_r[6] <= 300
 Overlap_12NotToRightOf_6: - 300 x_p[12,6] + xi_l[12] - xi_r[6] <= -1
 Overlap_12_Above_6: 155 z_p[12,6] - zi[12] + zi_hi[6] <= 155
 Overlap_12_ToRightOf_7: 300 x_p[12,7] - xi_l[12] + xi_r[7] <= 300
 Overlap_12NotToRightOf_7: - 300 x_p[12,7] + xi_l[12] - xi_r[7] <= -1
 Overlap_12_Above_7: 155 z_p[12,7] - zi[12] + zi_hi[7] <= 155
 Overlap_12_ToRightOf_8: 300 x_p[12,8] - xi_l[12] + xi_r[8] <= 300
 Overlap_12NotToRightOf_8: - 300 x_p[12,8] + xi_l[12] - xi_r[8] <= -1
 Overlap_12_Above_8: 155 z_p[12,8] - zi[12] + zi_hi[8] <= 155
 Overlap_12_ToRightOf_9: 300 x_p[12,9] - xi_l[12] + xi_r[9] <= 300
 Overlap_12NotToRightOf_9: - 300 x_p[12,9] + xi_l[12] - xi_r[9] <= -1
 Overlap_12_Above_9: 155 z_p[12,9] - zi[12] + zi_hi[9] <= 155
 Overlap_12_ToRightOf_10: 300 x_p[12,10] - xi_l[12] + xi_r[10] <= 300
 Overlap_12NotToRightOf_10: - 300 x_p[12,10] + xi_l[12] - xi_r[10] <= -1
 Overlap_12_Above_10: 155 z_p[12,10] - zi[12] + zi_hi[10] <= 155
 Overlap_12_ToRightOf_11: 300 x_p[12,11] - xi_l[12] + xi_r[11] <= 300
 Overlap_12NotToRightOf_11: - 300 x_p[12,11] + xi_l[12] - xi_r[11] <= -1
 Overlap_12_Above_11: 155 z_p[12,11] - zi[12] + zi_hi[11] <= 155
 OrientationLength_0: r[0,1,0] <= 1
 OrientationLength_1: r[1,1,0] <= 0
 OrientationLength_2: r[2,1,0] <= 1
 OrientationLength_3: r[3,1,0] <= 1
 OrientationLength_4: r[4,1,0] <= 1
 OrientationLength_5: r[5,1,0] <= 1
 OrientationLength_6: r[6,1,0] <= 1
 OrientationLength_7: r[7,1,0] <= 0
 OrientationLength_8: r[8,1,0] <= 1
 OrientationLength_9: r[9,1,0] <= 0
 OrientationLength_10: r[10,1,0] <= 1
 OrientationLength_11: r[11,1,0] <= 1
 OrientationLength_12: r[12,1,0] <= 1
 CutSupport_LowerBound_0_2: - 10000 p_ij[0,2] + 1.261904761904762 xi_l[0]
   + zi[0] >= -9947
 CutSupport_LowerBound_0_3: - 10000 p_ij[0,3] + 1.261904761904762 xi_l[0]
   + zi[0] >= -9947
 CutSupport_LowerBound_1_2: - 10000 p_ij[1,2] + 1.261904761904762 xi_l[1]
   + zi[1] >= -9947
 CutSupport_LowerBound_1_3: - 10000 p_ij[1,3] + 1.261904761904762 xi_l[1]
   + zi[1] >= -9947
 CutSupport_LowerBound_2_2: - 10000 p_ij[2,2] + 1.261904761904762 xi_l[2]
   + zi[2] >= -9947
 CutSupport_LowerBound_2_3: - 10000 p_ij[2,3] + 1.261904761904762 xi_l[2]
   + zi[2] >= -9947
 CutSupport_LowerBound_3_2: - 10000 p_ij[3,2] + 1.261904761904762 xi_l[3]
   + zi[3] >= -9947
 CutSupport_LowerBound_3_3: - 10000 p_ij[3,3] + 1.261904761904762 xi_l[3]
   + zi[3] >= -9947
 CutSupport_LowerBound_4_2: - 10000 p_ij[4,2] + 1.261904761904762 xi_l[4]
   + zi[4] >= -9947
 CutSupport_LowerBound_4_3: - 10000 p_ij[4,3] + 1.261904761904762 xi_l[4]
   + zi[4] >= -9947
 CutSupport_LowerBound_5_2: - 10000 p_ij[5,2] + 1.261904761904762 xi_l[5]
   + zi[5] >= -9947
 CutSupport_LowerBound_5_3: - 10000 p_ij[5,3] + 1.261904761904762 xi_l[5]
   + zi[5] >= -9947
 CutSupport_LowerBound_6_2: - 10000 p_ij[6,2] + 1.261904761904762 xi_l[6]
   + zi[6] >= -9947
 CutSupport_LowerBound_6_3: - 10000 p_ij[6,3] + 1.261904761904762 xi_l[6]
   + zi[6] >= -9947
 CutSupport_LowerBound_7_2: - 10000 p_ij[7,2] + 1.261904761904762 xi_l[7]
   + zi[7] >= -9947
 CutSupport_LowerBound_7_3: - 10000 p_ij[7,3] + 1.261904761904762 xi_l[7]
   + zi[7] >= -9947
 CutSupport_LowerBound_8_2: - 10000 p_ij[8,2] + 1.261904761904762 xi_l[8]
   + zi[8] >= -9947
 CutSupport_LowerBound_8_3: - 10000 p_ij[8,3] + 1.261904761904762 xi_l[8]
   + zi[8] >= -9947
 CutSupport_LowerBound_9_2: - 10000 p_ij[9,2] + 1.261904761904762 xi_l[9]
   + zi[9] >= -9947
 CutSupport_LowerBound_9_3: - 10000 p_ij[9,3] + 1.261904761904762 xi_l[9]
   + zi[9] >= -9947
 CutSupport_LowerBound_10_2: - 10000 p_ij[10,2]
   + 1.261904761904762 xi_l[10] + zi[10] >= -9947
 CutSupport_LowerBound_10_3: - 10000 p_ij[10,3]
   + 1.261904761904762 xi_l[10] + zi[10] >= -9947
 CutSupport_LowerBound_11_2: - 10000 p_ij[11,2]
   + 1.261904761904762 xi_l[11] + zi[11] >= -9947
 CutSupport_LowerBound_11_3: - 10000 p_ij[11,3]
   + 1.261904761904762 xi_l[11] + zi[11] >= -9947
 CutSupport_LowerBound_12_2: - 10000 p_ij[12,2]
   + 1.261904761904762 xi_l[12] + zi[12] >= -9947
 CutSupport_LowerBound_12_3: - 10000 p_ij[12,3]
   + 1.261904761904762 xi_l[12] + zi[12] >= -9947
 R1252: 2 g[0] + beta1[0,1] + beta1[0,2] + beta1[0,3] + beta1[0,4]
   + beta1[0,5] + beta1[0,6] + beta1[0,7] + beta1[0,8] + beta1[0,9]
   + beta1[0,10] + beta1[0,11] + beta1[0,12] + beta2[0,1] + beta2[0,2]
   + beta2[0,3] + beta2[0,4] + beta2[0,5] + beta2[0,6] + beta2[0,7]
   + beta2[0,8] + beta2[0,9] + beta2[0,10] + beta2[0,11] + beta2[0,12]
   + gamma[0] >= 2
 R1253: 2 g[1] + beta1[1,0] + beta1[1,2] + beta1[1,3] + beta1[1,4]
   + beta1[1,5] + beta1[1,6] + beta1[1,7] + beta1[1,8] + beta1[1,9]
   + beta1[1,10] + beta1[1,11] + beta1[1,12] + beta2[1,0] + beta2[1,2]
   + beta2[1,3] + beta2[1,4] + beta2[1,5] + beta2[1,6] + beta2[1,7]
   + beta2[1,8] + beta2[1,9] + beta2[1,10] + beta2[1,11] + beta2[1,12]
   + gamma[1] >= 2
 R1254: 2 g[2] + beta1[2,0] + beta1[2,1] + beta1[2,3] + beta1[2,4]
   + beta1[2,5] + beta1[2,6] + beta1[2,7] + beta1[2,8] + beta1[2,9]
   + beta1[2,10] + beta1[2,11] + beta1[2,12] + beta2[2,0] + beta2[2,1]
   + beta2[2,3] + beta2[2,4] + beta2[2,5] + beta2[2,6] + beta2[2,7]
   + beta2[2,8] + beta2[2,9] + beta2[2,10] + beta2[2,11] + beta2[2,12]
   + gamma[2] >= 2
 R1255: 2 g[3] + beta1[3,0] + beta1[3,1] + beta1[3,2] + beta1[3,4]
   + beta1[3,5] + beta1[3,6] + beta1[3,7] + beta1[3,8] + beta1[3,9]
   + beta1[3,10] + beta1[3,11] + beta1[3,12] + beta2[3,0] + beta2[3,1]
   + beta2[3,2] + beta2[3,4] + beta2[3,5] + beta2[3,6] + beta2[3,7]
   + beta2[3,8] + beta2[3,9] + beta2[3,10] + beta2[3,11] + beta2[3,12]
   + gamma[3] >= 2
 R1256: 2 g[4] + beta1[4,0] + beta1[4,1] + beta1[4,2] + beta1[4,3]
   + beta1[4,5] + beta1[4,6] + beta1[4,7] + beta1[4,8] + beta1[4,9]
   + beta1[4,10] + beta1[4,11] + beta1[4,12] + beta2[4,0] + beta2[4,1]
   + beta2[4,2] + beta2[4,3] + beta2[4,5] + beta2[4,6] + beta2[4,7]
   + beta2[4,8] + beta2[4,9] + beta2[4,10] + beta2[4,11] + beta2[4,12]
   + gamma[4] >= 2
 R1257: 2 g[5] + beta1[5,0] + beta1[5,1] + beta1[5,2] + beta1[5,3]
   + beta1[5,4] + beta1[5,6] + beta1[5,7] + beta1[5,8] + beta1[5,9]
   + beta1[5,10] + beta1[5,11] + beta1[5,12] + beta2[5,0] + beta2[5,1]
   + beta2[5,2] + beta2[5,3] + beta2[5,4] + beta2[5,6] + beta2[5,7]
   + beta2[5,8] + beta2[5,9] + beta2[5,10] + beta2[5,11] + beta2[5,12]
   + gamma[5] >= 2
 R1258: 2 g[6] + beta1[6,0] + beta1[6,1] + beta1[6,2] + beta1[6,3]
   + beta1[6,4] + beta1[6,5] + beta1[6,7] + beta1[6,8] + beta1[6,9]
   + beta1[6,10] + beta1[6,11] + beta1[6,12] + beta2[6,0] + beta2[6,1]
   + beta2[6,2] + beta2[6,3] + beta2[6,4] + beta2[6,5] + beta2[6,7]
   + beta2[6,8] + beta2[6,9] + beta2[6,10] + beta2[6,11] + beta2[6,12]
   + gamma[6] >= 2
 R1259: 2 g[7] + beta1[7,0] + beta1[7,1] + beta1[7,2] + beta1[7,3]
   + beta1[7,4] + beta1[7,5] + beta1[7,6] + beta1[7,8] + beta1[7,9]
   + beta1[7,10] + beta1[7,11] + beta1[7,12] + beta2[7,0] + beta2[7,1]
   + beta2[7,2] + beta2[7,3] + beta2[7,4] + beta2[7,5] + beta2[7,6]
   + beta2[7,8] + beta2[7,9] + beta2[7,10] + beta2[7,11] + beta2[7,12]
   + gamma[7] >= 2
 R1260: 2 g[8] + beta1[8,0] + beta1[8,1] + beta1[8,2] + beta1[8,3]
   + beta1[8,4] + beta1[8,5] + beta1[8,6] + beta1[8,7] + beta1[8,9]
   + beta1[8,10] + beta1[8,11] + beta1[8,12] + beta2[8,0] + beta2[8,1]
   + beta2[8,2] + beta2[8,3] + beta2[8,4] + beta2[8,5] + beta2[8,6]
   + beta2[8,7] + beta2[8,9] + beta2[8,10] + beta2[8,11] + beta2[8,12]
   + gamma[8] >= 2
 R1261: 2 g[9] + beta1[9,0] + beta1[9,1] + beta1[9,2] + beta1[9,3]
   + beta1[9,4] + beta1[9,5] + beta1[9,6] + beta1[9,7] + beta1[9,8]
   + beta1[9,10] + beta1[9,11] + beta1[9,12] + beta2[9,0] + beta2[9,1]
   + beta2[9,2] + beta2[9,3] + beta2[9,4] + beta2[9,5] + beta2[9,6]
   + beta2[9,7] + beta2[9,8] + beta2[9,10] + beta2[9,11] + beta2[9,12]
   + gamma[9] >= 2
 R1262: 2 g[10] + beta1[10,0] + beta1[10,1] + beta1[10,2] + beta1[10,3]
   + beta1[10,4] + beta1[10,5] + beta1[10,6] + beta1[10,7] + beta1[10,8]
   + beta1[10,9] + beta1[10,11] + beta1[10,12] + beta2[10,0] + beta2[10,1]
   + beta2[10,2] + beta2[10,3] + beta2[10,4] + beta2[10,5] + beta2[10,6]
   + beta2[10,7] + beta2[10,8] + beta2[10,9] + beta2[10,11] + beta2[10,12]
   + gamma[10] >= 2
 R1263: 2 g[11] + beta1[11,0] + beta1[11,1] + beta1[11,2] + beta1[11,3]
   + beta1[11,4] + beta1[11,5] + beta1[11,6] + beta1[11,7] + beta1[11,8]
   + beta1[11,9] + beta1[11,10] + beta1[11,12] + beta2[11,0] + beta2[11,1]
   + beta2[11,2] + beta2[11,3] + beta2[11,4] + beta2[11,5] + beta2[11,6]
   + beta2[11,7] + beta2[11,8] + beta2[11,9] + beta2[11,10] + beta2[11,12]
   + gamma[11] >= 2
 R1264: 2 g[12] + beta1[12,0] + beta1[12,1] + beta1[12,2] + beta1[12,3]
   + beta1[12,4] + beta1[12,5] + beta1[12,6] + beta1[12,7] + beta1[12,8]
   + beta1[12,9] + beta1[12,10] + beta1[12,11] + beta2[12,0] + beta2[12,1]
   + beta2[12,2] + beta2[12,3] + beta2[12,4] + beta2[12,5] + beta2[12,6]
   + beta2[12,7] + beta2[12,8] + beta2[12,9] + beta2[12,10] + beta2[12,11]
   + gamma[12] >= 2
 GroundConstraint_0: zi[0] + 155 g[0] <= 155
 0_OnCutIn2_0: 10000 p_ij[0,2] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 0_OnCutIn3_0: 10000 p_ij[0,3] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 AbsZ10_1: - zi[0] + zi_hi[1] - k[0,1] <= 0
 AbsZ20_1: zi[0] - zi_hi[1] - k[0,1] <= 0
 1_Below_0: zi[0] - zi_hi[1] + k[0,1] + 310 m[0,1] <= 310
 1_NotBelow_0: - zi[0] + zi_hi[1] + k[0,1] - 310 m[0,1] <= 0
 HeightCompat_0_1: h[0,1] - k[0,1] <= 0
 NoHeightCompat_0_1: - 155 h[0,1] + k[0,1] <= 0
 Overlap_0_1: - x_p[0,1] - x_p[1,0] + o[0,1] <= 0
 NoOverlap_0_1: x_p[0,1] + x_p[1,0] - 2 o[0,1] <= 0
 3_Supports_0_1: - h[0,1] - o[0,1] - s[0,1] <= -1
 3_DoesntSupport_0_1: h[0,1] + o[0,1] + 2 s[0,1] <= 2
 beta1_0_1: - s[0,1] + beta1[0,1] <= 0
 beta2_0_1: - s[0,1] + beta2[0,1] <= 0
 sBeta1Flag_0_1: - s[0,1] + beta1[0,1] <= 0
 sBeta2Flag_0_1: - s[0,1] + beta2[0,1] <= 0
 Eta1Overlap_0_1: eta1[0,1] + beta1[0,1] <= 1
 Eta2Overlap_0_1: eta2[0,1] + beta2[0,1] <= 1
 Eta1Flag_0_1: - xi_l[0] + xi_l[1] - 300 eta1[0,1] <= 0
 Eta2Flag_0_1: xi_r[0] - xi_r[1] - 300 eta2[0,1] <= 0
 1_Supports_0_InSameBin1_0: p_ij[0,0] - p_ij[1,0] - s[0,1] <= 0
 1_Supports_0_InSameBin2_0: - p_ij[0,0] + p_ij[1,0] - s[0,1] <= 0
 Gamma0ForItems_0_inBins_0_NoCut_1: p_ij[0,0] + gamma[0] <= 1
 1_Supports_0_InSameBin1_1: p_ij[0,1] - p_ij[1,1] - s[0,1] <= 0
 1_Supports_0_InSameBin2_1: - p_ij[0,1] + p_ij[1,1] - s[0,1] <= 0
 Gamma0ForItems_0_inBins_1_NoCut_1: p_ij[0,1] + gamma[0] <= 1
 1_Supports_0_InSameBin1_2: p_ij[0,2] - p_ij[1,2] - s[0,1] <= 0
 1_Supports_0_InSameBin2_2: - p_ij[0,2] + p_ij[1,2] - s[0,1] <= 0
 Gamma0ForItems_0_inBins_2_NoCut_1: p_ij[0,2] + gamma[0] <= 98
 1_Supports_0_InSameBin1_3: p_ij[0,3] - p_ij[1,3] - s[0,1] <= 0
 1_Supports_0_InSameBin2_3: - p_ij[0,3] + p_ij[1,3] - s[0,1] <= 0
 Gamma0ForItems_0_inBins_3_NoCut_1: p_ij[0,3] + gamma[0] <= 98
 0_OnCutIn2_1: 10000 p_ij[0,2] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 0_OnCutIn3_1: 10000 p_ij[0,3] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 AbsZ10_2: - zi[0] + zi_hi[2] - k[0,2] <= 0
 AbsZ20_2: zi[0] - zi_hi[2] - k[0,2] <= 0
 2_Below_0: zi[0] - zi_hi[2] + k[0,2] + 310 m[0,2] <= 310
 2_NotBelow_0: - zi[0] + zi_hi[2] + k[0,2] - 310 m[0,2] <= 0
 HeightCompat_0_2: h[0,2] - k[0,2] <= 0
 NoHeightCompat_0_2: - 155 h[0,2] + k[0,2] <= 0
 Overlap_0_2: - x_p[0,2] - x_p[2,0] + o[0,2] <= 0
 NoOverlap_0_2: x_p[0,2] + x_p[2,0] - 2 o[0,2] <= 0
 3_Supports_0_2: - h[0,2] - o[0,2] - s[0,2] <= -1
 3_DoesntSupport_0_2: h[0,2] + o[0,2] + 2 s[0,2] <= 2
 beta1_0_2: - s[0,2] + beta1[0,2] <= 0
 beta2_0_2: - s[0,2] + beta2[0,2] <= 0
 sBeta1Flag_0_2: - s[0,2] + beta1[0,2] <= 0
 sBeta2Flag_0_2: - s[0,2] + beta2[0,2] <= 0
 Eta1Overlap_0_2: eta1[0,2] + beta1[0,2] <= 1
 Eta2Overlap_0_2: eta2[0,2] + beta2[0,2] <= 1
 Eta1Flag_0_2: - xi_l[0] + xi_l[2] - 300 eta1[0,2] <= 0
 Eta2Flag_0_2: xi_r[0] - xi_r[2] - 300 eta2[0,2] <= 0
 2_Supports_0_InSameBin1_0: p_ij[0,0] - p_ij[2,0] - s[0,2] <= 0
 2_Supports_0_InSameBin2_0: - p_ij[0,0] + p_ij[2,0] - s[0,2] <= 0
 Gamma0ForItems_0_inBins_0_NoCut_2: p_ij[0,0] + gamma[0] <= 1
 2_Supports_0_InSameBin1_1: p_ij[0,1] - p_ij[2,1] - s[0,2] <= 0
 2_Supports_0_InSameBin2_1: - p_ij[0,1] + p_ij[2,1] - s[0,2] <= 0
 Gamma0ForItems_0_inBins_1_NoCut_2: p_ij[0,1] + gamma[0] <= 1
 2_Supports_0_InSameBin1_2: p_ij[0,2] - p_ij[2,2] - s[0,2] <= 0
 2_Supports_0_InSameBin2_2: - p_ij[0,2] + p_ij[2,2] - s[0,2] <= 0
 Gamma0ForItems_0_inBins_2_NoCut_2: p_ij[0,2] + gamma[0] <= 98
 2_Supports_0_InSameBin1_3: p_ij[0,3] - p_ij[2,3] - s[0,2] <= 0
 2_Supports_0_InSameBin2_3: - p_ij[0,3] + p_ij[2,3] - s[0,2] <= 0
 Gamma0ForItems_0_inBins_3_NoCut_2: p_ij[0,3] + gamma[0] <= 98
 0_OnCutIn2_2: 10000 p_ij[0,2] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 0_OnCutIn3_2: 10000 p_ij[0,3] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 AbsZ10_3: - zi[0] + zi_hi[3] - k[0,3] <= 0
 AbsZ20_3: zi[0] - zi_hi[3] - k[0,3] <= 0
 3_Below_0: zi[0] - zi_hi[3] + k[0,3] + 310 m[0,3] <= 310
 3_NotBelow_0: - zi[0] + zi_hi[3] + k[0,3] - 310 m[0,3] <= 0
 HeightCompat_0_3: h[0,3] - k[0,3] <= 0
 NoHeightCompat_0_3: - 155 h[0,3] + k[0,3] <= 0
 Overlap_0_3: - x_p[0,3] - x_p[3,0] + o[0,3] <= 0
 NoOverlap_0_3: x_p[0,3] + x_p[3,0] - 2 o[0,3] <= 0
 3_Supports_0_3: - h[0,3] - o[0,3] - s[0,3] <= -1
 3_DoesntSupport_0_3: h[0,3] + o[0,3] + 2 s[0,3] <= 2
 beta1_0_3: - s[0,3] + beta1[0,3] <= 0
 beta2_0_3: - s[0,3] + beta2[0,3] <= 0
 sBeta1Flag_0_3: - s[0,3] + beta1[0,3] <= 0
 sBeta2Flag_0_3: - s[0,3] + beta2[0,3] <= 0
 Eta1Overlap_0_3: eta1[0,3] + beta1[0,3] <= 1
 Eta2Overlap_0_3: eta2[0,3] + beta2[0,3] <= 1
 Eta1Flag_0_3: - xi_l[0] + xi_l[3] - 300 eta1[0,3] <= 0
 Eta2Flag_0_3: xi_r[0] - xi_r[3] - 300 eta2[0,3] <= 0
 3_Supports_0_InSameBin1_0: p_ij[0,0] - p_ij[3,0] - s[0,3] <= 0
 3_Supports_0_InSameBin2_0: - p_ij[0,0] + p_ij[3,0] - s[0,3] <= 0
 Gamma0ForItems_0_inBins_0_NoCut_3: p_ij[0,0] + gamma[0] <= 1
 3_Supports_0_InSameBin1_1: p_ij[0,1] - p_ij[3,1] - s[0,3] <= 0
 3_Supports_0_InSameBin2_1: - p_ij[0,1] + p_ij[3,1] - s[0,3] <= 0
 Gamma0ForItems_0_inBins_1_NoCut_3: p_ij[0,1] + gamma[0] <= 1
 3_Supports_0_InSameBin1_2: p_ij[0,2] - p_ij[3,2] - s[0,3] <= 0
 3_Supports_0_InSameBin2_2: - p_ij[0,2] + p_ij[3,2] - s[0,3] <= 0
 Gamma0ForItems_0_inBins_2_NoCut_3: p_ij[0,2] + gamma[0] <= 98
 3_Supports_0_InSameBin1_3: p_ij[0,3] - p_ij[3,3] - s[0,3] <= 0
 3_Supports_0_InSameBin2_3: - p_ij[0,3] + p_ij[3,3] - s[0,3] <= 0
 Gamma0ForItems_0_inBins_3_NoCut_3: p_ij[0,3] + gamma[0] <= 98
 0_OnCutIn2_3: 10000 p_ij[0,2] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 0_OnCutIn3_3: 10000 p_ij[0,3] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 AbsZ10_4: - zi[0] + zi_hi[4] - k[0,4] <= 0
 AbsZ20_4: zi[0] - zi_hi[4] - k[0,4] <= 0
 4_Below_0: zi[0] - zi_hi[4] + k[0,4] + 310 m[0,4] <= 310
 4_NotBelow_0: - zi[0] + zi_hi[4] + k[0,4] - 310 m[0,4] <= 0
 HeightCompat_0_4: h[0,4] - k[0,4] <= 0
 NoHeightCompat_0_4: - 155 h[0,4] + k[0,4] <= 0
 Overlap_0_4: - x_p[0,4] - x_p[4,0] + o[0,4] <= 0
 NoOverlap_0_4: x_p[0,4] + x_p[4,0] - 2 o[0,4] <= 0
 3_Supports_0_4: - h[0,4] - o[0,4] - s[0,4] <= -1
 3_DoesntSupport_0_4: h[0,4] + o[0,4] + 2 s[0,4] <= 2
 beta1_0_4: - s[0,4] + beta1[0,4] <= 0
 beta2_0_4: - s[0,4] + beta2[0,4] <= 0
 sBeta1Flag_0_4: - s[0,4] + beta1[0,4] <= 0
 sBeta2Flag_0_4: - s[0,4] + beta2[0,4] <= 0
 Eta1Overlap_0_4: eta1[0,4] + beta1[0,4] <= 1
 Eta2Overlap_0_4: eta2[0,4] + beta2[0,4] <= 1
 Eta1Flag_0_4: - xi_l[0] + xi_l[4] - 300 eta1[0,4] <= 0
 Eta2Flag_0_4: xi_r[0] - xi_r[4] - 300 eta2[0,4] <= 0
 4_Supports_0_InSameBin1_0: p_ij[0,0] - p_ij[4,0] - s[0,4] <= 0
 4_Supports_0_InSameBin2_0: - p_ij[0,0] + p_ij[4,0] - s[0,4] <= 0
 Gamma0ForItems_0_inBins_0_NoCut_4: p_ij[0,0] + gamma[0] <= 1
 4_Supports_0_InSameBin1_1: p_ij[0,1] - p_ij[4,1] - s[0,4] <= 0
 4_Supports_0_InSameBin2_1: - p_ij[0,1] + p_ij[4,1] - s[0,4] <= 0
 Gamma0ForItems_0_inBins_1_NoCut_4: p_ij[0,1] + gamma[0] <= 1
 4_Supports_0_InSameBin1_2: p_ij[0,2] - p_ij[4,2] - s[0,4] <= 0
 4_Supports_0_InSameBin2_2: - p_ij[0,2] + p_ij[4,2] - s[0,4] <= 0
 Gamma0ForItems_0_inBins_2_NoCut_4: p_ij[0,2] + gamma[0] <= 98
 4_Supports_0_InSameBin1_3: p_ij[0,3] - p_ij[4,3] - s[0,4] <= 0
 4_Supports_0_InSameBin2_3: - p_ij[0,3] + p_ij[4,3] - s[0,4] <= 0
 Gamma0ForItems_0_inBins_3_NoCut_4: p_ij[0,3] + gamma[0] <= 98
 0_OnCutIn2_4: 10000 p_ij[0,2] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 0_OnCutIn3_4: 10000 p_ij[0,3] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 AbsZ10_5: - zi[0] + zi_hi[5] - k[0,5] <= 0
 AbsZ20_5: zi[0] - zi_hi[5] - k[0,5] <= 0
 5_Below_0: zi[0] - zi_hi[5] + k[0,5] + 310 m[0,5] <= 310
 5_NotBelow_0: - zi[0] + zi_hi[5] + k[0,5] - 310 m[0,5] <= 0
 HeightCompat_0_5: h[0,5] - k[0,5] <= 0
 NoHeightCompat_0_5: - 155 h[0,5] + k[0,5] <= 0
 Overlap_0_5: - x_p[0,5] - x_p[5,0] + o[0,5] <= 0
 NoOverlap_0_5: x_p[0,5] + x_p[5,0] - 2 o[0,5] <= 0
 3_Supports_0_5: - h[0,5] - o[0,5] - s[0,5] <= -1
 3_DoesntSupport_0_5: h[0,5] + o[0,5] + 2 s[0,5] <= 2
 beta1_0_5: - s[0,5] + beta1[0,5] <= 0
 beta2_0_5: - s[0,5] + beta2[0,5] <= 0
 sBeta1Flag_0_5: - s[0,5] + beta1[0,5] <= 0
 sBeta2Flag_0_5: - s[0,5] + beta2[0,5] <= 0
 Eta1Overlap_0_5: eta1[0,5] + beta1[0,5] <= 1
 Eta2Overlap_0_5: eta2[0,5] + beta2[0,5] <= 1
 Eta1Flag_0_5: - xi_l[0] + xi_l[5] - 300 eta1[0,5] <= 0
 Eta2Flag_0_5: xi_r[0] - xi_r[5] - 300 eta2[0,5] <= 0
 5_Supports_0_InSameBin1_0: p_ij[0,0] - p_ij[5,0] - s[0,5] <= 0
 5_Supports_0_InSameBin2_0: - p_ij[0,0] + p_ij[5,0] - s[0,5] <= 0
 Gamma0ForItems_0_inBins_0_NoCut_5: p_ij[0,0] + gamma[0] <= 1
 5_Supports_0_InSameBin1_1: p_ij[0,1] - p_ij[5,1] - s[0,5] <= 0
 5_Supports_0_InSameBin2_1: - p_ij[0,1] + p_ij[5,1] - s[0,5] <= 0
 Gamma0ForItems_0_inBins_1_NoCut_5: p_ij[0,1] + gamma[0] <= 1
 5_Supports_0_InSameBin1_2: p_ij[0,2] - p_ij[5,2] - s[0,5] <= 0
 5_Supports_0_InSameBin2_2: - p_ij[0,2] + p_ij[5,2] - s[0,5] <= 0
 Gamma0ForItems_0_inBins_2_NoCut_5: p_ij[0,2] + gamma[0] <= 98
 5_Supports_0_InSameBin1_3: p_ij[0,3] - p_ij[5,3] - s[0,5] <= 0
 5_Supports_0_InSameBin2_3: - p_ij[0,3] + p_ij[5,3] - s[0,5] <= 0
 Gamma0ForItems_0_inBins_3_NoCut_5: p_ij[0,3] + gamma[0] <= 98
 0_OnCutIn2_5: 10000 p_ij[0,2] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 0_OnCutIn3_5: 10000 p_ij[0,3] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 AbsZ10_6: - zi[0] + zi_hi[6] - k[0,6] <= 0
 AbsZ20_6: zi[0] - zi_hi[6] - k[0,6] <= 0
 6_Below_0: zi[0] - zi_hi[6] + k[0,6] + 310 m[0,6] <= 310
 6_NotBelow_0: - zi[0] + zi_hi[6] + k[0,6] - 310 m[0,6] <= 0
 HeightCompat_0_6: h[0,6] - k[0,6] <= 0
 NoHeightCompat_0_6: - 155 h[0,6] + k[0,6] <= 0
 Overlap_0_6: - x_p[0,6] - x_p[6,0] + o[0,6] <= 0
 NoOverlap_0_6: x_p[0,6] + x_p[6,0] - 2 o[0,6] <= 0
 3_Supports_0_6: - h[0,6] - o[0,6] - s[0,6] <= -1
 3_DoesntSupport_0_6: h[0,6] + o[0,6] + 2 s[0,6] <= 2
 beta1_0_6: - s[0,6] + beta1[0,6] <= 0
 beta2_0_6: - s[0,6] + beta2[0,6] <= 0
 sBeta1Flag_0_6: - s[0,6] + beta1[0,6] <= 0
 sBeta2Flag_0_6: - s[0,6] + beta2[0,6] <= 0
 Eta1Overlap_0_6: eta1[0,6] + beta1[0,6] <= 1
 Eta2Overlap_0_6: eta2[0,6] + beta2[0,6] <= 1
 Eta1Flag_0_6: - xi_l[0] + xi_l[6] - 300 eta1[0,6] <= 0
 Eta2Flag_0_6: xi_r[0] - xi_r[6] - 300 eta2[0,6] <= 0
 6_Supports_0_InSameBin1_0: p_ij[0,0] - p_ij[6,0] - s[0,6] <= 0
 6_Supports_0_InSameBin2_0: - p_ij[0,0] + p_ij[6,0] - s[0,6] <= 0
 Gamma0ForItems_0_inBins_0_NoCut_6: p_ij[0,0] + gamma[0] <= 1
 6_Supports_0_InSameBin1_1: p_ij[0,1] - p_ij[6,1] - s[0,6] <= 0
 6_Supports_0_InSameBin2_1: - p_ij[0,1] + p_ij[6,1] - s[0,6] <= 0
 Gamma0ForItems_0_inBins_1_NoCut_6: p_ij[0,1] + gamma[0] <= 1
 6_Supports_0_InSameBin1_2: p_ij[0,2] - p_ij[6,2] - s[0,6] <= 0
 6_Supports_0_InSameBin2_2: - p_ij[0,2] + p_ij[6,2] - s[0,6] <= 0
 Gamma0ForItems_0_inBins_2_NoCut_6: p_ij[0,2] + gamma[0] <= 98
 6_Supports_0_InSameBin1_3: p_ij[0,3] - p_ij[6,3] - s[0,6] <= 0
 6_Supports_0_InSameBin2_3: - p_ij[0,3] + p_ij[6,3] - s[0,6] <= 0
 Gamma0ForItems_0_inBins_3_NoCut_6: p_ij[0,3] + gamma[0] <= 98
 0_OnCutIn2_6: 10000 p_ij[0,2] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 0_OnCutIn3_6: 10000 p_ij[0,3] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 AbsZ10_7: - zi[0] + zi_hi[7] - k[0,7] <= 0
 AbsZ20_7: zi[0] - zi_hi[7] - k[0,7] <= 0
 7_Below_0: zi[0] - zi_hi[7] + k[0,7] + 310 m[0,7] <= 310
 7_NotBelow_0: - zi[0] + zi_hi[7] + k[0,7] - 310 m[0,7] <= 0
 HeightCompat_0_7: h[0,7] - k[0,7] <= 0
 NoHeightCompat_0_7: - 155 h[0,7] + k[0,7] <= 0
 Overlap_0_7: - x_p[0,7] - x_p[7,0] + o[0,7] <= 0
 NoOverlap_0_7: x_p[0,7] + x_p[7,0] - 2 o[0,7] <= 0
 3_Supports_0_7: - h[0,7] - o[0,7] - s[0,7] <= -1
 3_DoesntSupport_0_7: h[0,7] + o[0,7] + 2 s[0,7] <= 2
 beta1_0_7: - s[0,7] + beta1[0,7] <= 0
 beta2_0_7: - s[0,7] + beta2[0,7] <= 0
 sBeta1Flag_0_7: - s[0,7] + beta1[0,7] <= 0
 sBeta2Flag_0_7: - s[0,7] + beta2[0,7] <= 0
 Eta1Overlap_0_7: eta1[0,7] + beta1[0,7] <= 1
 Eta2Overlap_0_7: eta2[0,7] + beta2[0,7] <= 1
 Eta1Flag_0_7: - xi_l[0] + xi_l[7] - 300 eta1[0,7] <= 0
 Eta2Flag_0_7: xi_r[0] - xi_r[7] - 300 eta2[0,7] <= 0
 7_Supports_0_InSameBin1_0: p_ij[0,0] - p_ij[7,0] - s[0,7] <= 0
 7_Supports_0_InSameBin2_0: - p_ij[0,0] + p_ij[7,0] - s[0,7] <= 0
 Gamma0ForItems_0_inBins_0_NoCut_7: p_ij[0,0] + gamma[0] <= 1
 7_Supports_0_InSameBin1_1: p_ij[0,1] - p_ij[7,1] - s[0,7] <= 0
 7_Supports_0_InSameBin2_1: - p_ij[0,1] + p_ij[7,1] - s[0,7] <= 0
 Gamma0ForItems_0_inBins_1_NoCut_7: p_ij[0,1] + gamma[0] <= 1
 7_Supports_0_InSameBin1_2: p_ij[0,2] - p_ij[7,2] - s[0,7] <= 0
 7_Supports_0_InSameBin2_2: - p_ij[0,2] + p_ij[7,2] - s[0,7] <= 0
 Gamma0ForItems_0_inBins_2_NoCut_7: p_ij[0,2] + gamma[0] <= 98
 7_Supports_0_InSameBin1_3: p_ij[0,3] - p_ij[7,3] - s[0,7] <= 0
 7_Supports_0_InSameBin2_3: - p_ij[0,3] + p_ij[7,3] - s[0,7] <= 0
 Gamma0ForItems_0_inBins_3_NoCut_7: p_ij[0,3] + gamma[0] <= 98
 0_OnCutIn2_7: 10000 p_ij[0,2] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 0_OnCutIn3_7: 10000 p_ij[0,3] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 AbsZ10_8: - zi[0] + zi_hi[8] - k[0,8] <= 0
 AbsZ20_8: zi[0] - zi_hi[8] - k[0,8] <= 0
 8_Below_0: zi[0] - zi_hi[8] + k[0,8] + 310 m[0,8] <= 310
 8_NotBelow_0: - zi[0] + zi_hi[8] + k[0,8] - 310 m[0,8] <= 0
 HeightCompat_0_8: h[0,8] - k[0,8] <= 0
 NoHeightCompat_0_8: - 155 h[0,8] + k[0,8] <= 0
 Overlap_0_8: - x_p[0,8] - x_p[8,0] + o[0,8] <= 0
 NoOverlap_0_8: x_p[0,8] + x_p[8,0] - 2 o[0,8] <= 0
 3_Supports_0_8: - h[0,8] - o[0,8] - s[0,8] <= -1
 3_DoesntSupport_0_8: h[0,8] + o[0,8] + 2 s[0,8] <= 2
 beta1_0_8: - s[0,8] + beta1[0,8] <= 0
 beta2_0_8: - s[0,8] + beta2[0,8] <= 0
 sBeta1Flag_0_8: - s[0,8] + beta1[0,8] <= 0
 sBeta2Flag_0_8: - s[0,8] + beta2[0,8] <= 0
 Eta1Overlap_0_8: eta1[0,8] + beta1[0,8] <= 1
 Eta2Overlap_0_8: eta2[0,8] + beta2[0,8] <= 1
 Eta1Flag_0_8: - xi_l[0] + xi_l[8] - 300 eta1[0,8] <= 0
 Eta2Flag_0_8: xi_r[0] - xi_r[8] - 300 eta2[0,8] <= 0
 8_Supports_0_InSameBin1_0: p_ij[0,0] - p_ij[8,0] - s[0,8] <= 0
 8_Supports_0_InSameBin2_0: - p_ij[0,0] + p_ij[8,0] - s[0,8] <= 0
 Gamma0ForItems_0_inBins_0_NoCut_8: p_ij[0,0] + gamma[0] <= 1
 8_Supports_0_InSameBin1_1: p_ij[0,1] - p_ij[8,1] - s[0,8] <= 0
 8_Supports_0_InSameBin2_1: - p_ij[0,1] + p_ij[8,1] - s[0,8] <= 0
 Gamma0ForItems_0_inBins_1_NoCut_8: p_ij[0,1] + gamma[0] <= 1
 8_Supports_0_InSameBin1_2: p_ij[0,2] - p_ij[8,2] - s[0,8] <= 0
 8_Supports_0_InSameBin2_2: - p_ij[0,2] + p_ij[8,2] - s[0,8] <= 0
 Gamma0ForItems_0_inBins_2_NoCut_8: p_ij[0,2] + gamma[0] <= 98
 8_Supports_0_InSameBin1_3: p_ij[0,3] - p_ij[8,3] - s[0,8] <= 0
 8_Supports_0_InSameBin2_3: - p_ij[0,3] + p_ij[8,3] - s[0,8] <= 0
 Gamma0ForItems_0_inBins_3_NoCut_8: p_ij[0,3] + gamma[0] <= 98
 0_OnCutIn2_8: 10000 p_ij[0,2] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 0_OnCutIn3_8: 10000 p_ij[0,3] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 AbsZ10_9: - zi[0] + zi_hi[9] - k[0,9] <= 0
 AbsZ20_9: zi[0] - zi_hi[9] - k[0,9] <= 0
 9_Below_0: zi[0] - zi_hi[9] + k[0,9] + 310 m[0,9] <= 310
 9_NotBelow_0: - zi[0] + zi_hi[9] + k[0,9] - 310 m[0,9] <= 0
 HeightCompat_0_9: h[0,9] - k[0,9] <= 0
 NoHeightCompat_0_9: - 155 h[0,9] + k[0,9] <= 0
 Overlap_0_9: - x_p[0,9] - x_p[9,0] + o[0,9] <= 0
 NoOverlap_0_9: x_p[0,9] + x_p[9,0] - 2 o[0,9] <= 0
 3_Supports_0_9: - h[0,9] - o[0,9] - s[0,9] <= -1
 3_DoesntSupport_0_9: h[0,9] + o[0,9] + 2 s[0,9] <= 2
 beta1_0_9: - s[0,9] + beta1[0,9] <= 0
 beta2_0_9: - s[0,9] + beta2[0,9] <= 0
 sBeta1Flag_0_9: - s[0,9] + beta1[0,9] <= 0
 sBeta2Flag_0_9: - s[0,9] + beta2[0,9] <= 0
 Eta1Overlap_0_9: eta1[0,9] + beta1[0,9] <= 1
 Eta2Overlap_0_9: eta2[0,9] + beta2[0,9] <= 1
 Eta1Flag_0_9: - xi_l[0] + xi_l[9] - 300 eta1[0,9] <= 0
 Eta2Flag_0_9: xi_r[0] - xi_r[9] - 300 eta2[0,9] <= 0
 9_Supports_0_InSameBin1_0: p_ij[0,0] - p_ij[9,0] - s[0,9] <= 0
 9_Supports_0_InSameBin2_0: - p_ij[0,0] + p_ij[9,0] - s[0,9] <= 0
 Gamma0ForItems_0_inBins_0_NoCut_9: p_ij[0,0] + gamma[0] <= 1
 9_Supports_0_InSameBin1_1: p_ij[0,1] - p_ij[9,1] - s[0,9] <= 0
 9_Supports_0_InSameBin2_1: - p_ij[0,1] + p_ij[9,1] - s[0,9] <= 0
 Gamma0ForItems_0_inBins_1_NoCut_9: p_ij[0,1] + gamma[0] <= 1
 9_Supports_0_InSameBin1_2: p_ij[0,2] - p_ij[9,2] - s[0,9] <= 0
 9_Supports_0_InSameBin2_2: - p_ij[0,2] + p_ij[9,2] - s[0,9] <= 0
 Gamma0ForItems_0_inBins_2_NoCut_9: p_ij[0,2] + gamma[0] <= 98
 9_Supports_0_InSameBin1_3: p_ij[0,3] - p_ij[9,3] - s[0,9] <= 0
 9_Supports_0_InSameBin2_3: - p_ij[0,3] + p_ij[9,3] - s[0,9] <= 0
 Gamma0ForItems_0_inBins_3_NoCut_9: p_ij[0,3] + gamma[0] <= 98
 0_OnCutIn2_9: 10000 p_ij[0,2] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 0_OnCutIn3_9: 10000 p_ij[0,3] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 AbsZ10_10: - zi[0] + zi_hi[10] - k[0,10] <= 0
 AbsZ20_10: zi[0] - zi_hi[10] - k[0,10] <= 0
 10_Below_0: zi[0] - zi_hi[10] + k[0,10] + 310 m[0,10] <= 310
 10_NotBelow_0: - zi[0] + zi_hi[10] + k[0,10] - 310 m[0,10] <= 0
 HeightCompat_0_10: h[0,10] - k[0,10] <= 0
 NoHeightCompat_0_10: - 155 h[0,10] + k[0,10] <= 0
 Overlap_0_10: - x_p[0,10] - x_p[10,0] + o[0,10] <= 0
 NoOverlap_0_10: x_p[0,10] + x_p[10,0] - 2 o[0,10] <= 0
 3_Supports_0_10: - h[0,10] - o[0,10] - s[0,10] <= -1
 3_DoesntSupport_0_10: h[0,10] + o[0,10] + 2 s[0,10] <= 2
 beta1_0_10: - s[0,10] + beta1[0,10] <= 0
 beta2_0_10: - s[0,10] + beta2[0,10] <= 0
 sBeta1Flag_0_10: - s[0,10] + beta1[0,10] <= 0
 sBeta2Flag_0_10: - s[0,10] + beta2[0,10] <= 0
 Eta1Overlap_0_10: eta1[0,10] + beta1[0,10] <= 1
 Eta2Overlap_0_10: eta2[0,10] + beta2[0,10] <= 1
 Eta1Flag_0_10: - xi_l[0] + xi_l[10] - 300 eta1[0,10] <= 0
 Eta2Flag_0_10: xi_r[0] - xi_r[10] - 300 eta2[0,10] <= 0
 10_Supports_0_InSameBin1_0: p_ij[0,0] - p_ij[10,0] - s[0,10] <= 0
 10_Supports_0_InSameBin2_0: - p_ij[0,0] + p_ij[10,0] - s[0,10] <= 0
 Gamma0ForItems_0_inBins_0_NoCut_10: p_ij[0,0] + gamma[0] <= 1
 10_Supports_0_InSameBin1_1: p_ij[0,1] - p_ij[10,1] - s[0,10] <= 0
 10_Supports_0_InSameBin2_1: - p_ij[0,1] + p_ij[10,1] - s[0,10] <= 0
 Gamma0ForItems_0_inBins_1_NoCut_10: p_ij[0,1] + gamma[0] <= 1
 10_Supports_0_InSameBin1_2: p_ij[0,2] - p_ij[10,2] - s[0,10] <= 0
 10_Supports_0_InSameBin2_2: - p_ij[0,2] + p_ij[10,2] - s[0,10] <= 0
 Gamma0ForItems_0_inBins_2_NoCut_10: p_ij[0,2] + gamma[0] <= 98
 10_Supports_0_InSameBin1_3: p_ij[0,3] - p_ij[10,3] - s[0,10] <= 0
 10_Supports_0_InSameBin2_3: - p_ij[0,3] + p_ij[10,3] - s[0,10] <= 0
 Gamma0ForItems_0_inBins_3_NoCut_10: p_ij[0,3] + gamma[0] <= 98
 0_OnCutIn2_10: 10000 p_ij[0,2] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 0_OnCutIn3_10: 10000 p_ij[0,3] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 AbsZ10_11: - zi[0] + zi_hi[11] - k[0,11] <= 0
 AbsZ20_11: zi[0] - zi_hi[11] - k[0,11] <= 0
 11_Below_0: zi[0] - zi_hi[11] + k[0,11] + 310 m[0,11] <= 310
 11_NotBelow_0: - zi[0] + zi_hi[11] + k[0,11] - 310 m[0,11] <= 0
 HeightCompat_0_11: h[0,11] - k[0,11] <= 0
 NoHeightCompat_0_11: - 155 h[0,11] + k[0,11] <= 0
 Overlap_0_11: - x_p[0,11] - x_p[11,0] + o[0,11] <= 0
 NoOverlap_0_11: x_p[0,11] + x_p[11,0] - 2 o[0,11] <= 0
 3_Supports_0_11: - h[0,11] - o[0,11] - s[0,11] <= -1
 3_DoesntSupport_0_11: h[0,11] + o[0,11] + 2 s[0,11] <= 2
 beta1_0_11: - s[0,11] + beta1[0,11] <= 0
 beta2_0_11: - s[0,11] + beta2[0,11] <= 0
 sBeta1Flag_0_11: - s[0,11] + beta1[0,11] <= 0
 sBeta2Flag_0_11: - s[0,11] + beta2[0,11] <= 0
 Eta1Overlap_0_11: eta1[0,11] + beta1[0,11] <= 1
 Eta2Overlap_0_11: eta2[0,11] + beta2[0,11] <= 1
 Eta1Flag_0_11: - xi_l[0] + xi_l[11] - 300 eta1[0,11] <= 0
 Eta2Flag_0_11: xi_r[0] - xi_r[11] - 300 eta2[0,11] <= 0
 11_Supports_0_InSameBin1_0: p_ij[0,0] - p_ij[11,0] - s[0,11] <= 0
 11_Supports_0_InSameBin2_0: - p_ij[0,0] + p_ij[11,0] - s[0,11] <= 0
 Gamma0ForItems_0_inBins_0_NoCut_11: p_ij[0,0] + gamma[0] <= 1
 11_Supports_0_InSameBin1_1: p_ij[0,1] - p_ij[11,1] - s[0,11] <= 0
 11_Supports_0_InSameBin2_1: - p_ij[0,1] + p_ij[11,1] - s[0,11] <= 0
 Gamma0ForItems_0_inBins_1_NoCut_11: p_ij[0,1] + gamma[0] <= 1
 11_Supports_0_InSameBin1_2: p_ij[0,2] - p_ij[11,2] - s[0,11] <= 0
 11_Supports_0_InSameBin2_2: - p_ij[0,2] + p_ij[11,2] - s[0,11] <= 0
 Gamma0ForItems_0_inBins_2_NoCut_11: p_ij[0,2] + gamma[0] <= 98
 11_Supports_0_InSameBin1_3: p_ij[0,3] - p_ij[11,3] - s[0,11] <= 0
 11_Supports_0_InSameBin2_3: - p_ij[0,3] + p_ij[11,3] - s[0,11] <= 0
 Gamma0ForItems_0_inBins_3_NoCut_11: p_ij[0,3] + gamma[0] <= 98
 0_OnCutIn2_11: 10000 p_ij[0,2] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 0_OnCutIn3_11: 10000 p_ij[0,3] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 AbsZ10_12: - zi[0] + zi_hi[12] - k[0,12] <= 0
 AbsZ20_12: zi[0] - zi_hi[12] - k[0,12] <= 0
 12_Below_0: zi[0] - zi_hi[12] + k[0,12] + 310 m[0,12] <= 310
 12_NotBelow_0: - zi[0] + zi_hi[12] + k[0,12] - 310 m[0,12] <= 0
 HeightCompat_0_12: h[0,12] - k[0,12] <= 0
 NoHeightCompat_0_12: - 155 h[0,12] + k[0,12] <= 0
 Overlap_0_12: - x_p[0,12] - x_p[12,0] + o[0,12] <= 0
 NoOverlap_0_12: x_p[0,12] + x_p[12,0] - 2 o[0,12] <= 0
 3_Supports_0_12: - h[0,12] - o[0,12] - s[0,12] <= -1
 3_DoesntSupport_0_12: h[0,12] + o[0,12] + 2 s[0,12] <= 2
 beta1_0_12: - s[0,12] + beta1[0,12] <= 0
 beta2_0_12: - s[0,12] + beta2[0,12] <= 0
 sBeta1Flag_0_12: - s[0,12] + beta1[0,12] <= 0
 sBeta2Flag_0_12: - s[0,12] + beta2[0,12] <= 0
 Eta1Overlap_0_12: eta1[0,12] + beta1[0,12] <= 1
 Eta2Overlap_0_12: eta2[0,12] + beta2[0,12] <= 1
 Eta1Flag_0_12: - xi_l[0] + xi_l[12] - 300 eta1[0,12] <= 0
 Eta2Flag_0_12: xi_r[0] - xi_r[12] - 300 eta2[0,12] <= 0
 12_Supports_0_InSameBin1_0: p_ij[0,0] - p_ij[12,0] - s[0,12] <= 0
 12_Supports_0_InSameBin2_0: - p_ij[0,0] + p_ij[12,0] - s[0,12] <= 0
 Gamma0ForItems_0_inBins_0_NoCut_12: p_ij[0,0] + gamma[0] <= 1
 12_Supports_0_InSameBin1_1: p_ij[0,1] - p_ij[12,1] - s[0,12] <= 0
 12_Supports_0_InSameBin2_1: - p_ij[0,1] + p_ij[12,1] - s[0,12] <= 0
 Gamma0ForItems_0_inBins_1_NoCut_12: p_ij[0,1] + gamma[0] <= 1
 12_Supports_0_InSameBin1_2: p_ij[0,2] - p_ij[12,2] - s[0,12] <= 0
 12_Supports_0_InSameBin2_2: - p_ij[0,2] + p_ij[12,2] - s[0,12] <= 0
 Gamma0ForItems_0_inBins_2_NoCut_12: p_ij[0,2] + gamma[0] <= 98
 12_Supports_0_InSameBin1_3: p_ij[0,3] - p_ij[12,3] - s[0,12] <= 0
 12_Supports_0_InSameBin2_3: - p_ij[0,3] + p_ij[12,3] - s[0,12] <= 0
 Gamma0ForItems_0_inBins_3_NoCut_12: p_ij[0,3] + gamma[0] <= 98
 0_OnCutIn2_12: 10000 p_ij[0,2] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 0_OnCutIn3_12: 10000 p_ij[0,3] + 1.261904761904762 xi_l[0] + zi[0]
   + 10000 gamma[0] <= 20053
 GroundConstraint_1: zi[1] + 155 g[1] <= 155
 AbsZ11_0: - zi[1] + zi_hi[0] - k[1,0] <= 0
 AbsZ21_0: zi[1] - zi_hi[0] - k[1,0] <= 0
 0_Below_1: zi[1] - zi_hi[0] + k[1,0] + 310 m[1,0] <= 310
 0_NotBelow_1: - zi[1] + zi_hi[0] + k[1,0] - 310 m[1,0] <= 0
 HeightCompat_1_0: h[1,0] - k[1,0] <= 0
 NoHeightCompat_1_0: - 155 h[1,0] + k[1,0] <= 0
 Overlap_1_0: - x_p[0,1] - x_p[1,0] + o[1,0] <= 0
 NoOverlap_1_0: x_p[0,1] + x_p[1,0] - 2 o[1,0] <= 0
 3_Supports_1_0: - h[1,0] - o[1,0] - s[1,0] <= -1
 3_DoesntSupport_1_0: h[1,0] + o[1,0] + 2 s[1,0] <= 2
 beta1_1_0: - s[1,0] + beta1[1,0] <= 0
 beta2_1_0: - s[1,0] + beta2[1,0] <= 0
 sBeta1Flag_1_0: - s[1,0] + beta1[1,0] <= 0
 sBeta2Flag_1_0: - s[1,0] + beta2[1,0] <= 0
 Eta1Overlap_1_0: eta1[1,0] + beta1[1,0] <= 1
 Eta2Overlap_1_0: eta2[1,0] + beta2[1,0] <= 1
 Eta1Flag_1_0: xi_l[0] - xi_l[1] - 300 eta1[1,0] <= 0
 Eta2Flag_1_0: - xi_r[0] + xi_r[1] - 300 eta2[1,0] <= 0
 0_Supports_1_InSameBin1_0: - p_ij[0,0] + p_ij[1,0] - s[1,0] <= 0
 0_Supports_1_InSameBin2_0: p_ij[0,0] - p_ij[1,0] - s[1,0] <= 0
 Gamma0ForItems_1_inBins_0_NoCut_0: p_ij[1,0] + gamma[1] <= 1
 0_Supports_1_InSameBin1_1: - p_ij[0,1] + p_ij[1,1] - s[1,0] <= 0
 0_Supports_1_InSameBin2_1: p_ij[0,1] - p_ij[1,1] - s[1,0] <= 0
 Gamma0ForItems_1_inBins_1_NoCut_0: p_ij[1,1] + gamma[1] <= 1
 0_Supports_1_InSameBin1_2: - p_ij[0,2] + p_ij[1,2] - s[1,0] <= 0
 0_Supports_1_InSameBin2_2: p_ij[0,2] - p_ij[1,2] - s[1,0] <= 0
 Gamma0ForItems_1_inBins_2_NoCut_0: p_ij[1,2] + gamma[1] <= 98
 0_Supports_1_InSameBin1_3: - p_ij[0,3] + p_ij[1,3] - s[1,0] <= 0
 0_Supports_1_InSameBin2_3: p_ij[0,3] - p_ij[1,3] - s[1,0] <= 0
 Gamma0ForItems_1_inBins_3_NoCut_0: p_ij[1,3] + gamma[1] <= 98
 1_OnCutIn2_0: 10000 p_ij[1,2] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 1_OnCutIn3_0: 10000 p_ij[1,3] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 1_OnCutIn2_1: 10000 p_ij[1,2] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 1_OnCutIn3_1: 10000 p_ij[1,3] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 AbsZ11_2: - zi[1] + zi_hi[2] - k[1,2] <= 0
 AbsZ21_2: zi[1] - zi_hi[2] - k[1,2] <= 0
 2_Below_1: zi[1] - zi_hi[2] + k[1,2] + 310 m[1,2] <= 310
 2_NotBelow_1: - zi[1] + zi_hi[2] + k[1,2] - 310 m[1,2] <= 0
 HeightCompat_1_2: h[1,2] - k[1,2] <= 0
 NoHeightCompat_1_2: - 155 h[1,2] + k[1,2] <= 0
 Overlap_1_2: - x_p[1,2] - x_p[2,1] + o[1,2] <= 0
 NoOverlap_1_2: x_p[1,2] + x_p[2,1] - 2 o[1,2] <= 0
 3_Supports_1_2: - h[1,2] - o[1,2] - s[1,2] <= -1
 3_DoesntSupport_1_2: h[1,2] + o[1,2] + 2 s[1,2] <= 2
 beta1_1_2: - s[1,2] + beta1[1,2] <= 0
 beta2_1_2: - s[1,2] + beta2[1,2] <= 0
 sBeta1Flag_1_2: - s[1,2] + beta1[1,2] <= 0
 sBeta2Flag_1_2: - s[1,2] + beta2[1,2] <= 0
 Eta1Overlap_1_2: eta1[1,2] + beta1[1,2] <= 1
 Eta2Overlap_1_2: eta2[1,2] + beta2[1,2] <= 1
 Eta1Flag_1_2: - xi_l[1] + xi_l[2] - 300 eta1[1,2] <= 0
 Eta2Flag_1_2: xi_r[1] - xi_r[2] - 300 eta2[1,2] <= 0
 2_Supports_1_InSameBin1_0: p_ij[1,0] - p_ij[2,0] - s[1,2] <= 0
 2_Supports_1_InSameBin2_0: - p_ij[1,0] + p_ij[2,0] - s[1,2] <= 0
 Gamma0ForItems_1_inBins_0_NoCut_2: p_ij[1,0] + gamma[1] <= 1
 2_Supports_1_InSameBin1_1: p_ij[1,1] - p_ij[2,1] - s[1,2] <= 0
 2_Supports_1_InSameBin2_1: - p_ij[1,1] + p_ij[2,1] - s[1,2] <= 0
 Gamma0ForItems_1_inBins_1_NoCut_2: p_ij[1,1] + gamma[1] <= 1
 2_Supports_1_InSameBin1_2: p_ij[1,2] - p_ij[2,2] - s[1,2] <= 0
 2_Supports_1_InSameBin2_2: - p_ij[1,2] + p_ij[2,2] - s[1,2] <= 0
 Gamma0ForItems_1_inBins_2_NoCut_2: p_ij[1,2] + gamma[1] <= 98
 2_Supports_1_InSameBin1_3: p_ij[1,3] - p_ij[2,3] - s[1,2] <= 0
 2_Supports_1_InSameBin2_3: - p_ij[1,3] + p_ij[2,3] - s[1,2] <= 0
 Gamma0ForItems_1_inBins_3_NoCut_2: p_ij[1,3] + gamma[1] <= 98
 1_OnCutIn2_2: 10000 p_ij[1,2] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 1_OnCutIn3_2: 10000 p_ij[1,3] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 AbsZ11_3: - zi[1] + zi_hi[3] - k[1,3] <= 0
 AbsZ21_3: zi[1] - zi_hi[3] - k[1,3] <= 0
 3_Below_1: zi[1] - zi_hi[3] + k[1,3] + 310 m[1,3] <= 310
 3_NotBelow_1: - zi[1] + zi_hi[3] + k[1,3] - 310 m[1,3] <= 0
 HeightCompat_1_3: h[1,3] - k[1,3] <= 0
 NoHeightCompat_1_3: - 155 h[1,3] + k[1,3] <= 0
 Overlap_1_3: - x_p[1,3] - x_p[3,1] + o[1,3] <= 0
 NoOverlap_1_3: x_p[1,3] + x_p[3,1] - 2 o[1,3] <= 0
 3_Supports_1_3: - h[1,3] - o[1,3] - s[1,3] <= -1
 3_DoesntSupport_1_3: h[1,3] + o[1,3] + 2 s[1,3] <= 2
 beta1_1_3: - s[1,3] + beta1[1,3] <= 0
 beta2_1_3: - s[1,3] + beta2[1,3] <= 0
 sBeta1Flag_1_3: - s[1,3] + beta1[1,3] <= 0
 sBeta2Flag_1_3: - s[1,3] + beta2[1,3] <= 0
 Eta1Overlap_1_3: eta1[1,3] + beta1[1,3] <= 1
 Eta2Overlap_1_3: eta2[1,3] + beta2[1,3] <= 1
 Eta1Flag_1_3: - xi_l[1] + xi_l[3] - 300 eta1[1,3] <= 0
 Eta2Flag_1_3: xi_r[1] - xi_r[3] - 300 eta2[1,3] <= 0
 3_Supports_1_InSameBin1_0: p_ij[1,0] - p_ij[3,0] - s[1,3] <= 0
 3_Supports_1_InSameBin2_0: - p_ij[1,0] + p_ij[3,0] - s[1,3] <= 0
 Gamma0ForItems_1_inBins_0_NoCut_3: p_ij[1,0] + gamma[1] <= 1
 3_Supports_1_InSameBin1_1: p_ij[1,1] - p_ij[3,1] - s[1,3] <= 0
 3_Supports_1_InSameBin2_1: - p_ij[1,1] + p_ij[3,1] - s[1,3] <= 0
 Gamma0ForItems_1_inBins_1_NoCut_3: p_ij[1,1] + gamma[1] <= 1
 3_Supports_1_InSameBin1_2: p_ij[1,2] - p_ij[3,2] - s[1,3] <= 0
 3_Supports_1_InSameBin2_2: - p_ij[1,2] + p_ij[3,2] - s[1,3] <= 0
 Gamma0ForItems_1_inBins_2_NoCut_3: p_ij[1,2] + gamma[1] <= 98
 3_Supports_1_InSameBin1_3: p_ij[1,3] - p_ij[3,3] - s[1,3] <= 0
 3_Supports_1_InSameBin2_3: - p_ij[1,3] + p_ij[3,3] - s[1,3] <= 0
 Gamma0ForItems_1_inBins_3_NoCut_3: p_ij[1,3] + gamma[1] <= 98
 1_OnCutIn2_3: 10000 p_ij[1,2] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 1_OnCutIn3_3: 10000 p_ij[1,3] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 AbsZ11_4: - zi[1] + zi_hi[4] - k[1,4] <= 0
 AbsZ21_4: zi[1] - zi_hi[4] - k[1,4] <= 0
 4_Below_1: zi[1] - zi_hi[4] + k[1,4] + 310 m[1,4] <= 310
 4_NotBelow_1: - zi[1] + zi_hi[4] + k[1,4] - 310 m[1,4] <= 0
 HeightCompat_1_4: h[1,4] - k[1,4] <= 0
 NoHeightCompat_1_4: - 155 h[1,4] + k[1,4] <= 0
 Overlap_1_4: - x_p[1,4] - x_p[4,1] + o[1,4] <= 0
 NoOverlap_1_4: x_p[1,4] + x_p[4,1] - 2 o[1,4] <= 0
 3_Supports_1_4: - h[1,4] - o[1,4] - s[1,4] <= -1
 3_DoesntSupport_1_4: h[1,4] + o[1,4] + 2 s[1,4] <= 2
 beta1_1_4: - s[1,4] + beta1[1,4] <= 0
 beta2_1_4: - s[1,4] + beta2[1,4] <= 0
 sBeta1Flag_1_4: - s[1,4] + beta1[1,4] <= 0
 sBeta2Flag_1_4: - s[1,4] + beta2[1,4] <= 0
 Eta1Overlap_1_4: eta1[1,4] + beta1[1,4] <= 1
 Eta2Overlap_1_4: eta2[1,4] + beta2[1,4] <= 1
 Eta1Flag_1_4: - xi_l[1] + xi_l[4] - 300 eta1[1,4] <= 0
 Eta2Flag_1_4: xi_r[1] - xi_r[4] - 300 eta2[1,4] <= 0
 4_Supports_1_InSameBin1_0: p_ij[1,0] - p_ij[4,0] - s[1,4] <= 0
 4_Supports_1_InSameBin2_0: - p_ij[1,0] + p_ij[4,0] - s[1,4] <= 0
 Gamma0ForItems_1_inBins_0_NoCut_4: p_ij[1,0] + gamma[1] <= 1
 4_Supports_1_InSameBin1_1: p_ij[1,1] - p_ij[4,1] - s[1,4] <= 0
 4_Supports_1_InSameBin2_1: - p_ij[1,1] + p_ij[4,1] - s[1,4] <= 0
 Gamma0ForItems_1_inBins_1_NoCut_4: p_ij[1,1] + gamma[1] <= 1
 4_Supports_1_InSameBin1_2: p_ij[1,2] - p_ij[4,2] - s[1,4] <= 0
 4_Supports_1_InSameBin2_2: - p_ij[1,2] + p_ij[4,2] - s[1,4] <= 0
 Gamma0ForItems_1_inBins_2_NoCut_4: p_ij[1,2] + gamma[1] <= 98
 4_Supports_1_InSameBin1_3: p_ij[1,3] - p_ij[4,3] - s[1,4] <= 0
 4_Supports_1_InSameBin2_3: - p_ij[1,3] + p_ij[4,3] - s[1,4] <= 0
 Gamma0ForItems_1_inBins_3_NoCut_4: p_ij[1,3] + gamma[1] <= 98
 1_OnCutIn2_4: 10000 p_ij[1,2] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 1_OnCutIn3_4: 10000 p_ij[1,3] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 AbsZ11_5: - zi[1] + zi_hi[5] - k[1,5] <= 0
 AbsZ21_5: zi[1] - zi_hi[5] - k[1,5] <= 0
 5_Below_1: zi[1] - zi_hi[5] + k[1,5] + 310 m[1,5] <= 310
 5_NotBelow_1: - zi[1] + zi_hi[5] + k[1,5] - 310 m[1,5] <= 0
 HeightCompat_1_5: h[1,5] - k[1,5] <= 0
 NoHeightCompat_1_5: - 155 h[1,5] + k[1,5] <= 0
 Overlap_1_5: - x_p[1,5] - x_p[5,1] + o[1,5] <= 0
 NoOverlap_1_5: x_p[1,5] + x_p[5,1] - 2 o[1,5] <= 0
 3_Supports_1_5: - h[1,5] - o[1,5] - s[1,5] <= -1
 3_DoesntSupport_1_5: h[1,5] + o[1,5] + 2 s[1,5] <= 2
 beta1_1_5: - s[1,5] + beta1[1,5] <= 0
 beta2_1_5: - s[1,5] + beta2[1,5] <= 0
 sBeta1Flag_1_5: - s[1,5] + beta1[1,5] <= 0
 sBeta2Flag_1_5: - s[1,5] + beta2[1,5] <= 0
 Eta1Overlap_1_5: eta1[1,5] + beta1[1,5] <= 1
 Eta2Overlap_1_5: eta2[1,5] + beta2[1,5] <= 1
 Eta1Flag_1_5: - xi_l[1] + xi_l[5] - 300 eta1[1,5] <= 0
 Eta2Flag_1_5: xi_r[1] - xi_r[5] - 300 eta2[1,5] <= 0
 5_Supports_1_InSameBin1_0: p_ij[1,0] - p_ij[5,0] - s[1,5] <= 0
 5_Supports_1_InSameBin2_0: - p_ij[1,0] + p_ij[5,0] - s[1,5] <= 0
 Gamma0ForItems_1_inBins_0_NoCut_5: p_ij[1,0] + gamma[1] <= 1
 5_Supports_1_InSameBin1_1: p_ij[1,1] - p_ij[5,1] - s[1,5] <= 0
 5_Supports_1_InSameBin2_1: - p_ij[1,1] + p_ij[5,1] - s[1,5] <= 0
 Gamma0ForItems_1_inBins_1_NoCut_5: p_ij[1,1] + gamma[1] <= 1
 5_Supports_1_InSameBin1_2: p_ij[1,2] - p_ij[5,2] - s[1,5] <= 0
 5_Supports_1_InSameBin2_2: - p_ij[1,2] + p_ij[5,2] - s[1,5] <= 0
 Gamma0ForItems_1_inBins_2_NoCut_5: p_ij[1,2] + gamma[1] <= 98
 5_Supports_1_InSameBin1_3: p_ij[1,3] - p_ij[5,3] - s[1,5] <= 0
 5_Supports_1_InSameBin2_3: - p_ij[1,3] + p_ij[5,3] - s[1,5] <= 0
 Gamma0ForItems_1_inBins_3_NoCut_5: p_ij[1,3] + gamma[1] <= 98
 1_OnCutIn2_5: 10000 p_ij[1,2] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 1_OnCutIn3_5: 10000 p_ij[1,3] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 AbsZ11_6: - zi[1] + zi_hi[6] - k[1,6] <= 0
 AbsZ21_6: zi[1] - zi_hi[6] - k[1,6] <= 0
 6_Below_1: zi[1] - zi_hi[6] + k[1,6] + 310 m[1,6] <= 310
 6_NotBelow_1: - zi[1] + zi_hi[6] + k[1,6] - 310 m[1,6] <= 0
 HeightCompat_1_6: h[1,6] - k[1,6] <= 0
 NoHeightCompat_1_6: - 155 h[1,6] + k[1,6] <= 0
 Overlap_1_6: - x_p[1,6] - x_p[6,1] + o[1,6] <= 0
 NoOverlap_1_6: x_p[1,6] + x_p[6,1] - 2 o[1,6] <= 0
 3_Supports_1_6: - h[1,6] - o[1,6] - s[1,6] <= -1
 3_DoesntSupport_1_6: h[1,6] + o[1,6] + 2 s[1,6] <= 2
 beta1_1_6: - s[1,6] + beta1[1,6] <= 0
 beta2_1_6: - s[1,6] + beta2[1,6] <= 0
 sBeta1Flag_1_6: - s[1,6] + beta1[1,6] <= 0
 sBeta2Flag_1_6: - s[1,6] + beta2[1,6] <= 0
 Eta1Overlap_1_6: eta1[1,6] + beta1[1,6] <= 1
 Eta2Overlap_1_6: eta2[1,6] + beta2[1,6] <= 1
 Eta1Flag_1_6: - xi_l[1] + xi_l[6] - 300 eta1[1,6] <= 0
 Eta2Flag_1_6: xi_r[1] - xi_r[6] - 300 eta2[1,6] <= 0
 6_Supports_1_InSameBin1_0: p_ij[1,0] - p_ij[6,0] - s[1,6] <= 0
 6_Supports_1_InSameBin2_0: - p_ij[1,0] + p_ij[6,0] - s[1,6] <= 0
 Gamma0ForItems_1_inBins_0_NoCut_6: p_ij[1,0] + gamma[1] <= 1
 6_Supports_1_InSameBin1_1: p_ij[1,1] - p_ij[6,1] - s[1,6] <= 0
 6_Supports_1_InSameBin2_1: - p_ij[1,1] + p_ij[6,1] - s[1,6] <= 0
 Gamma0ForItems_1_inBins_1_NoCut_6: p_ij[1,1] + gamma[1] <= 1
 6_Supports_1_InSameBin1_2: p_ij[1,2] - p_ij[6,2] - s[1,6] <= 0
 6_Supports_1_InSameBin2_2: - p_ij[1,2] + p_ij[6,2] - s[1,6] <= 0
 Gamma0ForItems_1_inBins_2_NoCut_6: p_ij[1,2] + gamma[1] <= 98
 6_Supports_1_InSameBin1_3: p_ij[1,3] - p_ij[6,3] - s[1,6] <= 0
 6_Supports_1_InSameBin2_3: - p_ij[1,3] + p_ij[6,3] - s[1,6] <= 0
 Gamma0ForItems_1_inBins_3_NoCut_6: p_ij[1,3] + gamma[1] <= 98
 1_OnCutIn2_6: 10000 p_ij[1,2] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 1_OnCutIn3_6: 10000 p_ij[1,3] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 AbsZ11_7: - zi[1] + zi_hi[7] - k[1,7] <= 0
 AbsZ21_7: zi[1] - zi_hi[7] - k[1,7] <= 0
 7_Below_1: zi[1] - zi_hi[7] + k[1,7] + 310 m[1,7] <= 310
 7_NotBelow_1: - zi[1] + zi_hi[7] + k[1,7] - 310 m[1,7] <= 0
 HeightCompat_1_7: h[1,7] - k[1,7] <= 0
 NoHeightCompat_1_7: - 155 h[1,7] + k[1,7] <= 0
 Overlap_1_7: - x_p[1,7] - x_p[7,1] + o[1,7] <= 0
 NoOverlap_1_7: x_p[1,7] + x_p[7,1] - 2 o[1,7] <= 0
 3_Supports_1_7: - h[1,7] - o[1,7] - s[1,7] <= -1
 3_DoesntSupport_1_7: h[1,7] + o[1,7] + 2 s[1,7] <= 2
 beta1_1_7: - s[1,7] + beta1[1,7] <= 0
 beta2_1_7: - s[1,7] + beta2[1,7] <= 0
 sBeta1Flag_1_7: - s[1,7] + beta1[1,7] <= 0
 sBeta2Flag_1_7: - s[1,7] + beta2[1,7] <= 0
 Eta1Overlap_1_7: eta1[1,7] + beta1[1,7] <= 1
 Eta2Overlap_1_7: eta2[1,7] + beta2[1,7] <= 1
 Eta1Flag_1_7: - xi_l[1] + xi_l[7] - 300 eta1[1,7] <= 0
 Eta2Flag_1_7: xi_r[1] - xi_r[7] - 300 eta2[1,7] <= 0
 7_Supports_1_InSameBin1_0: p_ij[1,0] - p_ij[7,0] - s[1,7] <= 0
 7_Supports_1_InSameBin2_0: - p_ij[1,0] + p_ij[7,0] - s[1,7] <= 0
 Gamma0ForItems_1_inBins_0_NoCut_7: p_ij[1,0] + gamma[1] <= 1
 7_Supports_1_InSameBin1_1: p_ij[1,1] - p_ij[7,1] - s[1,7] <= 0
 7_Supports_1_InSameBin2_1: - p_ij[1,1] + p_ij[7,1] - s[1,7] <= 0
 Gamma0ForItems_1_inBins_1_NoCut_7: p_ij[1,1] + gamma[1] <= 1
 7_Supports_1_InSameBin1_2: p_ij[1,2] - p_ij[7,2] - s[1,7] <= 0
 7_Supports_1_InSameBin2_2: - p_ij[1,2] + p_ij[7,2] - s[1,7] <= 0
 Gamma0ForItems_1_inBins_2_NoCut_7: p_ij[1,2] + gamma[1] <= 98
 7_Supports_1_InSameBin1_3: p_ij[1,3] - p_ij[7,3] - s[1,7] <= 0
 7_Supports_1_InSameBin2_3: - p_ij[1,3] + p_ij[7,3] - s[1,7] <= 0
 Gamma0ForItems_1_inBins_3_NoCut_7: p_ij[1,3] + gamma[1] <= 98
 1_OnCutIn2_7: 10000 p_ij[1,2] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 1_OnCutIn3_7: 10000 p_ij[1,3] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 AbsZ11_8: - zi[1] + zi_hi[8] - k[1,8] <= 0
 AbsZ21_8: zi[1] - zi_hi[8] - k[1,8] <= 0
 8_Below_1: zi[1] - zi_hi[8] + k[1,8] + 310 m[1,8] <= 310
 8_NotBelow_1: - zi[1] + zi_hi[8] + k[1,8] - 310 m[1,8] <= 0
 HeightCompat_1_8: h[1,8] - k[1,8] <= 0
 NoHeightCompat_1_8: - 155 h[1,8] + k[1,8] <= 0
 Overlap_1_8: - x_p[1,8] - x_p[8,1] + o[1,8] <= 0
 NoOverlap_1_8: x_p[1,8] + x_p[8,1] - 2 o[1,8] <= 0
 3_Supports_1_8: - h[1,8] - o[1,8] - s[1,8] <= -1
 3_DoesntSupport_1_8: h[1,8] + o[1,8] + 2 s[1,8] <= 2
 beta1_1_8: - s[1,8] + beta1[1,8] <= 0
 beta2_1_8: - s[1,8] + beta2[1,8] <= 0
 sBeta1Flag_1_8: - s[1,8] + beta1[1,8] <= 0
 sBeta2Flag_1_8: - s[1,8] + beta2[1,8] <= 0
 Eta1Overlap_1_8: eta1[1,8] + beta1[1,8] <= 1
 Eta2Overlap_1_8: eta2[1,8] + beta2[1,8] <= 1
 Eta1Flag_1_8: - xi_l[1] + xi_l[8] - 300 eta1[1,8] <= 0
 Eta2Flag_1_8: xi_r[1] - xi_r[8] - 300 eta2[1,8] <= 0
 8_Supports_1_InSameBin1_0: p_ij[1,0] - p_ij[8,0] - s[1,8] <= 0
 8_Supports_1_InSameBin2_0: - p_ij[1,0] + p_ij[8,0] - s[1,8] <= 0
 Gamma0ForItems_1_inBins_0_NoCut_8: p_ij[1,0] + gamma[1] <= 1
 8_Supports_1_InSameBin1_1: p_ij[1,1] - p_ij[8,1] - s[1,8] <= 0
 8_Supports_1_InSameBin2_1: - p_ij[1,1] + p_ij[8,1] - s[1,8] <= 0
 Gamma0ForItems_1_inBins_1_NoCut_8: p_ij[1,1] + gamma[1] <= 1
 8_Supports_1_InSameBin1_2: p_ij[1,2] - p_ij[8,2] - s[1,8] <= 0
 8_Supports_1_InSameBin2_2: - p_ij[1,2] + p_ij[8,2] - s[1,8] <= 0
 Gamma0ForItems_1_inBins_2_NoCut_8: p_ij[1,2] + gamma[1] <= 98
 8_Supports_1_InSameBin1_3: p_ij[1,3] - p_ij[8,3] - s[1,8] <= 0
 8_Supports_1_InSameBin2_3: - p_ij[1,3] + p_ij[8,3] - s[1,8] <= 0
 Gamma0ForItems_1_inBins_3_NoCut_8: p_ij[1,3] + gamma[1] <= 98
 1_OnCutIn2_8: 10000 p_ij[1,2] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 1_OnCutIn3_8: 10000 p_ij[1,3] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 AbsZ11_9: - zi[1] + zi_hi[9] - k[1,9] <= 0
 AbsZ21_9: zi[1] - zi_hi[9] - k[1,9] <= 0
 9_Below_1: zi[1] - zi_hi[9] + k[1,9] + 310 m[1,9] <= 310
 9_NotBelow_1: - zi[1] + zi_hi[9] + k[1,9] - 310 m[1,9] <= 0
 HeightCompat_1_9: h[1,9] - k[1,9] <= 0
 NoHeightCompat_1_9: - 155 h[1,9] + k[1,9] <= 0
 Overlap_1_9: - x_p[1,9] - x_p[9,1] + o[1,9] <= 0
 NoOverlap_1_9: x_p[1,9] + x_p[9,1] - 2 o[1,9] <= 0
 3_Supports_1_9: - h[1,9] - o[1,9] - s[1,9] <= -1
 3_DoesntSupport_1_9: h[1,9] + o[1,9] + 2 s[1,9] <= 2
 beta1_1_9: - s[1,9] + beta1[1,9] <= 0
 beta2_1_9: - s[1,9] + beta2[1,9] <= 0
 sBeta1Flag_1_9: - s[1,9] + beta1[1,9] <= 0
 sBeta2Flag_1_9: - s[1,9] + beta2[1,9] <= 0
 Eta1Overlap_1_9: eta1[1,9] + beta1[1,9] <= 1
 Eta2Overlap_1_9: eta2[1,9] + beta2[1,9] <= 1
 Eta1Flag_1_9: - xi_l[1] + xi_l[9] - 300 eta1[1,9] <= 0
 Eta2Flag_1_9: xi_r[1] - xi_r[9] - 300 eta2[1,9] <= 0
 9_Supports_1_InSameBin1_0: p_ij[1,0] - p_ij[9,0] - s[1,9] <= 0
 9_Supports_1_InSameBin2_0: - p_ij[1,0] + p_ij[9,0] - s[1,9] <= 0
 Gamma0ForItems_1_inBins_0_NoCut_9: p_ij[1,0] + gamma[1] <= 1
 9_Supports_1_InSameBin1_1: p_ij[1,1] - p_ij[9,1] - s[1,9] <= 0
 9_Supports_1_InSameBin2_1: - p_ij[1,1] + p_ij[9,1] - s[1,9] <= 0
 Gamma0ForItems_1_inBins_1_NoCut_9: p_ij[1,1] + gamma[1] <= 1
 9_Supports_1_InSameBin1_2: p_ij[1,2] - p_ij[9,2] - s[1,9] <= 0
 9_Supports_1_InSameBin2_2: - p_ij[1,2] + p_ij[9,2] - s[1,9] <= 0
 Gamma0ForItems_1_inBins_2_NoCut_9: p_ij[1,2] + gamma[1] <= 98
 9_Supports_1_InSameBin1_3: p_ij[1,3] - p_ij[9,3] - s[1,9] <= 0
 9_Supports_1_InSameBin2_3: - p_ij[1,3] + p_ij[9,3] - s[1,9] <= 0
 Gamma0ForItems_1_inBins_3_NoCut_9: p_ij[1,3] + gamma[1] <= 98
 1_OnCutIn2_9: 10000 p_ij[1,2] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 1_OnCutIn3_9: 10000 p_ij[1,3] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 AbsZ11_10: - zi[1] + zi_hi[10] - k[1,10] <= 0
 AbsZ21_10: zi[1] - zi_hi[10] - k[1,10] <= 0
 10_Below_1: zi[1] - zi_hi[10] + k[1,10] + 310 m[1,10] <= 310
 10_NotBelow_1: - zi[1] + zi_hi[10] + k[1,10] - 310 m[1,10] <= 0
 HeightCompat_1_10: h[1,10] - k[1,10] <= 0
 NoHeightCompat_1_10: - 155 h[1,10] + k[1,10] <= 0
 Overlap_1_10: - x_p[1,10] - x_p[10,1] + o[1,10] <= 0
 NoOverlap_1_10: x_p[1,10] + x_p[10,1] - 2 o[1,10] <= 0
 3_Supports_1_10: - h[1,10] - o[1,10] - s[1,10] <= -1
 3_DoesntSupport_1_10: h[1,10] + o[1,10] + 2 s[1,10] <= 2
 beta1_1_10: - s[1,10] + beta1[1,10] <= 0
 beta2_1_10: - s[1,10] + beta2[1,10] <= 0
 sBeta1Flag_1_10: - s[1,10] + beta1[1,10] <= 0
 sBeta2Flag_1_10: - s[1,10] + beta2[1,10] <= 0
 Eta1Overlap_1_10: eta1[1,10] + beta1[1,10] <= 1
 Eta2Overlap_1_10: eta2[1,10] + beta2[1,10] <= 1
 Eta1Flag_1_10: - xi_l[1] + xi_l[10] - 300 eta1[1,10] <= 0
 Eta2Flag_1_10: xi_r[1] - xi_r[10] - 300 eta2[1,10] <= 0
 10_Supports_1_InSameBin1_0: p_ij[1,0] - p_ij[10,0] - s[1,10] <= 0
 10_Supports_1_InSameBin2_0: - p_ij[1,0] + p_ij[10,0] - s[1,10] <= 0
 Gamma0ForItems_1_inBins_0_NoCut_10: p_ij[1,0] + gamma[1] <= 1
 10_Supports_1_InSameBin1_1: p_ij[1,1] - p_ij[10,1] - s[1,10] <= 0
 10_Supports_1_InSameBin2_1: - p_ij[1,1] + p_ij[10,1] - s[1,10] <= 0
 Gamma0ForItems_1_inBins_1_NoCut_10: p_ij[1,1] + gamma[1] <= 1
 10_Supports_1_InSameBin1_2: p_ij[1,2] - p_ij[10,2] - s[1,10] <= 0
 10_Supports_1_InSameBin2_2: - p_ij[1,2] + p_ij[10,2] - s[1,10] <= 0
 Gamma0ForItems_1_inBins_2_NoCut_10: p_ij[1,2] + gamma[1] <= 98
 10_Supports_1_InSameBin1_3: p_ij[1,3] - p_ij[10,3] - s[1,10] <= 0
 10_Supports_1_InSameBin2_3: - p_ij[1,3] + p_ij[10,3] - s[1,10] <= 0
 Gamma0ForItems_1_inBins_3_NoCut_10: p_ij[1,3] + gamma[1] <= 98
 1_OnCutIn2_10: 10000 p_ij[1,2] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 1_OnCutIn3_10: 10000 p_ij[1,3] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 AbsZ11_11: - zi[1] + zi_hi[11] - k[1,11] <= 0
 AbsZ21_11: zi[1] - zi_hi[11] - k[1,11] <= 0
 11_Below_1: zi[1] - zi_hi[11] + k[1,11] + 310 m[1,11] <= 310
 11_NotBelow_1: - zi[1] + zi_hi[11] + k[1,11] - 310 m[1,11] <= 0
 HeightCompat_1_11: h[1,11] - k[1,11] <= 0
 NoHeightCompat_1_11: - 155 h[1,11] + k[1,11] <= 0
 Overlap_1_11: - x_p[1,11] - x_p[11,1] + o[1,11] <= 0
 NoOverlap_1_11: x_p[1,11] + x_p[11,1] - 2 o[1,11] <= 0
 3_Supports_1_11: - h[1,11] - o[1,11] - s[1,11] <= -1
 3_DoesntSupport_1_11: h[1,11] + o[1,11] + 2 s[1,11] <= 2
 beta1_1_11: - s[1,11] + beta1[1,11] <= 0
 beta2_1_11: - s[1,11] + beta2[1,11] <= 0
 sBeta1Flag_1_11: - s[1,11] + beta1[1,11] <= 0
 sBeta2Flag_1_11: - s[1,11] + beta2[1,11] <= 0
 Eta1Overlap_1_11: eta1[1,11] + beta1[1,11] <= 1
 Eta2Overlap_1_11: eta2[1,11] + beta2[1,11] <= 1
 Eta1Flag_1_11: - xi_l[1] + xi_l[11] - 300 eta1[1,11] <= 0
 Eta2Flag_1_11: xi_r[1] - xi_r[11] - 300 eta2[1,11] <= 0
 11_Supports_1_InSameBin1_0: p_ij[1,0] - p_ij[11,0] - s[1,11] <= 0
 11_Supports_1_InSameBin2_0: - p_ij[1,0] + p_ij[11,0] - s[1,11] <= 0
 Gamma0ForItems_1_inBins_0_NoCut_11: p_ij[1,0] + gamma[1] <= 1
 11_Supports_1_InSameBin1_1: p_ij[1,1] - p_ij[11,1] - s[1,11] <= 0
 11_Supports_1_InSameBin2_1: - p_ij[1,1] + p_ij[11,1] - s[1,11] <= 0
 Gamma0ForItems_1_inBins_1_NoCut_11: p_ij[1,1] + gamma[1] <= 1
 11_Supports_1_InSameBin1_2: p_ij[1,2] - p_ij[11,2] - s[1,11] <= 0
 11_Supports_1_InSameBin2_2: - p_ij[1,2] + p_ij[11,2] - s[1,11] <= 0
 Gamma0ForItems_1_inBins_2_NoCut_11: p_ij[1,2] + gamma[1] <= 98
 11_Supports_1_InSameBin1_3: p_ij[1,3] - p_ij[11,3] - s[1,11] <= 0
 11_Supports_1_InSameBin2_3: - p_ij[1,3] + p_ij[11,3] - s[1,11] <= 0
 Gamma0ForItems_1_inBins_3_NoCut_11: p_ij[1,3] + gamma[1] <= 98
 1_OnCutIn2_11: 10000 p_ij[1,2] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 1_OnCutIn3_11: 10000 p_ij[1,3] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 AbsZ11_12: - zi[1] + zi_hi[12] - k[1,12] <= 0
 AbsZ21_12: zi[1] - zi_hi[12] - k[1,12] <= 0
 12_Below_1: zi[1] - zi_hi[12] + k[1,12] + 310 m[1,12] <= 310
 12_NotBelow_1: - zi[1] + zi_hi[12] + k[1,12] - 310 m[1,12] <= 0
 HeightCompat_1_12: h[1,12] - k[1,12] <= 0
 NoHeightCompat_1_12: - 155 h[1,12] + k[1,12] <= 0
 Overlap_1_12: - x_p[1,12] - x_p[12,1] + o[1,12] <= 0
 NoOverlap_1_12: x_p[1,12] + x_p[12,1] - 2 o[1,12] <= 0
 3_Supports_1_12: - h[1,12] - o[1,12] - s[1,12] <= -1
 3_DoesntSupport_1_12: h[1,12] + o[1,12] + 2 s[1,12] <= 2
 beta1_1_12: - s[1,12] + beta1[1,12] <= 0
 beta2_1_12: - s[1,12] + beta2[1,12] <= 0
 sBeta1Flag_1_12: - s[1,12] + beta1[1,12] <= 0
 sBeta2Flag_1_12: - s[1,12] + beta2[1,12] <= 0
 Eta1Overlap_1_12: eta1[1,12] + beta1[1,12] <= 1
 Eta2Overlap_1_12: eta2[1,12] + beta2[1,12] <= 1
 Eta1Flag_1_12: - xi_l[1] + xi_l[12] - 300 eta1[1,12] <= 0
 Eta2Flag_1_12: xi_r[1] - xi_r[12] - 300 eta2[1,12] <= 0
 12_Supports_1_InSameBin1_0: p_ij[1,0] - p_ij[12,0] - s[1,12] <= 0
 12_Supports_1_InSameBin2_0: - p_ij[1,0] + p_ij[12,0] - s[1,12] <= 0
 Gamma0ForItems_1_inBins_0_NoCut_12: p_ij[1,0] + gamma[1] <= 1
 12_Supports_1_InSameBin1_1: p_ij[1,1] - p_ij[12,1] - s[1,12] <= 0
 12_Supports_1_InSameBin2_1: - p_ij[1,1] + p_ij[12,1] - s[1,12] <= 0
 Gamma0ForItems_1_inBins_1_NoCut_12: p_ij[1,1] + gamma[1] <= 1
 12_Supports_1_InSameBin1_2: p_ij[1,2] - p_ij[12,2] - s[1,12] <= 0
 12_Supports_1_InSameBin2_2: - p_ij[1,2] + p_ij[12,2] - s[1,12] <= 0
 Gamma0ForItems_1_inBins_2_NoCut_12: p_ij[1,2] + gamma[1] <= 98
 12_Supports_1_InSameBin1_3: p_ij[1,3] - p_ij[12,3] - s[1,12] <= 0
 12_Supports_1_InSameBin2_3: - p_ij[1,3] + p_ij[12,3] - s[1,12] <= 0
 Gamma0ForItems_1_inBins_3_NoCut_12: p_ij[1,3] + gamma[1] <= 98
 1_OnCutIn2_12: 10000 p_ij[1,2] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 1_OnCutIn3_12: 10000 p_ij[1,3] + 1.261904761904762 xi_l[1] + zi[1]
   + 10000 gamma[1] <= 20053
 GroundConstraint_2: zi[2] + 155 g[2] <= 155
 AbsZ12_0: - zi[2] + zi_hi[0] - k[2,0] <= 0
 AbsZ22_0: zi[2] - zi_hi[0] - k[2,0] <= 0
 0_Below_2: zi[2] - zi_hi[0] + k[2,0] + 310 m[2,0] <= 310
 0_NotBelow_2: - zi[2] + zi_hi[0] + k[2,0] - 310 m[2,0] <= 0
 HeightCompat_2_0: h[2,0] - k[2,0] <= 0
 NoHeightCompat_2_0: - 155 h[2,0] + k[2,0] <= 0
 Overlap_2_0: - x_p[0,2] - x_p[2,0] + o[2,0] <= 0
 NoOverlap_2_0: x_p[0,2] + x_p[2,0] - 2 o[2,0] <= 0
 3_Supports_2_0: - h[2,0] - o[2,0] - s[2,0] <= -1
 3_DoesntSupport_2_0: h[2,0] + o[2,0] + 2 s[2,0] <= 2
 beta1_2_0: - s[2,0] + beta1[2,0] <= 0
 beta2_2_0: - s[2,0] + beta2[2,0] <= 0
 sBeta1Flag_2_0: - s[2,0] + beta1[2,0] <= 0
 sBeta2Flag_2_0: - s[2,0] + beta2[2,0] <= 0
 Eta1Overlap_2_0: eta1[2,0] + beta1[2,0] <= 1
 Eta2Overlap_2_0: eta2[2,0] + beta2[2,0] <= 1
 Eta1Flag_2_0: xi_l[0] - xi_l[2] - 300 eta1[2,0] <= 0
 Eta2Flag_2_0: - xi_r[0] + xi_r[2] - 300 eta2[2,0] <= 0
 0_Supports_2_InSameBin1_0: - p_ij[0,0] + p_ij[2,0] - s[2,0] <= 0
 0_Supports_2_InSameBin2_0: p_ij[0,0] - p_ij[2,0] - s[2,0] <= 0
 Gamma0ForItems_2_inBins_0_NoCut_0: p_ij[2,0] + gamma[2] <= 1
 0_Supports_2_InSameBin1_1: - p_ij[0,1] + p_ij[2,1] - s[2,0] <= 0
 0_Supports_2_InSameBin2_1: p_ij[0,1] - p_ij[2,1] - s[2,0] <= 0
 Gamma0ForItems_2_inBins_1_NoCut_0: p_ij[2,1] + gamma[2] <= 1
 0_Supports_2_InSameBin1_2: - p_ij[0,2] + p_ij[2,2] - s[2,0] <= 0
 0_Supports_2_InSameBin2_2: p_ij[0,2] - p_ij[2,2] - s[2,0] <= 0
 Gamma0ForItems_2_inBins_2_NoCut_0: p_ij[2,2] + gamma[2] <= 98
 0_Supports_2_InSameBin1_3: - p_ij[0,3] + p_ij[2,3] - s[2,0] <= 0
 0_Supports_2_InSameBin2_3: p_ij[0,3] - p_ij[2,3] - s[2,0] <= 0
 Gamma0ForItems_2_inBins_3_NoCut_0: p_ij[2,3] + gamma[2] <= 98
 2_OnCutIn2_0: 10000 p_ij[2,2] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 2_OnCutIn3_0: 10000 p_ij[2,3] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 AbsZ12_1: - zi[2] + zi_hi[1] - k[2,1] <= 0
 AbsZ22_1: zi[2] - zi_hi[1] - k[2,1] <= 0
 1_Below_2: zi[2] - zi_hi[1] + k[2,1] + 310 m[2,1] <= 310
 1_NotBelow_2: - zi[2] + zi_hi[1] + k[2,1] - 310 m[2,1] <= 0
 HeightCompat_2_1: h[2,1] - k[2,1] <= 0
 NoHeightCompat_2_1: - 155 h[2,1] + k[2,1] <= 0
 Overlap_2_1: - x_p[1,2] - x_p[2,1] + o[2,1] <= 0
 NoOverlap_2_1: x_p[1,2] + x_p[2,1] - 2 o[2,1] <= 0
 3_Supports_2_1: - h[2,1] - o[2,1] - s[2,1] <= -1
 3_DoesntSupport_2_1: h[2,1] + o[2,1] + 2 s[2,1] <= 2
 beta1_2_1: - s[2,1] + beta1[2,1] <= 0
 beta2_2_1: - s[2,1] + beta2[2,1] <= 0
 sBeta1Flag_2_1: - s[2,1] + beta1[2,1] <= 0
 sBeta2Flag_2_1: - s[2,1] + beta2[2,1] <= 0
 Eta1Overlap_2_1: eta1[2,1] + beta1[2,1] <= 1
 Eta2Overlap_2_1: eta2[2,1] + beta2[2,1] <= 1
 Eta1Flag_2_1: xi_l[1] - xi_l[2] - 300 eta1[2,1] <= 0
 Eta2Flag_2_1: - xi_r[1] + xi_r[2] - 300 eta2[2,1] <= 0
 1_Supports_2_InSameBin1_0: - p_ij[1,0] + p_ij[2,0] - s[2,1] <= 0
 1_Supports_2_InSameBin2_0: p_ij[1,0] - p_ij[2,0] - s[2,1] <= 0
 Gamma0ForItems_2_inBins_0_NoCut_1: p_ij[2,0] + gamma[2] <= 1
 1_Supports_2_InSameBin1_1: - p_ij[1,1] + p_ij[2,1] - s[2,1] <= 0
 1_Supports_2_InSameBin2_1: p_ij[1,1] - p_ij[2,1] - s[2,1] <= 0
 Gamma0ForItems_2_inBins_1_NoCut_1: p_ij[2,1] + gamma[2] <= 1
 1_Supports_2_InSameBin1_2: - p_ij[1,2] + p_ij[2,2] - s[2,1] <= 0
 1_Supports_2_InSameBin2_2: p_ij[1,2] - p_ij[2,2] - s[2,1] <= 0
 Gamma0ForItems_2_inBins_2_NoCut_1: p_ij[2,2] + gamma[2] <= 98
 1_Supports_2_InSameBin1_3: - p_ij[1,3] + p_ij[2,3] - s[2,1] <= 0
 1_Supports_2_InSameBin2_3: p_ij[1,3] - p_ij[2,3] - s[2,1] <= 0
 Gamma0ForItems_2_inBins_3_NoCut_1: p_ij[2,3] + gamma[2] <= 98
 2_OnCutIn2_1: 10000 p_ij[2,2] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 2_OnCutIn3_1: 10000 p_ij[2,3] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 2_OnCutIn2_2: 10000 p_ij[2,2] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 2_OnCutIn3_2: 10000 p_ij[2,3] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 AbsZ12_3: - zi[2] + zi_hi[3] - k[2,3] <= 0
 AbsZ22_3: zi[2] - zi_hi[3] - k[2,3] <= 0
 3_Below_2: zi[2] - zi_hi[3] + k[2,3] + 310 m[2,3] <= 310
 3_NotBelow_2: - zi[2] + zi_hi[3] + k[2,3] - 310 m[2,3] <= 0
 HeightCompat_2_3: h[2,3] - k[2,3] <= 0
 NoHeightCompat_2_3: - 155 h[2,3] + k[2,3] <= 0
 Overlap_2_3: - x_p[2,3] - x_p[3,2] + o[2,3] <= 0
 NoOverlap_2_3: x_p[2,3] + x_p[3,2] - 2 o[2,3] <= 0
 3_Supports_2_3: - h[2,3] - o[2,3] - s[2,3] <= -1
 3_DoesntSupport_2_3: h[2,3] + o[2,3] + 2 s[2,3] <= 2
 beta1_2_3: - s[2,3] + beta1[2,3] <= 0
 beta2_2_3: - s[2,3] + beta2[2,3] <= 0
 sBeta1Flag_2_3: - s[2,3] + beta1[2,3] <= 0
 sBeta2Flag_2_3: - s[2,3] + beta2[2,3] <= 0
 Eta1Overlap_2_3: eta1[2,3] + beta1[2,3] <= 1
 Eta2Overlap_2_3: eta2[2,3] + beta2[2,3] <= 1
 Eta1Flag_2_3: - xi_l[2] + xi_l[3] - 300 eta1[2,3] <= 0
 Eta2Flag_2_3: xi_r[2] - xi_r[3] - 300 eta2[2,3] <= 0
 3_Supports_2_InSameBin1_0: p_ij[2,0] - p_ij[3,0] - s[2,3] <= 0
 3_Supports_2_InSameBin2_0: - p_ij[2,0] + p_ij[3,0] - s[2,3] <= 0
 Gamma0ForItems_2_inBins_0_NoCut_3: p_ij[2,0] + gamma[2] <= 1
 3_Supports_2_InSameBin1_1: p_ij[2,1] - p_ij[3,1] - s[2,3] <= 0
 3_Supports_2_InSameBin2_1: - p_ij[2,1] + p_ij[3,1] - s[2,3] <= 0
 Gamma0ForItems_2_inBins_1_NoCut_3: p_ij[2,1] + gamma[2] <= 1
 3_Supports_2_InSameBin1_2: p_ij[2,2] - p_ij[3,2] - s[2,3] <= 0
 3_Supports_2_InSameBin2_2: - p_ij[2,2] + p_ij[3,2] - s[2,3] <= 0
 Gamma0ForItems_2_inBins_2_NoCut_3: p_ij[2,2] + gamma[2] <= 98
 3_Supports_2_InSameBin1_3: p_ij[2,3] - p_ij[3,3] - s[2,3] <= 0
 3_Supports_2_InSameBin2_3: - p_ij[2,3] + p_ij[3,3] - s[2,3] <= 0
 Gamma0ForItems_2_inBins_3_NoCut_3: p_ij[2,3] + gamma[2] <= 98
 2_OnCutIn2_3: 10000 p_ij[2,2] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 2_OnCutIn3_3: 10000 p_ij[2,3] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 AbsZ12_4: - zi[2] + zi_hi[4] - k[2,4] <= 0
 AbsZ22_4: zi[2] - zi_hi[4] - k[2,4] <= 0
 4_Below_2: zi[2] - zi_hi[4] + k[2,4] + 310 m[2,4] <= 310
 4_NotBelow_2: - zi[2] + zi_hi[4] + k[2,4] - 310 m[2,4] <= 0
 HeightCompat_2_4: h[2,4] - k[2,4] <= 0
 NoHeightCompat_2_4: - 155 h[2,4] + k[2,4] <= 0
 Overlap_2_4: - x_p[2,4] - x_p[4,2] + o[2,4] <= 0
 NoOverlap_2_4: x_p[2,4] + x_p[4,2] - 2 o[2,4] <= 0
 3_Supports_2_4: - h[2,4] - o[2,4] - s[2,4] <= -1
 3_DoesntSupport_2_4: h[2,4] + o[2,4] + 2 s[2,4] <= 2
 beta1_2_4: - s[2,4] + beta1[2,4] <= 0
 beta2_2_4: - s[2,4] + beta2[2,4] <= 0
 sBeta1Flag_2_4: - s[2,4] + beta1[2,4] <= 0
 sBeta2Flag_2_4: - s[2,4] + beta2[2,4] <= 0
 Eta1Overlap_2_4: eta1[2,4] + beta1[2,4] <= 1
 Eta2Overlap_2_4: eta2[2,4] + beta2[2,4] <= 1
 Eta1Flag_2_4: - xi_l[2] + xi_l[4] - 300 eta1[2,4] <= 0
 Eta2Flag_2_4: xi_r[2] - xi_r[4] - 300 eta2[2,4] <= 0
 4_Supports_2_InSameBin1_0: p_ij[2,0] - p_ij[4,0] - s[2,4] <= 0
 4_Supports_2_InSameBin2_0: - p_ij[2,0] + p_ij[4,0] - s[2,4] <= 0
 Gamma0ForItems_2_inBins_0_NoCut_4: p_ij[2,0] + gamma[2] <= 1
 4_Supports_2_InSameBin1_1: p_ij[2,1] - p_ij[4,1] - s[2,4] <= 0
 4_Supports_2_InSameBin2_1: - p_ij[2,1] + p_ij[4,1] - s[2,4] <= 0
 Gamma0ForItems_2_inBins_1_NoCut_4: p_ij[2,1] + gamma[2] <= 1
 4_Supports_2_InSameBin1_2: p_ij[2,2] - p_ij[4,2] - s[2,4] <= 0
 4_Supports_2_InSameBin2_2: - p_ij[2,2] + p_ij[4,2] - s[2,4] <= 0
 Gamma0ForItems_2_inBins_2_NoCut_4: p_ij[2,2] + gamma[2] <= 98
 4_Supports_2_InSameBin1_3: p_ij[2,3] - p_ij[4,3] - s[2,4] <= 0
 4_Supports_2_InSameBin2_3: - p_ij[2,3] + p_ij[4,3] - s[2,4] <= 0
 Gamma0ForItems_2_inBins_3_NoCut_4: p_ij[2,3] + gamma[2] <= 98
 2_OnCutIn2_4: 10000 p_ij[2,2] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 2_OnCutIn3_4: 10000 p_ij[2,3] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 AbsZ12_5: - zi[2] + zi_hi[5] - k[2,5] <= 0
 AbsZ22_5: zi[2] - zi_hi[5] - k[2,5] <= 0
 5_Below_2: zi[2] - zi_hi[5] + k[2,5] + 310 m[2,5] <= 310
 5_NotBelow_2: - zi[2] + zi_hi[5] + k[2,5] - 310 m[2,5] <= 0
 HeightCompat_2_5: h[2,5] - k[2,5] <= 0
 NoHeightCompat_2_5: - 155 h[2,5] + k[2,5] <= 0
 Overlap_2_5: - x_p[2,5] - x_p[5,2] + o[2,5] <= 0
 NoOverlap_2_5: x_p[2,5] + x_p[5,2] - 2 o[2,5] <= 0
 3_Supports_2_5: - h[2,5] - o[2,5] - s[2,5] <= -1
 3_DoesntSupport_2_5: h[2,5] + o[2,5] + 2 s[2,5] <= 2
 beta1_2_5: - s[2,5] + beta1[2,5] <= 0
 beta2_2_5: - s[2,5] + beta2[2,5] <= 0
 sBeta1Flag_2_5: - s[2,5] + beta1[2,5] <= 0
 sBeta2Flag_2_5: - s[2,5] + beta2[2,5] <= 0
 Eta1Overlap_2_5: eta1[2,5] + beta1[2,5] <= 1
 Eta2Overlap_2_5: eta2[2,5] + beta2[2,5] <= 1
 Eta1Flag_2_5: - xi_l[2] + xi_l[5] - 300 eta1[2,5] <= 0
 Eta2Flag_2_5: xi_r[2] - xi_r[5] - 300 eta2[2,5] <= 0
 5_Supports_2_InSameBin1_0: p_ij[2,0] - p_ij[5,0] - s[2,5] <= 0
 5_Supports_2_InSameBin2_0: - p_ij[2,0] + p_ij[5,0] - s[2,5] <= 0
 Gamma0ForItems_2_inBins_0_NoCut_5: p_ij[2,0] + gamma[2] <= 1
 5_Supports_2_InSameBin1_1: p_ij[2,1] - p_ij[5,1] - s[2,5] <= 0
 5_Supports_2_InSameBin2_1: - p_ij[2,1] + p_ij[5,1] - s[2,5] <= 0
 Gamma0ForItems_2_inBins_1_NoCut_5: p_ij[2,1] + gamma[2] <= 1
 5_Supports_2_InSameBin1_2: p_ij[2,2] - p_ij[5,2] - s[2,5] <= 0
 5_Supports_2_InSameBin2_2: - p_ij[2,2] + p_ij[5,2] - s[2,5] <= 0
 Gamma0ForItems_2_inBins_2_NoCut_5: p_ij[2,2] + gamma[2] <= 98
 5_Supports_2_InSameBin1_3: p_ij[2,3] - p_ij[5,3] - s[2,5] <= 0
 5_Supports_2_InSameBin2_3: - p_ij[2,3] + p_ij[5,3] - s[2,5] <= 0
 Gamma0ForItems_2_inBins_3_NoCut_5: p_ij[2,3] + gamma[2] <= 98
 2_OnCutIn2_5: 10000 p_ij[2,2] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 2_OnCutIn3_5: 10000 p_ij[2,3] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 AbsZ12_6: - zi[2] + zi_hi[6] - k[2,6] <= 0
 AbsZ22_6: zi[2] - zi_hi[6] - k[2,6] <= 0
 6_Below_2: zi[2] - zi_hi[6] + k[2,6] + 310 m[2,6] <= 310
 6_NotBelow_2: - zi[2] + zi_hi[6] + k[2,6] - 310 m[2,6] <= 0
 HeightCompat_2_6: h[2,6] - k[2,6] <= 0
 NoHeightCompat_2_6: - 155 h[2,6] + k[2,6] <= 0
 Overlap_2_6: - x_p[2,6] - x_p[6,2] + o[2,6] <= 0
 NoOverlap_2_6: x_p[2,6] + x_p[6,2] - 2 o[2,6] <= 0
 3_Supports_2_6: - h[2,6] - o[2,6] - s[2,6] <= -1
 3_DoesntSupport_2_6: h[2,6] + o[2,6] + 2 s[2,6] <= 2
 beta1_2_6: - s[2,6] + beta1[2,6] <= 0
 beta2_2_6: - s[2,6] + beta2[2,6] <= 0
 sBeta1Flag_2_6: - s[2,6] + beta1[2,6] <= 0
 sBeta2Flag_2_6: - s[2,6] + beta2[2,6] <= 0
 Eta1Overlap_2_6: eta1[2,6] + beta1[2,6] <= 1
 Eta2Overlap_2_6: eta2[2,6] + beta2[2,6] <= 1
 Eta1Flag_2_6: - xi_l[2] + xi_l[6] - 300 eta1[2,6] <= 0
 Eta2Flag_2_6: xi_r[2] - xi_r[6] - 300 eta2[2,6] <= 0
 6_Supports_2_InSameBin1_0: p_ij[2,0] - p_ij[6,0] - s[2,6] <= 0
 6_Supports_2_InSameBin2_0: - p_ij[2,0] + p_ij[6,0] - s[2,6] <= 0
 Gamma0ForItems_2_inBins_0_NoCut_6: p_ij[2,0] + gamma[2] <= 1
 6_Supports_2_InSameBin1_1: p_ij[2,1] - p_ij[6,1] - s[2,6] <= 0
 6_Supports_2_InSameBin2_1: - p_ij[2,1] + p_ij[6,1] - s[2,6] <= 0
 Gamma0ForItems_2_inBins_1_NoCut_6: p_ij[2,1] + gamma[2] <= 1
 6_Supports_2_InSameBin1_2: p_ij[2,2] - p_ij[6,2] - s[2,6] <= 0
 6_Supports_2_InSameBin2_2: - p_ij[2,2] + p_ij[6,2] - s[2,6] <= 0
 Gamma0ForItems_2_inBins_2_NoCut_6: p_ij[2,2] + gamma[2] <= 98
 6_Supports_2_InSameBin1_3: p_ij[2,3] - p_ij[6,3] - s[2,6] <= 0
 6_Supports_2_InSameBin2_3: - p_ij[2,3] + p_ij[6,3] - s[2,6] <= 0
 Gamma0ForItems_2_inBins_3_NoCut_6: p_ij[2,3] + gamma[2] <= 98
 2_OnCutIn2_6: 10000 p_ij[2,2] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 2_OnCutIn3_6: 10000 p_ij[2,3] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 AbsZ12_7: - zi[2] + zi_hi[7] - k[2,7] <= 0
 AbsZ22_7: zi[2] - zi_hi[7] - k[2,7] <= 0
 7_Below_2: zi[2] - zi_hi[7] + k[2,7] + 310 m[2,7] <= 310
 7_NotBelow_2: - zi[2] + zi_hi[7] + k[2,7] - 310 m[2,7] <= 0
 HeightCompat_2_7: h[2,7] - k[2,7] <= 0
 NoHeightCompat_2_7: - 155 h[2,7] + k[2,7] <= 0
 Overlap_2_7: - x_p[2,7] - x_p[7,2] + o[2,7] <= 0
 NoOverlap_2_7: x_p[2,7] + x_p[7,2] - 2 o[2,7] <= 0
 3_Supports_2_7: - h[2,7] - o[2,7] - s[2,7] <= -1
 3_DoesntSupport_2_7: h[2,7] + o[2,7] + 2 s[2,7] <= 2
 beta1_2_7: - s[2,7] + beta1[2,7] <= 0
 beta2_2_7: - s[2,7] + beta2[2,7] <= 0
 sBeta1Flag_2_7: - s[2,7] + beta1[2,7] <= 0
 sBeta2Flag_2_7: - s[2,7] + beta2[2,7] <= 0
 Eta1Overlap_2_7: eta1[2,7] + beta1[2,7] <= 1
 Eta2Overlap_2_7: eta2[2,7] + beta2[2,7] <= 1
 Eta1Flag_2_7: - xi_l[2] + xi_l[7] - 300 eta1[2,7] <= 0
 Eta2Flag_2_7: xi_r[2] - xi_r[7] - 300 eta2[2,7] <= 0
 7_Supports_2_InSameBin1_0: p_ij[2,0] - p_ij[7,0] - s[2,7] <= 0
 7_Supports_2_InSameBin2_0: - p_ij[2,0] + p_ij[7,0] - s[2,7] <= 0
 Gamma0ForItems_2_inBins_0_NoCut_7: p_ij[2,0] + gamma[2] <= 1
 7_Supports_2_InSameBin1_1: p_ij[2,1] - p_ij[7,1] - s[2,7] <= 0
 7_Supports_2_InSameBin2_1: - p_ij[2,1] + p_ij[7,1] - s[2,7] <= 0
 Gamma0ForItems_2_inBins_1_NoCut_7: p_ij[2,1] + gamma[2] <= 1
 7_Supports_2_InSameBin1_2: p_ij[2,2] - p_ij[7,2] - s[2,7] <= 0
 7_Supports_2_InSameBin2_2: - p_ij[2,2] + p_ij[7,2] - s[2,7] <= 0
 Gamma0ForItems_2_inBins_2_NoCut_7: p_ij[2,2] + gamma[2] <= 98
 7_Supports_2_InSameBin1_3: p_ij[2,3] - p_ij[7,3] - s[2,7] <= 0
 7_Supports_2_InSameBin2_3: - p_ij[2,3] + p_ij[7,3] - s[2,7] <= 0
 Gamma0ForItems_2_inBins_3_NoCut_7: p_ij[2,3] + gamma[2] <= 98
 2_OnCutIn2_7: 10000 p_ij[2,2] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 2_OnCutIn3_7: 10000 p_ij[2,3] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 AbsZ12_8: - zi[2] + zi_hi[8] - k[2,8] <= 0
 AbsZ22_8: zi[2] - zi_hi[8] - k[2,8] <= 0
 8_Below_2: zi[2] - zi_hi[8] + k[2,8] + 310 m[2,8] <= 310
 8_NotBelow_2: - zi[2] + zi_hi[8] + k[2,8] - 310 m[2,8] <= 0
 HeightCompat_2_8: h[2,8] - k[2,8] <= 0
 NoHeightCompat_2_8: - 155 h[2,8] + k[2,8] <= 0
 Overlap_2_8: - x_p[2,8] - x_p[8,2] + o[2,8] <= 0
 NoOverlap_2_8: x_p[2,8] + x_p[8,2] - 2 o[2,8] <= 0
 3_Supports_2_8: - h[2,8] - o[2,8] - s[2,8] <= -1
 3_DoesntSupport_2_8: h[2,8] + o[2,8] + 2 s[2,8] <= 2
 beta1_2_8: - s[2,8] + beta1[2,8] <= 0
 beta2_2_8: - s[2,8] + beta2[2,8] <= 0
 sBeta1Flag_2_8: - s[2,8] + beta1[2,8] <= 0
 sBeta2Flag_2_8: - s[2,8] + beta2[2,8] <= 0
 Eta1Overlap_2_8: eta1[2,8] + beta1[2,8] <= 1
 Eta2Overlap_2_8: eta2[2,8] + beta2[2,8] <= 1
 Eta1Flag_2_8: - xi_l[2] + xi_l[8] - 300 eta1[2,8] <= 0
 Eta2Flag_2_8: xi_r[2] - xi_r[8] - 300 eta2[2,8] <= 0
 8_Supports_2_InSameBin1_0: p_ij[2,0] - p_ij[8,0] - s[2,8] <= 0
 8_Supports_2_InSameBin2_0: - p_ij[2,0] + p_ij[8,0] - s[2,8] <= 0
 Gamma0ForItems_2_inBins_0_NoCut_8: p_ij[2,0] + gamma[2] <= 1
 8_Supports_2_InSameBin1_1: p_ij[2,1] - p_ij[8,1] - s[2,8] <= 0
 8_Supports_2_InSameBin2_1: - p_ij[2,1] + p_ij[8,1] - s[2,8] <= 0
 Gamma0ForItems_2_inBins_1_NoCut_8: p_ij[2,1] + gamma[2] <= 1
 8_Supports_2_InSameBin1_2: p_ij[2,2] - p_ij[8,2] - s[2,8] <= 0
 8_Supports_2_InSameBin2_2: - p_ij[2,2] + p_ij[8,2] - s[2,8] <= 0
 Gamma0ForItems_2_inBins_2_NoCut_8: p_ij[2,2] + gamma[2] <= 98
 8_Supports_2_InSameBin1_3: p_ij[2,3] - p_ij[8,3] - s[2,8] <= 0
 8_Supports_2_InSameBin2_3: - p_ij[2,3] + p_ij[8,3] - s[2,8] <= 0
 Gamma0ForItems_2_inBins_3_NoCut_8: p_ij[2,3] + gamma[2] <= 98
 2_OnCutIn2_8: 10000 p_ij[2,2] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 2_OnCutIn3_8: 10000 p_ij[2,3] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 AbsZ12_9: - zi[2] + zi_hi[9] - k[2,9] <= 0
 AbsZ22_9: zi[2] - zi_hi[9] - k[2,9] <= 0
 9_Below_2: zi[2] - zi_hi[9] + k[2,9] + 310 m[2,9] <= 310
 9_NotBelow_2: - zi[2] + zi_hi[9] + k[2,9] - 310 m[2,9] <= 0
 HeightCompat_2_9: h[2,9] - k[2,9] <= 0
 NoHeightCompat_2_9: - 155 h[2,9] + k[2,9] <= 0
 Overlap_2_9: - x_p[2,9] - x_p[9,2] + o[2,9] <= 0
 NoOverlap_2_9: x_p[2,9] + x_p[9,2] - 2 o[2,9] <= 0
 3_Supports_2_9: - h[2,9] - o[2,9] - s[2,9] <= -1
 3_DoesntSupport_2_9: h[2,9] + o[2,9] + 2 s[2,9] <= 2
 beta1_2_9: - s[2,9] + beta1[2,9] <= 0
 beta2_2_9: - s[2,9] + beta2[2,9] <= 0
 sBeta1Flag_2_9: - s[2,9] + beta1[2,9] <= 0
 sBeta2Flag_2_9: - s[2,9] + beta2[2,9] <= 0
 Eta1Overlap_2_9: eta1[2,9] + beta1[2,9] <= 1
 Eta2Overlap_2_9: eta2[2,9] + beta2[2,9] <= 1
 Eta1Flag_2_9: - xi_l[2] + xi_l[9] - 300 eta1[2,9] <= 0
 Eta2Flag_2_9: xi_r[2] - xi_r[9] - 300 eta2[2,9] <= 0
 9_Supports_2_InSameBin1_0: p_ij[2,0] - p_ij[9,0] - s[2,9] <= 0
 9_Supports_2_InSameBin2_0: - p_ij[2,0] + p_ij[9,0] - s[2,9] <= 0
 Gamma0ForItems_2_inBins_0_NoCut_9: p_ij[2,0] + gamma[2] <= 1
 9_Supports_2_InSameBin1_1: p_ij[2,1] - p_ij[9,1] - s[2,9] <= 0
 9_Supports_2_InSameBin2_1: - p_ij[2,1] + p_ij[9,1] - s[2,9] <= 0
 Gamma0ForItems_2_inBins_1_NoCut_9: p_ij[2,1] + gamma[2] <= 1
 9_Supports_2_InSameBin1_2: p_ij[2,2] - p_ij[9,2] - s[2,9] <= 0
 9_Supports_2_InSameBin2_2: - p_ij[2,2] + p_ij[9,2] - s[2,9] <= 0
 Gamma0ForItems_2_inBins_2_NoCut_9: p_ij[2,2] + gamma[2] <= 98
 9_Supports_2_InSameBin1_3: p_ij[2,3] - p_ij[9,3] - s[2,9] <= 0
 9_Supports_2_InSameBin2_3: - p_ij[2,3] + p_ij[9,3] - s[2,9] <= 0
 Gamma0ForItems_2_inBins_3_NoCut_9: p_ij[2,3] + gamma[2] <= 98
 2_OnCutIn2_9: 10000 p_ij[2,2] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 2_OnCutIn3_9: 10000 p_ij[2,3] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 AbsZ12_10: - zi[2] + zi_hi[10] - k[2,10] <= 0
 AbsZ22_10: zi[2] - zi_hi[10] - k[2,10] <= 0
 10_Below_2: zi[2] - zi_hi[10] + k[2,10] + 310 m[2,10] <= 310
 10_NotBelow_2: - zi[2] + zi_hi[10] + k[2,10] - 310 m[2,10] <= 0
 HeightCompat_2_10: h[2,10] - k[2,10] <= 0
 NoHeightCompat_2_10: - 155 h[2,10] + k[2,10] <= 0
 Overlap_2_10: - x_p[2,10] - x_p[10,2] + o[2,10] <= 0
 NoOverlap_2_10: x_p[2,10] + x_p[10,2] - 2 o[2,10] <= 0
 3_Supports_2_10: - h[2,10] - o[2,10] - s[2,10] <= -1
 3_DoesntSupport_2_10: h[2,10] + o[2,10] + 2 s[2,10] <= 2
 beta1_2_10: - s[2,10] + beta1[2,10] <= 0
 beta2_2_10: - s[2,10] + beta2[2,10] <= 0
 sBeta1Flag_2_10: - s[2,10] + beta1[2,10] <= 0
 sBeta2Flag_2_10: - s[2,10] + beta2[2,10] <= 0
 Eta1Overlap_2_10: eta1[2,10] + beta1[2,10] <= 1
 Eta2Overlap_2_10: eta2[2,10] + beta2[2,10] <= 1
 Eta1Flag_2_10: - xi_l[2] + xi_l[10] - 300 eta1[2,10] <= 0
 Eta2Flag_2_10: xi_r[2] - xi_r[10] - 300 eta2[2,10] <= 0
 10_Supports_2_InSameBin1_0: p_ij[2,0] - p_ij[10,0] - s[2,10] <= 0
 10_Supports_2_InSameBin2_0: - p_ij[2,0] + p_ij[10,0] - s[2,10] <= 0
 Gamma0ForItems_2_inBins_0_NoCut_10: p_ij[2,0] + gamma[2] <= 1
 10_Supports_2_InSameBin1_1: p_ij[2,1] - p_ij[10,1] - s[2,10] <= 0
 10_Supports_2_InSameBin2_1: - p_ij[2,1] + p_ij[10,1] - s[2,10] <= 0
 Gamma0ForItems_2_inBins_1_NoCut_10: p_ij[2,1] + gamma[2] <= 1
 10_Supports_2_InSameBin1_2: p_ij[2,2] - p_ij[10,2] - s[2,10] <= 0
 10_Supports_2_InSameBin2_2: - p_ij[2,2] + p_ij[10,2] - s[2,10] <= 0
 Gamma0ForItems_2_inBins_2_NoCut_10: p_ij[2,2] + gamma[2] <= 98
 10_Supports_2_InSameBin1_3: p_ij[2,3] - p_ij[10,3] - s[2,10] <= 0
 10_Supports_2_InSameBin2_3: - p_ij[2,3] + p_ij[10,3] - s[2,10] <= 0
 Gamma0ForItems_2_inBins_3_NoCut_10: p_ij[2,3] + gamma[2] <= 98
 2_OnCutIn2_10: 10000 p_ij[2,2] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 2_OnCutIn3_10: 10000 p_ij[2,3] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 AbsZ12_11: - zi[2] + zi_hi[11] - k[2,11] <= 0
 AbsZ22_11: zi[2] - zi_hi[11] - k[2,11] <= 0
 11_Below_2: zi[2] - zi_hi[11] + k[2,11] + 310 m[2,11] <= 310
 11_NotBelow_2: - zi[2] + zi_hi[11] + k[2,11] - 310 m[2,11] <= 0
 HeightCompat_2_11: h[2,11] - k[2,11] <= 0
 NoHeightCompat_2_11: - 155 h[2,11] + k[2,11] <= 0
 Overlap_2_11: - x_p[2,11] - x_p[11,2] + o[2,11] <= 0
 NoOverlap_2_11: x_p[2,11] + x_p[11,2] - 2 o[2,11] <= 0
 3_Supports_2_11: - h[2,11] - o[2,11] - s[2,11] <= -1
 3_DoesntSupport_2_11: h[2,11] + o[2,11] + 2 s[2,11] <= 2
 beta1_2_11: - s[2,11] + beta1[2,11] <= 0
 beta2_2_11: - s[2,11] + beta2[2,11] <= 0
 sBeta1Flag_2_11: - s[2,11] + beta1[2,11] <= 0
 sBeta2Flag_2_11: - s[2,11] + beta2[2,11] <= 0
 Eta1Overlap_2_11: eta1[2,11] + beta1[2,11] <= 1
 Eta2Overlap_2_11: eta2[2,11] + beta2[2,11] <= 1
 Eta1Flag_2_11: - xi_l[2] + xi_l[11] - 300 eta1[2,11] <= 0
 Eta2Flag_2_11: xi_r[2] - xi_r[11] - 300 eta2[2,11] <= 0
 11_Supports_2_InSameBin1_0: p_ij[2,0] - p_ij[11,0] - s[2,11] <= 0
 11_Supports_2_InSameBin2_0: - p_ij[2,0] + p_ij[11,0] - s[2,11] <= 0
 Gamma0ForItems_2_inBins_0_NoCut_11: p_ij[2,0] + gamma[2] <= 1
 11_Supports_2_InSameBin1_1: p_ij[2,1] - p_ij[11,1] - s[2,11] <= 0
 11_Supports_2_InSameBin2_1: - p_ij[2,1] + p_ij[11,1] - s[2,11] <= 0
 Gamma0ForItems_2_inBins_1_NoCut_11: p_ij[2,1] + gamma[2] <= 1
 11_Supports_2_InSameBin1_2: p_ij[2,2] - p_ij[11,2] - s[2,11] <= 0
 11_Supports_2_InSameBin2_2: - p_ij[2,2] + p_ij[11,2] - s[2,11] <= 0
 Gamma0ForItems_2_inBins_2_NoCut_11: p_ij[2,2] + gamma[2] <= 98
 11_Supports_2_InSameBin1_3: p_ij[2,3] - p_ij[11,3] - s[2,11] <= 0
 11_Supports_2_InSameBin2_3: - p_ij[2,3] + p_ij[11,3] - s[2,11] <= 0
 Gamma0ForItems_2_inBins_3_NoCut_11: p_ij[2,3] + gamma[2] <= 98
 2_OnCutIn2_11: 10000 p_ij[2,2] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 2_OnCutIn3_11: 10000 p_ij[2,3] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 AbsZ12_12: - zi[2] + zi_hi[12] - k[2,12] <= 0
 AbsZ22_12: zi[2] - zi_hi[12] - k[2,12] <= 0
 12_Below_2: zi[2] - zi_hi[12] + k[2,12] + 310 m[2,12] <= 310
 12_NotBelow_2: - zi[2] + zi_hi[12] + k[2,12] - 310 m[2,12] <= 0
 HeightCompat_2_12: h[2,12] - k[2,12] <= 0
 NoHeightCompat_2_12: - 155 h[2,12] + k[2,12] <= 0
 Overlap_2_12: - x_p[2,12] - x_p[12,2] + o[2,12] <= 0
 NoOverlap_2_12: x_p[2,12] + x_p[12,2] - 2 o[2,12] <= 0
 3_Supports_2_12: - h[2,12] - o[2,12] - s[2,12] <= -1
 3_DoesntSupport_2_12: h[2,12] + o[2,12] + 2 s[2,12] <= 2
 beta1_2_12: - s[2,12] + beta1[2,12] <= 0
 beta2_2_12: - s[2,12] + beta2[2,12] <= 0
 sBeta1Flag_2_12: - s[2,12] + beta1[2,12] <= 0
 sBeta2Flag_2_12: - s[2,12] + beta2[2,12] <= 0
 Eta1Overlap_2_12: eta1[2,12] + beta1[2,12] <= 1
 Eta2Overlap_2_12: eta2[2,12] + beta2[2,12] <= 1
 Eta1Flag_2_12: - xi_l[2] + xi_l[12] - 300 eta1[2,12] <= 0
 Eta2Flag_2_12: xi_r[2] - xi_r[12] - 300 eta2[2,12] <= 0
 12_Supports_2_InSameBin1_0: p_ij[2,0] - p_ij[12,0] - s[2,12] <= 0
 12_Supports_2_InSameBin2_0: - p_ij[2,0] + p_ij[12,0] - s[2,12] <= 0
 Gamma0ForItems_2_inBins_0_NoCut_12: p_ij[2,0] + gamma[2] <= 1
 12_Supports_2_InSameBin1_1: p_ij[2,1] - p_ij[12,1] - s[2,12] <= 0
 12_Supports_2_InSameBin2_1: - p_ij[2,1] + p_ij[12,1] - s[2,12] <= 0
 Gamma0ForItems_2_inBins_1_NoCut_12: p_ij[2,1] + gamma[2] <= 1
 12_Supports_2_InSameBin1_2: p_ij[2,2] - p_ij[12,2] - s[2,12] <= 0
 12_Supports_2_InSameBin2_2: - p_ij[2,2] + p_ij[12,2] - s[2,12] <= 0
 Gamma0ForItems_2_inBins_2_NoCut_12: p_ij[2,2] + gamma[2] <= 98
 12_Supports_2_InSameBin1_3: p_ij[2,3] - p_ij[12,3] - s[2,12] <= 0
 12_Supports_2_InSameBin2_3: - p_ij[2,3] + p_ij[12,3] - s[2,12] <= 0
 Gamma0ForItems_2_inBins_3_NoCut_12: p_ij[2,3] + gamma[2] <= 98
 2_OnCutIn2_12: 10000 p_ij[2,2] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 2_OnCutIn3_12: 10000 p_ij[2,3] + 1.261904761904762 xi_l[2] + zi[2]
   + 10000 gamma[2] <= 20053
 GroundConstraint_3: zi[3] + 155 g[3] <= 155
 AbsZ13_0: - zi[3] + zi_hi[0] - k[3,0] <= 0
 AbsZ23_0: zi[3] - zi_hi[0] - k[3,0] <= 0
 0_Below_3: zi[3] - zi_hi[0] + k[3,0] + 310 m[3,0] <= 310
 0_NotBelow_3: - zi[3] + zi_hi[0] + k[3,0] - 310 m[3,0] <= 0
 HeightCompat_3_0: h[3,0] - k[3,0] <= 0
 NoHeightCompat_3_0: - 155 h[3,0] + k[3,0] <= 0
 Overlap_3_0: - x_p[0,3] - x_p[3,0] + o[3,0] <= 0
 NoOverlap_3_0: x_p[0,3] + x_p[3,0] - 2 o[3,0] <= 0
 3_Supports_3_0: - h[3,0] - o[3,0] - s[3,0] <= -1
 3_DoesntSupport_3_0: h[3,0] + o[3,0] + 2 s[3,0] <= 2
 beta1_3_0: - s[3,0] + beta1[3,0] <= 0
 beta2_3_0: - s[3,0] + beta2[3,0] <= 0
 sBeta1Flag_3_0: - s[3,0] + beta1[3,0] <= 0
 sBeta2Flag_3_0: - s[3,0] + beta2[3,0] <= 0
 Eta1Overlap_3_0: eta1[3,0] + beta1[3,0] <= 1
 Eta2Overlap_3_0: eta2[3,0] + beta2[3,0] <= 1
 Eta1Flag_3_0: xi_l[0] - xi_l[3] - 300 eta1[3,0] <= 0
 Eta2Flag_3_0: - xi_r[0] + xi_r[3] - 300 eta2[3,0] <= 0
 0_Supports_3_InSameBin1_0: - p_ij[0,0] + p_ij[3,0] - s[3,0] <= 0
 0_Supports_3_InSameBin2_0: p_ij[0,0] - p_ij[3,0] - s[3,0] <= 0
 Gamma0ForItems_3_inBins_0_NoCut_0: p_ij[3,0] + gamma[3] <= 1
 0_Supports_3_InSameBin1_1: - p_ij[0,1] + p_ij[3,1] - s[3,0] <= 0
 0_Supports_3_InSameBin2_1: p_ij[0,1] - p_ij[3,1] - s[3,0] <= 0
 Gamma0ForItems_3_inBins_1_NoCut_0: p_ij[3,1] + gamma[3] <= 1
 0_Supports_3_InSameBin1_2: - p_ij[0,2] + p_ij[3,2] - s[3,0] <= 0
 0_Supports_3_InSameBin2_2: p_ij[0,2] - p_ij[3,2] - s[3,0] <= 0
 Gamma0ForItems_3_inBins_2_NoCut_0: p_ij[3,2] + gamma[3] <= 98
 0_Supports_3_InSameBin1_3: - p_ij[0,3] + p_ij[3,3] - s[3,0] <= 0
 0_Supports_3_InSameBin2_3: p_ij[0,3] - p_ij[3,3] - s[3,0] <= 0
 Gamma0ForItems_3_inBins_3_NoCut_0: p_ij[3,3] + gamma[3] <= 98
 3_OnCutIn2_0: 10000 p_ij[3,2] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 3_OnCutIn3_0: 10000 p_ij[3,3] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 AbsZ13_1: - zi[3] + zi_hi[1] - k[3,1] <= 0
 AbsZ23_1: zi[3] - zi_hi[1] - k[3,1] <= 0
 1_Below_3: zi[3] - zi_hi[1] + k[3,1] + 310 m[3,1] <= 310
 1_NotBelow_3: - zi[3] + zi_hi[1] + k[3,1] - 310 m[3,1] <= 0
 HeightCompat_3_1: h[3,1] - k[3,1] <= 0
 NoHeightCompat_3_1: - 155 h[3,1] + k[3,1] <= 0
 Overlap_3_1: - x_p[1,3] - x_p[3,1] + o[3,1] <= 0
 NoOverlap_3_1: x_p[1,3] + x_p[3,1] - 2 o[3,1] <= 0
 3_Supports_3_1: - h[3,1] - o[3,1] - s[3,1] <= -1
 3_DoesntSupport_3_1: h[3,1] + o[3,1] + 2 s[3,1] <= 2
 beta1_3_1: - s[3,1] + beta1[3,1] <= 0
 beta2_3_1: - s[3,1] + beta2[3,1] <= 0
 sBeta1Flag_3_1: - s[3,1] + beta1[3,1] <= 0
 sBeta2Flag_3_1: - s[3,1] + beta2[3,1] <= 0
 Eta1Overlap_3_1: eta1[3,1] + beta1[3,1] <= 1
 Eta2Overlap_3_1: eta2[3,1] + beta2[3,1] <= 1
 Eta1Flag_3_1: xi_l[1] - xi_l[3] - 300 eta1[3,1] <= 0
 Eta2Flag_3_1: - xi_r[1] + xi_r[3] - 300 eta2[3,1] <= 0
 1_Supports_3_InSameBin1_0: - p_ij[1,0] + p_ij[3,0] - s[3,1] <= 0
 1_Supports_3_InSameBin2_0: p_ij[1,0] - p_ij[3,0] - s[3,1] <= 0
 Gamma0ForItems_3_inBins_0_NoCut_1: p_ij[3,0] + gamma[3] <= 1
 1_Supports_3_InSameBin1_1: - p_ij[1,1] + p_ij[3,1] - s[3,1] <= 0
 1_Supports_3_InSameBin2_1: p_ij[1,1] - p_ij[3,1] - s[3,1] <= 0
 Gamma0ForItems_3_inBins_1_NoCut_1: p_ij[3,1] + gamma[3] <= 1
 1_Supports_3_InSameBin1_2: - p_ij[1,2] + p_ij[3,2] - s[3,1] <= 0
 1_Supports_3_InSameBin2_2: p_ij[1,2] - p_ij[3,2] - s[3,1] <= 0
 Gamma0ForItems_3_inBins_2_NoCut_1: p_ij[3,2] + gamma[3] <= 98
 1_Supports_3_InSameBin1_3: - p_ij[1,3] + p_ij[3,3] - s[3,1] <= 0
 1_Supports_3_InSameBin2_3: p_ij[1,3] - p_ij[3,3] - s[3,1] <= 0
 Gamma0ForItems_3_inBins_3_NoCut_1: p_ij[3,3] + gamma[3] <= 98
 3_OnCutIn2_1: 10000 p_ij[3,2] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 3_OnCutIn3_1: 10000 p_ij[3,3] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 AbsZ13_2: - zi[3] + zi_hi[2] - k[3,2] <= 0
 AbsZ23_2: zi[3] - zi_hi[2] - k[3,2] <= 0
 2_Below_3: zi[3] - zi_hi[2] + k[3,2] + 310 m[3,2] <= 310
 2_NotBelow_3: - zi[3] + zi_hi[2] + k[3,2] - 310 m[3,2] <= 0
 HeightCompat_3_2: h[3,2] - k[3,2] <= 0
 NoHeightCompat_3_2: - 155 h[3,2] + k[3,2] <= 0
 Overlap_3_2: - x_p[2,3] - x_p[3,2] + o[3,2] <= 0
 NoOverlap_3_2: x_p[2,3] + x_p[3,2] - 2 o[3,2] <= 0
 3_Supports_3_2: - h[3,2] - o[3,2] - s[3,2] <= -1
 3_DoesntSupport_3_2: h[3,2] + o[3,2] + 2 s[3,2] <= 2
 beta1_3_2: - s[3,2] + beta1[3,2] <= 0
 beta2_3_2: - s[3,2] + beta2[3,2] <= 0
 sBeta1Flag_3_2: - s[3,2] + beta1[3,2] <= 0
 sBeta2Flag_3_2: - s[3,2] + beta2[3,2] <= 0
 Eta1Overlap_3_2: eta1[3,2] + beta1[3,2] <= 1
 Eta2Overlap_3_2: eta2[3,2] + beta2[3,2] <= 1
 Eta1Flag_3_2: xi_l[2] - xi_l[3] - 300 eta1[3,2] <= 0
 Eta2Flag_3_2: - xi_r[2] + xi_r[3] - 300 eta2[3,2] <= 0
 2_Supports_3_InSameBin1_0: - p_ij[2,0] + p_ij[3,0] - s[3,2] <= 0
 2_Supports_3_InSameBin2_0: p_ij[2,0] - p_ij[3,0] - s[3,2] <= 0
 Gamma0ForItems_3_inBins_0_NoCut_2: p_ij[3,0] + gamma[3] <= 1
 2_Supports_3_InSameBin1_1: - p_ij[2,1] + p_ij[3,1] - s[3,2] <= 0
 2_Supports_3_InSameBin2_1: p_ij[2,1] - p_ij[3,1] - s[3,2] <= 0
 Gamma0ForItems_3_inBins_1_NoCut_2: p_ij[3,1] + gamma[3] <= 1
 2_Supports_3_InSameBin1_2: - p_ij[2,2] + p_ij[3,2] - s[3,2] <= 0
 2_Supports_3_InSameBin2_2: p_ij[2,2] - p_ij[3,2] - s[3,2] <= 0
 Gamma0ForItems_3_inBins_2_NoCut_2: p_ij[3,2] + gamma[3] <= 98
 2_Supports_3_InSameBin1_3: - p_ij[2,3] + p_ij[3,3] - s[3,2] <= 0
 2_Supports_3_InSameBin2_3: p_ij[2,3] - p_ij[3,3] - s[3,2] <= 0
 Gamma0ForItems_3_inBins_3_NoCut_2: p_ij[3,3] + gamma[3] <= 98
 3_OnCutIn2_2: 10000 p_ij[3,2] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 3_OnCutIn3_2: 10000 p_ij[3,3] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 3_OnCutIn2_3: 10000 p_ij[3,2] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 3_OnCutIn3_3: 10000 p_ij[3,3] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 AbsZ13_4: - zi[3] + zi_hi[4] - k[3,4] <= 0
 AbsZ23_4: zi[3] - zi_hi[4] - k[3,4] <= 0
 4_Below_3: zi[3] - zi_hi[4] + k[3,4] + 310 m[3,4] <= 310
 4_NotBelow_3: - zi[3] + zi_hi[4] + k[3,4] - 310 m[3,4] <= 0
 HeightCompat_3_4: h[3,4] - k[3,4] <= 0
 NoHeightCompat_3_4: - 155 h[3,4] + k[3,4] <= 0
 Overlap_3_4: - x_p[3,4] - x_p[4,3] + o[3,4] <= 0
 NoOverlap_3_4: x_p[3,4] + x_p[4,3] - 2 o[3,4] <= 0
 3_Supports_3_4: - h[3,4] - o[3,4] - s[3,4] <= -1
 3_DoesntSupport_3_4: h[3,4] + o[3,4] + 2 s[3,4] <= 2
 beta1_3_4: - s[3,4] + beta1[3,4] <= 0
 beta2_3_4: - s[3,4] + beta2[3,4] <= 0
 sBeta1Flag_3_4: - s[3,4] + beta1[3,4] <= 0
 sBeta2Flag_3_4: - s[3,4] + beta2[3,4] <= 0
 Eta1Overlap_3_4: eta1[3,4] + beta1[3,4] <= 1
 Eta2Overlap_3_4: eta2[3,4] + beta2[3,4] <= 1
 Eta1Flag_3_4: - xi_l[3] + xi_l[4] - 300 eta1[3,4] <= 0
 Eta2Flag_3_4: xi_r[3] - xi_r[4] - 300 eta2[3,4] <= 0
 4_Supports_3_InSameBin1_0: p_ij[3,0] - p_ij[4,0] - s[3,4] <= 0
 4_Supports_3_InSameBin2_0: - p_ij[3,0] + p_ij[4,0] - s[3,4] <= 0
 Gamma0ForItems_3_inBins_0_NoCut_4: p_ij[3,0] + gamma[3] <= 1
 4_Supports_3_InSameBin1_1: p_ij[3,1] - p_ij[4,1] - s[3,4] <= 0
 4_Supports_3_InSameBin2_1: - p_ij[3,1] + p_ij[4,1] - s[3,4] <= 0
 Gamma0ForItems_3_inBins_1_NoCut_4: p_ij[3,1] + gamma[3] <= 1
 4_Supports_3_InSameBin1_2: p_ij[3,2] - p_ij[4,2] - s[3,4] <= 0
 4_Supports_3_InSameBin2_2: - p_ij[3,2] + p_ij[4,2] - s[3,4] <= 0
 Gamma0ForItems_3_inBins_2_NoCut_4: p_ij[3,2] + gamma[3] <= 98
 4_Supports_3_InSameBin1_3: p_ij[3,3] - p_ij[4,3] - s[3,4] <= 0
 4_Supports_3_InSameBin2_3: - p_ij[3,3] + p_ij[4,3] - s[3,4] <= 0
 Gamma0ForItems_3_inBins_3_NoCut_4: p_ij[3,3] + gamma[3] <= 98
 3_OnCutIn2_4: 10000 p_ij[3,2] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 3_OnCutIn3_4: 10000 p_ij[3,3] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 AbsZ13_5: - zi[3] + zi_hi[5] - k[3,5] <= 0
 AbsZ23_5: zi[3] - zi_hi[5] - k[3,5] <= 0
 5_Below_3: zi[3] - zi_hi[5] + k[3,5] + 310 m[3,5] <= 310
 5_NotBelow_3: - zi[3] + zi_hi[5] + k[3,5] - 310 m[3,5] <= 0
 HeightCompat_3_5: h[3,5] - k[3,5] <= 0
 NoHeightCompat_3_5: - 155 h[3,5] + k[3,5] <= 0
 Overlap_3_5: - x_p[3,5] - x_p[5,3] + o[3,5] <= 0
 NoOverlap_3_5: x_p[3,5] + x_p[5,3] - 2 o[3,5] <= 0
 3_Supports_3_5: - h[3,5] - o[3,5] - s[3,5] <= -1
 3_DoesntSupport_3_5: h[3,5] + o[3,5] + 2 s[3,5] <= 2
 beta1_3_5: - s[3,5] + beta1[3,5] <= 0
 beta2_3_5: - s[3,5] + beta2[3,5] <= 0
 sBeta1Flag_3_5: - s[3,5] + beta1[3,5] <= 0
 sBeta2Flag_3_5: - s[3,5] + beta2[3,5] <= 0
 Eta1Overlap_3_5: eta1[3,5] + beta1[3,5] <= 1
 Eta2Overlap_3_5: eta2[3,5] + beta2[3,5] <= 1
 Eta1Flag_3_5: - xi_l[3] + xi_l[5] - 300 eta1[3,5] <= 0
 Eta2Flag_3_5: xi_r[3] - xi_r[5] - 300 eta2[3,5] <= 0
 5_Supports_3_InSameBin1_0: p_ij[3,0] - p_ij[5,0] - s[3,5] <= 0
 5_Supports_3_InSameBin2_0: - p_ij[3,0] + p_ij[5,0] - s[3,5] <= 0
 Gamma0ForItems_3_inBins_0_NoCut_5: p_ij[3,0] + gamma[3] <= 1
 5_Supports_3_InSameBin1_1: p_ij[3,1] - p_ij[5,1] - s[3,5] <= 0
 5_Supports_3_InSameBin2_1: - p_ij[3,1] + p_ij[5,1] - s[3,5] <= 0
 Gamma0ForItems_3_inBins_1_NoCut_5: p_ij[3,1] + gamma[3] <= 1
 5_Supports_3_InSameBin1_2: p_ij[3,2] - p_ij[5,2] - s[3,5] <= 0
 5_Supports_3_InSameBin2_2: - p_ij[3,2] + p_ij[5,2] - s[3,5] <= 0
 Gamma0ForItems_3_inBins_2_NoCut_5: p_ij[3,2] + gamma[3] <= 98
 5_Supports_3_InSameBin1_3: p_ij[3,3] - p_ij[5,3] - s[3,5] <= 0
 5_Supports_3_InSameBin2_3: - p_ij[3,3] + p_ij[5,3] - s[3,5] <= 0
 Gamma0ForItems_3_inBins_3_NoCut_5: p_ij[3,3] + gamma[3] <= 98
 3_OnCutIn2_5: 10000 p_ij[3,2] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 3_OnCutIn3_5: 10000 p_ij[3,3] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 AbsZ13_6: - zi[3] + zi_hi[6] - k[3,6] <= 0
 AbsZ23_6: zi[3] - zi_hi[6] - k[3,6] <= 0
 6_Below_3: zi[3] - zi_hi[6] + k[3,6] + 310 m[3,6] <= 310
 6_NotBelow_3: - zi[3] + zi_hi[6] + k[3,6] - 310 m[3,6] <= 0
 HeightCompat_3_6: h[3,6] - k[3,6] <= 0
 NoHeightCompat_3_6: - 155 h[3,6] + k[3,6] <= 0
 Overlap_3_6: - x_p[3,6] - x_p[6,3] + o[3,6] <= 0
 NoOverlap_3_6: x_p[3,6] + x_p[6,3] - 2 o[3,6] <= 0
 3_Supports_3_6: - h[3,6] - o[3,6] - s[3,6] <= -1
 3_DoesntSupport_3_6: h[3,6] + o[3,6] + 2 s[3,6] <= 2
 beta1_3_6: - s[3,6] + beta1[3,6] <= 0
 beta2_3_6: - s[3,6] + beta2[3,6] <= 0
 sBeta1Flag_3_6: - s[3,6] + beta1[3,6] <= 0
 sBeta2Flag_3_6: - s[3,6] + beta2[3,6] <= 0
 Eta1Overlap_3_6: eta1[3,6] + beta1[3,6] <= 1
 Eta2Overlap_3_6: eta2[3,6] + beta2[3,6] <= 1
 Eta1Flag_3_6: - xi_l[3] + xi_l[6] - 300 eta1[3,6] <= 0
 Eta2Flag_3_6: xi_r[3] - xi_r[6] - 300 eta2[3,6] <= 0
 6_Supports_3_InSameBin1_0: p_ij[3,0] - p_ij[6,0] - s[3,6] <= 0
 6_Supports_3_InSameBin2_0: - p_ij[3,0] + p_ij[6,0] - s[3,6] <= 0
 Gamma0ForItems_3_inBins_0_NoCut_6: p_ij[3,0] + gamma[3] <= 1
 6_Supports_3_InSameBin1_1: p_ij[3,1] - p_ij[6,1] - s[3,6] <= 0
 6_Supports_3_InSameBin2_1: - p_ij[3,1] + p_ij[6,1] - s[3,6] <= 0
 Gamma0ForItems_3_inBins_1_NoCut_6: p_ij[3,1] + gamma[3] <= 1
 6_Supports_3_InSameBin1_2: p_ij[3,2] - p_ij[6,2] - s[3,6] <= 0
 6_Supports_3_InSameBin2_2: - p_ij[3,2] + p_ij[6,2] - s[3,6] <= 0
 Gamma0ForItems_3_inBins_2_NoCut_6: p_ij[3,2] + gamma[3] <= 98
 6_Supports_3_InSameBin1_3: p_ij[3,3] - p_ij[6,3] - s[3,6] <= 0
 6_Supports_3_InSameBin2_3: - p_ij[3,3] + p_ij[6,3] - s[3,6] <= 0
 Gamma0ForItems_3_inBins_3_NoCut_6: p_ij[3,3] + gamma[3] <= 98
 3_OnCutIn2_6: 10000 p_ij[3,2] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 3_OnCutIn3_6: 10000 p_ij[3,3] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 AbsZ13_7: - zi[3] + zi_hi[7] - k[3,7] <= 0
 AbsZ23_7: zi[3] - zi_hi[7] - k[3,7] <= 0
 7_Below_3: zi[3] - zi_hi[7] + k[3,7] + 310 m[3,7] <= 310
 7_NotBelow_3: - zi[3] + zi_hi[7] + k[3,7] - 310 m[3,7] <= 0
 HeightCompat_3_7: h[3,7] - k[3,7] <= 0
 NoHeightCompat_3_7: - 155 h[3,7] + k[3,7] <= 0
 Overlap_3_7: - x_p[3,7] - x_p[7,3] + o[3,7] <= 0
 NoOverlap_3_7: x_p[3,7] + x_p[7,3] - 2 o[3,7] <= 0
 3_Supports_3_7: - h[3,7] - o[3,7] - s[3,7] <= -1
 3_DoesntSupport_3_7: h[3,7] + o[3,7] + 2 s[3,7] <= 2
 beta1_3_7: - s[3,7] + beta1[3,7] <= 0
 beta2_3_7: - s[3,7] + beta2[3,7] <= 0
 sBeta1Flag_3_7: - s[3,7] + beta1[3,7] <= 0
 sBeta2Flag_3_7: - s[3,7] + beta2[3,7] <= 0
 Eta1Overlap_3_7: eta1[3,7] + beta1[3,7] <= 1
 Eta2Overlap_3_7: eta2[3,7] + beta2[3,7] <= 1
 Eta1Flag_3_7: - xi_l[3] + xi_l[7] - 300 eta1[3,7] <= 0
 Eta2Flag_3_7: xi_r[3] - xi_r[7] - 300 eta2[3,7] <= 0
 7_Supports_3_InSameBin1_0: p_ij[3,0] - p_ij[7,0] - s[3,7] <= 0
 7_Supports_3_InSameBin2_0: - p_ij[3,0] + p_ij[7,0] - s[3,7] <= 0
 Gamma0ForItems_3_inBins_0_NoCut_7: p_ij[3,0] + gamma[3] <= 1
 7_Supports_3_InSameBin1_1: p_ij[3,1] - p_ij[7,1] - s[3,7] <= 0
 7_Supports_3_InSameBin2_1: - p_ij[3,1] + p_ij[7,1] - s[3,7] <= 0
 Gamma0ForItems_3_inBins_1_NoCut_7: p_ij[3,1] + gamma[3] <= 1
 7_Supports_3_InSameBin1_2: p_ij[3,2] - p_ij[7,2] - s[3,7] <= 0
 7_Supports_3_InSameBin2_2: - p_ij[3,2] + p_ij[7,2] - s[3,7] <= 0
 Gamma0ForItems_3_inBins_2_NoCut_7: p_ij[3,2] + gamma[3] <= 98
 7_Supports_3_InSameBin1_3: p_ij[3,3] - p_ij[7,3] - s[3,7] <= 0
 7_Supports_3_InSameBin2_3: - p_ij[3,3] + p_ij[7,3] - s[3,7] <= 0
 Gamma0ForItems_3_inBins_3_NoCut_7: p_ij[3,3] + gamma[3] <= 98
 3_OnCutIn2_7: 10000 p_ij[3,2] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 3_OnCutIn3_7: 10000 p_ij[3,3] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 AbsZ13_8: - zi[3] + zi_hi[8] - k[3,8] <= 0
 AbsZ23_8: zi[3] - zi_hi[8] - k[3,8] <= 0
 8_Below_3: zi[3] - zi_hi[8] + k[3,8] + 310 m[3,8] <= 310
 8_NotBelow_3: - zi[3] + zi_hi[8] + k[3,8] - 310 m[3,8] <= 0
 HeightCompat_3_8: h[3,8] - k[3,8] <= 0
 NoHeightCompat_3_8: - 155 h[3,8] + k[3,8] <= 0
 Overlap_3_8: - x_p[3,8] - x_p[8,3] + o[3,8] <= 0
 NoOverlap_3_8: x_p[3,8] + x_p[8,3] - 2 o[3,8] <= 0
 3_Supports_3_8: - h[3,8] - o[3,8] - s[3,8] <= -1
 3_DoesntSupport_3_8: h[3,8] + o[3,8] + 2 s[3,8] <= 2
 beta1_3_8: - s[3,8] + beta1[3,8] <= 0
 beta2_3_8: - s[3,8] + beta2[3,8] <= 0
 sBeta1Flag_3_8: - s[3,8] + beta1[3,8] <= 0
 sBeta2Flag_3_8: - s[3,8] + beta2[3,8] <= 0
 Eta1Overlap_3_8: eta1[3,8] + beta1[3,8] <= 1
 Eta2Overlap_3_8: eta2[3,8] + beta2[3,8] <= 1
 Eta1Flag_3_8: - xi_l[3] + xi_l[8] - 300 eta1[3,8] <= 0
 Eta2Flag_3_8: xi_r[3] - xi_r[8] - 300 eta2[3,8] <= 0
 8_Supports_3_InSameBin1_0: p_ij[3,0] - p_ij[8,0] - s[3,8] <= 0
 8_Supports_3_InSameBin2_0: - p_ij[3,0] + p_ij[8,0] - s[3,8] <= 0
 Gamma0ForItems_3_inBins_0_NoCut_8: p_ij[3,0] + gamma[3] <= 1
 8_Supports_3_InSameBin1_1: p_ij[3,1] - p_ij[8,1] - s[3,8] <= 0
 8_Supports_3_InSameBin2_1: - p_ij[3,1] + p_ij[8,1] - s[3,8] <= 0
 Gamma0ForItems_3_inBins_1_NoCut_8: p_ij[3,1] + gamma[3] <= 1
 8_Supports_3_InSameBin1_2: p_ij[3,2] - p_ij[8,2] - s[3,8] <= 0
 8_Supports_3_InSameBin2_2: - p_ij[3,2] + p_ij[8,2] - s[3,8] <= 0
 Gamma0ForItems_3_inBins_2_NoCut_8: p_ij[3,2] + gamma[3] <= 98
 8_Supports_3_InSameBin1_3: p_ij[3,3] - p_ij[8,3] - s[3,8] <= 0
 8_Supports_3_InSameBin2_3: - p_ij[3,3] + p_ij[8,3] - s[3,8] <= 0
 Gamma0ForItems_3_inBins_3_NoCut_8: p_ij[3,3] + gamma[3] <= 98
 3_OnCutIn2_8: 10000 p_ij[3,2] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 3_OnCutIn3_8: 10000 p_ij[3,3] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 AbsZ13_9: - zi[3] + zi_hi[9] - k[3,9] <= 0
 AbsZ23_9: zi[3] - zi_hi[9] - k[3,9] <= 0
 9_Below_3: zi[3] - zi_hi[9] + k[3,9] + 310 m[3,9] <= 310
 9_NotBelow_3: - zi[3] + zi_hi[9] + k[3,9] - 310 m[3,9] <= 0
 HeightCompat_3_9: h[3,9] - k[3,9] <= 0
 NoHeightCompat_3_9: - 155 h[3,9] + k[3,9] <= 0
 Overlap_3_9: - x_p[3,9] - x_p[9,3] + o[3,9] <= 0
 NoOverlap_3_9: x_p[3,9] + x_p[9,3] - 2 o[3,9] <= 0
 3_Supports_3_9: - h[3,9] - o[3,9] - s[3,9] <= -1
 3_DoesntSupport_3_9: h[3,9] + o[3,9] + 2 s[3,9] <= 2
 beta1_3_9: - s[3,9] + beta1[3,9] <= 0
 beta2_3_9: - s[3,9] + beta2[3,9] <= 0
 sBeta1Flag_3_9: - s[3,9] + beta1[3,9] <= 0
 sBeta2Flag_3_9: - s[3,9] + beta2[3,9] <= 0
 Eta1Overlap_3_9: eta1[3,9] + beta1[3,9] <= 1
 Eta2Overlap_3_9: eta2[3,9] + beta2[3,9] <= 1
 Eta1Flag_3_9: - xi_l[3] + xi_l[9] - 300 eta1[3,9] <= 0
 Eta2Flag_3_9: xi_r[3] - xi_r[9] - 300 eta2[3,9] <= 0
 9_Supports_3_InSameBin1_0: p_ij[3,0] - p_ij[9,0] - s[3,9] <= 0
 9_Supports_3_InSameBin2_0: - p_ij[3,0] + p_ij[9,0] - s[3,9] <= 0
 Gamma0ForItems_3_inBins_0_NoCut_9: p_ij[3,0] + gamma[3] <= 1
 9_Supports_3_InSameBin1_1: p_ij[3,1] - p_ij[9,1] - s[3,9] <= 0
 9_Supports_3_InSameBin2_1: - p_ij[3,1] + p_ij[9,1] - s[3,9] <= 0
 Gamma0ForItems_3_inBins_1_NoCut_9: p_ij[3,1] + gamma[3] <= 1
 9_Supports_3_InSameBin1_2: p_ij[3,2] - p_ij[9,2] - s[3,9] <= 0
 9_Supports_3_InSameBin2_2: - p_ij[3,2] + p_ij[9,2] - s[3,9] <= 0
 Gamma0ForItems_3_inBins_2_NoCut_9: p_ij[3,2] + gamma[3] <= 98
 9_Supports_3_InSameBin1_3: p_ij[3,3] - p_ij[9,3] - s[3,9] <= 0
 9_Supports_3_InSameBin2_3: - p_ij[3,3] + p_ij[9,3] - s[3,9] <= 0
 Gamma0ForItems_3_inBins_3_NoCut_9: p_ij[3,3] + gamma[3] <= 98
 3_OnCutIn2_9: 10000 p_ij[3,2] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 3_OnCutIn3_9: 10000 p_ij[3,3] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 AbsZ13_10: - zi[3] + zi_hi[10] - k[3,10] <= 0
 AbsZ23_10: zi[3] - zi_hi[10] - k[3,10] <= 0
 10_Below_3: zi[3] - zi_hi[10] + k[3,10] + 310 m[3,10] <= 310
 10_NotBelow_3: - zi[3] + zi_hi[10] + k[3,10] - 310 m[3,10] <= 0
 HeightCompat_3_10: h[3,10] - k[3,10] <= 0
 NoHeightCompat_3_10: - 155 h[3,10] + k[3,10] <= 0
 Overlap_3_10: - x_p[3,10] - x_p[10,3] + o[3,10] <= 0
 NoOverlap_3_10: x_p[3,10] + x_p[10,3] - 2 o[3,10] <= 0
 3_Supports_3_10: - h[3,10] - o[3,10] - s[3,10] <= -1
 3_DoesntSupport_3_10: h[3,10] + o[3,10] + 2 s[3,10] <= 2
 beta1_3_10: - s[3,10] + beta1[3,10] <= 0
 beta2_3_10: - s[3,10] + beta2[3,10] <= 0
 sBeta1Flag_3_10: - s[3,10] + beta1[3,10] <= 0
 sBeta2Flag_3_10: - s[3,10] + beta2[3,10] <= 0
 Eta1Overlap_3_10: eta1[3,10] + beta1[3,10] <= 1
 Eta2Overlap_3_10: eta2[3,10] + beta2[3,10] <= 1
 Eta1Flag_3_10: - xi_l[3] + xi_l[10] - 300 eta1[3,10] <= 0
 Eta2Flag_3_10: xi_r[3] - xi_r[10] - 300 eta2[3,10] <= 0
 10_Supports_3_InSameBin1_0: p_ij[3,0] - p_ij[10,0] - s[3,10] <= 0
 10_Supports_3_InSameBin2_0: - p_ij[3,0] + p_ij[10,0] - s[3,10] <= 0
 Gamma0ForItems_3_inBins_0_NoCut_10: p_ij[3,0] + gamma[3] <= 1
 10_Supports_3_InSameBin1_1: p_ij[3,1] - p_ij[10,1] - s[3,10] <= 0
 10_Supports_3_InSameBin2_1: - p_ij[3,1] + p_ij[10,1] - s[3,10] <= 0
 Gamma0ForItems_3_inBins_1_NoCut_10: p_ij[3,1] + gamma[3] <= 1
 10_Supports_3_InSameBin1_2: p_ij[3,2] - p_ij[10,2] - s[3,10] <= 0
 10_Supports_3_InSameBin2_2: - p_ij[3,2] + p_ij[10,2] - s[3,10] <= 0
 Gamma0ForItems_3_inBins_2_NoCut_10: p_ij[3,2] + gamma[3] <= 98
 10_Supports_3_InSameBin1_3: p_ij[3,3] - p_ij[10,3] - s[3,10] <= 0
 10_Supports_3_InSameBin2_3: - p_ij[3,3] + p_ij[10,3] - s[3,10] <= 0
 Gamma0ForItems_3_inBins_3_NoCut_10: p_ij[3,3] + gamma[3] <= 98
 3_OnCutIn2_10: 10000 p_ij[3,2] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 3_OnCutIn3_10: 10000 p_ij[3,3] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 AbsZ13_11: - zi[3] + zi_hi[11] - k[3,11] <= 0
 AbsZ23_11: zi[3] - zi_hi[11] - k[3,11] <= 0
 11_Below_3: zi[3] - zi_hi[11] + k[3,11] + 310 m[3,11] <= 310
 11_NotBelow_3: - zi[3] + zi_hi[11] + k[3,11] - 310 m[3,11] <= 0
 HeightCompat_3_11: h[3,11] - k[3,11] <= 0
 NoHeightCompat_3_11: - 155 h[3,11] + k[3,11] <= 0
 Overlap_3_11: - x_p[3,11] - x_p[11,3] + o[3,11] <= 0
 NoOverlap_3_11: x_p[3,11] + x_p[11,3] - 2 o[3,11] <= 0
 3_Supports_3_11: - h[3,11] - o[3,11] - s[3,11] <= -1
 3_DoesntSupport_3_11: h[3,11] + o[3,11] + 2 s[3,11] <= 2
 beta1_3_11: - s[3,11] + beta1[3,11] <= 0
 beta2_3_11: - s[3,11] + beta2[3,11] <= 0
 sBeta1Flag_3_11: - s[3,11] + beta1[3,11] <= 0
 sBeta2Flag_3_11: - s[3,11] + beta2[3,11] <= 0
 Eta1Overlap_3_11: eta1[3,11] + beta1[3,11] <= 1
 Eta2Overlap_3_11: eta2[3,11] + beta2[3,11] <= 1
 Eta1Flag_3_11: - xi_l[3] + xi_l[11] - 300 eta1[3,11] <= 0
 Eta2Flag_3_11: xi_r[3] - xi_r[11] - 300 eta2[3,11] <= 0
 11_Supports_3_InSameBin1_0: p_ij[3,0] - p_ij[11,0] - s[3,11] <= 0
 11_Supports_3_InSameBin2_0: - p_ij[3,0] + p_ij[11,0] - s[3,11] <= 0
 Gamma0ForItems_3_inBins_0_NoCut_11: p_ij[3,0] + gamma[3] <= 1
 11_Supports_3_InSameBin1_1: p_ij[3,1] - p_ij[11,1] - s[3,11] <= 0
 11_Supports_3_InSameBin2_1: - p_ij[3,1] + p_ij[11,1] - s[3,11] <= 0
 Gamma0ForItems_3_inBins_1_NoCut_11: p_ij[3,1] + gamma[3] <= 1
 11_Supports_3_InSameBin1_2: p_ij[3,2] - p_ij[11,2] - s[3,11] <= 0
 11_Supports_3_InSameBin2_2: - p_ij[3,2] + p_ij[11,2] - s[3,11] <= 0
 Gamma0ForItems_3_inBins_2_NoCut_11: p_ij[3,2] + gamma[3] <= 98
 11_Supports_3_InSameBin1_3: p_ij[3,3] - p_ij[11,3] - s[3,11] <= 0
 11_Supports_3_InSameBin2_3: - p_ij[3,3] + p_ij[11,3] - s[3,11] <= 0
 Gamma0ForItems_3_inBins_3_NoCut_11: p_ij[3,3] + gamma[3] <= 98
 3_OnCutIn2_11: 10000 p_ij[3,2] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 3_OnCutIn3_11: 10000 p_ij[3,3] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 AbsZ13_12: - zi[3] + zi_hi[12] - k[3,12] <= 0
 AbsZ23_12: zi[3] - zi_hi[12] - k[3,12] <= 0
 12_Below_3: zi[3] - zi_hi[12] + k[3,12] + 310 m[3,12] <= 310
 12_NotBelow_3: - zi[3] + zi_hi[12] + k[3,12] - 310 m[3,12] <= 0
 HeightCompat_3_12: h[3,12] - k[3,12] <= 0
 NoHeightCompat_3_12: - 155 h[3,12] + k[3,12] <= 0
 Overlap_3_12: - x_p[3,12] - x_p[12,3] + o[3,12] <= 0
 NoOverlap_3_12: x_p[3,12] + x_p[12,3] - 2 o[3,12] <= 0
 3_Supports_3_12: - h[3,12] - o[3,12] - s[3,12] <= -1
 3_DoesntSupport_3_12: h[3,12] + o[3,12] + 2 s[3,12] <= 2
 beta1_3_12: - s[3,12] + beta1[3,12] <= 0
 beta2_3_12: - s[3,12] + beta2[3,12] <= 0
 sBeta1Flag_3_12: - s[3,12] + beta1[3,12] <= 0
 sBeta2Flag_3_12: - s[3,12] + beta2[3,12] <= 0
 Eta1Overlap_3_12: eta1[3,12] + beta1[3,12] <= 1
 Eta2Overlap_3_12: eta2[3,12] + beta2[3,12] <= 1
 Eta1Flag_3_12: - xi_l[3] + xi_l[12] - 300 eta1[3,12] <= 0
 Eta2Flag_3_12: xi_r[3] - xi_r[12] - 300 eta2[3,12] <= 0
 12_Supports_3_InSameBin1_0: p_ij[3,0] - p_ij[12,0] - s[3,12] <= 0
 12_Supports_3_InSameBin2_0: - p_ij[3,0] + p_ij[12,0] - s[3,12] <= 0
 Gamma0ForItems_3_inBins_0_NoCut_12: p_ij[3,0] + gamma[3] <= 1
 12_Supports_3_InSameBin1_1: p_ij[3,1] - p_ij[12,1] - s[3,12] <= 0
 12_Supports_3_InSameBin2_1: - p_ij[3,1] + p_ij[12,1] - s[3,12] <= 0
 Gamma0ForItems_3_inBins_1_NoCut_12: p_ij[3,1] + gamma[3] <= 1
 12_Supports_3_InSameBin1_2: p_ij[3,2] - p_ij[12,2] - s[3,12] <= 0
 12_Supports_3_InSameBin2_2: - p_ij[3,2] + p_ij[12,2] - s[3,12] <= 0
 Gamma0ForItems_3_inBins_2_NoCut_12: p_ij[3,2] + gamma[3] <= 98
 12_Supports_3_InSameBin1_3: p_ij[3,3] - p_ij[12,3] - s[3,12] <= 0
 12_Supports_3_InSameBin2_3: - p_ij[3,3] + p_ij[12,3] - s[3,12] <= 0
 Gamma0ForItems_3_inBins_3_NoCut_12: p_ij[3,3] + gamma[3] <= 98
 3_OnCutIn2_12: 10000 p_ij[3,2] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 3_OnCutIn3_12: 10000 p_ij[3,3] + 1.261904761904762 xi_l[3] + zi[3]
   + 10000 gamma[3] <= 20053
 GroundConstraint_4: zi[4] + 155 g[4] <= 155
 AbsZ14_0: - zi[4] + zi_hi[0] - k[4,0] <= 0
 AbsZ24_0: zi[4] - zi_hi[0] - k[4,0] <= 0
 0_Below_4: zi[4] - zi_hi[0] + k[4,0] + 310 m[4,0] <= 310
 0_NotBelow_4: - zi[4] + zi_hi[0] + k[4,0] - 310 m[4,0] <= 0
 HeightCompat_4_0: h[4,0] - k[4,0] <= 0
 NoHeightCompat_4_0: - 155 h[4,0] + k[4,0] <= 0
 Overlap_4_0: - x_p[0,4] - x_p[4,0] + o[4,0] <= 0
 NoOverlap_4_0: x_p[0,4] + x_p[4,0] - 2 o[4,0] <= 0
 3_Supports_4_0: - h[4,0] - o[4,0] - s[4,0] <= -1
 3_DoesntSupport_4_0: h[4,0] + o[4,0] + 2 s[4,0] <= 2
 beta1_4_0: - s[4,0] + beta1[4,0] <= 0
 beta2_4_0: - s[4,0] + beta2[4,0] <= 0
 sBeta1Flag_4_0: - s[4,0] + beta1[4,0] <= 0
 sBeta2Flag_4_0: - s[4,0] + beta2[4,0] <= 0
 Eta1Overlap_4_0: eta1[4,0] + beta1[4,0] <= 1
 Eta2Overlap_4_0: eta2[4,0] + beta2[4,0] <= 1
 Eta1Flag_4_0: xi_l[0] - xi_l[4] - 300 eta1[4,0] <= 0
 Eta2Flag_4_0: - xi_r[0] + xi_r[4] - 300 eta2[4,0] <= 0
 0_Supports_4_InSameBin1_0: - p_ij[0,0] + p_ij[4,0] - s[4,0] <= 0
 0_Supports_4_InSameBin2_0: p_ij[0,0] - p_ij[4,0] - s[4,0] <= 0
 Gamma0ForItems_4_inBins_0_NoCut_0: p_ij[4,0] + gamma[4] <= 1
 0_Supports_4_InSameBin1_1: - p_ij[0,1] + p_ij[4,1] - s[4,0] <= 0
 0_Supports_4_InSameBin2_1: p_ij[0,1] - p_ij[4,1] - s[4,0] <= 0
 Gamma0ForItems_4_inBins_1_NoCut_0: p_ij[4,1] + gamma[4] <= 1
 0_Supports_4_InSameBin1_2: - p_ij[0,2] + p_ij[4,2] - s[4,0] <= 0
 0_Supports_4_InSameBin2_2: p_ij[0,2] - p_ij[4,2] - s[4,0] <= 0
 Gamma0ForItems_4_inBins_2_NoCut_0: p_ij[4,2] + gamma[4] <= 98
 0_Supports_4_InSameBin1_3: - p_ij[0,3] + p_ij[4,3] - s[4,0] <= 0
 0_Supports_4_InSameBin2_3: p_ij[0,3] - p_ij[4,3] - s[4,0] <= 0
 Gamma0ForItems_4_inBins_3_NoCut_0: p_ij[4,3] + gamma[4] <= 98
 4_OnCutIn2_0: 10000 p_ij[4,2] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 4_OnCutIn3_0: 10000 p_ij[4,3] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 AbsZ14_1: - zi[4] + zi_hi[1] - k[4,1] <= 0
 AbsZ24_1: zi[4] - zi_hi[1] - k[4,1] <= 0
 1_Below_4: zi[4] - zi_hi[1] + k[4,1] + 310 m[4,1] <= 310
 1_NotBelow_4: - zi[4] + zi_hi[1] + k[4,1] - 310 m[4,1] <= 0
 HeightCompat_4_1: h[4,1] - k[4,1] <= 0
 NoHeightCompat_4_1: - 155 h[4,1] + k[4,1] <= 0
 Overlap_4_1: - x_p[1,4] - x_p[4,1] + o[4,1] <= 0
 NoOverlap_4_1: x_p[1,4] + x_p[4,1] - 2 o[4,1] <= 0
 3_Supports_4_1: - h[4,1] - o[4,1] - s[4,1] <= -1
 3_DoesntSupport_4_1: h[4,1] + o[4,1] + 2 s[4,1] <= 2
 beta1_4_1: - s[4,1] + beta1[4,1] <= 0
 beta2_4_1: - s[4,1] + beta2[4,1] <= 0
 sBeta1Flag_4_1: - s[4,1] + beta1[4,1] <= 0
 sBeta2Flag_4_1: - s[4,1] + beta2[4,1] <= 0
 Eta1Overlap_4_1: eta1[4,1] + beta1[4,1] <= 1
 Eta2Overlap_4_1: eta2[4,1] + beta2[4,1] <= 1
 Eta1Flag_4_1: xi_l[1] - xi_l[4] - 300 eta1[4,1] <= 0
 Eta2Flag_4_1: - xi_r[1] + xi_r[4] - 300 eta2[4,1] <= 0
 1_Supports_4_InSameBin1_0: - p_ij[1,0] + p_ij[4,0] - s[4,1] <= 0
 1_Supports_4_InSameBin2_0: p_ij[1,0] - p_ij[4,0] - s[4,1] <= 0
 Gamma0ForItems_4_inBins_0_NoCut_1: p_ij[4,0] + gamma[4] <= 1
 1_Supports_4_InSameBin1_1: - p_ij[1,1] + p_ij[4,1] - s[4,1] <= 0
 1_Supports_4_InSameBin2_1: p_ij[1,1] - p_ij[4,1] - s[4,1] <= 0
 Gamma0ForItems_4_inBins_1_NoCut_1: p_ij[4,1] + gamma[4] <= 1
 1_Supports_4_InSameBin1_2: - p_ij[1,2] + p_ij[4,2] - s[4,1] <= 0
 1_Supports_4_InSameBin2_2: p_ij[1,2] - p_ij[4,2] - s[4,1] <= 0
 Gamma0ForItems_4_inBins_2_NoCut_1: p_ij[4,2] + gamma[4] <= 98
 1_Supports_4_InSameBin1_3: - p_ij[1,3] + p_ij[4,3] - s[4,1] <= 0
 1_Supports_4_InSameBin2_3: p_ij[1,3] - p_ij[4,3] - s[4,1] <= 0
 Gamma0ForItems_4_inBins_3_NoCut_1: p_ij[4,3] + gamma[4] <= 98
 4_OnCutIn2_1: 10000 p_ij[4,2] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 4_OnCutIn3_1: 10000 p_ij[4,3] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 AbsZ14_2: - zi[4] + zi_hi[2] - k[4,2] <= 0
 AbsZ24_2: zi[4] - zi_hi[2] - k[4,2] <= 0
 2_Below_4: zi[4] - zi_hi[2] + k[4,2] + 310 m[4,2] <= 310
 2_NotBelow_4: - zi[4] + zi_hi[2] + k[4,2] - 310 m[4,2] <= 0
 HeightCompat_4_2: h[4,2] - k[4,2] <= 0
 NoHeightCompat_4_2: - 155 h[4,2] + k[4,2] <= 0
 Overlap_4_2: - x_p[2,4] - x_p[4,2] + o[4,2] <= 0
 NoOverlap_4_2: x_p[2,4] + x_p[4,2] - 2 o[4,2] <= 0
 3_Supports_4_2: - h[4,2] - o[4,2] - s[4,2] <= -1
 3_DoesntSupport_4_2: h[4,2] + o[4,2] + 2 s[4,2] <= 2
 beta1_4_2: - s[4,2] + beta1[4,2] <= 0
 beta2_4_2: - s[4,2] + beta2[4,2] <= 0
 sBeta1Flag_4_2: - s[4,2] + beta1[4,2] <= 0
 sBeta2Flag_4_2: - s[4,2] + beta2[4,2] <= 0
 Eta1Overlap_4_2: eta1[4,2] + beta1[4,2] <= 1
 Eta2Overlap_4_2: eta2[4,2] + beta2[4,2] <= 1
 Eta1Flag_4_2: xi_l[2] - xi_l[4] - 300 eta1[4,2] <= 0
 Eta2Flag_4_2: - xi_r[2] + xi_r[4] - 300 eta2[4,2] <= 0
 2_Supports_4_InSameBin1_0: - p_ij[2,0] + p_ij[4,0] - s[4,2] <= 0
 2_Supports_4_InSameBin2_0: p_ij[2,0] - p_ij[4,0] - s[4,2] <= 0
 Gamma0ForItems_4_inBins_0_NoCut_2: p_ij[4,0] + gamma[4] <= 1
 2_Supports_4_InSameBin1_1: - p_ij[2,1] + p_ij[4,1] - s[4,2] <= 0
 2_Supports_4_InSameBin2_1: p_ij[2,1] - p_ij[4,1] - s[4,2] <= 0
 Gamma0ForItems_4_inBins_1_NoCut_2: p_ij[4,1] + gamma[4] <= 1
 2_Supports_4_InSameBin1_2: - p_ij[2,2] + p_ij[4,2] - s[4,2] <= 0
 2_Supports_4_InSameBin2_2: p_ij[2,2] - p_ij[4,2] - s[4,2] <= 0
 Gamma0ForItems_4_inBins_2_NoCut_2: p_ij[4,2] + gamma[4] <= 98
 2_Supports_4_InSameBin1_3: - p_ij[2,3] + p_ij[4,3] - s[4,2] <= 0
 2_Supports_4_InSameBin2_3: p_ij[2,3] - p_ij[4,3] - s[4,2] <= 0
 Gamma0ForItems_4_inBins_3_NoCut_2: p_ij[4,3] + gamma[4] <= 98
 4_OnCutIn2_2: 10000 p_ij[4,2] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 4_OnCutIn3_2: 10000 p_ij[4,3] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 AbsZ14_3: - zi[4] + zi_hi[3] - k[4,3] <= 0
 AbsZ24_3: zi[4] - zi_hi[3] - k[4,3] <= 0
 3_Below_4: zi[4] - zi_hi[3] + k[4,3] + 310 m[4,3] <= 310
 3_NotBelow_4: - zi[4] + zi_hi[3] + k[4,3] - 310 m[4,3] <= 0
 HeightCompat_4_3: h[4,3] - k[4,3] <= 0
 NoHeightCompat_4_3: - 155 h[4,3] + k[4,3] <= 0
 Overlap_4_3: - x_p[3,4] - x_p[4,3] + o[4,3] <= 0
 NoOverlap_4_3: x_p[3,4] + x_p[4,3] - 2 o[4,3] <= 0
 3_Supports_4_3: - h[4,3] - o[4,3] - s[4,3] <= -1
 3_DoesntSupport_4_3: h[4,3] + o[4,3] + 2 s[4,3] <= 2
 beta1_4_3: - s[4,3] + beta1[4,3] <= 0
 beta2_4_3: - s[4,3] + beta2[4,3] <= 0
 sBeta1Flag_4_3: - s[4,3] + beta1[4,3] <= 0
 sBeta2Flag_4_3: - s[4,3] + beta2[4,3] <= 0
 Eta1Overlap_4_3: eta1[4,3] + beta1[4,3] <= 1
 Eta2Overlap_4_3: eta2[4,3] + beta2[4,3] <= 1
 Eta1Flag_4_3: xi_l[3] - xi_l[4] - 300 eta1[4,3] <= 0
 Eta2Flag_4_3: - xi_r[3] + xi_r[4] - 300 eta2[4,3] <= 0
 3_Supports_4_InSameBin1_0: - p_ij[3,0] + p_ij[4,0] - s[4,3] <= 0
 3_Supports_4_InSameBin2_0: p_ij[3,0] - p_ij[4,0] - s[4,3] <= 0
 Gamma0ForItems_4_inBins_0_NoCut_3: p_ij[4,0] + gamma[4] <= 1
 3_Supports_4_InSameBin1_1: - p_ij[3,1] + p_ij[4,1] - s[4,3] <= 0
 3_Supports_4_InSameBin2_1: p_ij[3,1] - p_ij[4,1] - s[4,3] <= 0
 Gamma0ForItems_4_inBins_1_NoCut_3: p_ij[4,1] + gamma[4] <= 1
 3_Supports_4_InSameBin1_2: - p_ij[3,2] + p_ij[4,2] - s[4,3] <= 0
 3_Supports_4_InSameBin2_2: p_ij[3,2] - p_ij[4,2] - s[4,3] <= 0
 Gamma0ForItems_4_inBins_2_NoCut_3: p_ij[4,2] + gamma[4] <= 98
 3_Supports_4_InSameBin1_3: - p_ij[3,3] + p_ij[4,3] - s[4,3] <= 0
 3_Supports_4_InSameBin2_3: p_ij[3,3] - p_ij[4,3] - s[4,3] <= 0
 Gamma0ForItems_4_inBins_3_NoCut_3: p_ij[4,3] + gamma[4] <= 98
 4_OnCutIn2_3: 10000 p_ij[4,2] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 4_OnCutIn3_3: 10000 p_ij[4,3] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 4_OnCutIn2_4: 10000 p_ij[4,2] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 4_OnCutIn3_4: 10000 p_ij[4,3] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 AbsZ14_5: - zi[4] + zi_hi[5] - k[4,5] <= 0
 AbsZ24_5: zi[4] - zi_hi[5] - k[4,5] <= 0
 5_Below_4: zi[4] - zi_hi[5] + k[4,5] + 310 m[4,5] <= 310
 5_NotBelow_4: - zi[4] + zi_hi[5] + k[4,5] - 310 m[4,5] <= 0
 HeightCompat_4_5: h[4,5] - k[4,5] <= 0
 NoHeightCompat_4_5: - 155 h[4,5] + k[4,5] <= 0
 Overlap_4_5: - x_p[4,5] - x_p[5,4] + o[4,5] <= 0
 NoOverlap_4_5: x_p[4,5] + x_p[5,4] - 2 o[4,5] <= 0
 3_Supports_4_5: - h[4,5] - o[4,5] - s[4,5] <= -1
 3_DoesntSupport_4_5: h[4,5] + o[4,5] + 2 s[4,5] <= 2
 beta1_4_5: - s[4,5] + beta1[4,5] <= 0
 beta2_4_5: - s[4,5] + beta2[4,5] <= 0
 sBeta1Flag_4_5: - s[4,5] + beta1[4,5] <= 0
 sBeta2Flag_4_5: - s[4,5] + beta2[4,5] <= 0
 Eta1Overlap_4_5: eta1[4,5] + beta1[4,5] <= 1
 Eta2Overlap_4_5: eta2[4,5] + beta2[4,5] <= 1
 Eta1Flag_4_5: - xi_l[4] + xi_l[5] - 300 eta1[4,5] <= 0
 Eta2Flag_4_5: xi_r[4] - xi_r[5] - 300 eta2[4,5] <= 0
 5_Supports_4_InSameBin1_0: p_ij[4,0] - p_ij[5,0] - s[4,5] <= 0
 5_Supports_4_InSameBin2_0: - p_ij[4,0] + p_ij[5,0] - s[4,5] <= 0
 Gamma0ForItems_4_inBins_0_NoCut_5: p_ij[4,0] + gamma[4] <= 1
 5_Supports_4_InSameBin1_1: p_ij[4,1] - p_ij[5,1] - s[4,5] <= 0
 5_Supports_4_InSameBin2_1: - p_ij[4,1] + p_ij[5,1] - s[4,5] <= 0
 Gamma0ForItems_4_inBins_1_NoCut_5: p_ij[4,1] + gamma[4] <= 1
 5_Supports_4_InSameBin1_2: p_ij[4,2] - p_ij[5,2] - s[4,5] <= 0
 5_Supports_4_InSameBin2_2: - p_ij[4,2] + p_ij[5,2] - s[4,5] <= 0
 Gamma0ForItems_4_inBins_2_NoCut_5: p_ij[4,2] + gamma[4] <= 98
 5_Supports_4_InSameBin1_3: p_ij[4,3] - p_ij[5,3] - s[4,5] <= 0
 5_Supports_4_InSameBin2_3: - p_ij[4,3] + p_ij[5,3] - s[4,5] <= 0
 Gamma0ForItems_4_inBins_3_NoCut_5: p_ij[4,3] + gamma[4] <= 98
 4_OnCutIn2_5: 10000 p_ij[4,2] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 4_OnCutIn3_5: 10000 p_ij[4,3] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 AbsZ14_6: - zi[4] + zi_hi[6] - k[4,6] <= 0
 AbsZ24_6: zi[4] - zi_hi[6] - k[4,6] <= 0
 6_Below_4: zi[4] - zi_hi[6] + k[4,6] + 310 m[4,6] <= 310
 6_NotBelow_4: - zi[4] + zi_hi[6] + k[4,6] - 310 m[4,6] <= 0
 HeightCompat_4_6: h[4,6] - k[4,6] <= 0
 NoHeightCompat_4_6: - 155 h[4,6] + k[4,6] <= 0
 Overlap_4_6: - x_p[4,6] - x_p[6,4] + o[4,6] <= 0
 NoOverlap_4_6: x_p[4,6] + x_p[6,4] - 2 o[4,6] <= 0
 3_Supports_4_6: - h[4,6] - o[4,6] - s[4,6] <= -1
 3_DoesntSupport_4_6: h[4,6] + o[4,6] + 2 s[4,6] <= 2
 beta1_4_6: - s[4,6] + beta1[4,6] <= 0
 beta2_4_6: - s[4,6] + beta2[4,6] <= 0
 sBeta1Flag_4_6: - s[4,6] + beta1[4,6] <= 0
 sBeta2Flag_4_6: - s[4,6] + beta2[4,6] <= 0
 Eta1Overlap_4_6: eta1[4,6] + beta1[4,6] <= 1
 Eta2Overlap_4_6: eta2[4,6] + beta2[4,6] <= 1
 Eta1Flag_4_6: - xi_l[4] + xi_l[6] - 300 eta1[4,6] <= 0
 Eta2Flag_4_6: xi_r[4] - xi_r[6] - 300 eta2[4,6] <= 0
 6_Supports_4_InSameBin1_0: p_ij[4,0] - p_ij[6,0] - s[4,6] <= 0
 6_Supports_4_InSameBin2_0: - p_ij[4,0] + p_ij[6,0] - s[4,6] <= 0
 Gamma0ForItems_4_inBins_0_NoCut_6: p_ij[4,0] + gamma[4] <= 1
 6_Supports_4_InSameBin1_1: p_ij[4,1] - p_ij[6,1] - s[4,6] <= 0
 6_Supports_4_InSameBin2_1: - p_ij[4,1] + p_ij[6,1] - s[4,6] <= 0
 Gamma0ForItems_4_inBins_1_NoCut_6: p_ij[4,1] + gamma[4] <= 1
 6_Supports_4_InSameBin1_2: p_ij[4,2] - p_ij[6,2] - s[4,6] <= 0
 6_Supports_4_InSameBin2_2: - p_ij[4,2] + p_ij[6,2] - s[4,6] <= 0
 Gamma0ForItems_4_inBins_2_NoCut_6: p_ij[4,2] + gamma[4] <= 98
 6_Supports_4_InSameBin1_3: p_ij[4,3] - p_ij[6,3] - s[4,6] <= 0
 6_Supports_4_InSameBin2_3: - p_ij[4,3] + p_ij[6,3] - s[4,6] <= 0
 Gamma0ForItems_4_inBins_3_NoCut_6: p_ij[4,3] + gamma[4] <= 98
 4_OnCutIn2_6: 10000 p_ij[4,2] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 4_OnCutIn3_6: 10000 p_ij[4,3] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 AbsZ14_7: - zi[4] + zi_hi[7] - k[4,7] <= 0
 AbsZ24_7: zi[4] - zi_hi[7] - k[4,7] <= 0
 7_Below_4: zi[4] - zi_hi[7] + k[4,7] + 310 m[4,7] <= 310
 7_NotBelow_4: - zi[4] + zi_hi[7] + k[4,7] - 310 m[4,7] <= 0
 HeightCompat_4_7: h[4,7] - k[4,7] <= 0
 NoHeightCompat_4_7: - 155 h[4,7] + k[4,7] <= 0
 Overlap_4_7: - x_p[4,7] - x_p[7,4] + o[4,7] <= 0
 NoOverlap_4_7: x_p[4,7] + x_p[7,4] - 2 o[4,7] <= 0
 3_Supports_4_7: - h[4,7] - o[4,7] - s[4,7] <= -1
 3_DoesntSupport_4_7: h[4,7] + o[4,7] + 2 s[4,7] <= 2
 beta1_4_7: - s[4,7] + beta1[4,7] <= 0
 beta2_4_7: - s[4,7] + beta2[4,7] <= 0
 sBeta1Flag_4_7: - s[4,7] + beta1[4,7] <= 0
 sBeta2Flag_4_7: - s[4,7] + beta2[4,7] <= 0
 Eta1Overlap_4_7: eta1[4,7] + beta1[4,7] <= 1
 Eta2Overlap_4_7: eta2[4,7] + beta2[4,7] <= 1
 Eta1Flag_4_7: - xi_l[4] + xi_l[7] - 300 eta1[4,7] <= 0
 Eta2Flag_4_7: xi_r[4] - xi_r[7] - 300 eta2[4,7] <= 0
 7_Supports_4_InSameBin1_0: p_ij[4,0] - p_ij[7,0] - s[4,7] <= 0
 7_Supports_4_InSameBin2_0: - p_ij[4,0] + p_ij[7,0] - s[4,7] <= 0
 Gamma0ForItems_4_inBins_0_NoCut_7: p_ij[4,0] + gamma[4] <= 1
 7_Supports_4_InSameBin1_1: p_ij[4,1] - p_ij[7,1] - s[4,7] <= 0
 7_Supports_4_InSameBin2_1: - p_ij[4,1] + p_ij[7,1] - s[4,7] <= 0
 Gamma0ForItems_4_inBins_1_NoCut_7: p_ij[4,1] + gamma[4] <= 1
 7_Supports_4_InSameBin1_2: p_ij[4,2] - p_ij[7,2] - s[4,7] <= 0
 7_Supports_4_InSameBin2_2: - p_ij[4,2] + p_ij[7,2] - s[4,7] <= 0
 Gamma0ForItems_4_inBins_2_NoCut_7: p_ij[4,2] + gamma[4] <= 98
 7_Supports_4_InSameBin1_3: p_ij[4,3] - p_ij[7,3] - s[4,7] <= 0
 7_Supports_4_InSameBin2_3: - p_ij[4,3] + p_ij[7,3] - s[4,7] <= 0
 Gamma0ForItems_4_inBins_3_NoCut_7: p_ij[4,3] + gamma[4] <= 98
 4_OnCutIn2_7: 10000 p_ij[4,2] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 4_OnCutIn3_7: 10000 p_ij[4,3] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 AbsZ14_8: - zi[4] + zi_hi[8] - k[4,8] <= 0
 AbsZ24_8: zi[4] - zi_hi[8] - k[4,8] <= 0
 8_Below_4: zi[4] - zi_hi[8] + k[4,8] + 310 m[4,8] <= 310
 8_NotBelow_4: - zi[4] + zi_hi[8] + k[4,8] - 310 m[4,8] <= 0
 HeightCompat_4_8: h[4,8] - k[4,8] <= 0
 NoHeightCompat_4_8: - 155 h[4,8] + k[4,8] <= 0
 Overlap_4_8: - x_p[4,8] - x_p[8,4] + o[4,8] <= 0
 NoOverlap_4_8: x_p[4,8] + x_p[8,4] - 2 o[4,8] <= 0
 3_Supports_4_8: - h[4,8] - o[4,8] - s[4,8] <= -1
 3_DoesntSupport_4_8: h[4,8] + o[4,8] + 2 s[4,8] <= 2
 beta1_4_8: - s[4,8] + beta1[4,8] <= 0
 beta2_4_8: - s[4,8] + beta2[4,8] <= 0
 sBeta1Flag_4_8: - s[4,8] + beta1[4,8] <= 0
 sBeta2Flag_4_8: - s[4,8] + beta2[4,8] <= 0
 Eta1Overlap_4_8: eta1[4,8] + beta1[4,8] <= 1
 Eta2Overlap_4_8: eta2[4,8] + beta2[4,8] <= 1
 Eta1Flag_4_8: - xi_l[4] + xi_l[8] - 300 eta1[4,8] <= 0
 Eta2Flag_4_8: xi_r[4] - xi_r[8] - 300 eta2[4,8] <= 0
 8_Supports_4_InSameBin1_0: p_ij[4,0] - p_ij[8,0] - s[4,8] <= 0
 8_Supports_4_InSameBin2_0: - p_ij[4,0] + p_ij[8,0] - s[4,8] <= 0
 Gamma0ForItems_4_inBins_0_NoCut_8: p_ij[4,0] + gamma[4] <= 1
 8_Supports_4_InSameBin1_1: p_ij[4,1] - p_ij[8,1] - s[4,8] <= 0
 8_Supports_4_InSameBin2_1: - p_ij[4,1] + p_ij[8,1] - s[4,8] <= 0
 Gamma0ForItems_4_inBins_1_NoCut_8: p_ij[4,1] + gamma[4] <= 1
 8_Supports_4_InSameBin1_2: p_ij[4,2] - p_ij[8,2] - s[4,8] <= 0
 8_Supports_4_InSameBin2_2: - p_ij[4,2] + p_ij[8,2] - s[4,8] <= 0
 Gamma0ForItems_4_inBins_2_NoCut_8: p_ij[4,2] + gamma[4] <= 98
 8_Supports_4_InSameBin1_3: p_ij[4,3] - p_ij[8,3] - s[4,8] <= 0
 8_Supports_4_InSameBin2_3: - p_ij[4,3] + p_ij[8,3] - s[4,8] <= 0
 Gamma0ForItems_4_inBins_3_NoCut_8: p_ij[4,3] + gamma[4] <= 98
 4_OnCutIn2_8: 10000 p_ij[4,2] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 4_OnCutIn3_8: 10000 p_ij[4,3] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 AbsZ14_9: - zi[4] + zi_hi[9] - k[4,9] <= 0
 AbsZ24_9: zi[4] - zi_hi[9] - k[4,9] <= 0
 9_Below_4: zi[4] - zi_hi[9] + k[4,9] + 310 m[4,9] <= 310
 9_NotBelow_4: - zi[4] + zi_hi[9] + k[4,9] - 310 m[4,9] <= 0
 HeightCompat_4_9: h[4,9] - k[4,9] <= 0
 NoHeightCompat_4_9: - 155 h[4,9] + k[4,9] <= 0
 Overlap_4_9: - x_p[4,9] - x_p[9,4] + o[4,9] <= 0
 NoOverlap_4_9: x_p[4,9] + x_p[9,4] - 2 o[4,9] <= 0
 3_Supports_4_9: - h[4,9] - o[4,9] - s[4,9] <= -1
 3_DoesntSupport_4_9: h[4,9] + o[4,9] + 2 s[4,9] <= 2
 beta1_4_9: - s[4,9] + beta1[4,9] <= 0
 beta2_4_9: - s[4,9] + beta2[4,9] <= 0
 sBeta1Flag_4_9: - s[4,9] + beta1[4,9] <= 0
 sBeta2Flag_4_9: - s[4,9] + beta2[4,9] <= 0
 Eta1Overlap_4_9: eta1[4,9] + beta1[4,9] <= 1
 Eta2Overlap_4_9: eta2[4,9] + beta2[4,9] <= 1
 Eta1Flag_4_9: - xi_l[4] + xi_l[9] - 300 eta1[4,9] <= 0
 Eta2Flag_4_9: xi_r[4] - xi_r[9] - 300 eta2[4,9] <= 0
 9_Supports_4_InSameBin1_0: p_ij[4,0] - p_ij[9,0] - s[4,9] <= 0
 9_Supports_4_InSameBin2_0: - p_ij[4,0] + p_ij[9,0] - s[4,9] <= 0
 Gamma0ForItems_4_inBins_0_NoCut_9: p_ij[4,0] + gamma[4] <= 1
 9_Supports_4_InSameBin1_1: p_ij[4,1] - p_ij[9,1] - s[4,9] <= 0
 9_Supports_4_InSameBin2_1: - p_ij[4,1] + p_ij[9,1] - s[4,9] <= 0
 Gamma0ForItems_4_inBins_1_NoCut_9: p_ij[4,1] + gamma[4] <= 1
 9_Supports_4_InSameBin1_2: p_ij[4,2] - p_ij[9,2] - s[4,9] <= 0
 9_Supports_4_InSameBin2_2: - p_ij[4,2] + p_ij[9,2] - s[4,9] <= 0
 Gamma0ForItems_4_inBins_2_NoCut_9: p_ij[4,2] + gamma[4] <= 98
 9_Supports_4_InSameBin1_3: p_ij[4,3] - p_ij[9,3] - s[4,9] <= 0
 9_Supports_4_InSameBin2_3: - p_ij[4,3] + p_ij[9,3] - s[4,9] <= 0
 Gamma0ForItems_4_inBins_3_NoCut_9: p_ij[4,3] + gamma[4] <= 98
 4_OnCutIn2_9: 10000 p_ij[4,2] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 4_OnCutIn3_9: 10000 p_ij[4,3] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 AbsZ14_10: - zi[4] + zi_hi[10] - k[4,10] <= 0
 AbsZ24_10: zi[4] - zi_hi[10] - k[4,10] <= 0
 10_Below_4: zi[4] - zi_hi[10] + k[4,10] + 310 m[4,10] <= 310
 10_NotBelow_4: - zi[4] + zi_hi[10] + k[4,10] - 310 m[4,10] <= 0
 HeightCompat_4_10: h[4,10] - k[4,10] <= 0
 NoHeightCompat_4_10: - 155 h[4,10] + k[4,10] <= 0
 Overlap_4_10: - x_p[4,10] - x_p[10,4] + o[4,10] <= 0
 NoOverlap_4_10: x_p[4,10] + x_p[10,4] - 2 o[4,10] <= 0
 3_Supports_4_10: - h[4,10] - o[4,10] - s[4,10] <= -1
 3_DoesntSupport_4_10: h[4,10] + o[4,10] + 2 s[4,10] <= 2
 beta1_4_10: - s[4,10] + beta1[4,10] <= 0
 beta2_4_10: - s[4,10] + beta2[4,10] <= 0
 sBeta1Flag_4_10: - s[4,10] + beta1[4,10] <= 0
 sBeta2Flag_4_10: - s[4,10] + beta2[4,10] <= 0
 Eta1Overlap_4_10: eta1[4,10] + beta1[4,10] <= 1
 Eta2Overlap_4_10: eta2[4,10] + beta2[4,10] <= 1
 Eta1Flag_4_10: - xi_l[4] + xi_l[10] - 300 eta1[4,10] <= 0
 Eta2Flag_4_10: xi_r[4] - xi_r[10] - 300 eta2[4,10] <= 0
 10_Supports_4_InSameBin1_0: p_ij[4,0] - p_ij[10,0] - s[4,10] <= 0
 10_Supports_4_InSameBin2_0: - p_ij[4,0] + p_ij[10,0] - s[4,10] <= 0
 Gamma0ForItems_4_inBins_0_NoCut_10: p_ij[4,0] + gamma[4] <= 1
 10_Supports_4_InSameBin1_1: p_ij[4,1] - p_ij[10,1] - s[4,10] <= 0
 10_Supports_4_InSameBin2_1: - p_ij[4,1] + p_ij[10,1] - s[4,10] <= 0
 Gamma0ForItems_4_inBins_1_NoCut_10: p_ij[4,1] + gamma[4] <= 1
 10_Supports_4_InSameBin1_2: p_ij[4,2] - p_ij[10,2] - s[4,10] <= 0
 10_Supports_4_InSameBin2_2: - p_ij[4,2] + p_ij[10,2] - s[4,10] <= 0
 Gamma0ForItems_4_inBins_2_NoCut_10: p_ij[4,2] + gamma[4] <= 98
 10_Supports_4_InSameBin1_3: p_ij[4,3] - p_ij[10,3] - s[4,10] <= 0
 10_Supports_4_InSameBin2_3: - p_ij[4,3] + p_ij[10,3] - s[4,10] <= 0
 Gamma0ForItems_4_inBins_3_NoCut_10: p_ij[4,3] + gamma[4] <= 98
 4_OnCutIn2_10: 10000 p_ij[4,2] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 4_OnCutIn3_10: 10000 p_ij[4,3] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 AbsZ14_11: - zi[4] + zi_hi[11] - k[4,11] <= 0
 AbsZ24_11: zi[4] - zi_hi[11] - k[4,11] <= 0
 11_Below_4: zi[4] - zi_hi[11] + k[4,11] + 310 m[4,11] <= 310
 11_NotBelow_4: - zi[4] + zi_hi[11] + k[4,11] - 310 m[4,11] <= 0
 HeightCompat_4_11: h[4,11] - k[4,11] <= 0
 NoHeightCompat_4_11: - 155 h[4,11] + k[4,11] <= 0
 Overlap_4_11: - x_p[4,11] - x_p[11,4] + o[4,11] <= 0
 NoOverlap_4_11: x_p[4,11] + x_p[11,4] - 2 o[4,11] <= 0
 3_Supports_4_11: - h[4,11] - o[4,11] - s[4,11] <= -1
 3_DoesntSupport_4_11: h[4,11] + o[4,11] + 2 s[4,11] <= 2
 beta1_4_11: - s[4,11] + beta1[4,11] <= 0
 beta2_4_11: - s[4,11] + beta2[4,11] <= 0
 sBeta1Flag_4_11: - s[4,11] + beta1[4,11] <= 0
 sBeta2Flag_4_11: - s[4,11] + beta2[4,11] <= 0
 Eta1Overlap_4_11: eta1[4,11] + beta1[4,11] <= 1
 Eta2Overlap_4_11: eta2[4,11] + beta2[4,11] <= 1
 Eta1Flag_4_11: - xi_l[4] + xi_l[11] - 300 eta1[4,11] <= 0
 Eta2Flag_4_11: xi_r[4] - xi_r[11] - 300 eta2[4,11] <= 0
 11_Supports_4_InSameBin1_0: p_ij[4,0] - p_ij[11,0] - s[4,11] <= 0
 11_Supports_4_InSameBin2_0: - p_ij[4,0] + p_ij[11,0] - s[4,11] <= 0
 Gamma0ForItems_4_inBins_0_NoCut_11: p_ij[4,0] + gamma[4] <= 1
 11_Supports_4_InSameBin1_1: p_ij[4,1] - p_ij[11,1] - s[4,11] <= 0
 11_Supports_4_InSameBin2_1: - p_ij[4,1] + p_ij[11,1] - s[4,11] <= 0
 Gamma0ForItems_4_inBins_1_NoCut_11: p_ij[4,1] + gamma[4] <= 1
 11_Supports_4_InSameBin1_2: p_ij[4,2] - p_ij[11,2] - s[4,11] <= 0
 11_Supports_4_InSameBin2_2: - p_ij[4,2] + p_ij[11,2] - s[4,11] <= 0
 Gamma0ForItems_4_inBins_2_NoCut_11: p_ij[4,2] + gamma[4] <= 98
 11_Supports_4_InSameBin1_3: p_ij[4,3] - p_ij[11,3] - s[4,11] <= 0
 11_Supports_4_InSameBin2_3: - p_ij[4,3] + p_ij[11,3] - s[4,11] <= 0
 Gamma0ForItems_4_inBins_3_NoCut_11: p_ij[4,3] + gamma[4] <= 98
 4_OnCutIn2_11: 10000 p_ij[4,2] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 4_OnCutIn3_11: 10000 p_ij[4,3] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 AbsZ14_12: - zi[4] + zi_hi[12] - k[4,12] <= 0
 AbsZ24_12: zi[4] - zi_hi[12] - k[4,12] <= 0
 12_Below_4: zi[4] - zi_hi[12] + k[4,12] + 310 m[4,12] <= 310
 12_NotBelow_4: - zi[4] + zi_hi[12] + k[4,12] - 310 m[4,12] <= 0
 HeightCompat_4_12: h[4,12] - k[4,12] <= 0
 NoHeightCompat_4_12: - 155 h[4,12] + k[4,12] <= 0
 Overlap_4_12: - x_p[4,12] - x_p[12,4] + o[4,12] <= 0
 NoOverlap_4_12: x_p[4,12] + x_p[12,4] - 2 o[4,12] <= 0
 3_Supports_4_12: - h[4,12] - o[4,12] - s[4,12] <= -1
 3_DoesntSupport_4_12: h[4,12] + o[4,12] + 2 s[4,12] <= 2
 beta1_4_12: - s[4,12] + beta1[4,12] <= 0
 beta2_4_12: - s[4,12] + beta2[4,12] <= 0
 sBeta1Flag_4_12: - s[4,12] + beta1[4,12] <= 0
 sBeta2Flag_4_12: - s[4,12] + beta2[4,12] <= 0
 Eta1Overlap_4_12: eta1[4,12] + beta1[4,12] <= 1
 Eta2Overlap_4_12: eta2[4,12] + beta2[4,12] <= 1
 Eta1Flag_4_12: - xi_l[4] + xi_l[12] - 300 eta1[4,12] <= 0
 Eta2Flag_4_12: xi_r[4] - xi_r[12] - 300 eta2[4,12] <= 0
 12_Supports_4_InSameBin1_0: p_ij[4,0] - p_ij[12,0] - s[4,12] <= 0
 12_Supports_4_InSameBin2_0: - p_ij[4,0] + p_ij[12,0] - s[4,12] <= 0
 Gamma0ForItems_4_inBins_0_NoCut_12: p_ij[4,0] + gamma[4] <= 1
 12_Supports_4_InSameBin1_1: p_ij[4,1] - p_ij[12,1] - s[4,12] <= 0
 12_Supports_4_InSameBin2_1: - p_ij[4,1] + p_ij[12,1] - s[4,12] <= 0
 Gamma0ForItems_4_inBins_1_NoCut_12: p_ij[4,1] + gamma[4] <= 1
 12_Supports_4_InSameBin1_2: p_ij[4,2] - p_ij[12,2] - s[4,12] <= 0
 12_Supports_4_InSameBin2_2: - p_ij[4,2] + p_ij[12,2] - s[4,12] <= 0
 Gamma0ForItems_4_inBins_2_NoCut_12: p_ij[4,2] + gamma[4] <= 98
 12_Supports_4_InSameBin1_3: p_ij[4,3] - p_ij[12,3] - s[4,12] <= 0
 12_Supports_4_InSameBin2_3: - p_ij[4,3] + p_ij[12,3] - s[4,12] <= 0
 Gamma0ForItems_4_inBins_3_NoCut_12: p_ij[4,3] + gamma[4] <= 98
 4_OnCutIn2_12: 10000 p_ij[4,2] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 4_OnCutIn3_12: 10000 p_ij[4,3] + 1.261904761904762 xi_l[4] + zi[4]
   + 10000 gamma[4] <= 20053
 GroundConstraint_5: zi[5] + 155 g[5] <= 155
 AbsZ15_0: - zi[5] + zi_hi[0] - k[5,0] <= 0
 AbsZ25_0: zi[5] - zi_hi[0] - k[5,0] <= 0
 0_Below_5: zi[5] - zi_hi[0] + k[5,0] + 310 m[5,0] <= 310
 0_NotBelow_5: - zi[5] + zi_hi[0] + k[5,0] - 310 m[5,0] <= 0
 HeightCompat_5_0: h[5,0] - k[5,0] <= 0
 NoHeightCompat_5_0: - 155 h[5,0] + k[5,0] <= 0
 Overlap_5_0: - x_p[0,5] - x_p[5,0] + o[5,0] <= 0
 NoOverlap_5_0: x_p[0,5] + x_p[5,0] - 2 o[5,0] <= 0
 3_Supports_5_0: - h[5,0] - o[5,0] - s[5,0] <= -1
 3_DoesntSupport_5_0: h[5,0] + o[5,0] + 2 s[5,0] <= 2
 beta1_5_0: - s[5,0] + beta1[5,0] <= 0
 beta2_5_0: - s[5,0] + beta2[5,0] <= 0
 sBeta1Flag_5_0: - s[5,0] + beta1[5,0] <= 0
 sBeta2Flag_5_0: - s[5,0] + beta2[5,0] <= 0
 Eta1Overlap_5_0: eta1[5,0] + beta1[5,0] <= 1
 Eta2Overlap_5_0: eta2[5,0] + beta2[5,0] <= 1
 Eta1Flag_5_0: xi_l[0] - xi_l[5] - 300 eta1[5,0] <= 0
 Eta2Flag_5_0: - xi_r[0] + xi_r[5] - 300 eta2[5,0] <= 0
 0_Supports_5_InSameBin1_0: - p_ij[0,0] + p_ij[5,0] - s[5,0] <= 0
 0_Supports_5_InSameBin2_0: p_ij[0,0] - p_ij[5,0] - s[5,0] <= 0
 Gamma0ForItems_5_inBins_0_NoCut_0: p_ij[5,0] + gamma[5] <= 1
 0_Supports_5_InSameBin1_1: - p_ij[0,1] + p_ij[5,1] - s[5,0] <= 0
 0_Supports_5_InSameBin2_1: p_ij[0,1] - p_ij[5,1] - s[5,0] <= 0
 Gamma0ForItems_5_inBins_1_NoCut_0: p_ij[5,1] + gamma[5] <= 1
 0_Supports_5_InSameBin1_2: - p_ij[0,2] + p_ij[5,2] - s[5,0] <= 0
 0_Supports_5_InSameBin2_2: p_ij[0,2] - p_ij[5,2] - s[5,0] <= 0
 Gamma0ForItems_5_inBins_2_NoCut_0: p_ij[5,2] + gamma[5] <= 98
 0_Supports_5_InSameBin1_3: - p_ij[0,3] + p_ij[5,3] - s[5,0] <= 0
 0_Supports_5_InSameBin2_3: p_ij[0,3] - p_ij[5,3] - s[5,0] <= 0
 Gamma0ForItems_5_inBins_3_NoCut_0: p_ij[5,3] + gamma[5] <= 98
 5_OnCutIn2_0: 10000 p_ij[5,2] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 5_OnCutIn3_0: 10000 p_ij[5,3] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 AbsZ15_1: - zi[5] + zi_hi[1] - k[5,1] <= 0
 AbsZ25_1: zi[5] - zi_hi[1] - k[5,1] <= 0
 1_Below_5: zi[5] - zi_hi[1] + k[5,1] + 310 m[5,1] <= 310
 1_NotBelow_5: - zi[5] + zi_hi[1] + k[5,1] - 310 m[5,1] <= 0
 HeightCompat_5_1: h[5,1] - k[5,1] <= 0
 NoHeightCompat_5_1: - 155 h[5,1] + k[5,1] <= 0
 Overlap_5_1: - x_p[1,5] - x_p[5,1] + o[5,1] <= 0
 NoOverlap_5_1: x_p[1,5] + x_p[5,1] - 2 o[5,1] <= 0
 3_Supports_5_1: - h[5,1] - o[5,1] - s[5,1] <= -1
 3_DoesntSupport_5_1: h[5,1] + o[5,1] + 2 s[5,1] <= 2
 beta1_5_1: - s[5,1] + beta1[5,1] <= 0
 beta2_5_1: - s[5,1] + beta2[5,1] <= 0
 sBeta1Flag_5_1: - s[5,1] + beta1[5,1] <= 0
 sBeta2Flag_5_1: - s[5,1] + beta2[5,1] <= 0
 Eta1Overlap_5_1: eta1[5,1] + beta1[5,1] <= 1
 Eta2Overlap_5_1: eta2[5,1] + beta2[5,1] <= 1
 Eta1Flag_5_1: xi_l[1] - xi_l[5] - 300 eta1[5,1] <= 0
 Eta2Flag_5_1: - xi_r[1] + xi_r[5] - 300 eta2[5,1] <= 0
 1_Supports_5_InSameBin1_0: - p_ij[1,0] + p_ij[5,0] - s[5,1] <= 0
 1_Supports_5_InSameBin2_0: p_ij[1,0] - p_ij[5,0] - s[5,1] <= 0
 Gamma0ForItems_5_inBins_0_NoCut_1: p_ij[5,0] + gamma[5] <= 1
 1_Supports_5_InSameBin1_1: - p_ij[1,1] + p_ij[5,1] - s[5,1] <= 0
 1_Supports_5_InSameBin2_1: p_ij[1,1] - p_ij[5,1] - s[5,1] <= 0
 Gamma0ForItems_5_inBins_1_NoCut_1: p_ij[5,1] + gamma[5] <= 1
 1_Supports_5_InSameBin1_2: - p_ij[1,2] + p_ij[5,2] - s[5,1] <= 0
 1_Supports_5_InSameBin2_2: p_ij[1,2] - p_ij[5,2] - s[5,1] <= 0
 Gamma0ForItems_5_inBins_2_NoCut_1: p_ij[5,2] + gamma[5] <= 98
 1_Supports_5_InSameBin1_3: - p_ij[1,3] + p_ij[5,3] - s[5,1] <= 0
 1_Supports_5_InSameBin2_3: p_ij[1,3] - p_ij[5,3] - s[5,1] <= 0
 Gamma0ForItems_5_inBins_3_NoCut_1: p_ij[5,3] + gamma[5] <= 98
 5_OnCutIn2_1: 10000 p_ij[5,2] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 5_OnCutIn3_1: 10000 p_ij[5,3] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 AbsZ15_2: - zi[5] + zi_hi[2] - k[5,2] <= 0
 AbsZ25_2: zi[5] - zi_hi[2] - k[5,2] <= 0
 2_Below_5: zi[5] - zi_hi[2] + k[5,2] + 310 m[5,2] <= 310
 2_NotBelow_5: - zi[5] + zi_hi[2] + k[5,2] - 310 m[5,2] <= 0
 HeightCompat_5_2: h[5,2] - k[5,2] <= 0
 NoHeightCompat_5_2: - 155 h[5,2] + k[5,2] <= 0
 Overlap_5_2: - x_p[2,5] - x_p[5,2] + o[5,2] <= 0
 NoOverlap_5_2: x_p[2,5] + x_p[5,2] - 2 o[5,2] <= 0
 3_Supports_5_2: - h[5,2] - o[5,2] - s[5,2] <= -1
 3_DoesntSupport_5_2: h[5,2] + o[5,2] + 2 s[5,2] <= 2
 beta1_5_2: - s[5,2] + beta1[5,2] <= 0
 beta2_5_2: - s[5,2] + beta2[5,2] <= 0
 sBeta1Flag_5_2: - s[5,2] + beta1[5,2] <= 0
 sBeta2Flag_5_2: - s[5,2] + beta2[5,2] <= 0
 Eta1Overlap_5_2: eta1[5,2] + beta1[5,2] <= 1
 Eta2Overlap_5_2: eta2[5,2] + beta2[5,2] <= 1
 Eta1Flag_5_2: xi_l[2] - xi_l[5] - 300 eta1[5,2] <= 0
 Eta2Flag_5_2: - xi_r[2] + xi_r[5] - 300 eta2[5,2] <= 0
 2_Supports_5_InSameBin1_0: - p_ij[2,0] + p_ij[5,0] - s[5,2] <= 0
 2_Supports_5_InSameBin2_0: p_ij[2,0] - p_ij[5,0] - s[5,2] <= 0
 Gamma0ForItems_5_inBins_0_NoCut_2: p_ij[5,0] + gamma[5] <= 1
 2_Supports_5_InSameBin1_1: - p_ij[2,1] + p_ij[5,1] - s[5,2] <= 0
 2_Supports_5_InSameBin2_1: p_ij[2,1] - p_ij[5,1] - s[5,2] <= 0
 Gamma0ForItems_5_inBins_1_NoCut_2: p_ij[5,1] + gamma[5] <= 1
 2_Supports_5_InSameBin1_2: - p_ij[2,2] + p_ij[5,2] - s[5,2] <= 0
 2_Supports_5_InSameBin2_2: p_ij[2,2] - p_ij[5,2] - s[5,2] <= 0
 Gamma0ForItems_5_inBins_2_NoCut_2: p_ij[5,2] + gamma[5] <= 98
 2_Supports_5_InSameBin1_3: - p_ij[2,3] + p_ij[5,3] - s[5,2] <= 0
 2_Supports_5_InSameBin2_3: p_ij[2,3] - p_ij[5,3] - s[5,2] <= 0
 Gamma0ForItems_5_inBins_3_NoCut_2: p_ij[5,3] + gamma[5] <= 98
 5_OnCutIn2_2: 10000 p_ij[5,2] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 5_OnCutIn3_2: 10000 p_ij[5,3] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 AbsZ15_3: - zi[5] + zi_hi[3] - k[5,3] <= 0
 AbsZ25_3: zi[5] - zi_hi[3] - k[5,3] <= 0
 3_Below_5: zi[5] - zi_hi[3] + k[5,3] + 310 m[5,3] <= 310
 3_NotBelow_5: - zi[5] + zi_hi[3] + k[5,3] - 310 m[5,3] <= 0
 HeightCompat_5_3: h[5,3] - k[5,3] <= 0
 NoHeightCompat_5_3: - 155 h[5,3] + k[5,3] <= 0
 Overlap_5_3: - x_p[3,5] - x_p[5,3] + o[5,3] <= 0
 NoOverlap_5_3: x_p[3,5] + x_p[5,3] - 2 o[5,3] <= 0
 3_Supports_5_3: - h[5,3] - o[5,3] - s[5,3] <= -1
 3_DoesntSupport_5_3: h[5,3] + o[5,3] + 2 s[5,3] <= 2
 beta1_5_3: - s[5,3] + beta1[5,3] <= 0
 beta2_5_3: - s[5,3] + beta2[5,3] <= 0
 sBeta1Flag_5_3: - s[5,3] + beta1[5,3] <= 0
 sBeta2Flag_5_3: - s[5,3] + beta2[5,3] <= 0
 Eta1Overlap_5_3: eta1[5,3] + beta1[5,3] <= 1
 Eta2Overlap_5_3: eta2[5,3] + beta2[5,3] <= 1
 Eta1Flag_5_3: xi_l[3] - xi_l[5] - 300 eta1[5,3] <= 0
 Eta2Flag_5_3: - xi_r[3] + xi_r[5] - 300 eta2[5,3] <= 0
 3_Supports_5_InSameBin1_0: - p_ij[3,0] + p_ij[5,0] - s[5,3] <= 0
 3_Supports_5_InSameBin2_0: p_ij[3,0] - p_ij[5,0] - s[5,3] <= 0
 Gamma0ForItems_5_inBins_0_NoCut_3: p_ij[5,0] + gamma[5] <= 1
 3_Supports_5_InSameBin1_1: - p_ij[3,1] + p_ij[5,1] - s[5,3] <= 0
 3_Supports_5_InSameBin2_1: p_ij[3,1] - p_ij[5,1] - s[5,3] <= 0
 Gamma0ForItems_5_inBins_1_NoCut_3: p_ij[5,1] + gamma[5] <= 1
 3_Supports_5_InSameBin1_2: - p_ij[3,2] + p_ij[5,2] - s[5,3] <= 0
 3_Supports_5_InSameBin2_2: p_ij[3,2] - p_ij[5,2] - s[5,3] <= 0
 Gamma0ForItems_5_inBins_2_NoCut_3: p_ij[5,2] + gamma[5] <= 98
 3_Supports_5_InSameBin1_3: - p_ij[3,3] + p_ij[5,3] - s[5,3] <= 0
 3_Supports_5_InSameBin2_3: p_ij[3,3] - p_ij[5,3] - s[5,3] <= 0
 Gamma0ForItems_5_inBins_3_NoCut_3: p_ij[5,3] + gamma[5] <= 98
 5_OnCutIn2_3: 10000 p_ij[5,2] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 5_OnCutIn3_3: 10000 p_ij[5,3] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 AbsZ15_4: - zi[5] + zi_hi[4] - k[5,4] <= 0
 AbsZ25_4: zi[5] - zi_hi[4] - k[5,4] <= 0
 4_Below_5: zi[5] - zi_hi[4] + k[5,4] + 310 m[5,4] <= 310
 4_NotBelow_5: - zi[5] + zi_hi[4] + k[5,4] - 310 m[5,4] <= 0
 HeightCompat_5_4: h[5,4] - k[5,4] <= 0
 NoHeightCompat_5_4: - 155 h[5,4] + k[5,4] <= 0
 Overlap_5_4: - x_p[4,5] - x_p[5,4] + o[5,4] <= 0
 NoOverlap_5_4: x_p[4,5] + x_p[5,4] - 2 o[5,4] <= 0
 3_Supports_5_4: - h[5,4] - o[5,4] - s[5,4] <= -1
 3_DoesntSupport_5_4: h[5,4] + o[5,4] + 2 s[5,4] <= 2
 beta1_5_4: - s[5,4] + beta1[5,4] <= 0
 beta2_5_4: - s[5,4] + beta2[5,4] <= 0
 sBeta1Flag_5_4: - s[5,4] + beta1[5,4] <= 0
 sBeta2Flag_5_4: - s[5,4] + beta2[5,4] <= 0
 Eta1Overlap_5_4: eta1[5,4] + beta1[5,4] <= 1
 Eta2Overlap_5_4: eta2[5,4] + beta2[5,4] <= 1
 Eta1Flag_5_4: xi_l[4] - xi_l[5] - 300 eta1[5,4] <= 0
 Eta2Flag_5_4: - xi_r[4] + xi_r[5] - 300 eta2[5,4] <= 0
 4_Supports_5_InSameBin1_0: - p_ij[4,0] + p_ij[5,0] - s[5,4] <= 0
 4_Supports_5_InSameBin2_0: p_ij[4,0] - p_ij[5,0] - s[5,4] <= 0
 Gamma0ForItems_5_inBins_0_NoCut_4: p_ij[5,0] + gamma[5] <= 1
 4_Supports_5_InSameBin1_1: - p_ij[4,1] + p_ij[5,1] - s[5,4] <= 0
 4_Supports_5_InSameBin2_1: p_ij[4,1] - p_ij[5,1] - s[5,4] <= 0
 Gamma0ForItems_5_inBins_1_NoCut_4: p_ij[5,1] + gamma[5] <= 1
 4_Supports_5_InSameBin1_2: - p_ij[4,2] + p_ij[5,2] - s[5,4] <= 0
 4_Supports_5_InSameBin2_2: p_ij[4,2] - p_ij[5,2] - s[5,4] <= 0
 Gamma0ForItems_5_inBins_2_NoCut_4: p_ij[5,2] + gamma[5] <= 98
 4_Supports_5_InSameBin1_3: - p_ij[4,3] + p_ij[5,3] - s[5,4] <= 0
 4_Supports_5_InSameBin2_3: p_ij[4,3] - p_ij[5,3] - s[5,4] <= 0
 Gamma0ForItems_5_inBins_3_NoCut_4: p_ij[5,3] + gamma[5] <= 98
 5_OnCutIn2_4: 10000 p_ij[5,2] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 5_OnCutIn3_4: 10000 p_ij[5,3] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 5_OnCutIn2_5: 10000 p_ij[5,2] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 5_OnCutIn3_5: 10000 p_ij[5,3] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 AbsZ15_6: - zi[5] + zi_hi[6] - k[5,6] <= 0
 AbsZ25_6: zi[5] - zi_hi[6] - k[5,6] <= 0
 6_Below_5: zi[5] - zi_hi[6] + k[5,6] + 310 m[5,6] <= 310
 6_NotBelow_5: - zi[5] + zi_hi[6] + k[5,6] - 310 m[5,6] <= 0
 HeightCompat_5_6: h[5,6] - k[5,6] <= 0
 NoHeightCompat_5_6: - 155 h[5,6] + k[5,6] <= 0
 Overlap_5_6: - x_p[5,6] - x_p[6,5] + o[5,6] <= 0
 NoOverlap_5_6: x_p[5,6] + x_p[6,5] - 2 o[5,6] <= 0
 3_Supports_5_6: - h[5,6] - o[5,6] - s[5,6] <= -1
 3_DoesntSupport_5_6: h[5,6] + o[5,6] + 2 s[5,6] <= 2
 beta1_5_6: - s[5,6] + beta1[5,6] <= 0
 beta2_5_6: - s[5,6] + beta2[5,6] <= 0
 sBeta1Flag_5_6: - s[5,6] + beta1[5,6] <= 0
 sBeta2Flag_5_6: - s[5,6] + beta2[5,6] <= 0
 Eta1Overlap_5_6: eta1[5,6] + beta1[5,6] <= 1
 Eta2Overlap_5_6: eta2[5,6] + beta2[5,6] <= 1
 Eta1Flag_5_6: - xi_l[5] + xi_l[6] - 300 eta1[5,6] <= 0
 Eta2Flag_5_6: xi_r[5] - xi_r[6] - 300 eta2[5,6] <= 0
 6_Supports_5_InSameBin1_0: p_ij[5,0] - p_ij[6,0] - s[5,6] <= 0
 6_Supports_5_InSameBin2_0: - p_ij[5,0] + p_ij[6,0] - s[5,6] <= 0
 Gamma0ForItems_5_inBins_0_NoCut_6: p_ij[5,0] + gamma[5] <= 1
 6_Supports_5_InSameBin1_1: p_ij[5,1] - p_ij[6,1] - s[5,6] <= 0
 6_Supports_5_InSameBin2_1: - p_ij[5,1] + p_ij[6,1] - s[5,6] <= 0
 Gamma0ForItems_5_inBins_1_NoCut_6: p_ij[5,1] + gamma[5] <= 1
 6_Supports_5_InSameBin1_2: p_ij[5,2] - p_ij[6,2] - s[5,6] <= 0
 6_Supports_5_InSameBin2_2: - p_ij[5,2] + p_ij[6,2] - s[5,6] <= 0
 Gamma0ForItems_5_inBins_2_NoCut_6: p_ij[5,2] + gamma[5] <= 98
 6_Supports_5_InSameBin1_3: p_ij[5,3] - p_ij[6,3] - s[5,6] <= 0
 6_Supports_5_InSameBin2_3: - p_ij[5,3] + p_ij[6,3] - s[5,6] <= 0
 Gamma0ForItems_5_inBins_3_NoCut_6: p_ij[5,3] + gamma[5] <= 98
 5_OnCutIn2_6: 10000 p_ij[5,2] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 5_OnCutIn3_6: 10000 p_ij[5,3] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 AbsZ15_7: - zi[5] + zi_hi[7] - k[5,7] <= 0
 AbsZ25_7: zi[5] - zi_hi[7] - k[5,7] <= 0
 7_Below_5: zi[5] - zi_hi[7] + k[5,7] + 310 m[5,7] <= 310
 7_NotBelow_5: - zi[5] + zi_hi[7] + k[5,7] - 310 m[5,7] <= 0
 HeightCompat_5_7: h[5,7] - k[5,7] <= 0
 NoHeightCompat_5_7: - 155 h[5,7] + k[5,7] <= 0
 Overlap_5_7: - x_p[5,7] - x_p[7,5] + o[5,7] <= 0
 NoOverlap_5_7: x_p[5,7] + x_p[7,5] - 2 o[5,7] <= 0
 3_Supports_5_7: - h[5,7] - o[5,7] - s[5,7] <= -1
 3_DoesntSupport_5_7: h[5,7] + o[5,7] + 2 s[5,7] <= 2
 beta1_5_7: - s[5,7] + beta1[5,7] <= 0
 beta2_5_7: - s[5,7] + beta2[5,7] <= 0
 sBeta1Flag_5_7: - s[5,7] + beta1[5,7] <= 0
 sBeta2Flag_5_7: - s[5,7] + beta2[5,7] <= 0
 Eta1Overlap_5_7: eta1[5,7] + beta1[5,7] <= 1
 Eta2Overlap_5_7: eta2[5,7] + beta2[5,7] <= 1
 Eta1Flag_5_7: - xi_l[5] + xi_l[7] - 300 eta1[5,7] <= 0
 Eta2Flag_5_7: xi_r[5] - xi_r[7] - 300 eta2[5,7] <= 0
 7_Supports_5_InSameBin1_0: p_ij[5,0] - p_ij[7,0] - s[5,7] <= 0
 7_Supports_5_InSameBin2_0: - p_ij[5,0] + p_ij[7,0] - s[5,7] <= 0
 Gamma0ForItems_5_inBins_0_NoCut_7: p_ij[5,0] + gamma[5] <= 1
 7_Supports_5_InSameBin1_1: p_ij[5,1] - p_ij[7,1] - s[5,7] <= 0
 7_Supports_5_InSameBin2_1: - p_ij[5,1] + p_ij[7,1] - s[5,7] <= 0
 Gamma0ForItems_5_inBins_1_NoCut_7: p_ij[5,1] + gamma[5] <= 1
 7_Supports_5_InSameBin1_2: p_ij[5,2] - p_ij[7,2] - s[5,7] <= 0
 7_Supports_5_InSameBin2_2: - p_ij[5,2] + p_ij[7,2] - s[5,7] <= 0
 Gamma0ForItems_5_inBins_2_NoCut_7: p_ij[5,2] + gamma[5] <= 98
 7_Supports_5_InSameBin1_3: p_ij[5,3] - p_ij[7,3] - s[5,7] <= 0
 7_Supports_5_InSameBin2_3: - p_ij[5,3] + p_ij[7,3] - s[5,7] <= 0
 Gamma0ForItems_5_inBins_3_NoCut_7: p_ij[5,3] + gamma[5] <= 98
 5_OnCutIn2_7: 10000 p_ij[5,2] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 5_OnCutIn3_7: 10000 p_ij[5,3] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 AbsZ15_8: - zi[5] + zi_hi[8] - k[5,8] <= 0
 AbsZ25_8: zi[5] - zi_hi[8] - k[5,8] <= 0
 8_Below_5: zi[5] - zi_hi[8] + k[5,8] + 310 m[5,8] <= 310
 8_NotBelow_5: - zi[5] + zi_hi[8] + k[5,8] - 310 m[5,8] <= 0
 HeightCompat_5_8: h[5,8] - k[5,8] <= 0
 NoHeightCompat_5_8: - 155 h[5,8] + k[5,8] <= 0
 Overlap_5_8: - x_p[5,8] - x_p[8,5] + o[5,8] <= 0
 NoOverlap_5_8: x_p[5,8] + x_p[8,5] - 2 o[5,8] <= 0
 3_Supports_5_8: - h[5,8] - o[5,8] - s[5,8] <= -1
 3_DoesntSupport_5_8: h[5,8] + o[5,8] + 2 s[5,8] <= 2
 beta1_5_8: - s[5,8] + beta1[5,8] <= 0
 beta2_5_8: - s[5,8] + beta2[5,8] <= 0
 sBeta1Flag_5_8: - s[5,8] + beta1[5,8] <= 0
 sBeta2Flag_5_8: - s[5,8] + beta2[5,8] <= 0
 Eta1Overlap_5_8: eta1[5,8] + beta1[5,8] <= 1
 Eta2Overlap_5_8: eta2[5,8] + beta2[5,8] <= 1
 Eta1Flag_5_8: - xi_l[5] + xi_l[8] - 300 eta1[5,8] <= 0
 Eta2Flag_5_8: xi_r[5] - xi_r[8] - 300 eta2[5,8] <= 0
 8_Supports_5_InSameBin1_0: p_ij[5,0] - p_ij[8,0] - s[5,8] <= 0
 8_Supports_5_InSameBin2_0: - p_ij[5,0] + p_ij[8,0] - s[5,8] <= 0
 Gamma0ForItems_5_inBins_0_NoCut_8: p_ij[5,0] + gamma[5] <= 1
 8_Supports_5_InSameBin1_1: p_ij[5,1] - p_ij[8,1] - s[5,8] <= 0
 8_Supports_5_InSameBin2_1: - p_ij[5,1] + p_ij[8,1] - s[5,8] <= 0
 Gamma0ForItems_5_inBins_1_NoCut_8: p_ij[5,1] + gamma[5] <= 1
 8_Supports_5_InSameBin1_2: p_ij[5,2] - p_ij[8,2] - s[5,8] <= 0
 8_Supports_5_InSameBin2_2: - p_ij[5,2] + p_ij[8,2] - s[5,8] <= 0
 Gamma0ForItems_5_inBins_2_NoCut_8: p_ij[5,2] + gamma[5] <= 98
 8_Supports_5_InSameBin1_3: p_ij[5,3] - p_ij[8,3] - s[5,8] <= 0
 8_Supports_5_InSameBin2_3: - p_ij[5,3] + p_ij[8,3] - s[5,8] <= 0
 Gamma0ForItems_5_inBins_3_NoCut_8: p_ij[5,3] + gamma[5] <= 98
 5_OnCutIn2_8: 10000 p_ij[5,2] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 5_OnCutIn3_8: 10000 p_ij[5,3] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 AbsZ15_9: - zi[5] + zi_hi[9] - k[5,9] <= 0
 AbsZ25_9: zi[5] - zi_hi[9] - k[5,9] <= 0
 9_Below_5: zi[5] - zi_hi[9] + k[5,9] + 310 m[5,9] <= 310
 9_NotBelow_5: - zi[5] + zi_hi[9] + k[5,9] - 310 m[5,9] <= 0
 HeightCompat_5_9: h[5,9] - k[5,9] <= 0
 NoHeightCompat_5_9: - 155 h[5,9] + k[5,9] <= 0
 Overlap_5_9: - x_p[5,9] - x_p[9,5] + o[5,9] <= 0
 NoOverlap_5_9: x_p[5,9] + x_p[9,5] - 2 o[5,9] <= 0
 3_Supports_5_9: - h[5,9] - o[5,9] - s[5,9] <= -1
 3_DoesntSupport_5_9: h[5,9] + o[5,9] + 2 s[5,9] <= 2
 beta1_5_9: - s[5,9] + beta1[5,9] <= 0
 beta2_5_9: - s[5,9] + beta2[5,9] <= 0
 sBeta1Flag_5_9: - s[5,9] + beta1[5,9] <= 0
 sBeta2Flag_5_9: - s[5,9] + beta2[5,9] <= 0
 Eta1Overlap_5_9: eta1[5,9] + beta1[5,9] <= 1
 Eta2Overlap_5_9: eta2[5,9] + beta2[5,9] <= 1
 Eta1Flag_5_9: - xi_l[5] + xi_l[9] - 300 eta1[5,9] <= 0
 Eta2Flag_5_9: xi_r[5] - xi_r[9] - 300 eta2[5,9] <= 0
 9_Supports_5_InSameBin1_0: p_ij[5,0] - p_ij[9,0] - s[5,9] <= 0
 9_Supports_5_InSameBin2_0: - p_ij[5,0] + p_ij[9,0] - s[5,9] <= 0
 Gamma0ForItems_5_inBins_0_NoCut_9: p_ij[5,0] + gamma[5] <= 1
 9_Supports_5_InSameBin1_1: p_ij[5,1] - p_ij[9,1] - s[5,9] <= 0
 9_Supports_5_InSameBin2_1: - p_ij[5,1] + p_ij[9,1] - s[5,9] <= 0
 Gamma0ForItems_5_inBins_1_NoCut_9: p_ij[5,1] + gamma[5] <= 1
 9_Supports_5_InSameBin1_2: p_ij[5,2] - p_ij[9,2] - s[5,9] <= 0
 9_Supports_5_InSameBin2_2: - p_ij[5,2] + p_ij[9,2] - s[5,9] <= 0
 Gamma0ForItems_5_inBins_2_NoCut_9: p_ij[5,2] + gamma[5] <= 98
 9_Supports_5_InSameBin1_3: p_ij[5,3] - p_ij[9,3] - s[5,9] <= 0
 9_Supports_5_InSameBin2_3: - p_ij[5,3] + p_ij[9,3] - s[5,9] <= 0
 Gamma0ForItems_5_inBins_3_NoCut_9: p_ij[5,3] + gamma[5] <= 98
 5_OnCutIn2_9: 10000 p_ij[5,2] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 5_OnCutIn3_9: 10000 p_ij[5,3] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 AbsZ15_10: - zi[5] + zi_hi[10] - k[5,10] <= 0
 AbsZ25_10: zi[5] - zi_hi[10] - k[5,10] <= 0
 10_Below_5: zi[5] - zi_hi[10] + k[5,10] + 310 m[5,10] <= 310
 10_NotBelow_5: - zi[5] + zi_hi[10] + k[5,10] - 310 m[5,10] <= 0
 HeightCompat_5_10: h[5,10] - k[5,10] <= 0
 NoHeightCompat_5_10: - 155 h[5,10] + k[5,10] <= 0
 Overlap_5_10: - x_p[5,10] - x_p[10,5] + o[5,10] <= 0
 NoOverlap_5_10: x_p[5,10] + x_p[10,5] - 2 o[5,10] <= 0
 3_Supports_5_10: - h[5,10] - o[5,10] - s[5,10] <= -1
 3_DoesntSupport_5_10: h[5,10] + o[5,10] + 2 s[5,10] <= 2
 beta1_5_10: - s[5,10] + beta1[5,10] <= 0
 beta2_5_10: - s[5,10] + beta2[5,10] <= 0
 sBeta1Flag_5_10: - s[5,10] + beta1[5,10] <= 0
 sBeta2Flag_5_10: - s[5,10] + beta2[5,10] <= 0
 Eta1Overlap_5_10: eta1[5,10] + beta1[5,10] <= 1
 Eta2Overlap_5_10: eta2[5,10] + beta2[5,10] <= 1
 Eta1Flag_5_10: - xi_l[5] + xi_l[10] - 300 eta1[5,10] <= 0
 Eta2Flag_5_10: xi_r[5] - xi_r[10] - 300 eta2[5,10] <= 0
 10_Supports_5_InSameBin1_0: p_ij[5,0] - p_ij[10,0] - s[5,10] <= 0
 10_Supports_5_InSameBin2_0: - p_ij[5,0] + p_ij[10,0] - s[5,10] <= 0
 Gamma0ForItems_5_inBins_0_NoCut_10: p_ij[5,0] + gamma[5] <= 1
 10_Supports_5_InSameBin1_1: p_ij[5,1] - p_ij[10,1] - s[5,10] <= 0
 10_Supports_5_InSameBin2_1: - p_ij[5,1] + p_ij[10,1] - s[5,10] <= 0
 Gamma0ForItems_5_inBins_1_NoCut_10: p_ij[5,1] + gamma[5] <= 1
 10_Supports_5_InSameBin1_2: p_ij[5,2] - p_ij[10,2] - s[5,10] <= 0
 10_Supports_5_InSameBin2_2: - p_ij[5,2] + p_ij[10,2] - s[5,10] <= 0
 Gamma0ForItems_5_inBins_2_NoCut_10: p_ij[5,2] + gamma[5] <= 98
 10_Supports_5_InSameBin1_3: p_ij[5,3] - p_ij[10,3] - s[5,10] <= 0
 10_Supports_5_InSameBin2_3: - p_ij[5,3] + p_ij[10,3] - s[5,10] <= 0
 Gamma0ForItems_5_inBins_3_NoCut_10: p_ij[5,3] + gamma[5] <= 98
 5_OnCutIn2_10: 10000 p_ij[5,2] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 5_OnCutIn3_10: 10000 p_ij[5,3] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 AbsZ15_11: - zi[5] + zi_hi[11] - k[5,11] <= 0
 AbsZ25_11: zi[5] - zi_hi[11] - k[5,11] <= 0
 11_Below_5: zi[5] - zi_hi[11] + k[5,11] + 310 m[5,11] <= 310
 11_NotBelow_5: - zi[5] + zi_hi[11] + k[5,11] - 310 m[5,11] <= 0
 HeightCompat_5_11: h[5,11] - k[5,11] <= 0
 NoHeightCompat_5_11: - 155 h[5,11] + k[5,11] <= 0
 Overlap_5_11: - x_p[5,11] - x_p[11,5] + o[5,11] <= 0
 NoOverlap_5_11: x_p[5,11] + x_p[11,5] - 2 o[5,11] <= 0
 3_Supports_5_11: - h[5,11] - o[5,11] - s[5,11] <= -1
 3_DoesntSupport_5_11: h[5,11] + o[5,11] + 2 s[5,11] <= 2
 beta1_5_11: - s[5,11] + beta1[5,11] <= 0
 beta2_5_11: - s[5,11] + beta2[5,11] <= 0
 sBeta1Flag_5_11: - s[5,11] + beta1[5,11] <= 0
 sBeta2Flag_5_11: - s[5,11] + beta2[5,11] <= 0
 Eta1Overlap_5_11: eta1[5,11] + beta1[5,11] <= 1
 Eta2Overlap_5_11: eta2[5,11] + beta2[5,11] <= 1
 Eta1Flag_5_11: - xi_l[5] + xi_l[11] - 300 eta1[5,11] <= 0
 Eta2Flag_5_11: xi_r[5] - xi_r[11] - 300 eta2[5,11] <= 0
 11_Supports_5_InSameBin1_0: p_ij[5,0] - p_ij[11,0] - s[5,11] <= 0
 11_Supports_5_InSameBin2_0: - p_ij[5,0] + p_ij[11,0] - s[5,11] <= 0
 Gamma0ForItems_5_inBins_0_NoCut_11: p_ij[5,0] + gamma[5] <= 1
 11_Supports_5_InSameBin1_1: p_ij[5,1] - p_ij[11,1] - s[5,11] <= 0
 11_Supports_5_InSameBin2_1: - p_ij[5,1] + p_ij[11,1] - s[5,11] <= 0
 Gamma0ForItems_5_inBins_1_NoCut_11: p_ij[5,1] + gamma[5] <= 1
 11_Supports_5_InSameBin1_2: p_ij[5,2] - p_ij[11,2] - s[5,11] <= 0
 11_Supports_5_InSameBin2_2: - p_ij[5,2] + p_ij[11,2] - s[5,11] <= 0
 Gamma0ForItems_5_inBins_2_NoCut_11: p_ij[5,2] + gamma[5] <= 98
 11_Supports_5_InSameBin1_3: p_ij[5,3] - p_ij[11,3] - s[5,11] <= 0
 11_Supports_5_InSameBin2_3: - p_ij[5,3] + p_ij[11,3] - s[5,11] <= 0
 Gamma0ForItems_5_inBins_3_NoCut_11: p_ij[5,3] + gamma[5] <= 98
 5_OnCutIn2_11: 10000 p_ij[5,2] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 5_OnCutIn3_11: 10000 p_ij[5,3] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 AbsZ15_12: - zi[5] + zi_hi[12] - k[5,12] <= 0
 AbsZ25_12: zi[5] - zi_hi[12] - k[5,12] <= 0
 12_Below_5: zi[5] - zi_hi[12] + k[5,12] + 310 m[5,12] <= 310
 12_NotBelow_5: - zi[5] + zi_hi[12] + k[5,12] - 310 m[5,12] <= 0
 HeightCompat_5_12: h[5,12] - k[5,12] <= 0
 NoHeightCompat_5_12: - 155 h[5,12] + k[5,12] <= 0
 Overlap_5_12: - x_p[5,12] - x_p[12,5] + o[5,12] <= 0
 NoOverlap_5_12: x_p[5,12] + x_p[12,5] - 2 o[5,12] <= 0
 3_Supports_5_12: - h[5,12] - o[5,12] - s[5,12] <= -1
 3_DoesntSupport_5_12: h[5,12] + o[5,12] + 2 s[5,12] <= 2
 beta1_5_12: - s[5,12] + beta1[5,12] <= 0
 beta2_5_12: - s[5,12] + beta2[5,12] <= 0
 sBeta1Flag_5_12: - s[5,12] + beta1[5,12] <= 0
 sBeta2Flag_5_12: - s[5,12] + beta2[5,12] <= 0
 Eta1Overlap_5_12: eta1[5,12] + beta1[5,12] <= 1
 Eta2Overlap_5_12: eta2[5,12] + beta2[5,12] <= 1
 Eta1Flag_5_12: - xi_l[5] + xi_l[12] - 300 eta1[5,12] <= 0
 Eta2Flag_5_12: xi_r[5] - xi_r[12] - 300 eta2[5,12] <= 0
 12_Supports_5_InSameBin1_0: p_ij[5,0] - p_ij[12,0] - s[5,12] <= 0
 12_Supports_5_InSameBin2_0: - p_ij[5,0] + p_ij[12,0] - s[5,12] <= 0
 Gamma0ForItems_5_inBins_0_NoCut_12: p_ij[5,0] + gamma[5] <= 1
 12_Supports_5_InSameBin1_1: p_ij[5,1] - p_ij[12,1] - s[5,12] <= 0
 12_Supports_5_InSameBin2_1: - p_ij[5,1] + p_ij[12,1] - s[5,12] <= 0
 Gamma0ForItems_5_inBins_1_NoCut_12: p_ij[5,1] + gamma[5] <= 1
 12_Supports_5_InSameBin1_2: p_ij[5,2] - p_ij[12,2] - s[5,12] <= 0
 12_Supports_5_InSameBin2_2: - p_ij[5,2] + p_ij[12,2] - s[5,12] <= 0
 Gamma0ForItems_5_inBins_2_NoCut_12: p_ij[5,2] + gamma[5] <= 98
 12_Supports_5_InSameBin1_3: p_ij[5,3] - p_ij[12,3] - s[5,12] <= 0
 12_Supports_5_InSameBin2_3: - p_ij[5,3] + p_ij[12,3] - s[5,12] <= 0
 Gamma0ForItems_5_inBins_3_NoCut_12: p_ij[5,3] + gamma[5] <= 98
 5_OnCutIn2_12: 10000 p_ij[5,2] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 5_OnCutIn3_12: 10000 p_ij[5,3] + 1.261904761904762 xi_l[5] + zi[5]
   + 10000 gamma[5] <= 20053
 GroundConstraint_6: zi[6] + 155 g[6] <= 155
 AbsZ16_0: - zi[6] + zi_hi[0] - k[6,0] <= 0
 AbsZ26_0: zi[6] - zi_hi[0] - k[6,0] <= 0
 0_Below_6: zi[6] - zi_hi[0] + k[6,0] + 310 m[6,0] <= 310
 0_NotBelow_6: - zi[6] + zi_hi[0] + k[6,0] - 310 m[6,0] <= 0
 HeightCompat_6_0: h[6,0] - k[6,0] <= 0
 NoHeightCompat_6_0: - 155 h[6,0] + k[6,0] <= 0
 Overlap_6_0: - x_p[0,6] - x_p[6,0] + o[6,0] <= 0
 NoOverlap_6_0: x_p[0,6] + x_p[6,0] - 2 o[6,0] <= 0
 3_Supports_6_0: - h[6,0] - o[6,0] - s[6,0] <= -1
 3_DoesntSupport_6_0: h[6,0] + o[6,0] + 2 s[6,0] <= 2
 beta1_6_0: - s[6,0] + beta1[6,0] <= 0
 beta2_6_0: - s[6,0] + beta2[6,0] <= 0
 sBeta1Flag_6_0: - s[6,0] + beta1[6,0] <= 0
 sBeta2Flag_6_0: - s[6,0] + beta2[6,0] <= 0
 Eta1Overlap_6_0: eta1[6,0] + beta1[6,0] <= 1
 Eta2Overlap_6_0: eta2[6,0] + beta2[6,0] <= 1
 Eta1Flag_6_0: xi_l[0] - xi_l[6] - 300 eta1[6,0] <= 0
 Eta2Flag_6_0: - xi_r[0] + xi_r[6] - 300 eta2[6,0] <= 0
 0_Supports_6_InSameBin1_0: - p_ij[0,0] + p_ij[6,0] - s[6,0] <= 0
 0_Supports_6_InSameBin2_0: p_ij[0,0] - p_ij[6,0] - s[6,0] <= 0
 Gamma0ForItems_6_inBins_0_NoCut_0: p_ij[6,0] + gamma[6] <= 1
 0_Supports_6_InSameBin1_1: - p_ij[0,1] + p_ij[6,1] - s[6,0] <= 0
 0_Supports_6_InSameBin2_1: p_ij[0,1] - p_ij[6,1] - s[6,0] <= 0
 Gamma0ForItems_6_inBins_1_NoCut_0: p_ij[6,1] + gamma[6] <= 1
 0_Supports_6_InSameBin1_2: - p_ij[0,2] + p_ij[6,2] - s[6,0] <= 0
 0_Supports_6_InSameBin2_2: p_ij[0,2] - p_ij[6,2] - s[6,0] <= 0
 Gamma0ForItems_6_inBins_2_NoCut_0: p_ij[6,2] + gamma[6] <= 98
 0_Supports_6_InSameBin1_3: - p_ij[0,3] + p_ij[6,3] - s[6,0] <= 0
 0_Supports_6_InSameBin2_3: p_ij[0,3] - p_ij[6,3] - s[6,0] <= 0
 Gamma0ForItems_6_inBins_3_NoCut_0: p_ij[6,3] + gamma[6] <= 98
 6_OnCutIn2_0: 10000 p_ij[6,2] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 6_OnCutIn3_0: 10000 p_ij[6,3] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 AbsZ16_1: - zi[6] + zi_hi[1] - k[6,1] <= 0
 AbsZ26_1: zi[6] - zi_hi[1] - k[6,1] <= 0
 1_Below_6: zi[6] - zi_hi[1] + k[6,1] + 310 m[6,1] <= 310
 1_NotBelow_6: - zi[6] + zi_hi[1] + k[6,1] - 310 m[6,1] <= 0
 HeightCompat_6_1: h[6,1] - k[6,1] <= 0
 NoHeightCompat_6_1: - 155 h[6,1] + k[6,1] <= 0
 Overlap_6_1: - x_p[1,6] - x_p[6,1] + o[6,1] <= 0
 NoOverlap_6_1: x_p[1,6] + x_p[6,1] - 2 o[6,1] <= 0
 3_Supports_6_1: - h[6,1] - o[6,1] - s[6,1] <= -1
 3_DoesntSupport_6_1: h[6,1] + o[6,1] + 2 s[6,1] <= 2
 beta1_6_1: - s[6,1] + beta1[6,1] <= 0
 beta2_6_1: - s[6,1] + beta2[6,1] <= 0
 sBeta1Flag_6_1: - s[6,1] + beta1[6,1] <= 0
 sBeta2Flag_6_1: - s[6,1] + beta2[6,1] <= 0
 Eta1Overlap_6_1: eta1[6,1] + beta1[6,1] <= 1
 Eta2Overlap_6_1: eta2[6,1] + beta2[6,1] <= 1
 Eta1Flag_6_1: xi_l[1] - xi_l[6] - 300 eta1[6,1] <= 0
 Eta2Flag_6_1: - xi_r[1] + xi_r[6] - 300 eta2[6,1] <= 0
 1_Supports_6_InSameBin1_0: - p_ij[1,0] + p_ij[6,0] - s[6,1] <= 0
 1_Supports_6_InSameBin2_0: p_ij[1,0] - p_ij[6,0] - s[6,1] <= 0
 Gamma0ForItems_6_inBins_0_NoCut_1: p_ij[6,0] + gamma[6] <= 1
 1_Supports_6_InSameBin1_1: - p_ij[1,1] + p_ij[6,1] - s[6,1] <= 0
 1_Supports_6_InSameBin2_1: p_ij[1,1] - p_ij[6,1] - s[6,1] <= 0
 Gamma0ForItems_6_inBins_1_NoCut_1: p_ij[6,1] + gamma[6] <= 1
 1_Supports_6_InSameBin1_2: - p_ij[1,2] + p_ij[6,2] - s[6,1] <= 0
 1_Supports_6_InSameBin2_2: p_ij[1,2] - p_ij[6,2] - s[6,1] <= 0
 Gamma0ForItems_6_inBins_2_NoCut_1: p_ij[6,2] + gamma[6] <= 98
 1_Supports_6_InSameBin1_3: - p_ij[1,3] + p_ij[6,3] - s[6,1] <= 0
 1_Supports_6_InSameBin2_3: p_ij[1,3] - p_ij[6,3] - s[6,1] <= 0
 Gamma0ForItems_6_inBins_3_NoCut_1: p_ij[6,3] + gamma[6] <= 98
 6_OnCutIn2_1: 10000 p_ij[6,2] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 6_OnCutIn3_1: 10000 p_ij[6,3] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 AbsZ16_2: - zi[6] + zi_hi[2] - k[6,2] <= 0
 AbsZ26_2: zi[6] - zi_hi[2] - k[6,2] <= 0
 2_Below_6: zi[6] - zi_hi[2] + k[6,2] + 310 m[6,2] <= 310
 2_NotBelow_6: - zi[6] + zi_hi[2] + k[6,2] - 310 m[6,2] <= 0
 HeightCompat_6_2: h[6,2] - k[6,2] <= 0
 NoHeightCompat_6_2: - 155 h[6,2] + k[6,2] <= 0
 Overlap_6_2: - x_p[2,6] - x_p[6,2] + o[6,2] <= 0
 NoOverlap_6_2: x_p[2,6] + x_p[6,2] - 2 o[6,2] <= 0
 3_Supports_6_2: - h[6,2] - o[6,2] - s[6,2] <= -1
 3_DoesntSupport_6_2: h[6,2] + o[6,2] + 2 s[6,2] <= 2
 beta1_6_2: - s[6,2] + beta1[6,2] <= 0
 beta2_6_2: - s[6,2] + beta2[6,2] <= 0
 sBeta1Flag_6_2: - s[6,2] + beta1[6,2] <= 0
 sBeta2Flag_6_2: - s[6,2] + beta2[6,2] <= 0
 Eta1Overlap_6_2: eta1[6,2] + beta1[6,2] <= 1
 Eta2Overlap_6_2: eta2[6,2] + beta2[6,2] <= 1
 Eta1Flag_6_2: xi_l[2] - xi_l[6] - 300 eta1[6,2] <= 0
 Eta2Flag_6_2: - xi_r[2] + xi_r[6] - 300 eta2[6,2] <= 0
 2_Supports_6_InSameBin1_0: - p_ij[2,0] + p_ij[6,0] - s[6,2] <= 0
 2_Supports_6_InSameBin2_0: p_ij[2,0] - p_ij[6,0] - s[6,2] <= 0
 Gamma0ForItems_6_inBins_0_NoCut_2: p_ij[6,0] + gamma[6] <= 1
 2_Supports_6_InSameBin1_1: - p_ij[2,1] + p_ij[6,1] - s[6,2] <= 0
 2_Supports_6_InSameBin2_1: p_ij[2,1] - p_ij[6,1] - s[6,2] <= 0
 Gamma0ForItems_6_inBins_1_NoCut_2: p_ij[6,1] + gamma[6] <= 1
 2_Supports_6_InSameBin1_2: - p_ij[2,2] + p_ij[6,2] - s[6,2] <= 0
 2_Supports_6_InSameBin2_2: p_ij[2,2] - p_ij[6,2] - s[6,2] <= 0
 Gamma0ForItems_6_inBins_2_NoCut_2: p_ij[6,2] + gamma[6] <= 98
 2_Supports_6_InSameBin1_3: - p_ij[2,3] + p_ij[6,3] - s[6,2] <= 0
 2_Supports_6_InSameBin2_3: p_ij[2,3] - p_ij[6,3] - s[6,2] <= 0
 Gamma0ForItems_6_inBins_3_NoCut_2: p_ij[6,3] + gamma[6] <= 98
 6_OnCutIn2_2: 10000 p_ij[6,2] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 6_OnCutIn3_2: 10000 p_ij[6,3] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 AbsZ16_3: - zi[6] + zi_hi[3] - k[6,3] <= 0
 AbsZ26_3: zi[6] - zi_hi[3] - k[6,3] <= 0
 3_Below_6: zi[6] - zi_hi[3] + k[6,3] + 310 m[6,3] <= 310
 3_NotBelow_6: - zi[6] + zi_hi[3] + k[6,3] - 310 m[6,3] <= 0
 HeightCompat_6_3: h[6,3] - k[6,3] <= 0
 NoHeightCompat_6_3: - 155 h[6,3] + k[6,3] <= 0
 Overlap_6_3: - x_p[3,6] - x_p[6,3] + o[6,3] <= 0
 NoOverlap_6_3: x_p[3,6] + x_p[6,3] - 2 o[6,3] <= 0
 3_Supports_6_3: - h[6,3] - o[6,3] - s[6,3] <= -1
 3_DoesntSupport_6_3: h[6,3] + o[6,3] + 2 s[6,3] <= 2
 beta1_6_3: - s[6,3] + beta1[6,3] <= 0
 beta2_6_3: - s[6,3] + beta2[6,3] <= 0
 sBeta1Flag_6_3: - s[6,3] + beta1[6,3] <= 0
 sBeta2Flag_6_3: - s[6,3] + beta2[6,3] <= 0
 Eta1Overlap_6_3: eta1[6,3] + beta1[6,3] <= 1
 Eta2Overlap_6_3: eta2[6,3] + beta2[6,3] <= 1
 Eta1Flag_6_3: xi_l[3] - xi_l[6] - 300 eta1[6,3] <= 0
 Eta2Flag_6_3: - xi_r[3] + xi_r[6] - 300 eta2[6,3] <= 0
 3_Supports_6_InSameBin1_0: - p_ij[3,0] + p_ij[6,0] - s[6,3] <= 0
 3_Supports_6_InSameBin2_0: p_ij[3,0] - p_ij[6,0] - s[6,3] <= 0
 Gamma0ForItems_6_inBins_0_NoCut_3: p_ij[6,0] + gamma[6] <= 1
 3_Supports_6_InSameBin1_1: - p_ij[3,1] + p_ij[6,1] - s[6,3] <= 0
 3_Supports_6_InSameBin2_1: p_ij[3,1] - p_ij[6,1] - s[6,3] <= 0
 Gamma0ForItems_6_inBins_1_NoCut_3: p_ij[6,1] + gamma[6] <= 1
 3_Supports_6_InSameBin1_2: - p_ij[3,2] + p_ij[6,2] - s[6,3] <= 0
 3_Supports_6_InSameBin2_2: p_ij[3,2] - p_ij[6,2] - s[6,3] <= 0
 Gamma0ForItems_6_inBins_2_NoCut_3: p_ij[6,2] + gamma[6] <= 98
 3_Supports_6_InSameBin1_3: - p_ij[3,3] + p_ij[6,3] - s[6,3] <= 0
 3_Supports_6_InSameBin2_3: p_ij[3,3] - p_ij[6,3] - s[6,3] <= 0
 Gamma0ForItems_6_inBins_3_NoCut_3: p_ij[6,3] + gamma[6] <= 98
 6_OnCutIn2_3: 10000 p_ij[6,2] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 6_OnCutIn3_3: 10000 p_ij[6,3] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 AbsZ16_4: - zi[6] + zi_hi[4] - k[6,4] <= 0
 AbsZ26_4: zi[6] - zi_hi[4] - k[6,4] <= 0
 4_Below_6: zi[6] - zi_hi[4] + k[6,4] + 310 m[6,4] <= 310
 4_NotBelow_6: - zi[6] + zi_hi[4] + k[6,4] - 310 m[6,4] <= 0
 HeightCompat_6_4: h[6,4] - k[6,4] <= 0
 NoHeightCompat_6_4: - 155 h[6,4] + k[6,4] <= 0
 Overlap_6_4: - x_p[4,6] - x_p[6,4] + o[6,4] <= 0
 NoOverlap_6_4: x_p[4,6] + x_p[6,4] - 2 o[6,4] <= 0
 3_Supports_6_4: - h[6,4] - o[6,4] - s[6,4] <= -1
 3_DoesntSupport_6_4: h[6,4] + o[6,4] + 2 s[6,4] <= 2
 beta1_6_4: - s[6,4] + beta1[6,4] <= 0
 beta2_6_4: - s[6,4] + beta2[6,4] <= 0
 sBeta1Flag_6_4: - s[6,4] + beta1[6,4] <= 0
 sBeta2Flag_6_4: - s[6,4] + beta2[6,4] <= 0
 Eta1Overlap_6_4: eta1[6,4] + beta1[6,4] <= 1
 Eta2Overlap_6_4: eta2[6,4] + beta2[6,4] <= 1
 Eta1Flag_6_4: xi_l[4] - xi_l[6] - 300 eta1[6,4] <= 0
 Eta2Flag_6_4: - xi_r[4] + xi_r[6] - 300 eta2[6,4] <= 0
 4_Supports_6_InSameBin1_0: - p_ij[4,0] + p_ij[6,0] - s[6,4] <= 0
 4_Supports_6_InSameBin2_0: p_ij[4,0] - p_ij[6,0] - s[6,4] <= 0
 Gamma0ForItems_6_inBins_0_NoCut_4: p_ij[6,0] + gamma[6] <= 1
 4_Supports_6_InSameBin1_1: - p_ij[4,1] + p_ij[6,1] - s[6,4] <= 0
 4_Supports_6_InSameBin2_1: p_ij[4,1] - p_ij[6,1] - s[6,4] <= 0
 Gamma0ForItems_6_inBins_1_NoCut_4: p_ij[6,1] + gamma[6] <= 1
 4_Supports_6_InSameBin1_2: - p_ij[4,2] + p_ij[6,2] - s[6,4] <= 0
 4_Supports_6_InSameBin2_2: p_ij[4,2] - p_ij[6,2] - s[6,4] <= 0
 Gamma0ForItems_6_inBins_2_NoCut_4: p_ij[6,2] + gamma[6] <= 98
 4_Supports_6_InSameBin1_3: - p_ij[4,3] + p_ij[6,3] - s[6,4] <= 0
 4_Supports_6_InSameBin2_3: p_ij[4,3] - p_ij[6,3] - s[6,4] <= 0
 Gamma0ForItems_6_inBins_3_NoCut_4: p_ij[6,3] + gamma[6] <= 98
 6_OnCutIn2_4: 10000 p_ij[6,2] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 6_OnCutIn3_4: 10000 p_ij[6,3] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 AbsZ16_5: - zi[6] + zi_hi[5] - k[6,5] <= 0
 AbsZ26_5: zi[6] - zi_hi[5] - k[6,5] <= 0
 5_Below_6: zi[6] - zi_hi[5] + k[6,5] + 310 m[6,5] <= 310
 5_NotBelow_6: - zi[6] + zi_hi[5] + k[6,5] - 310 m[6,5] <= 0
 HeightCompat_6_5: h[6,5] - k[6,5] <= 0
 NoHeightCompat_6_5: - 155 h[6,5] + k[6,5] <= 0
 Overlap_6_5: - x_p[5,6] - x_p[6,5] + o[6,5] <= 0
 NoOverlap_6_5: x_p[5,6] + x_p[6,5] - 2 o[6,5] <= 0
 3_Supports_6_5: - h[6,5] - o[6,5] - s[6,5] <= -1
 3_DoesntSupport_6_5: h[6,5] + o[6,5] + 2 s[6,5] <= 2
 beta1_6_5: - s[6,5] + beta1[6,5] <= 0
 beta2_6_5: - s[6,5] + beta2[6,5] <= 0
 sBeta1Flag_6_5: - s[6,5] + beta1[6,5] <= 0
 sBeta2Flag_6_5: - s[6,5] + beta2[6,5] <= 0
 Eta1Overlap_6_5: eta1[6,5] + beta1[6,5] <= 1
 Eta2Overlap_6_5: eta2[6,5] + beta2[6,5] <= 1
 Eta1Flag_6_5: xi_l[5] - xi_l[6] - 300 eta1[6,5] <= 0
 Eta2Flag_6_5: - xi_r[5] + xi_r[6] - 300 eta2[6,5] <= 0
 5_Supports_6_InSameBin1_0: - p_ij[5,0] + p_ij[6,0] - s[6,5] <= 0
 5_Supports_6_InSameBin2_0: p_ij[5,0] - p_ij[6,0] - s[6,5] <= 0
 Gamma0ForItems_6_inBins_0_NoCut_5: p_ij[6,0] + gamma[6] <= 1
 5_Supports_6_InSameBin1_1: - p_ij[5,1] + p_ij[6,1] - s[6,5] <= 0
 5_Supports_6_InSameBin2_1: p_ij[5,1] - p_ij[6,1] - s[6,5] <= 0
 Gamma0ForItems_6_inBins_1_NoCut_5: p_ij[6,1] + gamma[6] <= 1
 5_Supports_6_InSameBin1_2: - p_ij[5,2] + p_ij[6,2] - s[6,5] <= 0
 5_Supports_6_InSameBin2_2: p_ij[5,2] - p_ij[6,2] - s[6,5] <= 0
 Gamma0ForItems_6_inBins_2_NoCut_5: p_ij[6,2] + gamma[6] <= 98
 5_Supports_6_InSameBin1_3: - p_ij[5,3] + p_ij[6,3] - s[6,5] <= 0
 5_Supports_6_InSameBin2_3: p_ij[5,3] - p_ij[6,3] - s[6,5] <= 0
 Gamma0ForItems_6_inBins_3_NoCut_5: p_ij[6,3] + gamma[6] <= 98
 6_OnCutIn2_5: 10000 p_ij[6,2] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 6_OnCutIn3_5: 10000 p_ij[6,3] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 6_OnCutIn2_6: 10000 p_ij[6,2] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 6_OnCutIn3_6: 10000 p_ij[6,3] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 AbsZ16_7: - zi[6] + zi_hi[7] - k[6,7] <= 0
 AbsZ26_7: zi[6] - zi_hi[7] - k[6,7] <= 0
 7_Below_6: zi[6] - zi_hi[7] + k[6,7] + 310 m[6,7] <= 310
 7_NotBelow_6: - zi[6] + zi_hi[7] + k[6,7] - 310 m[6,7] <= 0
 HeightCompat_6_7: h[6,7] - k[6,7] <= 0
 NoHeightCompat_6_7: - 155 h[6,7] + k[6,7] <= 0
 Overlap_6_7: - x_p[6,7] - x_p[7,6] + o[6,7] <= 0
 NoOverlap_6_7: x_p[6,7] + x_p[7,6] - 2 o[6,7] <= 0
 3_Supports_6_7: - h[6,7] - o[6,7] - s[6,7] <= -1
 3_DoesntSupport_6_7: h[6,7] + o[6,7] + 2 s[6,7] <= 2
 beta1_6_7: - s[6,7] + beta1[6,7] <= 0
 beta2_6_7: - s[6,7] + beta2[6,7] <= 0
 sBeta1Flag_6_7: - s[6,7] + beta1[6,7] <= 0
 sBeta2Flag_6_7: - s[6,7] + beta2[6,7] <= 0
 Eta1Overlap_6_7: eta1[6,7] + beta1[6,7] <= 1
 Eta2Overlap_6_7: eta2[6,7] + beta2[6,7] <= 1
 Eta1Flag_6_7: - xi_l[6] + xi_l[7] - 300 eta1[6,7] <= 0
 Eta2Flag_6_7: xi_r[6] - xi_r[7] - 300 eta2[6,7] <= 0
 7_Supports_6_InSameBin1_0: p_ij[6,0] - p_ij[7,0] - s[6,7] <= 0
 7_Supports_6_InSameBin2_0: - p_ij[6,0] + p_ij[7,0] - s[6,7] <= 0
 Gamma0ForItems_6_inBins_0_NoCut_7: p_ij[6,0] + gamma[6] <= 1
 7_Supports_6_InSameBin1_1: p_ij[6,1] - p_ij[7,1] - s[6,7] <= 0
 7_Supports_6_InSameBin2_1: - p_ij[6,1] + p_ij[7,1] - s[6,7] <= 0
 Gamma0ForItems_6_inBins_1_NoCut_7: p_ij[6,1] + gamma[6] <= 1
 7_Supports_6_InSameBin1_2: p_ij[6,2] - p_ij[7,2] - s[6,7] <= 0
 7_Supports_6_InSameBin2_2: - p_ij[6,2] + p_ij[7,2] - s[6,7] <= 0
 Gamma0ForItems_6_inBins_2_NoCut_7: p_ij[6,2] + gamma[6] <= 98
 7_Supports_6_InSameBin1_3: p_ij[6,3] - p_ij[7,3] - s[6,7] <= 0
 7_Supports_6_InSameBin2_3: - p_ij[6,3] + p_ij[7,3] - s[6,7] <= 0
 Gamma0ForItems_6_inBins_3_NoCut_7: p_ij[6,3] + gamma[6] <= 98
 6_OnCutIn2_7: 10000 p_ij[6,2] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 6_OnCutIn3_7: 10000 p_ij[6,3] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 AbsZ16_8: - zi[6] + zi_hi[8] - k[6,8] <= 0
 AbsZ26_8: zi[6] - zi_hi[8] - k[6,8] <= 0
 8_Below_6: zi[6] - zi_hi[8] + k[6,8] + 310 m[6,8] <= 310
 8_NotBelow_6: - zi[6] + zi_hi[8] + k[6,8] - 310 m[6,8] <= 0
 HeightCompat_6_8: h[6,8] - k[6,8] <= 0
 NoHeightCompat_6_8: - 155 h[6,8] + k[6,8] <= 0
 Overlap_6_8: - x_p[6,8] - x_p[8,6] + o[6,8] <= 0
 NoOverlap_6_8: x_p[6,8] + x_p[8,6] - 2 o[6,8] <= 0
 3_Supports_6_8: - h[6,8] - o[6,8] - s[6,8] <= -1
 3_DoesntSupport_6_8: h[6,8] + o[6,8] + 2 s[6,8] <= 2
 beta1_6_8: - s[6,8] + beta1[6,8] <= 0
 beta2_6_8: - s[6,8] + beta2[6,8] <= 0
 sBeta1Flag_6_8: - s[6,8] + beta1[6,8] <= 0
 sBeta2Flag_6_8: - s[6,8] + beta2[6,8] <= 0
 Eta1Overlap_6_8: eta1[6,8] + beta1[6,8] <= 1
 Eta2Overlap_6_8: eta2[6,8] + beta2[6,8] <= 1
 Eta1Flag_6_8: - xi_l[6] + xi_l[8] - 300 eta1[6,8] <= 0
 Eta2Flag_6_8: xi_r[6] - xi_r[8] - 300 eta2[6,8] <= 0
 8_Supports_6_InSameBin1_0: p_ij[6,0] - p_ij[8,0] - s[6,8] <= 0
 8_Supports_6_InSameBin2_0: - p_ij[6,0] + p_ij[8,0] - s[6,8] <= 0
 Gamma0ForItems_6_inBins_0_NoCut_8: p_ij[6,0] + gamma[6] <= 1
 8_Supports_6_InSameBin1_1: p_ij[6,1] - p_ij[8,1] - s[6,8] <= 0
 8_Supports_6_InSameBin2_1: - p_ij[6,1] + p_ij[8,1] - s[6,8] <= 0
 Gamma0ForItems_6_inBins_1_NoCut_8: p_ij[6,1] + gamma[6] <= 1
 8_Supports_6_InSameBin1_2: p_ij[6,2] - p_ij[8,2] - s[6,8] <= 0
 8_Supports_6_InSameBin2_2: - p_ij[6,2] + p_ij[8,2] - s[6,8] <= 0
 Gamma0ForItems_6_inBins_2_NoCut_8: p_ij[6,2] + gamma[6] <= 98
 8_Supports_6_InSameBin1_3: p_ij[6,3] - p_ij[8,3] - s[6,8] <= 0
 8_Supports_6_InSameBin2_3: - p_ij[6,3] + p_ij[8,3] - s[6,8] <= 0
 Gamma0ForItems_6_inBins_3_NoCut_8: p_ij[6,3] + gamma[6] <= 98
 6_OnCutIn2_8: 10000 p_ij[6,2] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 6_OnCutIn3_8: 10000 p_ij[6,3] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 AbsZ16_9: - zi[6] + zi_hi[9] - k[6,9] <= 0
 AbsZ26_9: zi[6] - zi_hi[9] - k[6,9] <= 0
 9_Below_6: zi[6] - zi_hi[9] + k[6,9] + 310 m[6,9] <= 310
 9_NotBelow_6: - zi[6] + zi_hi[9] + k[6,9] - 310 m[6,9] <= 0
 HeightCompat_6_9: h[6,9] - k[6,9] <= 0
 NoHeightCompat_6_9: - 155 h[6,9] + k[6,9] <= 0
 Overlap_6_9: - x_p[6,9] - x_p[9,6] + o[6,9] <= 0
 NoOverlap_6_9: x_p[6,9] + x_p[9,6] - 2 o[6,9] <= 0
 3_Supports_6_9: - h[6,9] - o[6,9] - s[6,9] <= -1
 3_DoesntSupport_6_9: h[6,9] + o[6,9] + 2 s[6,9] <= 2
 beta1_6_9: - s[6,9] + beta1[6,9] <= 0
 beta2_6_9: - s[6,9] + beta2[6,9] <= 0
 sBeta1Flag_6_9: - s[6,9] + beta1[6,9] <= 0
 sBeta2Flag_6_9: - s[6,9] + beta2[6,9] <= 0
 Eta1Overlap_6_9: eta1[6,9] + beta1[6,9] <= 1
 Eta2Overlap_6_9: eta2[6,9] + beta2[6,9] <= 1
 Eta1Flag_6_9: - xi_l[6] + xi_l[9] - 300 eta1[6,9] <= 0
 Eta2Flag_6_9: xi_r[6] - xi_r[9] - 300 eta2[6,9] <= 0
 9_Supports_6_InSameBin1_0: p_ij[6,0] - p_ij[9,0] - s[6,9] <= 0
 9_Supports_6_InSameBin2_0: - p_ij[6,0] + p_ij[9,0] - s[6,9] <= 0
 Gamma0ForItems_6_inBins_0_NoCut_9: p_ij[6,0] + gamma[6] <= 1
 9_Supports_6_InSameBin1_1: p_ij[6,1] - p_ij[9,1] - s[6,9] <= 0
 9_Supports_6_InSameBin2_1: - p_ij[6,1] + p_ij[9,1] - s[6,9] <= 0
 Gamma0ForItems_6_inBins_1_NoCut_9: p_ij[6,1] + gamma[6] <= 1
 9_Supports_6_InSameBin1_2: p_ij[6,2] - p_ij[9,2] - s[6,9] <= 0
 9_Supports_6_InSameBin2_2: - p_ij[6,2] + p_ij[9,2] - s[6,9] <= 0
 Gamma0ForItems_6_inBins_2_NoCut_9: p_ij[6,2] + gamma[6] <= 98
 9_Supports_6_InSameBin1_3: p_ij[6,3] - p_ij[9,3] - s[6,9] <= 0
 9_Supports_6_InSameBin2_3: - p_ij[6,3] + p_ij[9,3] - s[6,9] <= 0
 Gamma0ForItems_6_inBins_3_NoCut_9: p_ij[6,3] + gamma[6] <= 98
 6_OnCutIn2_9: 10000 p_ij[6,2] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 6_OnCutIn3_9: 10000 p_ij[6,3] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 AbsZ16_10: - zi[6] + zi_hi[10] - k[6,10] <= 0
 AbsZ26_10: zi[6] - zi_hi[10] - k[6,10] <= 0
 10_Below_6: zi[6] - zi_hi[10] + k[6,10] + 310 m[6,10] <= 310
 10_NotBelow_6: - zi[6] + zi_hi[10] + k[6,10] - 310 m[6,10] <= 0
 HeightCompat_6_10: h[6,10] - k[6,10] <= 0
 NoHeightCompat_6_10: - 155 h[6,10] + k[6,10] <= 0
 Overlap_6_10: - x_p[6,10] - x_p[10,6] + o[6,10] <= 0
 NoOverlap_6_10: x_p[6,10] + x_p[10,6] - 2 o[6,10] <= 0
 3_Supports_6_10: - h[6,10] - o[6,10] - s[6,10] <= -1
 3_DoesntSupport_6_10: h[6,10] + o[6,10] + 2 s[6,10] <= 2
 beta1_6_10: - s[6,10] + beta1[6,10] <= 0
 beta2_6_10: - s[6,10] + beta2[6,10] <= 0
 sBeta1Flag_6_10: - s[6,10] + beta1[6,10] <= 0
 sBeta2Flag_6_10: - s[6,10] + beta2[6,10] <= 0
 Eta1Overlap_6_10: eta1[6,10] + beta1[6,10] <= 1
 Eta2Overlap_6_10: eta2[6,10] + beta2[6,10] <= 1
 Eta1Flag_6_10: - xi_l[6] + xi_l[10] - 300 eta1[6,10] <= 0
 Eta2Flag_6_10: xi_r[6] - xi_r[10] - 300 eta2[6,10] <= 0
 10_Supports_6_InSameBin1_0: p_ij[6,0] - p_ij[10,0] - s[6,10] <= 0
 10_Supports_6_InSameBin2_0: - p_ij[6,0] + p_ij[10,0] - s[6,10] <= 0
 Gamma0ForItems_6_inBins_0_NoCut_10: p_ij[6,0] + gamma[6] <= 1
 10_Supports_6_InSameBin1_1: p_ij[6,1] - p_ij[10,1] - s[6,10] <= 0
 10_Supports_6_InSameBin2_1: - p_ij[6,1] + p_ij[10,1] - s[6,10] <= 0
 Gamma0ForItems_6_inBins_1_NoCut_10: p_ij[6,1] + gamma[6] <= 1
 10_Supports_6_InSameBin1_2: p_ij[6,2] - p_ij[10,2] - s[6,10] <= 0
 10_Supports_6_InSameBin2_2: - p_ij[6,2] + p_ij[10,2] - s[6,10] <= 0
 Gamma0ForItems_6_inBins_2_NoCut_10: p_ij[6,2] + gamma[6] <= 98
 10_Supports_6_InSameBin1_3: p_ij[6,3] - p_ij[10,3] - s[6,10] <= 0
 10_Supports_6_InSameBin2_3: - p_ij[6,3] + p_ij[10,3] - s[6,10] <= 0
 Gamma0ForItems_6_inBins_3_NoCut_10: p_ij[6,3] + gamma[6] <= 98
 6_OnCutIn2_10: 10000 p_ij[6,2] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 6_OnCutIn3_10: 10000 p_ij[6,3] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 AbsZ16_11: - zi[6] + zi_hi[11] - k[6,11] <= 0
 AbsZ26_11: zi[6] - zi_hi[11] - k[6,11] <= 0
 11_Below_6: zi[6] - zi_hi[11] + k[6,11] + 310 m[6,11] <= 310
 11_NotBelow_6: - zi[6] + zi_hi[11] + k[6,11] - 310 m[6,11] <= 0
 HeightCompat_6_11: h[6,11] - k[6,11] <= 0
 NoHeightCompat_6_11: - 155 h[6,11] + k[6,11] <= 0
 Overlap_6_11: - x_p[6,11] - x_p[11,6] + o[6,11] <= 0
 NoOverlap_6_11: x_p[6,11] + x_p[11,6] - 2 o[6,11] <= 0
 3_Supports_6_11: - h[6,11] - o[6,11] - s[6,11] <= -1
 3_DoesntSupport_6_11: h[6,11] + o[6,11] + 2 s[6,11] <= 2
 beta1_6_11: - s[6,11] + beta1[6,11] <= 0
 beta2_6_11: - s[6,11] + beta2[6,11] <= 0
 sBeta1Flag_6_11: - s[6,11] + beta1[6,11] <= 0
 sBeta2Flag_6_11: - s[6,11] + beta2[6,11] <= 0
 Eta1Overlap_6_11: eta1[6,11] + beta1[6,11] <= 1
 Eta2Overlap_6_11: eta2[6,11] + beta2[6,11] <= 1
 Eta1Flag_6_11: - xi_l[6] + xi_l[11] - 300 eta1[6,11] <= 0
 Eta2Flag_6_11: xi_r[6] - xi_r[11] - 300 eta2[6,11] <= 0
 11_Supports_6_InSameBin1_0: p_ij[6,0] - p_ij[11,0] - s[6,11] <= 0
 11_Supports_6_InSameBin2_0: - p_ij[6,0] + p_ij[11,0] - s[6,11] <= 0
 Gamma0ForItems_6_inBins_0_NoCut_11: p_ij[6,0] + gamma[6] <= 1
 11_Supports_6_InSameBin1_1: p_ij[6,1] - p_ij[11,1] - s[6,11] <= 0
 11_Supports_6_InSameBin2_1: - p_ij[6,1] + p_ij[11,1] - s[6,11] <= 0
 Gamma0ForItems_6_inBins_1_NoCut_11: p_ij[6,1] + gamma[6] <= 1
 11_Supports_6_InSameBin1_2: p_ij[6,2] - p_ij[11,2] - s[6,11] <= 0
 11_Supports_6_InSameBin2_2: - p_ij[6,2] + p_ij[11,2] - s[6,11] <= 0
 Gamma0ForItems_6_inBins_2_NoCut_11: p_ij[6,2] + gamma[6] <= 98
 11_Supports_6_InSameBin1_3: p_ij[6,3] - p_ij[11,3] - s[6,11] <= 0
 11_Supports_6_InSameBin2_3: - p_ij[6,3] + p_ij[11,3] - s[6,11] <= 0
 Gamma0ForItems_6_inBins_3_NoCut_11: p_ij[6,3] + gamma[6] <= 98
 6_OnCutIn2_11: 10000 p_ij[6,2] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 6_OnCutIn3_11: 10000 p_ij[6,3] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 AbsZ16_12: - zi[6] + zi_hi[12] - k[6,12] <= 0
 AbsZ26_12: zi[6] - zi_hi[12] - k[6,12] <= 0
 12_Below_6: zi[6] - zi_hi[12] + k[6,12] + 310 m[6,12] <= 310
 12_NotBelow_6: - zi[6] + zi_hi[12] + k[6,12] - 310 m[6,12] <= 0
 HeightCompat_6_12: h[6,12] - k[6,12] <= 0
 NoHeightCompat_6_12: - 155 h[6,12] + k[6,12] <= 0
 Overlap_6_12: - x_p[6,12] - x_p[12,6] + o[6,12] <= 0
 NoOverlap_6_12: x_p[6,12] + x_p[12,6] - 2 o[6,12] <= 0
 3_Supports_6_12: - h[6,12] - o[6,12] - s[6,12] <= -1
 3_DoesntSupport_6_12: h[6,12] + o[6,12] + 2 s[6,12] <= 2
 beta1_6_12: - s[6,12] + beta1[6,12] <= 0
 beta2_6_12: - s[6,12] + beta2[6,12] <= 0
 sBeta1Flag_6_12: - s[6,12] + beta1[6,12] <= 0
 sBeta2Flag_6_12: - s[6,12] + beta2[6,12] <= 0
 Eta1Overlap_6_12: eta1[6,12] + beta1[6,12] <= 1
 Eta2Overlap_6_12: eta2[6,12] + beta2[6,12] <= 1
 Eta1Flag_6_12: - xi_l[6] + xi_l[12] - 300 eta1[6,12] <= 0
 Eta2Flag_6_12: xi_r[6] - xi_r[12] - 300 eta2[6,12] <= 0
 12_Supports_6_InSameBin1_0: p_ij[6,0] - p_ij[12,0] - s[6,12] <= 0
 12_Supports_6_InSameBin2_0: - p_ij[6,0] + p_ij[12,0] - s[6,12] <= 0
 Gamma0ForItems_6_inBins_0_NoCut_12: p_ij[6,0] + gamma[6] <= 1
 12_Supports_6_InSameBin1_1: p_ij[6,1] - p_ij[12,1] - s[6,12] <= 0
 12_Supports_6_InSameBin2_1: - p_ij[6,1] + p_ij[12,1] - s[6,12] <= 0
 Gamma0ForItems_6_inBins_1_NoCut_12: p_ij[6,1] + gamma[6] <= 1
 12_Supports_6_InSameBin1_2: p_ij[6,2] - p_ij[12,2] - s[6,12] <= 0
 12_Supports_6_InSameBin2_2: - p_ij[6,2] + p_ij[12,2] - s[6,12] <= 0
 Gamma0ForItems_6_inBins_2_NoCut_12: p_ij[6,2] + gamma[6] <= 98
 12_Supports_6_InSameBin1_3: p_ij[6,3] - p_ij[12,3] - s[6,12] <= 0
 12_Supports_6_InSameBin2_3: - p_ij[6,3] + p_ij[12,3] - s[6,12] <= 0
 Gamma0ForItems_6_inBins_3_NoCut_12: p_ij[6,3] + gamma[6] <= 98
 6_OnCutIn2_12: 10000 p_ij[6,2] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 6_OnCutIn3_12: 10000 p_ij[6,3] + 1.261904761904762 xi_l[6] + zi[6]
   + 10000 gamma[6] <= 20053
 GroundConstraint_7: zi[7] + 155 g[7] <= 155
 AbsZ17_0: - zi[7] + zi_hi[0] - k[7,0] <= 0
 AbsZ27_0: zi[7] - zi_hi[0] - k[7,0] <= 0
 0_Below_7: zi[7] - zi_hi[0] + k[7,0] + 310 m[7,0] <= 310
 0_NotBelow_7: - zi[7] + zi_hi[0] + k[7,0] - 310 m[7,0] <= 0
 HeightCompat_7_0: h[7,0] - k[7,0] <= 0
 NoHeightCompat_7_0: - 155 h[7,0] + k[7,0] <= 0
 Overlap_7_0: - x_p[0,7] - x_p[7,0] + o[7,0] <= 0
 NoOverlap_7_0: x_p[0,7] + x_p[7,0] - 2 o[7,0] <= 0
 3_Supports_7_0: - h[7,0] - o[7,0] - s[7,0] <= -1
 3_DoesntSupport_7_0: h[7,0] + o[7,0] + 2 s[7,0] <= 2
 beta1_7_0: - s[7,0] + beta1[7,0] <= 0
 beta2_7_0: - s[7,0] + beta2[7,0] <= 0
 sBeta1Flag_7_0: - s[7,0] + beta1[7,0] <= 0
 sBeta2Flag_7_0: - s[7,0] + beta2[7,0] <= 0
 Eta1Overlap_7_0: eta1[7,0] + beta1[7,0] <= 1
 Eta2Overlap_7_0: eta2[7,0] + beta2[7,0] <= 1
 Eta1Flag_7_0: xi_l[0] - xi_l[7] - 300 eta1[7,0] <= 0
 Eta2Flag_7_0: - xi_r[0] + xi_r[7] - 300 eta2[7,0] <= 0
 0_Supports_7_InSameBin1_0: - p_ij[0,0] + p_ij[7,0] - s[7,0] <= 0
 0_Supports_7_InSameBin2_0: p_ij[0,0] - p_ij[7,0] - s[7,0] <= 0
 Gamma0ForItems_7_inBins_0_NoCut_0: p_ij[7,0] + gamma[7] <= 1
 0_Supports_7_InSameBin1_1: - p_ij[0,1] + p_ij[7,1] - s[7,0] <= 0
 0_Supports_7_InSameBin2_1: p_ij[0,1] - p_ij[7,1] - s[7,0] <= 0
 Gamma0ForItems_7_inBins_1_NoCut_0: p_ij[7,1] + gamma[7] <= 1
 0_Supports_7_InSameBin1_2: - p_ij[0,2] + p_ij[7,2] - s[7,0] <= 0
 0_Supports_7_InSameBin2_2: p_ij[0,2] - p_ij[7,2] - s[7,0] <= 0
 Gamma0ForItems_7_inBins_2_NoCut_0: p_ij[7,2] + gamma[7] <= 98
 0_Supports_7_InSameBin1_3: - p_ij[0,3] + p_ij[7,3] - s[7,0] <= 0
 0_Supports_7_InSameBin2_3: p_ij[0,3] - p_ij[7,3] - s[7,0] <= 0
 Gamma0ForItems_7_inBins_3_NoCut_0: p_ij[7,3] + gamma[7] <= 98
 7_OnCutIn2_0: 10000 p_ij[7,2] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 7_OnCutIn3_0: 10000 p_ij[7,3] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 AbsZ17_1: - zi[7] + zi_hi[1] - k[7,1] <= 0
 AbsZ27_1: zi[7] - zi_hi[1] - k[7,1] <= 0
 1_Below_7: zi[7] - zi_hi[1] + k[7,1] + 310 m[7,1] <= 310
 1_NotBelow_7: - zi[7] + zi_hi[1] + k[7,1] - 310 m[7,1] <= 0
 HeightCompat_7_1: h[7,1] - k[7,1] <= 0
 NoHeightCompat_7_1: - 155 h[7,1] + k[7,1] <= 0
 Overlap_7_1: - x_p[1,7] - x_p[7,1] + o[7,1] <= 0
 NoOverlap_7_1: x_p[1,7] + x_p[7,1] - 2 o[7,1] <= 0
 3_Supports_7_1: - h[7,1] - o[7,1] - s[7,1] <= -1
 3_DoesntSupport_7_1: h[7,1] + o[7,1] + 2 s[7,1] <= 2
 beta1_7_1: - s[7,1] + beta1[7,1] <= 0
 beta2_7_1: - s[7,1] + beta2[7,1] <= 0
 sBeta1Flag_7_1: - s[7,1] + beta1[7,1] <= 0
 sBeta2Flag_7_1: - s[7,1] + beta2[7,1] <= 0
 Eta1Overlap_7_1: eta1[7,1] + beta1[7,1] <= 1
 Eta2Overlap_7_1: eta2[7,1] + beta2[7,1] <= 1
 Eta1Flag_7_1: xi_l[1] - xi_l[7] - 300 eta1[7,1] <= 0
 Eta2Flag_7_1: - xi_r[1] + xi_r[7] - 300 eta2[7,1] <= 0
 1_Supports_7_InSameBin1_0: - p_ij[1,0] + p_ij[7,0] - s[7,1] <= 0
 1_Supports_7_InSameBin2_0: p_ij[1,0] - p_ij[7,0] - s[7,1] <= 0
 Gamma0ForItems_7_inBins_0_NoCut_1: p_ij[7,0] + gamma[7] <= 1
 1_Supports_7_InSameBin1_1: - p_ij[1,1] + p_ij[7,1] - s[7,1] <= 0
 1_Supports_7_InSameBin2_1: p_ij[1,1] - p_ij[7,1] - s[7,1] <= 0
 Gamma0ForItems_7_inBins_1_NoCut_1: p_ij[7,1] + gamma[7] <= 1
 1_Supports_7_InSameBin1_2: - p_ij[1,2] + p_ij[7,2] - s[7,1] <= 0
 1_Supports_7_InSameBin2_2: p_ij[1,2] - p_ij[7,2] - s[7,1] <= 0
 Gamma0ForItems_7_inBins_2_NoCut_1: p_ij[7,2] + gamma[7] <= 98
 1_Supports_7_InSameBin1_3: - p_ij[1,3] + p_ij[7,3] - s[7,1] <= 0
 1_Supports_7_InSameBin2_3: p_ij[1,3] - p_ij[7,3] - s[7,1] <= 0
 Gamma0ForItems_7_inBins_3_NoCut_1: p_ij[7,3] + gamma[7] <= 98
 7_OnCutIn2_1: 10000 p_ij[7,2] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 7_OnCutIn3_1: 10000 p_ij[7,3] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 AbsZ17_2: - zi[7] + zi_hi[2] - k[7,2] <= 0
 AbsZ27_2: zi[7] - zi_hi[2] - k[7,2] <= 0
 2_Below_7: zi[7] - zi_hi[2] + k[7,2] + 310 m[7,2] <= 310
 2_NotBelow_7: - zi[7] + zi_hi[2] + k[7,2] - 310 m[7,2] <= 0
 HeightCompat_7_2: h[7,2] - k[7,2] <= 0
 NoHeightCompat_7_2: - 155 h[7,2] + k[7,2] <= 0
 Overlap_7_2: - x_p[2,7] - x_p[7,2] + o[7,2] <= 0
 NoOverlap_7_2: x_p[2,7] + x_p[7,2] - 2 o[7,2] <= 0
 3_Supports_7_2: - h[7,2] - o[7,2] - s[7,2] <= -1
 3_DoesntSupport_7_2: h[7,2] + o[7,2] + 2 s[7,2] <= 2
 beta1_7_2: - s[7,2] + beta1[7,2] <= 0
 beta2_7_2: - s[7,2] + beta2[7,2] <= 0
 sBeta1Flag_7_2: - s[7,2] + beta1[7,2] <= 0
 sBeta2Flag_7_2: - s[7,2] + beta2[7,2] <= 0
 Eta1Overlap_7_2: eta1[7,2] + beta1[7,2] <= 1
 Eta2Overlap_7_2: eta2[7,2] + beta2[7,2] <= 1
 Eta1Flag_7_2: xi_l[2] - xi_l[7] - 300 eta1[7,2] <= 0
 Eta2Flag_7_2: - xi_r[2] + xi_r[7] - 300 eta2[7,2] <= 0
 2_Supports_7_InSameBin1_0: - p_ij[2,0] + p_ij[7,0] - s[7,2] <= 0
 2_Supports_7_InSameBin2_0: p_ij[2,0] - p_ij[7,0] - s[7,2] <= 0
 Gamma0ForItems_7_inBins_0_NoCut_2: p_ij[7,0] + gamma[7] <= 1
 2_Supports_7_InSameBin1_1: - p_ij[2,1] + p_ij[7,1] - s[7,2] <= 0
 2_Supports_7_InSameBin2_1: p_ij[2,1] - p_ij[7,1] - s[7,2] <= 0
 Gamma0ForItems_7_inBins_1_NoCut_2: p_ij[7,1] + gamma[7] <= 1
 2_Supports_7_InSameBin1_2: - p_ij[2,2] + p_ij[7,2] - s[7,2] <= 0
 2_Supports_7_InSameBin2_2: p_ij[2,2] - p_ij[7,2] - s[7,2] <= 0
 Gamma0ForItems_7_inBins_2_NoCut_2: p_ij[7,2] + gamma[7] <= 98
 2_Supports_7_InSameBin1_3: - p_ij[2,3] + p_ij[7,3] - s[7,2] <= 0
 2_Supports_7_InSameBin2_3: p_ij[2,3] - p_ij[7,3] - s[7,2] <= 0
 Gamma0ForItems_7_inBins_3_NoCut_2: p_ij[7,3] + gamma[7] <= 98
 7_OnCutIn2_2: 10000 p_ij[7,2] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 7_OnCutIn3_2: 10000 p_ij[7,3] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 AbsZ17_3: - zi[7] + zi_hi[3] - k[7,3] <= 0
 AbsZ27_3: zi[7] - zi_hi[3] - k[7,3] <= 0
 3_Below_7: zi[7] - zi_hi[3] + k[7,3] + 310 m[7,3] <= 310
 3_NotBelow_7: - zi[7] + zi_hi[3] + k[7,3] - 310 m[7,3] <= 0
 HeightCompat_7_3: h[7,3] - k[7,3] <= 0
 NoHeightCompat_7_3: - 155 h[7,3] + k[7,3] <= 0
 Overlap_7_3: - x_p[3,7] - x_p[7,3] + o[7,3] <= 0
 NoOverlap_7_3: x_p[3,7] + x_p[7,3] - 2 o[7,3] <= 0
 3_Supports_7_3: - h[7,3] - o[7,3] - s[7,3] <= -1
 3_DoesntSupport_7_3: h[7,3] + o[7,3] + 2 s[7,3] <= 2
 beta1_7_3: - s[7,3] + beta1[7,3] <= 0
 beta2_7_3: - s[7,3] + beta2[7,3] <= 0
 sBeta1Flag_7_3: - s[7,3] + beta1[7,3] <= 0
 sBeta2Flag_7_3: - s[7,3] + beta2[7,3] <= 0
 Eta1Overlap_7_3: eta1[7,3] + beta1[7,3] <= 1
 Eta2Overlap_7_3: eta2[7,3] + beta2[7,3] <= 1
 Eta1Flag_7_3: xi_l[3] - xi_l[7] - 300 eta1[7,3] <= 0
 Eta2Flag_7_3: - xi_r[3] + xi_r[7] - 300 eta2[7,3] <= 0
 3_Supports_7_InSameBin1_0: - p_ij[3,0] + p_ij[7,0] - s[7,3] <= 0
 3_Supports_7_InSameBin2_0: p_ij[3,0] - p_ij[7,0] - s[7,3] <= 0
 Gamma0ForItems_7_inBins_0_NoCut_3: p_ij[7,0] + gamma[7] <= 1
 3_Supports_7_InSameBin1_1: - p_ij[3,1] + p_ij[7,1] - s[7,3] <= 0
 3_Supports_7_InSameBin2_1: p_ij[3,1] - p_ij[7,1] - s[7,3] <= 0
 Gamma0ForItems_7_inBins_1_NoCut_3: p_ij[7,1] + gamma[7] <= 1
 3_Supports_7_InSameBin1_2: - p_ij[3,2] + p_ij[7,2] - s[7,3] <= 0
 3_Supports_7_InSameBin2_2: p_ij[3,2] - p_ij[7,2] - s[7,3] <= 0
 Gamma0ForItems_7_inBins_2_NoCut_3: p_ij[7,2] + gamma[7] <= 98
 3_Supports_7_InSameBin1_3: - p_ij[3,3] + p_ij[7,3] - s[7,3] <= 0
 3_Supports_7_InSameBin2_3: p_ij[3,3] - p_ij[7,3] - s[7,3] <= 0
 Gamma0ForItems_7_inBins_3_NoCut_3: p_ij[7,3] + gamma[7] <= 98
 7_OnCutIn2_3: 10000 p_ij[7,2] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 7_OnCutIn3_3: 10000 p_ij[7,3] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 AbsZ17_4: - zi[7] + zi_hi[4] - k[7,4] <= 0
 AbsZ27_4: zi[7] - zi_hi[4] - k[7,4] <= 0
 4_Below_7: zi[7] - zi_hi[4] + k[7,4] + 310 m[7,4] <= 310
 4_NotBelow_7: - zi[7] + zi_hi[4] + k[7,4] - 310 m[7,4] <= 0
 HeightCompat_7_4: h[7,4] - k[7,4] <= 0
 NoHeightCompat_7_4: - 155 h[7,4] + k[7,4] <= 0
 Overlap_7_4: - x_p[4,7] - x_p[7,4] + o[7,4] <= 0
 NoOverlap_7_4: x_p[4,7] + x_p[7,4] - 2 o[7,4] <= 0
 3_Supports_7_4: - h[7,4] - o[7,4] - s[7,4] <= -1
 3_DoesntSupport_7_4: h[7,4] + o[7,4] + 2 s[7,4] <= 2
 beta1_7_4: - s[7,4] + beta1[7,4] <= 0
 beta2_7_4: - s[7,4] + beta2[7,4] <= 0
 sBeta1Flag_7_4: - s[7,4] + beta1[7,4] <= 0
 sBeta2Flag_7_4: - s[7,4] + beta2[7,4] <= 0
 Eta1Overlap_7_4: eta1[7,4] + beta1[7,4] <= 1
 Eta2Overlap_7_4: eta2[7,4] + beta2[7,4] <= 1
 Eta1Flag_7_4: xi_l[4] - xi_l[7] - 300 eta1[7,4] <= 0
 Eta2Flag_7_4: - xi_r[4] + xi_r[7] - 300 eta2[7,4] <= 0
 4_Supports_7_InSameBin1_0: - p_ij[4,0] + p_ij[7,0] - s[7,4] <= 0
 4_Supports_7_InSameBin2_0: p_ij[4,0] - p_ij[7,0] - s[7,4] <= 0
 Gamma0ForItems_7_inBins_0_NoCut_4: p_ij[7,0] + gamma[7] <= 1
 4_Supports_7_InSameBin1_1: - p_ij[4,1] + p_ij[7,1] - s[7,4] <= 0
 4_Supports_7_InSameBin2_1: p_ij[4,1] - p_ij[7,1] - s[7,4] <= 0
 Gamma0ForItems_7_inBins_1_NoCut_4: p_ij[7,1] + gamma[7] <= 1
 4_Supports_7_InSameBin1_2: - p_ij[4,2] + p_ij[7,2] - s[7,4] <= 0
 4_Supports_7_InSameBin2_2: p_ij[4,2] - p_ij[7,2] - s[7,4] <= 0
 Gamma0ForItems_7_inBins_2_NoCut_4: p_ij[7,2] + gamma[7] <= 98
 4_Supports_7_InSameBin1_3: - p_ij[4,3] + p_ij[7,3] - s[7,4] <= 0
 4_Supports_7_InSameBin2_3: p_ij[4,3] - p_ij[7,3] - s[7,4] <= 0
 Gamma0ForItems_7_inBins_3_NoCut_4: p_ij[7,3] + gamma[7] <= 98
 7_OnCutIn2_4: 10000 p_ij[7,2] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 7_OnCutIn3_4: 10000 p_ij[7,3] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 AbsZ17_5: - zi[7] + zi_hi[5] - k[7,5] <= 0
 AbsZ27_5: zi[7] - zi_hi[5] - k[7,5] <= 0
 5_Below_7: zi[7] - zi_hi[5] + k[7,5] + 310 m[7,5] <= 310
 5_NotBelow_7: - zi[7] + zi_hi[5] + k[7,5] - 310 m[7,5] <= 0
 HeightCompat_7_5: h[7,5] - k[7,5] <= 0
 NoHeightCompat_7_5: - 155 h[7,5] + k[7,5] <= 0
 Overlap_7_5: - x_p[5,7] - x_p[7,5] + o[7,5] <= 0
 NoOverlap_7_5: x_p[5,7] + x_p[7,5] - 2 o[7,5] <= 0
 3_Supports_7_5: - h[7,5] - o[7,5] - s[7,5] <= -1
 3_DoesntSupport_7_5: h[7,5] + o[7,5] + 2 s[7,5] <= 2
 beta1_7_5: - s[7,5] + beta1[7,5] <= 0
 beta2_7_5: - s[7,5] + beta2[7,5] <= 0
 sBeta1Flag_7_5: - s[7,5] + beta1[7,5] <= 0
 sBeta2Flag_7_5: - s[7,5] + beta2[7,5] <= 0
 Eta1Overlap_7_5: eta1[7,5] + beta1[7,5] <= 1
 Eta2Overlap_7_5: eta2[7,5] + beta2[7,5] <= 1
 Eta1Flag_7_5: xi_l[5] - xi_l[7] - 300 eta1[7,5] <= 0
 Eta2Flag_7_5: - xi_r[5] + xi_r[7] - 300 eta2[7,5] <= 0
 5_Supports_7_InSameBin1_0: - p_ij[5,0] + p_ij[7,0] - s[7,5] <= 0
 5_Supports_7_InSameBin2_0: p_ij[5,0] - p_ij[7,0] - s[7,5] <= 0
 Gamma0ForItems_7_inBins_0_NoCut_5: p_ij[7,0] + gamma[7] <= 1
 5_Supports_7_InSameBin1_1: - p_ij[5,1] + p_ij[7,1] - s[7,5] <= 0
 5_Supports_7_InSameBin2_1: p_ij[5,1] - p_ij[7,1] - s[7,5] <= 0
 Gamma0ForItems_7_inBins_1_NoCut_5: p_ij[7,1] + gamma[7] <= 1
 5_Supports_7_InSameBin1_2: - p_ij[5,2] + p_ij[7,2] - s[7,5] <= 0
 5_Supports_7_InSameBin2_2: p_ij[5,2] - p_ij[7,2] - s[7,5] <= 0
 Gamma0ForItems_7_inBins_2_NoCut_5: p_ij[7,2] + gamma[7] <= 98
 5_Supports_7_InSameBin1_3: - p_ij[5,3] + p_ij[7,3] - s[7,5] <= 0
 5_Supports_7_InSameBin2_3: p_ij[5,3] - p_ij[7,3] - s[7,5] <= 0
 Gamma0ForItems_7_inBins_3_NoCut_5: p_ij[7,3] + gamma[7] <= 98
 7_OnCutIn2_5: 10000 p_ij[7,2] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 7_OnCutIn3_5: 10000 p_ij[7,3] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 AbsZ17_6: - zi[7] + zi_hi[6] - k[7,6] <= 0
 AbsZ27_6: zi[7] - zi_hi[6] - k[7,6] <= 0
 6_Below_7: zi[7] - zi_hi[6] + k[7,6] + 310 m[7,6] <= 310
 6_NotBelow_7: - zi[7] + zi_hi[6] + k[7,6] - 310 m[7,6] <= 0
 HeightCompat_7_6: h[7,6] - k[7,6] <= 0
 NoHeightCompat_7_6: - 155 h[7,6] + k[7,6] <= 0
 Overlap_7_6: - x_p[6,7] - x_p[7,6] + o[7,6] <= 0
 NoOverlap_7_6: x_p[6,7] + x_p[7,6] - 2 o[7,6] <= 0
 3_Supports_7_6: - h[7,6] - o[7,6] - s[7,6] <= -1
 3_DoesntSupport_7_6: h[7,6] + o[7,6] + 2 s[7,6] <= 2
 beta1_7_6: - s[7,6] + beta1[7,6] <= 0
 beta2_7_6: - s[7,6] + beta2[7,6] <= 0
 sBeta1Flag_7_6: - s[7,6] + beta1[7,6] <= 0
 sBeta2Flag_7_6: - s[7,6] + beta2[7,6] <= 0
 Eta1Overlap_7_6: eta1[7,6] + beta1[7,6] <= 1
 Eta2Overlap_7_6: eta2[7,6] + beta2[7,6] <= 1
 Eta1Flag_7_6: xi_l[6] - xi_l[7] - 300 eta1[7,6] <= 0
 Eta2Flag_7_6: - xi_r[6] + xi_r[7] - 300 eta2[7,6] <= 0
 6_Supports_7_InSameBin1_0: - p_ij[6,0] + p_ij[7,0] - s[7,6] <= 0
 6_Supports_7_InSameBin2_0: p_ij[6,0] - p_ij[7,0] - s[7,6] <= 0
 Gamma0ForItems_7_inBins_0_NoCut_6: p_ij[7,0] + gamma[7] <= 1
 6_Supports_7_InSameBin1_1: - p_ij[6,1] + p_ij[7,1] - s[7,6] <= 0
 6_Supports_7_InSameBin2_1: p_ij[6,1] - p_ij[7,1] - s[7,6] <= 0
 Gamma0ForItems_7_inBins_1_NoCut_6: p_ij[7,1] + gamma[7] <= 1
 6_Supports_7_InSameBin1_2: - p_ij[6,2] + p_ij[7,2] - s[7,6] <= 0
 6_Supports_7_InSameBin2_2: p_ij[6,2] - p_ij[7,2] - s[7,6] <= 0
 Gamma0ForItems_7_inBins_2_NoCut_6: p_ij[7,2] + gamma[7] <= 98
 6_Supports_7_InSameBin1_3: - p_ij[6,3] + p_ij[7,3] - s[7,6] <= 0
 6_Supports_7_InSameBin2_3: p_ij[6,3] - p_ij[7,3] - s[7,6] <= 0
 Gamma0ForItems_7_inBins_3_NoCut_6: p_ij[7,3] + gamma[7] <= 98
 7_OnCutIn2_6: 10000 p_ij[7,2] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 7_OnCutIn3_6: 10000 p_ij[7,3] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 7_OnCutIn2_7: 10000 p_ij[7,2] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 7_OnCutIn3_7: 10000 p_ij[7,3] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 AbsZ17_8: - zi[7] + zi_hi[8] - k[7,8] <= 0
 AbsZ27_8: zi[7] - zi_hi[8] - k[7,8] <= 0
 8_Below_7: zi[7] - zi_hi[8] + k[7,8] + 310 m[7,8] <= 310
 8_NotBelow_7: - zi[7] + zi_hi[8] + k[7,8] - 310 m[7,8] <= 0
 HeightCompat_7_8: h[7,8] - k[7,8] <= 0
 NoHeightCompat_7_8: - 155 h[7,8] + k[7,8] <= 0
 Overlap_7_8: - x_p[7,8] - x_p[8,7] + o[7,8] <= 0
 NoOverlap_7_8: x_p[7,8] + x_p[8,7] - 2 o[7,8] <= 0
 3_Supports_7_8: - h[7,8] - o[7,8] - s[7,8] <= -1
 3_DoesntSupport_7_8: h[7,8] + o[7,8] + 2 s[7,8] <= 2
 beta1_7_8: - s[7,8] + beta1[7,8] <= 0
 beta2_7_8: - s[7,8] + beta2[7,8] <= 0
 sBeta1Flag_7_8: - s[7,8] + beta1[7,8] <= 0
 sBeta2Flag_7_8: - s[7,8] + beta2[7,8] <= 0
 Eta1Overlap_7_8: eta1[7,8] + beta1[7,8] <= 1
 Eta2Overlap_7_8: eta2[7,8] + beta2[7,8] <= 1
 Eta1Flag_7_8: - xi_l[7] + xi_l[8] - 300 eta1[7,8] <= 0
 Eta2Flag_7_8: xi_r[7] - xi_r[8] - 300 eta2[7,8] <= 0
 8_Supports_7_InSameBin1_0: p_ij[7,0] - p_ij[8,0] - s[7,8] <= 0
 8_Supports_7_InSameBin2_0: - p_ij[7,0] + p_ij[8,0] - s[7,8] <= 0
 Gamma0ForItems_7_inBins_0_NoCut_8: p_ij[7,0] + gamma[7] <= 1
 8_Supports_7_InSameBin1_1: p_ij[7,1] - p_ij[8,1] - s[7,8] <= 0
 8_Supports_7_InSameBin2_1: - p_ij[7,1] + p_ij[8,1] - s[7,8] <= 0
 Gamma0ForItems_7_inBins_1_NoCut_8: p_ij[7,1] + gamma[7] <= 1
 8_Supports_7_InSameBin1_2: p_ij[7,2] - p_ij[8,2] - s[7,8] <= 0
 8_Supports_7_InSameBin2_2: - p_ij[7,2] + p_ij[8,2] - s[7,8] <= 0
 Gamma0ForItems_7_inBins_2_NoCut_8: p_ij[7,2] + gamma[7] <= 98
 8_Supports_7_InSameBin1_3: p_ij[7,3] - p_ij[8,3] - s[7,8] <= 0
 8_Supports_7_InSameBin2_3: - p_ij[7,3] + p_ij[8,3] - s[7,8] <= 0
 Gamma0ForItems_7_inBins_3_NoCut_8: p_ij[7,3] + gamma[7] <= 98
 7_OnCutIn2_8: 10000 p_ij[7,2] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 7_OnCutIn3_8: 10000 p_ij[7,3] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 AbsZ17_9: - zi[7] + zi_hi[9] - k[7,9] <= 0
 AbsZ27_9: zi[7] - zi_hi[9] - k[7,9] <= 0
 9_Below_7: zi[7] - zi_hi[9] + k[7,9] + 310 m[7,9] <= 310
 9_NotBelow_7: - zi[7] + zi_hi[9] + k[7,9] - 310 m[7,9] <= 0
 HeightCompat_7_9: h[7,9] - k[7,9] <= 0
 NoHeightCompat_7_9: - 155 h[7,9] + k[7,9] <= 0
 Overlap_7_9: - x_p[7,9] - x_p[9,7] + o[7,9] <= 0
 NoOverlap_7_9: x_p[7,9] + x_p[9,7] - 2 o[7,9] <= 0
 3_Supports_7_9: - h[7,9] - o[7,9] - s[7,9] <= -1
 3_DoesntSupport_7_9: h[7,9] + o[7,9] + 2 s[7,9] <= 2
 beta1_7_9: - s[7,9] + beta1[7,9] <= 0
 beta2_7_9: - s[7,9] + beta2[7,9] <= 0
 sBeta1Flag_7_9: - s[7,9] + beta1[7,9] <= 0
 sBeta2Flag_7_9: - s[7,9] + beta2[7,9] <= 0
 Eta1Overlap_7_9: eta1[7,9] + beta1[7,9] <= 1
 Eta2Overlap_7_9: eta2[7,9] + beta2[7,9] <= 1
 Eta1Flag_7_9: - xi_l[7] + xi_l[9] - 300 eta1[7,9] <= 0
 Eta2Flag_7_9: xi_r[7] - xi_r[9] - 300 eta2[7,9] <= 0
 9_Supports_7_InSameBin1_0: p_ij[7,0] - p_ij[9,0] - s[7,9] <= 0
 9_Supports_7_InSameBin2_0: - p_ij[7,0] + p_ij[9,0] - s[7,9] <= 0
 Gamma0ForItems_7_inBins_0_NoCut_9: p_ij[7,0] + gamma[7] <= 1
 9_Supports_7_InSameBin1_1: p_ij[7,1] - p_ij[9,1] - s[7,9] <= 0
 9_Supports_7_InSameBin2_1: - p_ij[7,1] + p_ij[9,1] - s[7,9] <= 0
 Gamma0ForItems_7_inBins_1_NoCut_9: p_ij[7,1] + gamma[7] <= 1
 9_Supports_7_InSameBin1_2: p_ij[7,2] - p_ij[9,2] - s[7,9] <= 0
 9_Supports_7_InSameBin2_2: - p_ij[7,2] + p_ij[9,2] - s[7,9] <= 0
 Gamma0ForItems_7_inBins_2_NoCut_9: p_ij[7,2] + gamma[7] <= 98
 9_Supports_7_InSameBin1_3: p_ij[7,3] - p_ij[9,3] - s[7,9] <= 0
 9_Supports_7_InSameBin2_3: - p_ij[7,3] + p_ij[9,3] - s[7,9] <= 0
 Gamma0ForItems_7_inBins_3_NoCut_9: p_ij[7,3] + gamma[7] <= 98
 7_OnCutIn2_9: 10000 p_ij[7,2] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 7_OnCutIn3_9: 10000 p_ij[7,3] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 AbsZ17_10: - zi[7] + zi_hi[10] - k[7,10] <= 0
 AbsZ27_10: zi[7] - zi_hi[10] - k[7,10] <= 0
 10_Below_7: zi[7] - zi_hi[10] + k[7,10] + 310 m[7,10] <= 310
 10_NotBelow_7: - zi[7] + zi_hi[10] + k[7,10] - 310 m[7,10] <= 0
 HeightCompat_7_10: h[7,10] - k[7,10] <= 0
 NoHeightCompat_7_10: - 155 h[7,10] + k[7,10] <= 0
 Overlap_7_10: - x_p[7,10] - x_p[10,7] + o[7,10] <= 0
 NoOverlap_7_10: x_p[7,10] + x_p[10,7] - 2 o[7,10] <= 0
 3_Supports_7_10: - h[7,10] - o[7,10] - s[7,10] <= -1
 3_DoesntSupport_7_10: h[7,10] + o[7,10] + 2 s[7,10] <= 2
 beta1_7_10: - s[7,10] + beta1[7,10] <= 0
 beta2_7_10: - s[7,10] + beta2[7,10] <= 0
 sBeta1Flag_7_10: - s[7,10] + beta1[7,10] <= 0
 sBeta2Flag_7_10: - s[7,10] + beta2[7,10] <= 0
 Eta1Overlap_7_10: eta1[7,10] + beta1[7,10] <= 1
 Eta2Overlap_7_10: eta2[7,10] + beta2[7,10] <= 1
 Eta1Flag_7_10: - xi_l[7] + xi_l[10] - 300 eta1[7,10] <= 0
 Eta2Flag_7_10: xi_r[7] - xi_r[10] - 300 eta2[7,10] <= 0
 10_Supports_7_InSameBin1_0: p_ij[7,0] - p_ij[10,0] - s[7,10] <= 0
 10_Supports_7_InSameBin2_0: - p_ij[7,0] + p_ij[10,0] - s[7,10] <= 0
 Gamma0ForItems_7_inBins_0_NoCut_10: p_ij[7,0] + gamma[7] <= 1
 10_Supports_7_InSameBin1_1: p_ij[7,1] - p_ij[10,1] - s[7,10] <= 0
 10_Supports_7_InSameBin2_1: - p_ij[7,1] + p_ij[10,1] - s[7,10] <= 0
 Gamma0ForItems_7_inBins_1_NoCut_10: p_ij[7,1] + gamma[7] <= 1
 10_Supports_7_InSameBin1_2: p_ij[7,2] - p_ij[10,2] - s[7,10] <= 0
 10_Supports_7_InSameBin2_2: - p_ij[7,2] + p_ij[10,2] - s[7,10] <= 0
 Gamma0ForItems_7_inBins_2_NoCut_10: p_ij[7,2] + gamma[7] <= 98
 10_Supports_7_InSameBin1_3: p_ij[7,3] - p_ij[10,3] - s[7,10] <= 0
 10_Supports_7_InSameBin2_3: - p_ij[7,3] + p_ij[10,3] - s[7,10] <= 0
 Gamma0ForItems_7_inBins_3_NoCut_10: p_ij[7,3] + gamma[7] <= 98
 7_OnCutIn2_10: 10000 p_ij[7,2] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 7_OnCutIn3_10: 10000 p_ij[7,3] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 AbsZ17_11: - zi[7] + zi_hi[11] - k[7,11] <= 0
 AbsZ27_11: zi[7] - zi_hi[11] - k[7,11] <= 0
 11_Below_7: zi[7] - zi_hi[11] + k[7,11] + 310 m[7,11] <= 310
 11_NotBelow_7: - zi[7] + zi_hi[11] + k[7,11] - 310 m[7,11] <= 0
 HeightCompat_7_11: h[7,11] - k[7,11] <= 0
 NoHeightCompat_7_11: - 155 h[7,11] + k[7,11] <= 0
 Overlap_7_11: - x_p[7,11] - x_p[11,7] + o[7,11] <= 0
 NoOverlap_7_11: x_p[7,11] + x_p[11,7] - 2 o[7,11] <= 0
 3_Supports_7_11: - h[7,11] - o[7,11] - s[7,11] <= -1
 3_DoesntSupport_7_11: h[7,11] + o[7,11] + 2 s[7,11] <= 2
 beta1_7_11: - s[7,11] + beta1[7,11] <= 0
 beta2_7_11: - s[7,11] + beta2[7,11] <= 0
 sBeta1Flag_7_11: - s[7,11] + beta1[7,11] <= 0
 sBeta2Flag_7_11: - s[7,11] + beta2[7,11] <= 0
 Eta1Overlap_7_11: eta1[7,11] + beta1[7,11] <= 1
 Eta2Overlap_7_11: eta2[7,11] + beta2[7,11] <= 1
 Eta1Flag_7_11: - xi_l[7] + xi_l[11] - 300 eta1[7,11] <= 0
 Eta2Flag_7_11: xi_r[7] - xi_r[11] - 300 eta2[7,11] <= 0
 11_Supports_7_InSameBin1_0: p_ij[7,0] - p_ij[11,0] - s[7,11] <= 0
 11_Supports_7_InSameBin2_0: - p_ij[7,0] + p_ij[11,0] - s[7,11] <= 0
 Gamma0ForItems_7_inBins_0_NoCut_11: p_ij[7,0] + gamma[7] <= 1
 11_Supports_7_InSameBin1_1: p_ij[7,1] - p_ij[11,1] - s[7,11] <= 0
 11_Supports_7_InSameBin2_1: - p_ij[7,1] + p_ij[11,1] - s[7,11] <= 0
 Gamma0ForItems_7_inBins_1_NoCut_11: p_ij[7,1] + gamma[7] <= 1
 11_Supports_7_InSameBin1_2: p_ij[7,2] - p_ij[11,2] - s[7,11] <= 0
 11_Supports_7_InSameBin2_2: - p_ij[7,2] + p_ij[11,2] - s[7,11] <= 0
 Gamma0ForItems_7_inBins_2_NoCut_11: p_ij[7,2] + gamma[7] <= 98
 11_Supports_7_InSameBin1_3: p_ij[7,3] - p_ij[11,3] - s[7,11] <= 0
 11_Supports_7_InSameBin2_3: - p_ij[7,3] + p_ij[11,3] - s[7,11] <= 0
 Gamma0ForItems_7_inBins_3_NoCut_11: p_ij[7,3] + gamma[7] <= 98
 7_OnCutIn2_11: 10000 p_ij[7,2] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 7_OnCutIn3_11: 10000 p_ij[7,3] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 AbsZ17_12: - zi[7] + zi_hi[12] - k[7,12] <= 0
 AbsZ27_12: zi[7] - zi_hi[12] - k[7,12] <= 0
 12_Below_7: zi[7] - zi_hi[12] + k[7,12] + 310 m[7,12] <= 310
 12_NotBelow_7: - zi[7] + zi_hi[12] + k[7,12] - 310 m[7,12] <= 0
 HeightCompat_7_12: h[7,12] - k[7,12] <= 0
 NoHeightCompat_7_12: - 155 h[7,12] + k[7,12] <= 0
 Overlap_7_12: - x_p[7,12] - x_p[12,7] + o[7,12] <= 0
 NoOverlap_7_12: x_p[7,12] + x_p[12,7] - 2 o[7,12] <= 0
 3_Supports_7_12: - h[7,12] - o[7,12] - s[7,12] <= -1
 3_DoesntSupport_7_12: h[7,12] + o[7,12] + 2 s[7,12] <= 2
 beta1_7_12: - s[7,12] + beta1[7,12] <= 0
 beta2_7_12: - s[7,12] + beta2[7,12] <= 0
 sBeta1Flag_7_12: - s[7,12] + beta1[7,12] <= 0
 sBeta2Flag_7_12: - s[7,12] + beta2[7,12] <= 0
 Eta1Overlap_7_12: eta1[7,12] + beta1[7,12] <= 1
 Eta2Overlap_7_12: eta2[7,12] + beta2[7,12] <= 1
 Eta1Flag_7_12: - xi_l[7] + xi_l[12] - 300 eta1[7,12] <= 0
 Eta2Flag_7_12: xi_r[7] - xi_r[12] - 300 eta2[7,12] <= 0
 12_Supports_7_InSameBin1_0: p_ij[7,0] - p_ij[12,0] - s[7,12] <= 0
 12_Supports_7_InSameBin2_0: - p_ij[7,0] + p_ij[12,0] - s[7,12] <= 0
 Gamma0ForItems_7_inBins_0_NoCut_12: p_ij[7,0] + gamma[7] <= 1
 12_Supports_7_InSameBin1_1: p_ij[7,1] - p_ij[12,1] - s[7,12] <= 0
 12_Supports_7_InSameBin2_1: - p_ij[7,1] + p_ij[12,1] - s[7,12] <= 0
 Gamma0ForItems_7_inBins_1_NoCut_12: p_ij[7,1] + gamma[7] <= 1
 12_Supports_7_InSameBin1_2: p_ij[7,2] - p_ij[12,2] - s[7,12] <= 0
 12_Supports_7_InSameBin2_2: - p_ij[7,2] + p_ij[12,2] - s[7,12] <= 0
 Gamma0ForItems_7_inBins_2_NoCut_12: p_ij[7,2] + gamma[7] <= 98
 12_Supports_7_InSameBin1_3: p_ij[7,3] - p_ij[12,3] - s[7,12] <= 0
 12_Supports_7_InSameBin2_3: - p_ij[7,3] + p_ij[12,3] - s[7,12] <= 0
 Gamma0ForItems_7_inBins_3_NoCut_12: p_ij[7,3] + gamma[7] <= 98
 7_OnCutIn2_12: 10000 p_ij[7,2] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 7_OnCutIn3_12: 10000 p_ij[7,3] + 1.261904761904762 xi_l[7] + zi[7]
   + 10000 gamma[7] <= 20053
 GroundConstraint_8: zi[8] + 155 g[8] <= 155
 AbsZ18_0: - zi[8] + zi_hi[0] - k[8,0] <= 0
 AbsZ28_0: zi[8] - zi_hi[0] - k[8,0] <= 0
 0_Below_8: zi[8] - zi_hi[0] + k[8,0] + 310 m[8,0] <= 310
 0_NotBelow_8: - zi[8] + zi_hi[0] + k[8,0] - 310 m[8,0] <= 0
 HeightCompat_8_0: h[8,0] - k[8,0] <= 0
 NoHeightCompat_8_0: - 155 h[8,0] + k[8,0] <= 0
 Overlap_8_0: - x_p[0,8] - x_p[8,0] + o[8,0] <= 0
 NoOverlap_8_0: x_p[0,8] + x_p[8,0] - 2 o[8,0] <= 0
 3_Supports_8_0: - h[8,0] - o[8,0] - s[8,0] <= -1
 3_DoesntSupport_8_0: h[8,0] + o[8,0] + 2 s[8,0] <= 2
 beta1_8_0: - s[8,0] + beta1[8,0] <= 0
 beta2_8_0: - s[8,0] + beta2[8,0] <= 0
 sBeta1Flag_8_0: - s[8,0] + beta1[8,0] <= 0
 sBeta2Flag_8_0: - s[8,0] + beta2[8,0] <= 0
 Eta1Overlap_8_0: eta1[8,0] + beta1[8,0] <= 1
 Eta2Overlap_8_0: eta2[8,0] + beta2[8,0] <= 1
 Eta1Flag_8_0: xi_l[0] - xi_l[8] - 300 eta1[8,0] <= 0
 Eta2Flag_8_0: - xi_r[0] + xi_r[8] - 300 eta2[8,0] <= 0
 0_Supports_8_InSameBin1_0: - p_ij[0,0] + p_ij[8,0] - s[8,0] <= 0
 0_Supports_8_InSameBin2_0: p_ij[0,0] - p_ij[8,0] - s[8,0] <= 0
 Gamma0ForItems_8_inBins_0_NoCut_0: p_ij[8,0] + gamma[8] <= 1
 0_Supports_8_InSameBin1_1: - p_ij[0,1] + p_ij[8,1] - s[8,0] <= 0
 0_Supports_8_InSameBin2_1: p_ij[0,1] - p_ij[8,1] - s[8,0] <= 0
 Gamma0ForItems_8_inBins_1_NoCut_0: p_ij[8,1] + gamma[8] <= 1
 0_Supports_8_InSameBin1_2: - p_ij[0,2] + p_ij[8,2] - s[8,0] <= 0
 0_Supports_8_InSameBin2_2: p_ij[0,2] - p_ij[8,2] - s[8,0] <= 0
 Gamma0ForItems_8_inBins_2_NoCut_0: p_ij[8,2] + gamma[8] <= 98
 0_Supports_8_InSameBin1_3: - p_ij[0,3] + p_ij[8,3] - s[8,0] <= 0
 0_Supports_8_InSameBin2_3: p_ij[0,3] - p_ij[8,3] - s[8,0] <= 0
 Gamma0ForItems_8_inBins_3_NoCut_0: p_ij[8,3] + gamma[8] <= 98
 8_OnCutIn2_0: 10000 p_ij[8,2] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 8_OnCutIn3_0: 10000 p_ij[8,3] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 AbsZ18_1: - zi[8] + zi_hi[1] - k[8,1] <= 0
 AbsZ28_1: zi[8] - zi_hi[1] - k[8,1] <= 0
 1_Below_8: zi[8] - zi_hi[1] + k[8,1] + 310 m[8,1] <= 310
 1_NotBelow_8: - zi[8] + zi_hi[1] + k[8,1] - 310 m[8,1] <= 0
 HeightCompat_8_1: h[8,1] - k[8,1] <= 0
 NoHeightCompat_8_1: - 155 h[8,1] + k[8,1] <= 0
 Overlap_8_1: - x_p[1,8] - x_p[8,1] + o[8,1] <= 0
 NoOverlap_8_1: x_p[1,8] + x_p[8,1] - 2 o[8,1] <= 0
 3_Supports_8_1: - h[8,1] - o[8,1] - s[8,1] <= -1
 3_DoesntSupport_8_1: h[8,1] + o[8,1] + 2 s[8,1] <= 2
 beta1_8_1: - s[8,1] + beta1[8,1] <= 0
 beta2_8_1: - s[8,1] + beta2[8,1] <= 0
 sBeta1Flag_8_1: - s[8,1] + beta1[8,1] <= 0
 sBeta2Flag_8_1: - s[8,1] + beta2[8,1] <= 0
 Eta1Overlap_8_1: eta1[8,1] + beta1[8,1] <= 1
 Eta2Overlap_8_1: eta2[8,1] + beta2[8,1] <= 1
 Eta1Flag_8_1: xi_l[1] - xi_l[8] - 300 eta1[8,1] <= 0
 Eta2Flag_8_1: - xi_r[1] + xi_r[8] - 300 eta2[8,1] <= 0
 1_Supports_8_InSameBin1_0: - p_ij[1,0] + p_ij[8,0] - s[8,1] <= 0
 1_Supports_8_InSameBin2_0: p_ij[1,0] - p_ij[8,0] - s[8,1] <= 0
 Gamma0ForItems_8_inBins_0_NoCut_1: p_ij[8,0] + gamma[8] <= 1
 1_Supports_8_InSameBin1_1: - p_ij[1,1] + p_ij[8,1] - s[8,1] <= 0
 1_Supports_8_InSameBin2_1: p_ij[1,1] - p_ij[8,1] - s[8,1] <= 0
 Gamma0ForItems_8_inBins_1_NoCut_1: p_ij[8,1] + gamma[8] <= 1
 1_Supports_8_InSameBin1_2: - p_ij[1,2] + p_ij[8,2] - s[8,1] <= 0
 1_Supports_8_InSameBin2_2: p_ij[1,2] - p_ij[8,2] - s[8,1] <= 0
 Gamma0ForItems_8_inBins_2_NoCut_1: p_ij[8,2] + gamma[8] <= 98
 1_Supports_8_InSameBin1_3: - p_ij[1,3] + p_ij[8,3] - s[8,1] <= 0
 1_Supports_8_InSameBin2_3: p_ij[1,3] - p_ij[8,3] - s[8,1] <= 0
 Gamma0ForItems_8_inBins_3_NoCut_1: p_ij[8,3] + gamma[8] <= 98
 8_OnCutIn2_1: 10000 p_ij[8,2] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 8_OnCutIn3_1: 10000 p_ij[8,3] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 AbsZ18_2: - zi[8] + zi_hi[2] - k[8,2] <= 0
 AbsZ28_2: zi[8] - zi_hi[2] - k[8,2] <= 0
 2_Below_8: zi[8] - zi_hi[2] + k[8,2] + 310 m[8,2] <= 310
 2_NotBelow_8: - zi[8] + zi_hi[2] + k[8,2] - 310 m[8,2] <= 0
 HeightCompat_8_2: h[8,2] - k[8,2] <= 0
 NoHeightCompat_8_2: - 155 h[8,2] + k[8,2] <= 0
 Overlap_8_2: - x_p[2,8] - x_p[8,2] + o[8,2] <= 0
 NoOverlap_8_2: x_p[2,8] + x_p[8,2] - 2 o[8,2] <= 0
 3_Supports_8_2: - h[8,2] - o[8,2] - s[8,2] <= -1
 3_DoesntSupport_8_2: h[8,2] + o[8,2] + 2 s[8,2] <= 2
 beta1_8_2: - s[8,2] + beta1[8,2] <= 0
 beta2_8_2: - s[8,2] + beta2[8,2] <= 0
 sBeta1Flag_8_2: - s[8,2] + beta1[8,2] <= 0
 sBeta2Flag_8_2: - s[8,2] + beta2[8,2] <= 0
 Eta1Overlap_8_2: eta1[8,2] + beta1[8,2] <= 1
 Eta2Overlap_8_2: eta2[8,2] + beta2[8,2] <= 1
 Eta1Flag_8_2: xi_l[2] - xi_l[8] - 300 eta1[8,2] <= 0
 Eta2Flag_8_2: - xi_r[2] + xi_r[8] - 300 eta2[8,2] <= 0
 2_Supports_8_InSameBin1_0: - p_ij[2,0] + p_ij[8,0] - s[8,2] <= 0
 2_Supports_8_InSameBin2_0: p_ij[2,0] - p_ij[8,0] - s[8,2] <= 0
 Gamma0ForItems_8_inBins_0_NoCut_2: p_ij[8,0] + gamma[8] <= 1
 2_Supports_8_InSameBin1_1: - p_ij[2,1] + p_ij[8,1] - s[8,2] <= 0
 2_Supports_8_InSameBin2_1: p_ij[2,1] - p_ij[8,1] - s[8,2] <= 0
 Gamma0ForItems_8_inBins_1_NoCut_2: p_ij[8,1] + gamma[8] <= 1
 2_Supports_8_InSameBin1_2: - p_ij[2,2] + p_ij[8,2] - s[8,2] <= 0
 2_Supports_8_InSameBin2_2: p_ij[2,2] - p_ij[8,2] - s[8,2] <= 0
 Gamma0ForItems_8_inBins_2_NoCut_2: p_ij[8,2] + gamma[8] <= 98
 2_Supports_8_InSameBin1_3: - p_ij[2,3] + p_ij[8,3] - s[8,2] <= 0
 2_Supports_8_InSameBin2_3: p_ij[2,3] - p_ij[8,3] - s[8,2] <= 0
 Gamma0ForItems_8_inBins_3_NoCut_2: p_ij[8,3] + gamma[8] <= 98
 8_OnCutIn2_2: 10000 p_ij[8,2] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 8_OnCutIn3_2: 10000 p_ij[8,3] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 AbsZ18_3: - zi[8] + zi_hi[3] - k[8,3] <= 0
 AbsZ28_3: zi[8] - zi_hi[3] - k[8,3] <= 0
 3_Below_8: zi[8] - zi_hi[3] + k[8,3] + 310 m[8,3] <= 310
 3_NotBelow_8: - zi[8] + zi_hi[3] + k[8,3] - 310 m[8,3] <= 0
 HeightCompat_8_3: h[8,3] - k[8,3] <= 0
 NoHeightCompat_8_3: - 155 h[8,3] + k[8,3] <= 0
 Overlap_8_3: - x_p[3,8] - x_p[8,3] + o[8,3] <= 0
 NoOverlap_8_3: x_p[3,8] + x_p[8,3] - 2 o[8,3] <= 0
 3_Supports_8_3: - h[8,3] - o[8,3] - s[8,3] <= -1
 3_DoesntSupport_8_3: h[8,3] + o[8,3] + 2 s[8,3] <= 2
 beta1_8_3: - s[8,3] + beta1[8,3] <= 0
 beta2_8_3: - s[8,3] + beta2[8,3] <= 0
 sBeta1Flag_8_3: - s[8,3] + beta1[8,3] <= 0
 sBeta2Flag_8_3: - s[8,3] + beta2[8,3] <= 0
 Eta1Overlap_8_3: eta1[8,3] + beta1[8,3] <= 1
 Eta2Overlap_8_3: eta2[8,3] + beta2[8,3] <= 1
 Eta1Flag_8_3: xi_l[3] - xi_l[8] - 300 eta1[8,3] <= 0
 Eta2Flag_8_3: - xi_r[3] + xi_r[8] - 300 eta2[8,3] <= 0
 3_Supports_8_InSameBin1_0: - p_ij[3,0] + p_ij[8,0] - s[8,3] <= 0
 3_Supports_8_InSameBin2_0: p_ij[3,0] - p_ij[8,0] - s[8,3] <= 0
 Gamma0ForItems_8_inBins_0_NoCut_3: p_ij[8,0] + gamma[8] <= 1
 3_Supports_8_InSameBin1_1: - p_ij[3,1] + p_ij[8,1] - s[8,3] <= 0
 3_Supports_8_InSameBin2_1: p_ij[3,1] - p_ij[8,1] - s[8,3] <= 0
 Gamma0ForItems_8_inBins_1_NoCut_3: p_ij[8,1] + gamma[8] <= 1
 3_Supports_8_InSameBin1_2: - p_ij[3,2] + p_ij[8,2] - s[8,3] <= 0
 3_Supports_8_InSameBin2_2: p_ij[3,2] - p_ij[8,2] - s[8,3] <= 0
 Gamma0ForItems_8_inBins_2_NoCut_3: p_ij[8,2] + gamma[8] <= 98
 3_Supports_8_InSameBin1_3: - p_ij[3,3] + p_ij[8,3] - s[8,3] <= 0
 3_Supports_8_InSameBin2_3: p_ij[3,3] - p_ij[8,3] - s[8,3] <= 0
 Gamma0ForItems_8_inBins_3_NoCut_3: p_ij[8,3] + gamma[8] <= 98
 8_OnCutIn2_3: 10000 p_ij[8,2] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 8_OnCutIn3_3: 10000 p_ij[8,3] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 AbsZ18_4: - zi[8] + zi_hi[4] - k[8,4] <= 0
 AbsZ28_4: zi[8] - zi_hi[4] - k[8,4] <= 0
 4_Below_8: zi[8] - zi_hi[4] + k[8,4] + 310 m[8,4] <= 310
 4_NotBelow_8: - zi[8] + zi_hi[4] + k[8,4] - 310 m[8,4] <= 0
 HeightCompat_8_4: h[8,4] - k[8,4] <= 0
 NoHeightCompat_8_4: - 155 h[8,4] + k[8,4] <= 0
 Overlap_8_4: - x_p[4,8] - x_p[8,4] + o[8,4] <= 0
 NoOverlap_8_4: x_p[4,8] + x_p[8,4] - 2 o[8,4] <= 0
 3_Supports_8_4: - h[8,4] - o[8,4] - s[8,4] <= -1
 3_DoesntSupport_8_4: h[8,4] + o[8,4] + 2 s[8,4] <= 2
 beta1_8_4: - s[8,4] + beta1[8,4] <= 0
 beta2_8_4: - s[8,4] + beta2[8,4] <= 0
 sBeta1Flag_8_4: - s[8,4] + beta1[8,4] <= 0
 sBeta2Flag_8_4: - s[8,4] + beta2[8,4] <= 0
 Eta1Overlap_8_4: eta1[8,4] + beta1[8,4] <= 1
 Eta2Overlap_8_4: eta2[8,4] + beta2[8,4] <= 1
 Eta1Flag_8_4: xi_l[4] - xi_l[8] - 300 eta1[8,4] <= 0
 Eta2Flag_8_4: - xi_r[4] + xi_r[8] - 300 eta2[8,4] <= 0
 4_Supports_8_InSameBin1_0: - p_ij[4,0] + p_ij[8,0] - s[8,4] <= 0
 4_Supports_8_InSameBin2_0: p_ij[4,0] - p_ij[8,0] - s[8,4] <= 0
 Gamma0ForItems_8_inBins_0_NoCut_4: p_ij[8,0] + gamma[8] <= 1
 4_Supports_8_InSameBin1_1: - p_ij[4,1] + p_ij[8,1] - s[8,4] <= 0
 4_Supports_8_InSameBin2_1: p_ij[4,1] - p_ij[8,1] - s[8,4] <= 0
 Gamma0ForItems_8_inBins_1_NoCut_4: p_ij[8,1] + gamma[8] <= 1
 4_Supports_8_InSameBin1_2: - p_ij[4,2] + p_ij[8,2] - s[8,4] <= 0
 4_Supports_8_InSameBin2_2: p_ij[4,2] - p_ij[8,2] - s[8,4] <= 0
 Gamma0ForItems_8_inBins_2_NoCut_4: p_ij[8,2] + gamma[8] <= 98
 4_Supports_8_InSameBin1_3: - p_ij[4,3] + p_ij[8,3] - s[8,4] <= 0
 4_Supports_8_InSameBin2_3: p_ij[4,3] - p_ij[8,3] - s[8,4] <= 0
 Gamma0ForItems_8_inBins_3_NoCut_4: p_ij[8,3] + gamma[8] <= 98
 8_OnCutIn2_4: 10000 p_ij[8,2] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 8_OnCutIn3_4: 10000 p_ij[8,3] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 AbsZ18_5: - zi[8] + zi_hi[5] - k[8,5] <= 0
 AbsZ28_5: zi[8] - zi_hi[5] - k[8,5] <= 0
 5_Below_8: zi[8] - zi_hi[5] + k[8,5] + 310 m[8,5] <= 310
 5_NotBelow_8: - zi[8] + zi_hi[5] + k[8,5] - 310 m[8,5] <= 0
 HeightCompat_8_5: h[8,5] - k[8,5] <= 0
 NoHeightCompat_8_5: - 155 h[8,5] + k[8,5] <= 0
 Overlap_8_5: - x_p[5,8] - x_p[8,5] + o[8,5] <= 0
 NoOverlap_8_5: x_p[5,8] + x_p[8,5] - 2 o[8,5] <= 0
 3_Supports_8_5: - h[8,5] - o[8,5] - s[8,5] <= -1
 3_DoesntSupport_8_5: h[8,5] + o[8,5] + 2 s[8,5] <= 2
 beta1_8_5: - s[8,5] + beta1[8,5] <= 0
 beta2_8_5: - s[8,5] + beta2[8,5] <= 0
 sBeta1Flag_8_5: - s[8,5] + beta1[8,5] <= 0
 sBeta2Flag_8_5: - s[8,5] + beta2[8,5] <= 0
 Eta1Overlap_8_5: eta1[8,5] + beta1[8,5] <= 1
 Eta2Overlap_8_5: eta2[8,5] + beta2[8,5] <= 1
 Eta1Flag_8_5: xi_l[5] - xi_l[8] - 300 eta1[8,5] <= 0
 Eta2Flag_8_5: - xi_r[5] + xi_r[8] - 300 eta2[8,5] <= 0
 5_Supports_8_InSameBin1_0: - p_ij[5,0] + p_ij[8,0] - s[8,5] <= 0
 5_Supports_8_InSameBin2_0: p_ij[5,0] - p_ij[8,0] - s[8,5] <= 0
 Gamma0ForItems_8_inBins_0_NoCut_5: p_ij[8,0] + gamma[8] <= 1
 5_Supports_8_InSameBin1_1: - p_ij[5,1] + p_ij[8,1] - s[8,5] <= 0
 5_Supports_8_InSameBin2_1: p_ij[5,1] - p_ij[8,1] - s[8,5] <= 0
 Gamma0ForItems_8_inBins_1_NoCut_5: p_ij[8,1] + gamma[8] <= 1
 5_Supports_8_InSameBin1_2: - p_ij[5,2] + p_ij[8,2] - s[8,5] <= 0
 5_Supports_8_InSameBin2_2: p_ij[5,2] - p_ij[8,2] - s[8,5] <= 0
 Gamma0ForItems_8_inBins_2_NoCut_5: p_ij[8,2] + gamma[8] <= 98
 5_Supports_8_InSameBin1_3: - p_ij[5,3] + p_ij[8,3] - s[8,5] <= 0
 5_Supports_8_InSameBin2_3: p_ij[5,3] - p_ij[8,3] - s[8,5] <= 0
 Gamma0ForItems_8_inBins_3_NoCut_5: p_ij[8,3] + gamma[8] <= 98
 8_OnCutIn2_5: 10000 p_ij[8,2] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 8_OnCutIn3_5: 10000 p_ij[8,3] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 AbsZ18_6: - zi[8] + zi_hi[6] - k[8,6] <= 0
 AbsZ28_6: zi[8] - zi_hi[6] - k[8,6] <= 0
 6_Below_8: zi[8] - zi_hi[6] + k[8,6] + 310 m[8,6] <= 310
 6_NotBelow_8: - zi[8] + zi_hi[6] + k[8,6] - 310 m[8,6] <= 0
 HeightCompat_8_6: h[8,6] - k[8,6] <= 0
 NoHeightCompat_8_6: - 155 h[8,6] + k[8,6] <= 0
 Overlap_8_6: - x_p[6,8] - x_p[8,6] + o[8,6] <= 0
 NoOverlap_8_6: x_p[6,8] + x_p[8,6] - 2 o[8,6] <= 0
 3_Supports_8_6: - h[8,6] - o[8,6] - s[8,6] <= -1
 3_DoesntSupport_8_6: h[8,6] + o[8,6] + 2 s[8,6] <= 2
 beta1_8_6: - s[8,6] + beta1[8,6] <= 0
 beta2_8_6: - s[8,6] + beta2[8,6] <= 0
 sBeta1Flag_8_6: - s[8,6] + beta1[8,6] <= 0
 sBeta2Flag_8_6: - s[8,6] + beta2[8,6] <= 0
 Eta1Overlap_8_6: eta1[8,6] + beta1[8,6] <= 1
 Eta2Overlap_8_6: eta2[8,6] + beta2[8,6] <= 1
 Eta1Flag_8_6: xi_l[6] - xi_l[8] - 300 eta1[8,6] <= 0
 Eta2Flag_8_6: - xi_r[6] + xi_r[8] - 300 eta2[8,6] <= 0
 6_Supports_8_InSameBin1_0: - p_ij[6,0] + p_ij[8,0] - s[8,6] <= 0
 6_Supports_8_InSameBin2_0: p_ij[6,0] - p_ij[8,0] - s[8,6] <= 0
 Gamma0ForItems_8_inBins_0_NoCut_6: p_ij[8,0] + gamma[8] <= 1
 6_Supports_8_InSameBin1_1: - p_ij[6,1] + p_ij[8,1] - s[8,6] <= 0
 6_Supports_8_InSameBin2_1: p_ij[6,1] - p_ij[8,1] - s[8,6] <= 0
 Gamma0ForItems_8_inBins_1_NoCut_6: p_ij[8,1] + gamma[8] <= 1
 6_Supports_8_InSameBin1_2: - p_ij[6,2] + p_ij[8,2] - s[8,6] <= 0
 6_Supports_8_InSameBin2_2: p_ij[6,2] - p_ij[8,2] - s[8,6] <= 0
 Gamma0ForItems_8_inBins_2_NoCut_6: p_ij[8,2] + gamma[8] <= 98
 6_Supports_8_InSameBin1_3: - p_ij[6,3] + p_ij[8,3] - s[8,6] <= 0
 6_Supports_8_InSameBin2_3: p_ij[6,3] - p_ij[8,3] - s[8,6] <= 0
 Gamma0ForItems_8_inBins_3_NoCut_6: p_ij[8,3] + gamma[8] <= 98
 8_OnCutIn2_6: 10000 p_ij[8,2] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 8_OnCutIn3_6: 10000 p_ij[8,3] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 AbsZ18_7: - zi[8] + zi_hi[7] - k[8,7] <= 0
 AbsZ28_7: zi[8] - zi_hi[7] - k[8,7] <= 0
 7_Below_8: zi[8] - zi_hi[7] + k[8,7] + 310 m[8,7] <= 310
 7_NotBelow_8: - zi[8] + zi_hi[7] + k[8,7] - 310 m[8,7] <= 0
 HeightCompat_8_7: h[8,7] - k[8,7] <= 0
 NoHeightCompat_8_7: - 155 h[8,7] + k[8,7] <= 0
 Overlap_8_7: - x_p[7,8] - x_p[8,7] + o[8,7] <= 0
 NoOverlap_8_7: x_p[7,8] + x_p[8,7] - 2 o[8,7] <= 0
 3_Supports_8_7: - h[8,7] - o[8,7] - s[8,7] <= -1
 3_DoesntSupport_8_7: h[8,7] + o[8,7] + 2 s[8,7] <= 2
 beta1_8_7: - s[8,7] + beta1[8,7] <= 0
 beta2_8_7: - s[8,7] + beta2[8,7] <= 0
 sBeta1Flag_8_7: - s[8,7] + beta1[8,7] <= 0
 sBeta2Flag_8_7: - s[8,7] + beta2[8,7] <= 0
 Eta1Overlap_8_7: eta1[8,7] + beta1[8,7] <= 1
 Eta2Overlap_8_7: eta2[8,7] + beta2[8,7] <= 1
 Eta1Flag_8_7: xi_l[7] - xi_l[8] - 300 eta1[8,7] <= 0
 Eta2Flag_8_7: - xi_r[7] + xi_r[8] - 300 eta2[8,7] <= 0
 7_Supports_8_InSameBin1_0: - p_ij[7,0] + p_ij[8,0] - s[8,7] <= 0
 7_Supports_8_InSameBin2_0: p_ij[7,0] - p_ij[8,0] - s[8,7] <= 0
 Gamma0ForItems_8_inBins_0_NoCut_7: p_ij[8,0] + gamma[8] <= 1
 7_Supports_8_InSameBin1_1: - p_ij[7,1] + p_ij[8,1] - s[8,7] <= 0
 7_Supports_8_InSameBin2_1: p_ij[7,1] - p_ij[8,1] - s[8,7] <= 0
 Gamma0ForItems_8_inBins_1_NoCut_7: p_ij[8,1] + gamma[8] <= 1
 7_Supports_8_InSameBin1_2: - p_ij[7,2] + p_ij[8,2] - s[8,7] <= 0
 7_Supports_8_InSameBin2_2: p_ij[7,2] - p_ij[8,2] - s[8,7] <= 0
 Gamma0ForItems_8_inBins_2_NoCut_7: p_ij[8,2] + gamma[8] <= 98
 7_Supports_8_InSameBin1_3: - p_ij[7,3] + p_ij[8,3] - s[8,7] <= 0
 7_Supports_8_InSameBin2_3: p_ij[7,3] - p_ij[8,3] - s[8,7] <= 0
 Gamma0ForItems_8_inBins_3_NoCut_7: p_ij[8,3] + gamma[8] <= 98
 8_OnCutIn2_7: 10000 p_ij[8,2] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 8_OnCutIn3_7: 10000 p_ij[8,3] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 8_OnCutIn2_8: 10000 p_ij[8,2] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 8_OnCutIn3_8: 10000 p_ij[8,3] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 AbsZ18_9: - zi[8] + zi_hi[9] - k[8,9] <= 0
 AbsZ28_9: zi[8] - zi_hi[9] - k[8,9] <= 0
 9_Below_8: zi[8] - zi_hi[9] + k[8,9] + 310 m[8,9] <= 310
 9_NotBelow_8: - zi[8] + zi_hi[9] + k[8,9] - 310 m[8,9] <= 0
 HeightCompat_8_9: h[8,9] - k[8,9] <= 0
 NoHeightCompat_8_9: - 155 h[8,9] + k[8,9] <= 0
 Overlap_8_9: - x_p[8,9] - x_p[9,8] + o[8,9] <= 0
 NoOverlap_8_9: x_p[8,9] + x_p[9,8] - 2 o[8,9] <= 0
 3_Supports_8_9: - h[8,9] - o[8,9] - s[8,9] <= -1
 3_DoesntSupport_8_9: h[8,9] + o[8,9] + 2 s[8,9] <= 2
 beta1_8_9: - s[8,9] + beta1[8,9] <= 0
 beta2_8_9: - s[8,9] + beta2[8,9] <= 0
 sBeta1Flag_8_9: - s[8,9] + beta1[8,9] <= 0
 sBeta2Flag_8_9: - s[8,9] + beta2[8,9] <= 0
 Eta1Overlap_8_9: eta1[8,9] + beta1[8,9] <= 1
 Eta2Overlap_8_9: eta2[8,9] + beta2[8,9] <= 1
 Eta1Flag_8_9: - xi_l[8] + xi_l[9] - 300 eta1[8,9] <= 0
 Eta2Flag_8_9: xi_r[8] - xi_r[9] - 300 eta2[8,9] <= 0
 9_Supports_8_InSameBin1_0: p_ij[8,0] - p_ij[9,0] - s[8,9] <= 0
 9_Supports_8_InSameBin2_0: - p_ij[8,0] + p_ij[9,0] - s[8,9] <= 0
 Gamma0ForItems_8_inBins_0_NoCut_9: p_ij[8,0] + gamma[8] <= 1
 9_Supports_8_InSameBin1_1: p_ij[8,1] - p_ij[9,1] - s[8,9] <= 0
 9_Supports_8_InSameBin2_1: - p_ij[8,1] + p_ij[9,1] - s[8,9] <= 0
 Gamma0ForItems_8_inBins_1_NoCut_9: p_ij[8,1] + gamma[8] <= 1
 9_Supports_8_InSameBin1_2: p_ij[8,2] - p_ij[9,2] - s[8,9] <= 0
 9_Supports_8_InSameBin2_2: - p_ij[8,2] + p_ij[9,2] - s[8,9] <= 0
 Gamma0ForItems_8_inBins_2_NoCut_9: p_ij[8,2] + gamma[8] <= 98
 9_Supports_8_InSameBin1_3: p_ij[8,3] - p_ij[9,3] - s[8,9] <= 0
 9_Supports_8_InSameBin2_3: - p_ij[8,3] + p_ij[9,3] - s[8,9] <= 0
 Gamma0ForItems_8_inBins_3_NoCut_9: p_ij[8,3] + gamma[8] <= 98
 8_OnCutIn2_9: 10000 p_ij[8,2] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 8_OnCutIn3_9: 10000 p_ij[8,3] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 AbsZ18_10: - zi[8] + zi_hi[10] - k[8,10] <= 0
 AbsZ28_10: zi[8] - zi_hi[10] - k[8,10] <= 0
 10_Below_8: zi[8] - zi_hi[10] + k[8,10] + 310 m[8,10] <= 310
 10_NotBelow_8: - zi[8] + zi_hi[10] + k[8,10] - 310 m[8,10] <= 0
 HeightCompat_8_10: h[8,10] - k[8,10] <= 0
 NoHeightCompat_8_10: - 155 h[8,10] + k[8,10] <= 0
 Overlap_8_10: - x_p[8,10] - x_p[10,8] + o[8,10] <= 0
 NoOverlap_8_10: x_p[8,10] + x_p[10,8] - 2 o[8,10] <= 0
 3_Supports_8_10: - h[8,10] - o[8,10] - s[8,10] <= -1
 3_DoesntSupport_8_10: h[8,10] + o[8,10] + 2 s[8,10] <= 2
 beta1_8_10: - s[8,10] + beta1[8,10] <= 0
 beta2_8_10: - s[8,10] + beta2[8,10] <= 0
 sBeta1Flag_8_10: - s[8,10] + beta1[8,10] <= 0
 sBeta2Flag_8_10: - s[8,10] + beta2[8,10] <= 0
 Eta1Overlap_8_10: eta1[8,10] + beta1[8,10] <= 1
 Eta2Overlap_8_10: eta2[8,10] + beta2[8,10] <= 1
 Eta1Flag_8_10: - xi_l[8] + xi_l[10] - 300 eta1[8,10] <= 0
 Eta2Flag_8_10: xi_r[8] - xi_r[10] - 300 eta2[8,10] <= 0
 10_Supports_8_InSameBin1_0: p_ij[8,0] - p_ij[10,0] - s[8,10] <= 0
 10_Supports_8_InSameBin2_0: - p_ij[8,0] + p_ij[10,0] - s[8,10] <= 0
 Gamma0ForItems_8_inBins_0_NoCut_10: p_ij[8,0] + gamma[8] <= 1
 10_Supports_8_InSameBin1_1: p_ij[8,1] - p_ij[10,1] - s[8,10] <= 0
 10_Supports_8_InSameBin2_1: - p_ij[8,1] + p_ij[10,1] - s[8,10] <= 0
 Gamma0ForItems_8_inBins_1_NoCut_10: p_ij[8,1] + gamma[8] <= 1
 10_Supports_8_InSameBin1_2: p_ij[8,2] - p_ij[10,2] - s[8,10] <= 0
 10_Supports_8_InSameBin2_2: - p_ij[8,2] + p_ij[10,2] - s[8,10] <= 0
 Gamma0ForItems_8_inBins_2_NoCut_10: p_ij[8,2] + gamma[8] <= 98
 10_Supports_8_InSameBin1_3: p_ij[8,3] - p_ij[10,3] - s[8,10] <= 0
 10_Supports_8_InSameBin2_3: - p_ij[8,3] + p_ij[10,3] - s[8,10] <= 0
 Gamma0ForItems_8_inBins_3_NoCut_10: p_ij[8,3] + gamma[8] <= 98
 8_OnCutIn2_10: 10000 p_ij[8,2] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 8_OnCutIn3_10: 10000 p_ij[8,3] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 AbsZ18_11: - zi[8] + zi_hi[11] - k[8,11] <= 0
 AbsZ28_11: zi[8] - zi_hi[11] - k[8,11] <= 0
 11_Below_8: zi[8] - zi_hi[11] + k[8,11] + 310 m[8,11] <= 310
 11_NotBelow_8: - zi[8] + zi_hi[11] + k[8,11] - 310 m[8,11] <= 0
 HeightCompat_8_11: h[8,11] - k[8,11] <= 0
 NoHeightCompat_8_11: - 155 h[8,11] + k[8,11] <= 0
 Overlap_8_11: - x_p[8,11] - x_p[11,8] + o[8,11] <= 0
 NoOverlap_8_11: x_p[8,11] + x_p[11,8] - 2 o[8,11] <= 0
 3_Supports_8_11: - h[8,11] - o[8,11] - s[8,11] <= -1
 3_DoesntSupport_8_11: h[8,11] + o[8,11] + 2 s[8,11] <= 2
 beta1_8_11: - s[8,11] + beta1[8,11] <= 0
 beta2_8_11: - s[8,11] + beta2[8,11] <= 0
 sBeta1Flag_8_11: - s[8,11] + beta1[8,11] <= 0
 sBeta2Flag_8_11: - s[8,11] + beta2[8,11] <= 0
 Eta1Overlap_8_11: eta1[8,11] + beta1[8,11] <= 1
 Eta2Overlap_8_11: eta2[8,11] + beta2[8,11] <= 1
 Eta1Flag_8_11: - xi_l[8] + xi_l[11] - 300 eta1[8,11] <= 0
 Eta2Flag_8_11: xi_r[8] - xi_r[11] - 300 eta2[8,11] <= 0
 11_Supports_8_InSameBin1_0: p_ij[8,0] - p_ij[11,0] - s[8,11] <= 0
 11_Supports_8_InSameBin2_0: - p_ij[8,0] + p_ij[11,0] - s[8,11] <= 0
 Gamma0ForItems_8_inBins_0_NoCut_11: p_ij[8,0] + gamma[8] <= 1
 11_Supports_8_InSameBin1_1: p_ij[8,1] - p_ij[11,1] - s[8,11] <= 0
 11_Supports_8_InSameBin2_1: - p_ij[8,1] + p_ij[11,1] - s[8,11] <= 0
 Gamma0ForItems_8_inBins_1_NoCut_11: p_ij[8,1] + gamma[8] <= 1
 11_Supports_8_InSameBin1_2: p_ij[8,2] - p_ij[11,2] - s[8,11] <= 0
 11_Supports_8_InSameBin2_2: - p_ij[8,2] + p_ij[11,2] - s[8,11] <= 0
 Gamma0ForItems_8_inBins_2_NoCut_11: p_ij[8,2] + gamma[8] <= 98
 11_Supports_8_InSameBin1_3: p_ij[8,3] - p_ij[11,3] - s[8,11] <= 0
 11_Supports_8_InSameBin2_3: - p_ij[8,3] + p_ij[11,3] - s[8,11] <= 0
 Gamma0ForItems_8_inBins_3_NoCut_11: p_ij[8,3] + gamma[8] <= 98
 8_OnCutIn2_11: 10000 p_ij[8,2] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 8_OnCutIn3_11: 10000 p_ij[8,3] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 AbsZ18_12: - zi[8] + zi_hi[12] - k[8,12] <= 0
 AbsZ28_12: zi[8] - zi_hi[12] - k[8,12] <= 0
 12_Below_8: zi[8] - zi_hi[12] + k[8,12] + 310 m[8,12] <= 310
 12_NotBelow_8: - zi[8] + zi_hi[12] + k[8,12] - 310 m[8,12] <= 0
 HeightCompat_8_12: h[8,12] - k[8,12] <= 0
 NoHeightCompat_8_12: - 155 h[8,12] + k[8,12] <= 0
 Overlap_8_12: - x_p[8,12] - x_p[12,8] + o[8,12] <= 0
 NoOverlap_8_12: x_p[8,12] + x_p[12,8] - 2 o[8,12] <= 0
 3_Supports_8_12: - h[8,12] - o[8,12] - s[8,12] <= -1
 3_DoesntSupport_8_12: h[8,12] + o[8,12] + 2 s[8,12] <= 2
 beta1_8_12: - s[8,12] + beta1[8,12] <= 0
 beta2_8_12: - s[8,12] + beta2[8,12] <= 0
 sBeta1Flag_8_12: - s[8,12] + beta1[8,12] <= 0
 sBeta2Flag_8_12: - s[8,12] + beta2[8,12] <= 0
 Eta1Overlap_8_12: eta1[8,12] + beta1[8,12] <= 1
 Eta2Overlap_8_12: eta2[8,12] + beta2[8,12] <= 1
 Eta1Flag_8_12: - xi_l[8] + xi_l[12] - 300 eta1[8,12] <= 0
 Eta2Flag_8_12: xi_r[8] - xi_r[12] - 300 eta2[8,12] <= 0
 12_Supports_8_InSameBin1_0: p_ij[8,0] - p_ij[12,0] - s[8,12] <= 0
 12_Supports_8_InSameBin2_0: - p_ij[8,0] + p_ij[12,0] - s[8,12] <= 0
 Gamma0ForItems_8_inBins_0_NoCut_12: p_ij[8,0] + gamma[8] <= 1
 12_Supports_8_InSameBin1_1: p_ij[8,1] - p_ij[12,1] - s[8,12] <= 0
 12_Supports_8_InSameBin2_1: - p_ij[8,1] + p_ij[12,1] - s[8,12] <= 0
 Gamma0ForItems_8_inBins_1_NoCut_12: p_ij[8,1] + gamma[8] <= 1
 12_Supports_8_InSameBin1_2: p_ij[8,2] - p_ij[12,2] - s[8,12] <= 0
 12_Supports_8_InSameBin2_2: - p_ij[8,2] + p_ij[12,2] - s[8,12] <= 0
 Gamma0ForItems_8_inBins_2_NoCut_12: p_ij[8,2] + gamma[8] <= 98
 12_Supports_8_InSameBin1_3: p_ij[8,3] - p_ij[12,3] - s[8,12] <= 0
 12_Supports_8_InSameBin2_3: - p_ij[8,3] + p_ij[12,3] - s[8,12] <= 0
 Gamma0ForItems_8_inBins_3_NoCut_12: p_ij[8,3] + gamma[8] <= 98
 8_OnCutIn2_12: 10000 p_ij[8,2] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 8_OnCutIn3_12: 10000 p_ij[8,3] + 1.261904761904762 xi_l[8] + zi[8]
   + 10000 gamma[8] <= 20053
 GroundConstraint_9: zi[9] + 155 g[9] <= 155
 AbsZ19_0: - zi[9] + zi_hi[0] - k[9,0] <= 0
 AbsZ29_0: zi[9] - zi_hi[0] - k[9,0] <= 0
 0_Below_9: zi[9] - zi_hi[0] + k[9,0] + 310 m[9,0] <= 310
 0_NotBelow_9: - zi[9] + zi_hi[0] + k[9,0] - 310 m[9,0] <= 0
 HeightCompat_9_0: h[9,0] - k[9,0] <= 0
 NoHeightCompat_9_0: - 155 h[9,0] + k[9,0] <= 0
 Overlap_9_0: - x_p[0,9] - x_p[9,0] + o[9,0] <= 0
 NoOverlap_9_0: x_p[0,9] + x_p[9,0] - 2 o[9,0] <= 0
 3_Supports_9_0: - h[9,0] - o[9,0] - s[9,0] <= -1
 3_DoesntSupport_9_0: h[9,0] + o[9,0] + 2 s[9,0] <= 2
 beta1_9_0: - s[9,0] + beta1[9,0] <= 0
 beta2_9_0: - s[9,0] + beta2[9,0] <= 0
 sBeta1Flag_9_0: - s[9,0] + beta1[9,0] <= 0
 sBeta2Flag_9_0: - s[9,0] + beta2[9,0] <= 0
 Eta1Overlap_9_0: eta1[9,0] + beta1[9,0] <= 1
 Eta2Overlap_9_0: eta2[9,0] + beta2[9,0] <= 1
 Eta1Flag_9_0: xi_l[0] - xi_l[9] - 300 eta1[9,0] <= 0
 Eta2Flag_9_0: - xi_r[0] + xi_r[9] - 300 eta2[9,0] <= 0
 0_Supports_9_InSameBin1_0: - p_ij[0,0] + p_ij[9,0] - s[9,0] <= 0
 0_Supports_9_InSameBin2_0: p_ij[0,0] - p_ij[9,0] - s[9,0] <= 0
 Gamma0ForItems_9_inBins_0_NoCut_0: p_ij[9,0] + gamma[9] <= 1
 0_Supports_9_InSameBin1_1: - p_ij[0,1] + p_ij[9,1] - s[9,0] <= 0
 0_Supports_9_InSameBin2_1: p_ij[0,1] - p_ij[9,1] - s[9,0] <= 0
 Gamma0ForItems_9_inBins_1_NoCut_0: p_ij[9,1] + gamma[9] <= 1
 0_Supports_9_InSameBin1_2: - p_ij[0,2] + p_ij[9,2] - s[9,0] <= 0
 0_Supports_9_InSameBin2_2: p_ij[0,2] - p_ij[9,2] - s[9,0] <= 0
 Gamma0ForItems_9_inBins_2_NoCut_0: p_ij[9,2] + gamma[9] <= 98
 0_Supports_9_InSameBin1_3: - p_ij[0,3] + p_ij[9,3] - s[9,0] <= 0
 0_Supports_9_InSameBin2_3: p_ij[0,3] - p_ij[9,3] - s[9,0] <= 0
 Gamma0ForItems_9_inBins_3_NoCut_0: p_ij[9,3] + gamma[9] <= 98
 9_OnCutIn2_0: 10000 p_ij[9,2] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 9_OnCutIn3_0: 10000 p_ij[9,3] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 AbsZ19_1: - zi[9] + zi_hi[1] - k[9,1] <= 0
 AbsZ29_1: zi[9] - zi_hi[1] - k[9,1] <= 0
 1_Below_9: zi[9] - zi_hi[1] + k[9,1] + 310 m[9,1] <= 310
 1_NotBelow_9: - zi[9] + zi_hi[1] + k[9,1] - 310 m[9,1] <= 0
 HeightCompat_9_1: h[9,1] - k[9,1] <= 0
 NoHeightCompat_9_1: - 155 h[9,1] + k[9,1] <= 0
 Overlap_9_1: - x_p[1,9] - x_p[9,1] + o[9,1] <= 0
 NoOverlap_9_1: x_p[1,9] + x_p[9,1] - 2 o[9,1] <= 0
 3_Supports_9_1: - h[9,1] - o[9,1] - s[9,1] <= -1
 3_DoesntSupport_9_1: h[9,1] + o[9,1] + 2 s[9,1] <= 2
 beta1_9_1: - s[9,1] + beta1[9,1] <= 0
 beta2_9_1: - s[9,1] + beta2[9,1] <= 0
 sBeta1Flag_9_1: - s[9,1] + beta1[9,1] <= 0
 sBeta2Flag_9_1: - s[9,1] + beta2[9,1] <= 0
 Eta1Overlap_9_1: eta1[9,1] + beta1[9,1] <= 1
 Eta2Overlap_9_1: eta2[9,1] + beta2[9,1] <= 1
 Eta1Flag_9_1: xi_l[1] - xi_l[9] - 300 eta1[9,1] <= 0
 Eta2Flag_9_1: - xi_r[1] + xi_r[9] - 300 eta2[9,1] <= 0
 1_Supports_9_InSameBin1_0: - p_ij[1,0] + p_ij[9,0] - s[9,1] <= 0
 1_Supports_9_InSameBin2_0: p_ij[1,0] - p_ij[9,0] - s[9,1] <= 0
 Gamma0ForItems_9_inBins_0_NoCut_1: p_ij[9,0] + gamma[9] <= 1
 1_Supports_9_InSameBin1_1: - p_ij[1,1] + p_ij[9,1] - s[9,1] <= 0
 1_Supports_9_InSameBin2_1: p_ij[1,1] - p_ij[9,1] - s[9,1] <= 0
 Gamma0ForItems_9_inBins_1_NoCut_1: p_ij[9,1] + gamma[9] <= 1
 1_Supports_9_InSameBin1_2: - p_ij[1,2] + p_ij[9,2] - s[9,1] <= 0
 1_Supports_9_InSameBin2_2: p_ij[1,2] - p_ij[9,2] - s[9,1] <= 0
 Gamma0ForItems_9_inBins_2_NoCut_1: p_ij[9,2] + gamma[9] <= 98
 1_Supports_9_InSameBin1_3: - p_ij[1,3] + p_ij[9,3] - s[9,1] <= 0
 1_Supports_9_InSameBin2_3: p_ij[1,3] - p_ij[9,3] - s[9,1] <= 0
 Gamma0ForItems_9_inBins_3_NoCut_1: p_ij[9,3] + gamma[9] <= 98
 9_OnCutIn2_1: 10000 p_ij[9,2] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 9_OnCutIn3_1: 10000 p_ij[9,3] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 AbsZ19_2: - zi[9] + zi_hi[2] - k[9,2] <= 0
 AbsZ29_2: zi[9] - zi_hi[2] - k[9,2] <= 0
 2_Below_9: zi[9] - zi_hi[2] + k[9,2] + 310 m[9,2] <= 310
 2_NotBelow_9: - zi[9] + zi_hi[2] + k[9,2] - 310 m[9,2] <= 0
 HeightCompat_9_2: h[9,2] - k[9,2] <= 0
 NoHeightCompat_9_2: - 155 h[9,2] + k[9,2] <= 0
 Overlap_9_2: - x_p[2,9] - x_p[9,2] + o[9,2] <= 0
 NoOverlap_9_2: x_p[2,9] + x_p[9,2] - 2 o[9,2] <= 0
 3_Supports_9_2: - h[9,2] - o[9,2] - s[9,2] <= -1
 3_DoesntSupport_9_2: h[9,2] + o[9,2] + 2 s[9,2] <= 2
 beta1_9_2: - s[9,2] + beta1[9,2] <= 0
 beta2_9_2: - s[9,2] + beta2[9,2] <= 0
 sBeta1Flag_9_2: - s[9,2] + beta1[9,2] <= 0
 sBeta2Flag_9_2: - s[9,2] + beta2[9,2] <= 0
 Eta1Overlap_9_2: eta1[9,2] + beta1[9,2] <= 1
 Eta2Overlap_9_2: eta2[9,2] + beta2[9,2] <= 1
 Eta1Flag_9_2: xi_l[2] - xi_l[9] - 300 eta1[9,2] <= 0
 Eta2Flag_9_2: - xi_r[2] + xi_r[9] - 300 eta2[9,2] <= 0
 2_Supports_9_InSameBin1_0: - p_ij[2,0] + p_ij[9,0] - s[9,2] <= 0
 2_Supports_9_InSameBin2_0: p_ij[2,0] - p_ij[9,0] - s[9,2] <= 0
 Gamma0ForItems_9_inBins_0_NoCut_2: p_ij[9,0] + gamma[9] <= 1
 2_Supports_9_InSameBin1_1: - p_ij[2,1] + p_ij[9,1] - s[9,2] <= 0
 2_Supports_9_InSameBin2_1: p_ij[2,1] - p_ij[9,1] - s[9,2] <= 0
 Gamma0ForItems_9_inBins_1_NoCut_2: p_ij[9,1] + gamma[9] <= 1
 2_Supports_9_InSameBin1_2: - p_ij[2,2] + p_ij[9,2] - s[9,2] <= 0
 2_Supports_9_InSameBin2_2: p_ij[2,2] - p_ij[9,2] - s[9,2] <= 0
 Gamma0ForItems_9_inBins_2_NoCut_2: p_ij[9,2] + gamma[9] <= 98
 2_Supports_9_InSameBin1_3: - p_ij[2,3] + p_ij[9,3] - s[9,2] <= 0
 2_Supports_9_InSameBin2_3: p_ij[2,3] - p_ij[9,3] - s[9,2] <= 0
 Gamma0ForItems_9_inBins_3_NoCut_2: p_ij[9,3] + gamma[9] <= 98
 9_OnCutIn2_2: 10000 p_ij[9,2] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 9_OnCutIn3_2: 10000 p_ij[9,3] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 AbsZ19_3: - zi[9] + zi_hi[3] - k[9,3] <= 0
 AbsZ29_3: zi[9] - zi_hi[3] - k[9,3] <= 0
 3_Below_9: zi[9] - zi_hi[3] + k[9,3] + 310 m[9,3] <= 310
 3_NotBelow_9: - zi[9] + zi_hi[3] + k[9,3] - 310 m[9,3] <= 0
 HeightCompat_9_3: h[9,3] - k[9,3] <= 0
 NoHeightCompat_9_3: - 155 h[9,3] + k[9,3] <= 0
 Overlap_9_3: - x_p[3,9] - x_p[9,3] + o[9,3] <= 0
 NoOverlap_9_3: x_p[3,9] + x_p[9,3] - 2 o[9,3] <= 0
 3_Supports_9_3: - h[9,3] - o[9,3] - s[9,3] <= -1
 3_DoesntSupport_9_3: h[9,3] + o[9,3] + 2 s[9,3] <= 2
 beta1_9_3: - s[9,3] + beta1[9,3] <= 0
 beta2_9_3: - s[9,3] + beta2[9,3] <= 0
 sBeta1Flag_9_3: - s[9,3] + beta1[9,3] <= 0
 sBeta2Flag_9_3: - s[9,3] + beta2[9,3] <= 0
 Eta1Overlap_9_3: eta1[9,3] + beta1[9,3] <= 1
 Eta2Overlap_9_3: eta2[9,3] + beta2[9,3] <= 1
 Eta1Flag_9_3: xi_l[3] - xi_l[9] - 300 eta1[9,3] <= 0
 Eta2Flag_9_3: - xi_r[3] + xi_r[9] - 300 eta2[9,3] <= 0
 3_Supports_9_InSameBin1_0: - p_ij[3,0] + p_ij[9,0] - s[9,3] <= 0
 3_Supports_9_InSameBin2_0: p_ij[3,0] - p_ij[9,0] - s[9,3] <= 0
 Gamma0ForItems_9_inBins_0_NoCut_3: p_ij[9,0] + gamma[9] <= 1
 3_Supports_9_InSameBin1_1: - p_ij[3,1] + p_ij[9,1] - s[9,3] <= 0
 3_Supports_9_InSameBin2_1: p_ij[3,1] - p_ij[9,1] - s[9,3] <= 0
 Gamma0ForItems_9_inBins_1_NoCut_3: p_ij[9,1] + gamma[9] <= 1
 3_Supports_9_InSameBin1_2: - p_ij[3,2] + p_ij[9,2] - s[9,3] <= 0
 3_Supports_9_InSameBin2_2: p_ij[3,2] - p_ij[9,2] - s[9,3] <= 0
 Gamma0ForItems_9_inBins_2_NoCut_3: p_ij[9,2] + gamma[9] <= 98
 3_Supports_9_InSameBin1_3: - p_ij[3,3] + p_ij[9,3] - s[9,3] <= 0
 3_Supports_9_InSameBin2_3: p_ij[3,3] - p_ij[9,3] - s[9,3] <= 0
 Gamma0ForItems_9_inBins_3_NoCut_3: p_ij[9,3] + gamma[9] <= 98
 9_OnCutIn2_3: 10000 p_ij[9,2] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 9_OnCutIn3_3: 10000 p_ij[9,3] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 AbsZ19_4: - zi[9] + zi_hi[4] - k[9,4] <= 0
 AbsZ29_4: zi[9] - zi_hi[4] - k[9,4] <= 0
 4_Below_9: zi[9] - zi_hi[4] + k[9,4] + 310 m[9,4] <= 310
 4_NotBelow_9: - zi[9] + zi_hi[4] + k[9,4] - 310 m[9,4] <= 0
 HeightCompat_9_4: h[9,4] - k[9,4] <= 0
 NoHeightCompat_9_4: - 155 h[9,4] + k[9,4] <= 0
 Overlap_9_4: - x_p[4,9] - x_p[9,4] + o[9,4] <= 0
 NoOverlap_9_4: x_p[4,9] + x_p[9,4] - 2 o[9,4] <= 0
 3_Supports_9_4: - h[9,4] - o[9,4] - s[9,4] <= -1
 3_DoesntSupport_9_4: h[9,4] + o[9,4] + 2 s[9,4] <= 2
 beta1_9_4: - s[9,4] + beta1[9,4] <= 0
 beta2_9_4: - s[9,4] + beta2[9,4] <= 0
 sBeta1Flag_9_4: - s[9,4] + beta1[9,4] <= 0
 sBeta2Flag_9_4: - s[9,4] + beta2[9,4] <= 0
 Eta1Overlap_9_4: eta1[9,4] + beta1[9,4] <= 1
 Eta2Overlap_9_4: eta2[9,4] + beta2[9,4] <= 1
 Eta1Flag_9_4: xi_l[4] - xi_l[9] - 300 eta1[9,4] <= 0
 Eta2Flag_9_4: - xi_r[4] + xi_r[9] - 300 eta2[9,4] <= 0
 4_Supports_9_InSameBin1_0: - p_ij[4,0] + p_ij[9,0] - s[9,4] <= 0
 4_Supports_9_InSameBin2_0: p_ij[4,0] - p_ij[9,0] - s[9,4] <= 0
 Gamma0ForItems_9_inBins_0_NoCut_4: p_ij[9,0] + gamma[9] <= 1
 4_Supports_9_InSameBin1_1: - p_ij[4,1] + p_ij[9,1] - s[9,4] <= 0
 4_Supports_9_InSameBin2_1: p_ij[4,1] - p_ij[9,1] - s[9,4] <= 0
 Gamma0ForItems_9_inBins_1_NoCut_4: p_ij[9,1] + gamma[9] <= 1
 4_Supports_9_InSameBin1_2: - p_ij[4,2] + p_ij[9,2] - s[9,4] <= 0
 4_Supports_9_InSameBin2_2: p_ij[4,2] - p_ij[9,2] - s[9,4] <= 0
 Gamma0ForItems_9_inBins_2_NoCut_4: p_ij[9,2] + gamma[9] <= 98
 4_Supports_9_InSameBin1_3: - p_ij[4,3] + p_ij[9,3] - s[9,4] <= 0
 4_Supports_9_InSameBin2_3: p_ij[4,3] - p_ij[9,3] - s[9,4] <= 0
 Gamma0ForItems_9_inBins_3_NoCut_4: p_ij[9,3] + gamma[9] <= 98
 9_OnCutIn2_4: 10000 p_ij[9,2] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 9_OnCutIn3_4: 10000 p_ij[9,3] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 AbsZ19_5: - zi[9] + zi_hi[5] - k[9,5] <= 0
 AbsZ29_5: zi[9] - zi_hi[5] - k[9,5] <= 0
 5_Below_9: zi[9] - zi_hi[5] + k[9,5] + 310 m[9,5] <= 310
 5_NotBelow_9: - zi[9] + zi_hi[5] + k[9,5] - 310 m[9,5] <= 0
 HeightCompat_9_5: h[9,5] - k[9,5] <= 0
 NoHeightCompat_9_5: - 155 h[9,5] + k[9,5] <= 0
 Overlap_9_5: - x_p[5,9] - x_p[9,5] + o[9,5] <= 0
 NoOverlap_9_5: x_p[5,9] + x_p[9,5] - 2 o[9,5] <= 0
 3_Supports_9_5: - h[9,5] - o[9,5] - s[9,5] <= -1
 3_DoesntSupport_9_5: h[9,5] + o[9,5] + 2 s[9,5] <= 2
 beta1_9_5: - s[9,5] + beta1[9,5] <= 0
 beta2_9_5: - s[9,5] + beta2[9,5] <= 0
 sBeta1Flag_9_5: - s[9,5] + beta1[9,5] <= 0
 sBeta2Flag_9_5: - s[9,5] + beta2[9,5] <= 0
 Eta1Overlap_9_5: eta1[9,5] + beta1[9,5] <= 1
 Eta2Overlap_9_5: eta2[9,5] + beta2[9,5] <= 1
 Eta1Flag_9_5: xi_l[5] - xi_l[9] - 300 eta1[9,5] <= 0
 Eta2Flag_9_5: - xi_r[5] + xi_r[9] - 300 eta2[9,5] <= 0
 5_Supports_9_InSameBin1_0: - p_ij[5,0] + p_ij[9,0] - s[9,5] <= 0
 5_Supports_9_InSameBin2_0: p_ij[5,0] - p_ij[9,0] - s[9,5] <= 0
 Gamma0ForItems_9_inBins_0_NoCut_5: p_ij[9,0] + gamma[9] <= 1
 5_Supports_9_InSameBin1_1: - p_ij[5,1] + p_ij[9,1] - s[9,5] <= 0
 5_Supports_9_InSameBin2_1: p_ij[5,1] - p_ij[9,1] - s[9,5] <= 0
 Gamma0ForItems_9_inBins_1_NoCut_5: p_ij[9,1] + gamma[9] <= 1
 5_Supports_9_InSameBin1_2: - p_ij[5,2] + p_ij[9,2] - s[9,5] <= 0
 5_Supports_9_InSameBin2_2: p_ij[5,2] - p_ij[9,2] - s[9,5] <= 0
 Gamma0ForItems_9_inBins_2_NoCut_5: p_ij[9,2] + gamma[9] <= 98
 5_Supports_9_InSameBin1_3: - p_ij[5,3] + p_ij[9,3] - s[9,5] <= 0
 5_Supports_9_InSameBin2_3: p_ij[5,3] - p_ij[9,3] - s[9,5] <= 0
 Gamma0ForItems_9_inBins_3_NoCut_5: p_ij[9,3] + gamma[9] <= 98
 9_OnCutIn2_5: 10000 p_ij[9,2] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 9_OnCutIn3_5: 10000 p_ij[9,3] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 AbsZ19_6: - zi[9] + zi_hi[6] - k[9,6] <= 0
 AbsZ29_6: zi[9] - zi_hi[6] - k[9,6] <= 0
 6_Below_9: zi[9] - zi_hi[6] + k[9,6] + 310 m[9,6] <= 310
 6_NotBelow_9: - zi[9] + zi_hi[6] + k[9,6] - 310 m[9,6] <= 0
 HeightCompat_9_6: h[9,6] - k[9,6] <= 0
 NoHeightCompat_9_6: - 155 h[9,6] + k[9,6] <= 0
 Overlap_9_6: - x_p[6,9] - x_p[9,6] + o[9,6] <= 0
 NoOverlap_9_6: x_p[6,9] + x_p[9,6] - 2 o[9,6] <= 0
 3_Supports_9_6: - h[9,6] - o[9,6] - s[9,6] <= -1
 3_DoesntSupport_9_6: h[9,6] + o[9,6] + 2 s[9,6] <= 2
 beta1_9_6: - s[9,6] + beta1[9,6] <= 0
 beta2_9_6: - s[9,6] + beta2[9,6] <= 0
 sBeta1Flag_9_6: - s[9,6] + beta1[9,6] <= 0
 sBeta2Flag_9_6: - s[9,6] + beta2[9,6] <= 0
 Eta1Overlap_9_6: eta1[9,6] + beta1[9,6] <= 1
 Eta2Overlap_9_6: eta2[9,6] + beta2[9,6] <= 1
 Eta1Flag_9_6: xi_l[6] - xi_l[9] - 300 eta1[9,6] <= 0
 Eta2Flag_9_6: - xi_r[6] + xi_r[9] - 300 eta2[9,6] <= 0
 6_Supports_9_InSameBin1_0: - p_ij[6,0] + p_ij[9,0] - s[9,6] <= 0
 6_Supports_9_InSameBin2_0: p_ij[6,0] - p_ij[9,0] - s[9,6] <= 0
 Gamma0ForItems_9_inBins_0_NoCut_6: p_ij[9,0] + gamma[9] <= 1
 6_Supports_9_InSameBin1_1: - p_ij[6,1] + p_ij[9,1] - s[9,6] <= 0
 6_Supports_9_InSameBin2_1: p_ij[6,1] - p_ij[9,1] - s[9,6] <= 0
 Gamma0ForItems_9_inBins_1_NoCut_6: p_ij[9,1] + gamma[9] <= 1
 6_Supports_9_InSameBin1_2: - p_ij[6,2] + p_ij[9,2] - s[9,6] <= 0
 6_Supports_9_InSameBin2_2: p_ij[6,2] - p_ij[9,2] - s[9,6] <= 0
 Gamma0ForItems_9_inBins_2_NoCut_6: p_ij[9,2] + gamma[9] <= 98
 6_Supports_9_InSameBin1_3: - p_ij[6,3] + p_ij[9,3] - s[9,6] <= 0
 6_Supports_9_InSameBin2_3: p_ij[6,3] - p_ij[9,3] - s[9,6] <= 0
 Gamma0ForItems_9_inBins_3_NoCut_6: p_ij[9,3] + gamma[9] <= 98
 9_OnCutIn2_6: 10000 p_ij[9,2] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 9_OnCutIn3_6: 10000 p_ij[9,3] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 AbsZ19_7: - zi[9] + zi_hi[7] - k[9,7] <= 0
 AbsZ29_7: zi[9] - zi_hi[7] - k[9,7] <= 0
 7_Below_9: zi[9] - zi_hi[7] + k[9,7] + 310 m[9,7] <= 310
 7_NotBelow_9: - zi[9] + zi_hi[7] + k[9,7] - 310 m[9,7] <= 0
 HeightCompat_9_7: h[9,7] - k[9,7] <= 0
 NoHeightCompat_9_7: - 155 h[9,7] + k[9,7] <= 0
 Overlap_9_7: - x_p[7,9] - x_p[9,7] + o[9,7] <= 0
 NoOverlap_9_7: x_p[7,9] + x_p[9,7] - 2 o[9,7] <= 0
 3_Supports_9_7: - h[9,7] - o[9,7] - s[9,7] <= -1
 3_DoesntSupport_9_7: h[9,7] + o[9,7] + 2 s[9,7] <= 2
 beta1_9_7: - s[9,7] + beta1[9,7] <= 0
 beta2_9_7: - s[9,7] + beta2[9,7] <= 0
 sBeta1Flag_9_7: - s[9,7] + beta1[9,7] <= 0
 sBeta2Flag_9_7: - s[9,7] + beta2[9,7] <= 0
 Eta1Overlap_9_7: eta1[9,7] + beta1[9,7] <= 1
 Eta2Overlap_9_7: eta2[9,7] + beta2[9,7] <= 1
 Eta1Flag_9_7: xi_l[7] - xi_l[9] - 300 eta1[9,7] <= 0
 Eta2Flag_9_7: - xi_r[7] + xi_r[9] - 300 eta2[9,7] <= 0
 7_Supports_9_InSameBin1_0: - p_ij[7,0] + p_ij[9,0] - s[9,7] <= 0
 7_Supports_9_InSameBin2_0: p_ij[7,0] - p_ij[9,0] - s[9,7] <= 0
 Gamma0ForItems_9_inBins_0_NoCut_7: p_ij[9,0] + gamma[9] <= 1
 7_Supports_9_InSameBin1_1: - p_ij[7,1] + p_ij[9,1] - s[9,7] <= 0
 7_Supports_9_InSameBin2_1: p_ij[7,1] - p_ij[9,1] - s[9,7] <= 0
 Gamma0ForItems_9_inBins_1_NoCut_7: p_ij[9,1] + gamma[9] <= 1
 7_Supports_9_InSameBin1_2: - p_ij[7,2] + p_ij[9,2] - s[9,7] <= 0
 7_Supports_9_InSameBin2_2: p_ij[7,2] - p_ij[9,2] - s[9,7] <= 0
 Gamma0ForItems_9_inBins_2_NoCut_7: p_ij[9,2] + gamma[9] <= 98
 7_Supports_9_InSameBin1_3: - p_ij[7,3] + p_ij[9,3] - s[9,7] <= 0
 7_Supports_9_InSameBin2_3: p_ij[7,3] - p_ij[9,3] - s[9,7] <= 0
 Gamma0ForItems_9_inBins_3_NoCut_7: p_ij[9,3] + gamma[9] <= 98
 9_OnCutIn2_7: 10000 p_ij[9,2] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 9_OnCutIn3_7: 10000 p_ij[9,3] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 AbsZ19_8: - zi[9] + zi_hi[8] - k[9,8] <= 0
 AbsZ29_8: zi[9] - zi_hi[8] - k[9,8] <= 0
 8_Below_9: zi[9] - zi_hi[8] + k[9,8] + 310 m[9,8] <= 310
 8_NotBelow_9: - zi[9] + zi_hi[8] + k[9,8] - 310 m[9,8] <= 0
 HeightCompat_9_8: h[9,8] - k[9,8] <= 0
 NoHeightCompat_9_8: - 155 h[9,8] + k[9,8] <= 0
 Overlap_9_8: - x_p[8,9] - x_p[9,8] + o[9,8] <= 0
 NoOverlap_9_8: x_p[8,9] + x_p[9,8] - 2 o[9,8] <= 0
 3_Supports_9_8: - h[9,8] - o[9,8] - s[9,8] <= -1
 3_DoesntSupport_9_8: h[9,8] + o[9,8] + 2 s[9,8] <= 2
 beta1_9_8: - s[9,8] + beta1[9,8] <= 0
 beta2_9_8: - s[9,8] + beta2[9,8] <= 0
 sBeta1Flag_9_8: - s[9,8] + beta1[9,8] <= 0
 sBeta2Flag_9_8: - s[9,8] + beta2[9,8] <= 0
 Eta1Overlap_9_8: eta1[9,8] + beta1[9,8] <= 1
 Eta2Overlap_9_8: eta2[9,8] + beta2[9,8] <= 1
 Eta1Flag_9_8: xi_l[8] - xi_l[9] - 300 eta1[9,8] <= 0
 Eta2Flag_9_8: - xi_r[8] + xi_r[9] - 300 eta2[9,8] <= 0
 8_Supports_9_InSameBin1_0: - p_ij[8,0] + p_ij[9,0] - s[9,8] <= 0
 8_Supports_9_InSameBin2_0: p_ij[8,0] - p_ij[9,0] - s[9,8] <= 0
 Gamma0ForItems_9_inBins_0_NoCut_8: p_ij[9,0] + gamma[9] <= 1
 8_Supports_9_InSameBin1_1: - p_ij[8,1] + p_ij[9,1] - s[9,8] <= 0
 8_Supports_9_InSameBin2_1: p_ij[8,1] - p_ij[9,1] - s[9,8] <= 0
 Gamma0ForItems_9_inBins_1_NoCut_8: p_ij[9,1] + gamma[9] <= 1
 8_Supports_9_InSameBin1_2: - p_ij[8,2] + p_ij[9,2] - s[9,8] <= 0
 8_Supports_9_InSameBin2_2: p_ij[8,2] - p_ij[9,2] - s[9,8] <= 0
 Gamma0ForItems_9_inBins_2_NoCut_8: p_ij[9,2] + gamma[9] <= 98
 8_Supports_9_InSameBin1_3: - p_ij[8,3] + p_ij[9,3] - s[9,8] <= 0
 8_Supports_9_InSameBin2_3: p_ij[8,3] - p_ij[9,3] - s[9,8] <= 0
 Gamma0ForItems_9_inBins_3_NoCut_8: p_ij[9,3] + gamma[9] <= 98
 9_OnCutIn2_8: 10000 p_ij[9,2] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 9_OnCutIn3_8: 10000 p_ij[9,3] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 9_OnCutIn2_9: 10000 p_ij[9,2] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 9_OnCutIn3_9: 10000 p_ij[9,3] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 AbsZ19_10: - zi[9] + zi_hi[10] - k[9,10] <= 0
 AbsZ29_10: zi[9] - zi_hi[10] - k[9,10] <= 0
 10_Below_9: zi[9] - zi_hi[10] + k[9,10] + 310 m[9,10] <= 310
 10_NotBelow_9: - zi[9] + zi_hi[10] + k[9,10] - 310 m[9,10] <= 0
 HeightCompat_9_10: h[9,10] - k[9,10] <= 0
 NoHeightCompat_9_10: - 155 h[9,10] + k[9,10] <= 0
 Overlap_9_10: - x_p[9,10] - x_p[10,9] + o[9,10] <= 0
 NoOverlap_9_10: x_p[9,10] + x_p[10,9] - 2 o[9,10] <= 0
 3_Supports_9_10: - h[9,10] - o[9,10] - s[9,10] <= -1
 3_DoesntSupport_9_10: h[9,10] + o[9,10] + 2 s[9,10] <= 2
 beta1_9_10: - s[9,10] + beta1[9,10] <= 0
 beta2_9_10: - s[9,10] + beta2[9,10] <= 0
 sBeta1Flag_9_10: - s[9,10] + beta1[9,10] <= 0
 sBeta2Flag_9_10: - s[9,10] + beta2[9,10] <= 0
 Eta1Overlap_9_10: eta1[9,10] + beta1[9,10] <= 1
 Eta2Overlap_9_10: eta2[9,10] + beta2[9,10] <= 1
 Eta1Flag_9_10: - xi_l[9] + xi_l[10] - 300 eta1[9,10] <= 0
 Eta2Flag_9_10: xi_r[9] - xi_r[10] - 300 eta2[9,10] <= 0
 10_Supports_9_InSameBin1_0: p_ij[9,0] - p_ij[10,0] - s[9,10] <= 0
 10_Supports_9_InSameBin2_0: - p_ij[9,0] + p_ij[10,0] - s[9,10] <= 0
 Gamma0ForItems_9_inBins_0_NoCut_10: p_ij[9,0] + gamma[9] <= 1
 10_Supports_9_InSameBin1_1: p_ij[9,1] - p_ij[10,1] - s[9,10] <= 0
 10_Supports_9_InSameBin2_1: - p_ij[9,1] + p_ij[10,1] - s[9,10] <= 0
 Gamma0ForItems_9_inBins_1_NoCut_10: p_ij[9,1] + gamma[9] <= 1
 10_Supports_9_InSameBin1_2: p_ij[9,2] - p_ij[10,2] - s[9,10] <= 0
 10_Supports_9_InSameBin2_2: - p_ij[9,2] + p_ij[10,2] - s[9,10] <= 0
 Gamma0ForItems_9_inBins_2_NoCut_10: p_ij[9,2] + gamma[9] <= 98
 10_Supports_9_InSameBin1_3: p_ij[9,3] - p_ij[10,3] - s[9,10] <= 0
 10_Supports_9_InSameBin2_3: - p_ij[9,3] + p_ij[10,3] - s[9,10] <= 0
 Gamma0ForItems_9_inBins_3_NoCut_10: p_ij[9,3] + gamma[9] <= 98
 9_OnCutIn2_10: 10000 p_ij[9,2] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 9_OnCutIn3_10: 10000 p_ij[9,3] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 AbsZ19_11: - zi[9] + zi_hi[11] - k[9,11] <= 0
 AbsZ29_11: zi[9] - zi_hi[11] - k[9,11] <= 0
 11_Below_9: zi[9] - zi_hi[11] + k[9,11] + 310 m[9,11] <= 310
 11_NotBelow_9: - zi[9] + zi_hi[11] + k[9,11] - 310 m[9,11] <= 0
 HeightCompat_9_11: h[9,11] - k[9,11] <= 0
 NoHeightCompat_9_11: - 155 h[9,11] + k[9,11] <= 0
 Overlap_9_11: - x_p[9,11] - x_p[11,9] + o[9,11] <= 0
 NoOverlap_9_11: x_p[9,11] + x_p[11,9] - 2 o[9,11] <= 0
 3_Supports_9_11: - h[9,11] - o[9,11] - s[9,11] <= -1
 3_DoesntSupport_9_11: h[9,11] + o[9,11] + 2 s[9,11] <= 2
 beta1_9_11: - s[9,11] + beta1[9,11] <= 0
 beta2_9_11: - s[9,11] + beta2[9,11] <= 0
 sBeta1Flag_9_11: - s[9,11] + beta1[9,11] <= 0
 sBeta2Flag_9_11: - s[9,11] + beta2[9,11] <= 0
 Eta1Overlap_9_11: eta1[9,11] + beta1[9,11] <= 1
 Eta2Overlap_9_11: eta2[9,11] + beta2[9,11] <= 1
 Eta1Flag_9_11: - xi_l[9] + xi_l[11] - 300 eta1[9,11] <= 0
 Eta2Flag_9_11: xi_r[9] - xi_r[11] - 300 eta2[9,11] <= 0
 11_Supports_9_InSameBin1_0: p_ij[9,0] - p_ij[11,0] - s[9,11] <= 0
 11_Supports_9_InSameBin2_0: - p_ij[9,0] + p_ij[11,0] - s[9,11] <= 0
 Gamma0ForItems_9_inBins_0_NoCut_11: p_ij[9,0] + gamma[9] <= 1
 11_Supports_9_InSameBin1_1: p_ij[9,1] - p_ij[11,1] - s[9,11] <= 0
 11_Supports_9_InSameBin2_1: - p_ij[9,1] + p_ij[11,1] - s[9,11] <= 0
 Gamma0ForItems_9_inBins_1_NoCut_11: p_ij[9,1] + gamma[9] <= 1
 11_Supports_9_InSameBin1_2: p_ij[9,2] - p_ij[11,2] - s[9,11] <= 0
 11_Supports_9_InSameBin2_2: - p_ij[9,2] + p_ij[11,2] - s[9,11] <= 0
 Gamma0ForItems_9_inBins_2_NoCut_11: p_ij[9,2] + gamma[9] <= 98
 11_Supports_9_InSameBin1_3: p_ij[9,3] - p_ij[11,3] - s[9,11] <= 0
 11_Supports_9_InSameBin2_3: - p_ij[9,3] + p_ij[11,3] - s[9,11] <= 0
 Gamma0ForItems_9_inBins_3_NoCut_11: p_ij[9,3] + gamma[9] <= 98
 9_OnCutIn2_11: 10000 p_ij[9,2] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 9_OnCutIn3_11: 10000 p_ij[9,3] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 AbsZ19_12: - zi[9] + zi_hi[12] - k[9,12] <= 0
 AbsZ29_12: zi[9] - zi_hi[12] - k[9,12] <= 0
 12_Below_9: zi[9] - zi_hi[12] + k[9,12] + 310 m[9,12] <= 310
 12_NotBelow_9: - zi[9] + zi_hi[12] + k[9,12] - 310 m[9,12] <= 0
 HeightCompat_9_12: h[9,12] - k[9,12] <= 0
 NoHeightCompat_9_12: - 155 h[9,12] + k[9,12] <= 0
 Overlap_9_12: - x_p[9,12] - x_p[12,9] + o[9,12] <= 0
 NoOverlap_9_12: x_p[9,12] + x_p[12,9] - 2 o[9,12] <= 0
 3_Supports_9_12: - h[9,12] - o[9,12] - s[9,12] <= -1
 3_DoesntSupport_9_12: h[9,12] + o[9,12] + 2 s[9,12] <= 2
 beta1_9_12: - s[9,12] + beta1[9,12] <= 0
 beta2_9_12: - s[9,12] + beta2[9,12] <= 0
 sBeta1Flag_9_12: - s[9,12] + beta1[9,12] <= 0
 sBeta2Flag_9_12: - s[9,12] + beta2[9,12] <= 0
 Eta1Overlap_9_12: eta1[9,12] + beta1[9,12] <= 1
 Eta2Overlap_9_12: eta2[9,12] + beta2[9,12] <= 1
 Eta1Flag_9_12: - xi_l[9] + xi_l[12] - 300 eta1[9,12] <= 0
 Eta2Flag_9_12: xi_r[9] - xi_r[12] - 300 eta2[9,12] <= 0
 12_Supports_9_InSameBin1_0: p_ij[9,0] - p_ij[12,0] - s[9,12] <= 0
 12_Supports_9_InSameBin2_0: - p_ij[9,0] + p_ij[12,0] - s[9,12] <= 0
 Gamma0ForItems_9_inBins_0_NoCut_12: p_ij[9,0] + gamma[9] <= 1
 12_Supports_9_InSameBin1_1: p_ij[9,1] - p_ij[12,1] - s[9,12] <= 0
 12_Supports_9_InSameBin2_1: - p_ij[9,1] + p_ij[12,1] - s[9,12] <= 0
 Gamma0ForItems_9_inBins_1_NoCut_12: p_ij[9,1] + gamma[9] <= 1
 12_Supports_9_InSameBin1_2: p_ij[9,2] - p_ij[12,2] - s[9,12] <= 0
 12_Supports_9_InSameBin2_2: - p_ij[9,2] + p_ij[12,2] - s[9,12] <= 0
 Gamma0ForItems_9_inBins_2_NoCut_12: p_ij[9,2] + gamma[9] <= 98
 12_Supports_9_InSameBin1_3: p_ij[9,3] - p_ij[12,3] - s[9,12] <= 0
 12_Supports_9_InSameBin2_3: - p_ij[9,3] + p_ij[12,3] - s[9,12] <= 0
 Gamma0ForItems_9_inBins_3_NoCut_12: p_ij[9,3] + gamma[9] <= 98
 9_OnCutIn2_12: 10000 p_ij[9,2] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 9_OnCutIn3_12: 10000 p_ij[9,3] + 1.261904761904762 xi_l[9] + zi[9]
   + 10000 gamma[9] <= 20053
 GroundConstraint_10: zi[10] + 155 g[10] <= 155
 AbsZ110_0: - zi[10] + zi_hi[0] - k[10,0] <= 0
 AbsZ210_0: zi[10] - zi_hi[0] - k[10,0] <= 0
 0_Below_10: zi[10] - zi_hi[0] + k[10,0] + 310 m[10,0] <= 310
 0_NotBelow_10: - zi[10] + zi_hi[0] + k[10,0] - 310 m[10,0] <= 0
 HeightCompat_10_0: h[10,0] - k[10,0] <= 0
 NoHeightCompat_10_0: - 155 h[10,0] + k[10,0] <= 0
 Overlap_10_0: - x_p[0,10] - x_p[10,0] + o[10,0] <= 0
 NoOverlap_10_0: x_p[0,10] + x_p[10,0] - 2 o[10,0] <= 0
 3_Supports_10_0: - h[10,0] - o[10,0] - s[10,0] <= -1
 3_DoesntSupport_10_0: h[10,0] + o[10,0] + 2 s[10,0] <= 2
 beta1_10_0: - s[10,0] + beta1[10,0] <= 0
 beta2_10_0: - s[10,0] + beta2[10,0] <= 0
 sBeta1Flag_10_0: - s[10,0] + beta1[10,0] <= 0
 sBeta2Flag_10_0: - s[10,0] + beta2[10,0] <= 0
 Eta1Overlap_10_0: eta1[10,0] + beta1[10,0] <= 1
 Eta2Overlap_10_0: eta2[10,0] + beta2[10,0] <= 1
 Eta1Flag_10_0: xi_l[0] - xi_l[10] - 300 eta1[10,0] <= 0
 Eta2Flag_10_0: - xi_r[0] + xi_r[10] - 300 eta2[10,0] <= 0
 0_Supports_10_InSameBin1_0: - p_ij[0,0] + p_ij[10,0] - s[10,0] <= 0
 0_Supports_10_InSameBin2_0: p_ij[0,0] - p_ij[10,0] - s[10,0] <= 0
 Gamma0ForItems_10_inBins_0_NoCut_0: p_ij[10,0] + gamma[10] <= 1
 0_Supports_10_InSameBin1_1: - p_ij[0,1] + p_ij[10,1] - s[10,0] <= 0
 0_Supports_10_InSameBin2_1: p_ij[0,1] - p_ij[10,1] - s[10,0] <= 0
 Gamma0ForItems_10_inBins_1_NoCut_0: p_ij[10,1] + gamma[10] <= 1
 0_Supports_10_InSameBin1_2: - p_ij[0,2] + p_ij[10,2] - s[10,0] <= 0
 0_Supports_10_InSameBin2_2: p_ij[0,2] - p_ij[10,2] - s[10,0] <= 0
 Gamma0ForItems_10_inBins_2_NoCut_0: p_ij[10,2] + gamma[10] <= 98
 0_Supports_10_InSameBin1_3: - p_ij[0,3] + p_ij[10,3] - s[10,0] <= 0
 0_Supports_10_InSameBin2_3: p_ij[0,3] - p_ij[10,3] - s[10,0] <= 0
 Gamma0ForItems_10_inBins_3_NoCut_0: p_ij[10,3] + gamma[10] <= 98
 10_OnCutIn2_0: 10000 p_ij[10,2] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 10_OnCutIn3_0: 10000 p_ij[10,3] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 AbsZ110_1: - zi[10] + zi_hi[1] - k[10,1] <= 0
 AbsZ210_1: zi[10] - zi_hi[1] - k[10,1] <= 0
 1_Below_10: zi[10] - zi_hi[1] + k[10,1] + 310 m[10,1] <= 310
 1_NotBelow_10: - zi[10] + zi_hi[1] + k[10,1] - 310 m[10,1] <= 0
 HeightCompat_10_1: h[10,1] - k[10,1] <= 0
 NoHeightCompat_10_1: - 155 h[10,1] + k[10,1] <= 0
 Overlap_10_1: - x_p[1,10] - x_p[10,1] + o[10,1] <= 0
 NoOverlap_10_1: x_p[1,10] + x_p[10,1] - 2 o[10,1] <= 0
 3_Supports_10_1: - h[10,1] - o[10,1] - s[10,1] <= -1
 3_DoesntSupport_10_1: h[10,1] + o[10,1] + 2 s[10,1] <= 2
 beta1_10_1: - s[10,1] + beta1[10,1] <= 0
 beta2_10_1: - s[10,1] + beta2[10,1] <= 0
 sBeta1Flag_10_1: - s[10,1] + beta1[10,1] <= 0
 sBeta2Flag_10_1: - s[10,1] + beta2[10,1] <= 0
 Eta1Overlap_10_1: eta1[10,1] + beta1[10,1] <= 1
 Eta2Overlap_10_1: eta2[10,1] + beta2[10,1] <= 1
 Eta1Flag_10_1: xi_l[1] - xi_l[10] - 300 eta1[10,1] <= 0
 Eta2Flag_10_1: - xi_r[1] + xi_r[10] - 300 eta2[10,1] <= 0
 1_Supports_10_InSameBin1_0: - p_ij[1,0] + p_ij[10,0] - s[10,1] <= 0
 1_Supports_10_InSameBin2_0: p_ij[1,0] - p_ij[10,0] - s[10,1] <= 0
 Gamma0ForItems_10_inBins_0_NoCut_1: p_ij[10,0] + gamma[10] <= 1
 1_Supports_10_InSameBin1_1: - p_ij[1,1] + p_ij[10,1] - s[10,1] <= 0
 1_Supports_10_InSameBin2_1: p_ij[1,1] - p_ij[10,1] - s[10,1] <= 0
 Gamma0ForItems_10_inBins_1_NoCut_1: p_ij[10,1] + gamma[10] <= 1
 1_Supports_10_InSameBin1_2: - p_ij[1,2] + p_ij[10,2] - s[10,1] <= 0
 1_Supports_10_InSameBin2_2: p_ij[1,2] - p_ij[10,2] - s[10,1] <= 0
 Gamma0ForItems_10_inBins_2_NoCut_1: p_ij[10,2] + gamma[10] <= 98
 1_Supports_10_InSameBin1_3: - p_ij[1,3] + p_ij[10,3] - s[10,1] <= 0
 1_Supports_10_InSameBin2_3: p_ij[1,3] - p_ij[10,3] - s[10,1] <= 0
 Gamma0ForItems_10_inBins_3_NoCut_1: p_ij[10,3] + gamma[10] <= 98
 10_OnCutIn2_1: 10000 p_ij[10,2] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 10_OnCutIn3_1: 10000 p_ij[10,3] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 AbsZ110_2: - zi[10] + zi_hi[2] - k[10,2] <= 0
 AbsZ210_2: zi[10] - zi_hi[2] - k[10,2] <= 0
 2_Below_10: zi[10] - zi_hi[2] + k[10,2] + 310 m[10,2] <= 310
 2_NotBelow_10: - zi[10] + zi_hi[2] + k[10,2] - 310 m[10,2] <= 0
 HeightCompat_10_2: h[10,2] - k[10,2] <= 0
 NoHeightCompat_10_2: - 155 h[10,2] + k[10,2] <= 0
 Overlap_10_2: - x_p[2,10] - x_p[10,2] + o[10,2] <= 0
 NoOverlap_10_2: x_p[2,10] + x_p[10,2] - 2 o[10,2] <= 0
 3_Supports_10_2: - h[10,2] - o[10,2] - s[10,2] <= -1
 3_DoesntSupport_10_2: h[10,2] + o[10,2] + 2 s[10,2] <= 2
 beta1_10_2: - s[10,2] + beta1[10,2] <= 0
 beta2_10_2: - s[10,2] + beta2[10,2] <= 0
 sBeta1Flag_10_2: - s[10,2] + beta1[10,2] <= 0
 sBeta2Flag_10_2: - s[10,2] + beta2[10,2] <= 0
 Eta1Overlap_10_2: eta1[10,2] + beta1[10,2] <= 1
 Eta2Overlap_10_2: eta2[10,2] + beta2[10,2] <= 1
 Eta1Flag_10_2: xi_l[2] - xi_l[10] - 300 eta1[10,2] <= 0
 Eta2Flag_10_2: - xi_r[2] + xi_r[10] - 300 eta2[10,2] <= 0
 2_Supports_10_InSameBin1_0: - p_ij[2,0] + p_ij[10,0] - s[10,2] <= 0
 2_Supports_10_InSameBin2_0: p_ij[2,0] - p_ij[10,0] - s[10,2] <= 0
 Gamma0ForItems_10_inBins_0_NoCut_2: p_ij[10,0] + gamma[10] <= 1
 2_Supports_10_InSameBin1_1: - p_ij[2,1] + p_ij[10,1] - s[10,2] <= 0
 2_Supports_10_InSameBin2_1: p_ij[2,1] - p_ij[10,1] - s[10,2] <= 0
 Gamma0ForItems_10_inBins_1_NoCut_2: p_ij[10,1] + gamma[10] <= 1
 2_Supports_10_InSameBin1_2: - p_ij[2,2] + p_ij[10,2] - s[10,2] <= 0
 2_Supports_10_InSameBin2_2: p_ij[2,2] - p_ij[10,2] - s[10,2] <= 0
 Gamma0ForItems_10_inBins_2_NoCut_2: p_ij[10,2] + gamma[10] <= 98
 2_Supports_10_InSameBin1_3: - p_ij[2,3] + p_ij[10,3] - s[10,2] <= 0
 2_Supports_10_InSameBin2_3: p_ij[2,3] - p_ij[10,3] - s[10,2] <= 0
 Gamma0ForItems_10_inBins_3_NoCut_2: p_ij[10,3] + gamma[10] <= 98
 10_OnCutIn2_2: 10000 p_ij[10,2] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 10_OnCutIn3_2: 10000 p_ij[10,3] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 AbsZ110_3: - zi[10] + zi_hi[3] - k[10,3] <= 0
 AbsZ210_3: zi[10] - zi_hi[3] - k[10,3] <= 0
 3_Below_10: zi[10] - zi_hi[3] + k[10,3] + 310 m[10,3] <= 310
 3_NotBelow_10: - zi[10] + zi_hi[3] + k[10,3] - 310 m[10,3] <= 0
 HeightCompat_10_3: h[10,3] - k[10,3] <= 0
 NoHeightCompat_10_3: - 155 h[10,3] + k[10,3] <= 0
 Overlap_10_3: - x_p[3,10] - x_p[10,3] + o[10,3] <= 0
 NoOverlap_10_3: x_p[3,10] + x_p[10,3] - 2 o[10,3] <= 0
 3_Supports_10_3: - h[10,3] - o[10,3] - s[10,3] <= -1
 3_DoesntSupport_10_3: h[10,3] + o[10,3] + 2 s[10,3] <= 2
 beta1_10_3: - s[10,3] + beta1[10,3] <= 0
 beta2_10_3: - s[10,3] + beta2[10,3] <= 0
 sBeta1Flag_10_3: - s[10,3] + beta1[10,3] <= 0
 sBeta2Flag_10_3: - s[10,3] + beta2[10,3] <= 0
 Eta1Overlap_10_3: eta1[10,3] + beta1[10,3] <= 1
 Eta2Overlap_10_3: eta2[10,3] + beta2[10,3] <= 1
 Eta1Flag_10_3: xi_l[3] - xi_l[10] - 300 eta1[10,3] <= 0
 Eta2Flag_10_3: - xi_r[3] + xi_r[10] - 300 eta2[10,3] <= 0
 3_Supports_10_InSameBin1_0: - p_ij[3,0] + p_ij[10,0] - s[10,3] <= 0
 3_Supports_10_InSameBin2_0: p_ij[3,0] - p_ij[10,0] - s[10,3] <= 0
 Gamma0ForItems_10_inBins_0_NoCut_3: p_ij[10,0] + gamma[10] <= 1
 3_Supports_10_InSameBin1_1: - p_ij[3,1] + p_ij[10,1] - s[10,3] <= 0
 3_Supports_10_InSameBin2_1: p_ij[3,1] - p_ij[10,1] - s[10,3] <= 0
 Gamma0ForItems_10_inBins_1_NoCut_3: p_ij[10,1] + gamma[10] <= 1
 3_Supports_10_InSameBin1_2: - p_ij[3,2] + p_ij[10,2] - s[10,3] <= 0
 3_Supports_10_InSameBin2_2: p_ij[3,2] - p_ij[10,2] - s[10,3] <= 0
 Gamma0ForItems_10_inBins_2_NoCut_3: p_ij[10,2] + gamma[10] <= 98
 3_Supports_10_InSameBin1_3: - p_ij[3,3] + p_ij[10,3] - s[10,3] <= 0
 3_Supports_10_InSameBin2_3: p_ij[3,3] - p_ij[10,3] - s[10,3] <= 0
 Gamma0ForItems_10_inBins_3_NoCut_3: p_ij[10,3] + gamma[10] <= 98
 10_OnCutIn2_3: 10000 p_ij[10,2] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 10_OnCutIn3_3: 10000 p_ij[10,3] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 AbsZ110_4: - zi[10] + zi_hi[4] - k[10,4] <= 0
 AbsZ210_4: zi[10] - zi_hi[4] - k[10,4] <= 0
 4_Below_10: zi[10] - zi_hi[4] + k[10,4] + 310 m[10,4] <= 310
 4_NotBelow_10: - zi[10] + zi_hi[4] + k[10,4] - 310 m[10,4] <= 0
 HeightCompat_10_4: h[10,4] - k[10,4] <= 0
 NoHeightCompat_10_4: - 155 h[10,4] + k[10,4] <= 0
 Overlap_10_4: - x_p[4,10] - x_p[10,4] + o[10,4] <= 0
 NoOverlap_10_4: x_p[4,10] + x_p[10,4] - 2 o[10,4] <= 0
 3_Supports_10_4: - h[10,4] - o[10,4] - s[10,4] <= -1
 3_DoesntSupport_10_4: h[10,4] + o[10,4] + 2 s[10,4] <= 2
 beta1_10_4: - s[10,4] + beta1[10,4] <= 0
 beta2_10_4: - s[10,4] + beta2[10,4] <= 0
 sBeta1Flag_10_4: - s[10,4] + beta1[10,4] <= 0
 sBeta2Flag_10_4: - s[10,4] + beta2[10,4] <= 0
 Eta1Overlap_10_4: eta1[10,4] + beta1[10,4] <= 1
 Eta2Overlap_10_4: eta2[10,4] + beta2[10,4] <= 1
 Eta1Flag_10_4: xi_l[4] - xi_l[10] - 300 eta1[10,4] <= 0
 Eta2Flag_10_4: - xi_r[4] + xi_r[10] - 300 eta2[10,4] <= 0
 4_Supports_10_InSameBin1_0: - p_ij[4,0] + p_ij[10,0] - s[10,4] <= 0
 4_Supports_10_InSameBin2_0: p_ij[4,0] - p_ij[10,0] - s[10,4] <= 0
 Gamma0ForItems_10_inBins_0_NoCut_4: p_ij[10,0] + gamma[10] <= 1
 4_Supports_10_InSameBin1_1: - p_ij[4,1] + p_ij[10,1] - s[10,4] <= 0
 4_Supports_10_InSameBin2_1: p_ij[4,1] - p_ij[10,1] - s[10,4] <= 0
 Gamma0ForItems_10_inBins_1_NoCut_4: p_ij[10,1] + gamma[10] <= 1
 4_Supports_10_InSameBin1_2: - p_ij[4,2] + p_ij[10,2] - s[10,4] <= 0
 4_Supports_10_InSameBin2_2: p_ij[4,2] - p_ij[10,2] - s[10,4] <= 0
 Gamma0ForItems_10_inBins_2_NoCut_4: p_ij[10,2] + gamma[10] <= 98
 4_Supports_10_InSameBin1_3: - p_ij[4,3] + p_ij[10,3] - s[10,4] <= 0
 4_Supports_10_InSameBin2_3: p_ij[4,3] - p_ij[10,3] - s[10,4] <= 0
 Gamma0ForItems_10_inBins_3_NoCut_4: p_ij[10,3] + gamma[10] <= 98
 10_OnCutIn2_4: 10000 p_ij[10,2] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 10_OnCutIn3_4: 10000 p_ij[10,3] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 AbsZ110_5: - zi[10] + zi_hi[5] - k[10,5] <= 0
 AbsZ210_5: zi[10] - zi_hi[5] - k[10,5] <= 0
 5_Below_10: zi[10] - zi_hi[5] + k[10,5] + 310 m[10,5] <= 310
 5_NotBelow_10: - zi[10] + zi_hi[5] + k[10,5] - 310 m[10,5] <= 0
 HeightCompat_10_5: h[10,5] - k[10,5] <= 0
 NoHeightCompat_10_5: - 155 h[10,5] + k[10,5] <= 0
 Overlap_10_5: - x_p[5,10] - x_p[10,5] + o[10,5] <= 0
 NoOverlap_10_5: x_p[5,10] + x_p[10,5] - 2 o[10,5] <= 0
 3_Supports_10_5: - h[10,5] - o[10,5] - s[10,5] <= -1
 3_DoesntSupport_10_5: h[10,5] + o[10,5] + 2 s[10,5] <= 2
 beta1_10_5: - s[10,5] + beta1[10,5] <= 0
 beta2_10_5: - s[10,5] + beta2[10,5] <= 0
 sBeta1Flag_10_5: - s[10,5] + beta1[10,5] <= 0
 sBeta2Flag_10_5: - s[10,5] + beta2[10,5] <= 0
 Eta1Overlap_10_5: eta1[10,5] + beta1[10,5] <= 1
 Eta2Overlap_10_5: eta2[10,5] + beta2[10,5] <= 1
 Eta1Flag_10_5: xi_l[5] - xi_l[10] - 300 eta1[10,5] <= 0
 Eta2Flag_10_5: - xi_r[5] + xi_r[10] - 300 eta2[10,5] <= 0
 5_Supports_10_InSameBin1_0: - p_ij[5,0] + p_ij[10,0] - s[10,5] <= 0
 5_Supports_10_InSameBin2_0: p_ij[5,0] - p_ij[10,0] - s[10,5] <= 0
 Gamma0ForItems_10_inBins_0_NoCut_5: p_ij[10,0] + gamma[10] <= 1
 5_Supports_10_InSameBin1_1: - p_ij[5,1] + p_ij[10,1] - s[10,5] <= 0
 5_Supports_10_InSameBin2_1: p_ij[5,1] - p_ij[10,1] - s[10,5] <= 0
 Gamma0ForItems_10_inBins_1_NoCut_5: p_ij[10,1] + gamma[10] <= 1
 5_Supports_10_InSameBin1_2: - p_ij[5,2] + p_ij[10,2] - s[10,5] <= 0
 5_Supports_10_InSameBin2_2: p_ij[5,2] - p_ij[10,2] - s[10,5] <= 0
 Gamma0ForItems_10_inBins_2_NoCut_5: p_ij[10,2] + gamma[10] <= 98
 5_Supports_10_InSameBin1_3: - p_ij[5,3] + p_ij[10,3] - s[10,5] <= 0
 5_Supports_10_InSameBin2_3: p_ij[5,3] - p_ij[10,3] - s[10,5] <= 0
 Gamma0ForItems_10_inBins_3_NoCut_5: p_ij[10,3] + gamma[10] <= 98
 10_OnCutIn2_5: 10000 p_ij[10,2] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 10_OnCutIn3_5: 10000 p_ij[10,3] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 AbsZ110_6: - zi[10] + zi_hi[6] - k[10,6] <= 0
 AbsZ210_6: zi[10] - zi_hi[6] - k[10,6] <= 0
 6_Below_10: zi[10] - zi_hi[6] + k[10,6] + 310 m[10,6] <= 310
 6_NotBelow_10: - zi[10] + zi_hi[6] + k[10,6] - 310 m[10,6] <= 0
 HeightCompat_10_6: h[10,6] - k[10,6] <= 0
 NoHeightCompat_10_6: - 155 h[10,6] + k[10,6] <= 0
 Overlap_10_6: - x_p[6,10] - x_p[10,6] + o[10,6] <= 0
 NoOverlap_10_6: x_p[6,10] + x_p[10,6] - 2 o[10,6] <= 0
 3_Supports_10_6: - h[10,6] - o[10,6] - s[10,6] <= -1
 3_DoesntSupport_10_6: h[10,6] + o[10,6] + 2 s[10,6] <= 2
 beta1_10_6: - s[10,6] + beta1[10,6] <= 0
 beta2_10_6: - s[10,6] + beta2[10,6] <= 0
 sBeta1Flag_10_6: - s[10,6] + beta1[10,6] <= 0
 sBeta2Flag_10_6: - s[10,6] + beta2[10,6] <= 0
 Eta1Overlap_10_6: eta1[10,6] + beta1[10,6] <= 1
 Eta2Overlap_10_6: eta2[10,6] + beta2[10,6] <= 1
 Eta1Flag_10_6: xi_l[6] - xi_l[10] - 300 eta1[10,6] <= 0
 Eta2Flag_10_6: - xi_r[6] + xi_r[10] - 300 eta2[10,6] <= 0
 6_Supports_10_InSameBin1_0: - p_ij[6,0] + p_ij[10,0] - s[10,6] <= 0
 6_Supports_10_InSameBin2_0: p_ij[6,0] - p_ij[10,0] - s[10,6] <= 0
 Gamma0ForItems_10_inBins_0_NoCut_6: p_ij[10,0] + gamma[10] <= 1
 6_Supports_10_InSameBin1_1: - p_ij[6,1] + p_ij[10,1] - s[10,6] <= 0
 6_Supports_10_InSameBin2_1: p_ij[6,1] - p_ij[10,1] - s[10,6] <= 0
 Gamma0ForItems_10_inBins_1_NoCut_6: p_ij[10,1] + gamma[10] <= 1
 6_Supports_10_InSameBin1_2: - p_ij[6,2] + p_ij[10,2] - s[10,6] <= 0
 6_Supports_10_InSameBin2_2: p_ij[6,2] - p_ij[10,2] - s[10,6] <= 0
 Gamma0ForItems_10_inBins_2_NoCut_6: p_ij[10,2] + gamma[10] <= 98
 6_Supports_10_InSameBin1_3: - p_ij[6,3] + p_ij[10,3] - s[10,6] <= 0
 6_Supports_10_InSameBin2_3: p_ij[6,3] - p_ij[10,3] - s[10,6] <= 0
 Gamma0ForItems_10_inBins_3_NoCut_6: p_ij[10,3] + gamma[10] <= 98
 10_OnCutIn2_6: 10000 p_ij[10,2] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 10_OnCutIn3_6: 10000 p_ij[10,3] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 AbsZ110_7: - zi[10] + zi_hi[7] - k[10,7] <= 0
 AbsZ210_7: zi[10] - zi_hi[7] - k[10,7] <= 0
 7_Below_10: zi[10] - zi_hi[7] + k[10,7] + 310 m[10,7] <= 310
 7_NotBelow_10: - zi[10] + zi_hi[7] + k[10,7] - 310 m[10,7] <= 0
 HeightCompat_10_7: h[10,7] - k[10,7] <= 0
 NoHeightCompat_10_7: - 155 h[10,7] + k[10,7] <= 0
 Overlap_10_7: - x_p[7,10] - x_p[10,7] + o[10,7] <= 0
 NoOverlap_10_7: x_p[7,10] + x_p[10,7] - 2 o[10,7] <= 0
 3_Supports_10_7: - h[10,7] - o[10,7] - s[10,7] <= -1
 3_DoesntSupport_10_7: h[10,7] + o[10,7] + 2 s[10,7] <= 2
 beta1_10_7: - s[10,7] + beta1[10,7] <= 0
 beta2_10_7: - s[10,7] + beta2[10,7] <= 0
 sBeta1Flag_10_7: - s[10,7] + beta1[10,7] <= 0
 sBeta2Flag_10_7: - s[10,7] + beta2[10,7] <= 0
 Eta1Overlap_10_7: eta1[10,7] + beta1[10,7] <= 1
 Eta2Overlap_10_7: eta2[10,7] + beta2[10,7] <= 1
 Eta1Flag_10_7: xi_l[7] - xi_l[10] - 300 eta1[10,7] <= 0
 Eta2Flag_10_7: - xi_r[7] + xi_r[10] - 300 eta2[10,7] <= 0
 7_Supports_10_InSameBin1_0: - p_ij[7,0] + p_ij[10,0] - s[10,7] <= 0
 7_Supports_10_InSameBin2_0: p_ij[7,0] - p_ij[10,0] - s[10,7] <= 0
 Gamma0ForItems_10_inBins_0_NoCut_7: p_ij[10,0] + gamma[10] <= 1
 7_Supports_10_InSameBin1_1: - p_ij[7,1] + p_ij[10,1] - s[10,7] <= 0
 7_Supports_10_InSameBin2_1: p_ij[7,1] - p_ij[10,1] - s[10,7] <= 0
 Gamma0ForItems_10_inBins_1_NoCut_7: p_ij[10,1] + gamma[10] <= 1
 7_Supports_10_InSameBin1_2: - p_ij[7,2] + p_ij[10,2] - s[10,7] <= 0
 7_Supports_10_InSameBin2_2: p_ij[7,2] - p_ij[10,2] - s[10,7] <= 0
 Gamma0ForItems_10_inBins_2_NoCut_7: p_ij[10,2] + gamma[10] <= 98
 7_Supports_10_InSameBin1_3: - p_ij[7,3] + p_ij[10,3] - s[10,7] <= 0
 7_Supports_10_InSameBin2_3: p_ij[7,3] - p_ij[10,3] - s[10,7] <= 0
 Gamma0ForItems_10_inBins_3_NoCut_7: p_ij[10,3] + gamma[10] <= 98
 10_OnCutIn2_7: 10000 p_ij[10,2] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 10_OnCutIn3_7: 10000 p_ij[10,3] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 AbsZ110_8: - zi[10] + zi_hi[8] - k[10,8] <= 0
 AbsZ210_8: zi[10] - zi_hi[8] - k[10,8] <= 0
 8_Below_10: zi[10] - zi_hi[8] + k[10,8] + 310 m[10,8] <= 310
 8_NotBelow_10: - zi[10] + zi_hi[8] + k[10,8] - 310 m[10,8] <= 0
 HeightCompat_10_8: h[10,8] - k[10,8] <= 0
 NoHeightCompat_10_8: - 155 h[10,8] + k[10,8] <= 0
 Overlap_10_8: - x_p[8,10] - x_p[10,8] + o[10,8] <= 0
 NoOverlap_10_8: x_p[8,10] + x_p[10,8] - 2 o[10,8] <= 0
 3_Supports_10_8: - h[10,8] - o[10,8] - s[10,8] <= -1
 3_DoesntSupport_10_8: h[10,8] + o[10,8] + 2 s[10,8] <= 2
 beta1_10_8: - s[10,8] + beta1[10,8] <= 0
 beta2_10_8: - s[10,8] + beta2[10,8] <= 0
 sBeta1Flag_10_8: - s[10,8] + beta1[10,8] <= 0
 sBeta2Flag_10_8: - s[10,8] + beta2[10,8] <= 0
 Eta1Overlap_10_8: eta1[10,8] + beta1[10,8] <= 1
 Eta2Overlap_10_8: eta2[10,8] + beta2[10,8] <= 1
 Eta1Flag_10_8: xi_l[8] - xi_l[10] - 300 eta1[10,8] <= 0
 Eta2Flag_10_8: - xi_r[8] + xi_r[10] - 300 eta2[10,8] <= 0
 8_Supports_10_InSameBin1_0: - p_ij[8,0] + p_ij[10,0] - s[10,8] <= 0
 8_Supports_10_InSameBin2_0: p_ij[8,0] - p_ij[10,0] - s[10,8] <= 0
 Gamma0ForItems_10_inBins_0_NoCut_8: p_ij[10,0] + gamma[10] <= 1
 8_Supports_10_InSameBin1_1: - p_ij[8,1] + p_ij[10,1] - s[10,8] <= 0
 8_Supports_10_InSameBin2_1: p_ij[8,1] - p_ij[10,1] - s[10,8] <= 0
 Gamma0ForItems_10_inBins_1_NoCut_8: p_ij[10,1] + gamma[10] <= 1
 8_Supports_10_InSameBin1_2: - p_ij[8,2] + p_ij[10,2] - s[10,8] <= 0
 8_Supports_10_InSameBin2_2: p_ij[8,2] - p_ij[10,2] - s[10,8] <= 0
 Gamma0ForItems_10_inBins_2_NoCut_8: p_ij[10,2] + gamma[10] <= 98
 8_Supports_10_InSameBin1_3: - p_ij[8,3] + p_ij[10,3] - s[10,8] <= 0
 8_Supports_10_InSameBin2_3: p_ij[8,3] - p_ij[10,3] - s[10,8] <= 0
 Gamma0ForItems_10_inBins_3_NoCut_8: p_ij[10,3] + gamma[10] <= 98
 10_OnCutIn2_8: 10000 p_ij[10,2] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 10_OnCutIn3_8: 10000 p_ij[10,3] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 AbsZ110_9: - zi[10] + zi_hi[9] - k[10,9] <= 0
 AbsZ210_9: zi[10] - zi_hi[9] - k[10,9] <= 0
 9_Below_10: zi[10] - zi_hi[9] + k[10,9] + 310 m[10,9] <= 310
 9_NotBelow_10: - zi[10] + zi_hi[9] + k[10,9] - 310 m[10,9] <= 0
 HeightCompat_10_9: h[10,9] - k[10,9] <= 0
 NoHeightCompat_10_9: - 155 h[10,9] + k[10,9] <= 0
 Overlap_10_9: - x_p[9,10] - x_p[10,9] + o[10,9] <= 0
 NoOverlap_10_9: x_p[9,10] + x_p[10,9] - 2 o[10,9] <= 0
 3_Supports_10_9: - h[10,9] - o[10,9] - s[10,9] <= -1
 3_DoesntSupport_10_9: h[10,9] + o[10,9] + 2 s[10,9] <= 2
 beta1_10_9: - s[10,9] + beta1[10,9] <= 0
 beta2_10_9: - s[10,9] + beta2[10,9] <= 0
 sBeta1Flag_10_9: - s[10,9] + beta1[10,9] <= 0
 sBeta2Flag_10_9: - s[10,9] + beta2[10,9] <= 0
 Eta1Overlap_10_9: eta1[10,9] + beta1[10,9] <= 1
 Eta2Overlap_10_9: eta2[10,9] + beta2[10,9] <= 1
 Eta1Flag_10_9: xi_l[9] - xi_l[10] - 300 eta1[10,9] <= 0
 Eta2Flag_10_9: - xi_r[9] + xi_r[10] - 300 eta2[10,9] <= 0
 9_Supports_10_InSameBin1_0: - p_ij[9,0] + p_ij[10,0] - s[10,9] <= 0
 9_Supports_10_InSameBin2_0: p_ij[9,0] - p_ij[10,0] - s[10,9] <= 0
 Gamma0ForItems_10_inBins_0_NoCut_9: p_ij[10,0] + gamma[10] <= 1
 9_Supports_10_InSameBin1_1: - p_ij[9,1] + p_ij[10,1] - s[10,9] <= 0
 9_Supports_10_InSameBin2_1: p_ij[9,1] - p_ij[10,1] - s[10,9] <= 0
 Gamma0ForItems_10_inBins_1_NoCut_9: p_ij[10,1] + gamma[10] <= 1
 9_Supports_10_InSameBin1_2: - p_ij[9,2] + p_ij[10,2] - s[10,9] <= 0
 9_Supports_10_InSameBin2_2: p_ij[9,2] - p_ij[10,2] - s[10,9] <= 0
 Gamma0ForItems_10_inBins_2_NoCut_9: p_ij[10,2] + gamma[10] <= 98
 9_Supports_10_InSameBin1_3: - p_ij[9,3] + p_ij[10,3] - s[10,9] <= 0
 9_Supports_10_InSameBin2_3: p_ij[9,3] - p_ij[10,3] - s[10,9] <= 0
 Gamma0ForItems_10_inBins_3_NoCut_9: p_ij[10,3] + gamma[10] <= 98
 10_OnCutIn2_9: 10000 p_ij[10,2] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 10_OnCutIn3_9: 10000 p_ij[10,3] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 10_OnCutIn2_10: 10000 p_ij[10,2] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 10_OnCutIn3_10: 10000 p_ij[10,3] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 AbsZ110_11: - zi[10] + zi_hi[11] - k[10,11] <= 0
 AbsZ210_11: zi[10] - zi_hi[11] - k[10,11] <= 0
 11_Below_10: zi[10] - zi_hi[11] + k[10,11] + 310 m[10,11] <= 310
 11_NotBelow_10: - zi[10] + zi_hi[11] + k[10,11] - 310 m[10,11] <= 0
 HeightCompat_10_11: h[10,11] - k[10,11] <= 0
 NoHeightCompat_10_11: - 155 h[10,11] + k[10,11] <= 0
 Overlap_10_11: - x_p[10,11] - x_p[11,10] + o[10,11] <= 0
 NoOverlap_10_11: x_p[10,11] + x_p[11,10] - 2 o[10,11] <= 0
 3_Supports_10_11: - h[10,11] - o[10,11] - s[10,11] <= -1
 3_DoesntSupport_10_11: h[10,11] + o[10,11] + 2 s[10,11] <= 2
 beta1_10_11: - s[10,11] + beta1[10,11] <= 0
 beta2_10_11: - s[10,11] + beta2[10,11] <= 0
 sBeta1Flag_10_11: - s[10,11] + beta1[10,11] <= 0
 sBeta2Flag_10_11: - s[10,11] + beta2[10,11] <= 0
 Eta1Overlap_10_11: eta1[10,11] + beta1[10,11] <= 1
 Eta2Overlap_10_11: eta2[10,11] + beta2[10,11] <= 1
 Eta1Flag_10_11: - xi_l[10] + xi_l[11] - 300 eta1[10,11] <= 0
 Eta2Flag_10_11: xi_r[10] - xi_r[11] - 300 eta2[10,11] <= 0
 11_Supports_10_InSameBin1_0: p_ij[10,0] - p_ij[11,0] - s[10,11] <= 0
 11_Supports_10_InSameBin2_0: - p_ij[10,0] + p_ij[11,0] - s[10,11] <= 0
 Gamma0ForItems_10_inBins_0_NoCut_11: p_ij[10,0] + gamma[10] <= 1
 11_Supports_10_InSameBin1_1: p_ij[10,1] - p_ij[11,1] - s[10,11] <= 0
 11_Supports_10_InSameBin2_1: - p_ij[10,1] + p_ij[11,1] - s[10,11] <= 0
 Gamma0ForItems_10_inBins_1_NoCut_11: p_ij[10,1] + gamma[10] <= 1
 11_Supports_10_InSameBin1_2: p_ij[10,2] - p_ij[11,2] - s[10,11] <= 0
 11_Supports_10_InSameBin2_2: - p_ij[10,2] + p_ij[11,2] - s[10,11] <= 0
 Gamma0ForItems_10_inBins_2_NoCut_11: p_ij[10,2] + gamma[10] <= 98
 11_Supports_10_InSameBin1_3: p_ij[10,3] - p_ij[11,3] - s[10,11] <= 0
 11_Supports_10_InSameBin2_3: - p_ij[10,3] + p_ij[11,3] - s[10,11] <= 0
 Gamma0ForItems_10_inBins_3_NoCut_11: p_ij[10,3] + gamma[10] <= 98
 10_OnCutIn2_11: 10000 p_ij[10,2] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 10_OnCutIn3_11: 10000 p_ij[10,3] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 AbsZ110_12: - zi[10] + zi_hi[12] - k[10,12] <= 0
 AbsZ210_12: zi[10] - zi_hi[12] - k[10,12] <= 0
 12_Below_10: zi[10] - zi_hi[12] + k[10,12] + 310 m[10,12] <= 310
 12_NotBelow_10: - zi[10] + zi_hi[12] + k[10,12] - 310 m[10,12] <= 0
 HeightCompat_10_12: h[10,12] - k[10,12] <= 0
 NoHeightCompat_10_12: - 155 h[10,12] + k[10,12] <= 0
 Overlap_10_12: - x_p[10,12] - x_p[12,10] + o[10,12] <= 0
 NoOverlap_10_12: x_p[10,12] + x_p[12,10] - 2 o[10,12] <= 0
 3_Supports_10_12: - h[10,12] - o[10,12] - s[10,12] <= -1
 3_DoesntSupport_10_12: h[10,12] + o[10,12] + 2 s[10,12] <= 2
 beta1_10_12: - s[10,12] + beta1[10,12] <= 0
 beta2_10_12: - s[10,12] + beta2[10,12] <= 0
 sBeta1Flag_10_12: - s[10,12] + beta1[10,12] <= 0
 sBeta2Flag_10_12: - s[10,12] + beta2[10,12] <= 0
 Eta1Overlap_10_12: eta1[10,12] + beta1[10,12] <= 1
 Eta2Overlap_10_12: eta2[10,12] + beta2[10,12] <= 1
 Eta1Flag_10_12: - xi_l[10] + xi_l[12] - 300 eta1[10,12] <= 0
 Eta2Flag_10_12: xi_r[10] - xi_r[12] - 300 eta2[10,12] <= 0
 12_Supports_10_InSameBin1_0: p_ij[10,0] - p_ij[12,0] - s[10,12] <= 0
 12_Supports_10_InSameBin2_0: - p_ij[10,0] + p_ij[12,0] - s[10,12] <= 0
 Gamma0ForItems_10_inBins_0_NoCut_12: p_ij[10,0] + gamma[10] <= 1
 12_Supports_10_InSameBin1_1: p_ij[10,1] - p_ij[12,1] - s[10,12] <= 0
 12_Supports_10_InSameBin2_1: - p_ij[10,1] + p_ij[12,1] - s[10,12] <= 0
 Gamma0ForItems_10_inBins_1_NoCut_12: p_ij[10,1] + gamma[10] <= 1
 12_Supports_10_InSameBin1_2: p_ij[10,2] - p_ij[12,2] - s[10,12] <= 0
 12_Supports_10_InSameBin2_2: - p_ij[10,2] + p_ij[12,2] - s[10,12] <= 0
 Gamma0ForItems_10_inBins_2_NoCut_12: p_ij[10,2] + gamma[10] <= 98
 12_Supports_10_InSameBin1_3: p_ij[10,3] - p_ij[12,3] - s[10,12] <= 0
 12_Supports_10_InSameBin2_3: - p_ij[10,3] + p_ij[12,3] - s[10,12] <= 0
 Gamma0ForItems_10_inBins_3_NoCut_12: p_ij[10,3] + gamma[10] <= 98
 10_OnCutIn2_12: 10000 p_ij[10,2] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 10_OnCutIn3_12: 10000 p_ij[10,3] + 1.261904761904762 xi_l[10] + zi[10]
   + 10000 gamma[10] <= 20053
 GroundConstraint_11: zi[11] + 155 g[11] <= 155
 AbsZ111_0: - zi[11] + zi_hi[0] - k[11,0] <= 0
 AbsZ211_0: zi[11] - zi_hi[0] - k[11,0] <= 0
 0_Below_11: zi[11] - zi_hi[0] + k[11,0] + 310 m[11,0] <= 310
 0_NotBelow_11: - zi[11] + zi_hi[0] + k[11,0] - 310 m[11,0] <= 0
 HeightCompat_11_0: h[11,0] - k[11,0] <= 0
 NoHeightCompat_11_0: - 155 h[11,0] + k[11,0] <= 0
 Overlap_11_0: - x_p[0,11] - x_p[11,0] + o[11,0] <= 0
 NoOverlap_11_0: x_p[0,11] + x_p[11,0] - 2 o[11,0] <= 0
 3_Supports_11_0: - h[11,0] - o[11,0] - s[11,0] <= -1
 3_DoesntSupport_11_0: h[11,0] + o[11,0] + 2 s[11,0] <= 2
 beta1_11_0: - s[11,0] + beta1[11,0] <= 0
 beta2_11_0: - s[11,0] + beta2[11,0] <= 0
 sBeta1Flag_11_0: - s[11,0] + beta1[11,0] <= 0
 sBeta2Flag_11_0: - s[11,0] + beta2[11,0] <= 0
 Eta1Overlap_11_0: eta1[11,0] + beta1[11,0] <= 1
 Eta2Overlap_11_0: eta2[11,0] + beta2[11,0] <= 1
 Eta1Flag_11_0: xi_l[0] - xi_l[11] - 300 eta1[11,0] <= 0
 Eta2Flag_11_0: - xi_r[0] + xi_r[11] - 300 eta2[11,0] <= 0
 0_Supports_11_InSameBin1_0: - p_ij[0,0] + p_ij[11,0] - s[11,0] <= 0
 0_Supports_11_InSameBin2_0: p_ij[0,0] - p_ij[11,0] - s[11,0] <= 0
 Gamma0ForItems_11_inBins_0_NoCut_0: p_ij[11,0] + gamma[11] <= 1
 0_Supports_11_InSameBin1_1: - p_ij[0,1] + p_ij[11,1] - s[11,0] <= 0
 0_Supports_11_InSameBin2_1: p_ij[0,1] - p_ij[11,1] - s[11,0] <= 0
 Gamma0ForItems_11_inBins_1_NoCut_0: p_ij[11,1] + gamma[11] <= 1
 0_Supports_11_InSameBin1_2: - p_ij[0,2] + p_ij[11,2] - s[11,0] <= 0
 0_Supports_11_InSameBin2_2: p_ij[0,2] - p_ij[11,2] - s[11,0] <= 0
 Gamma0ForItems_11_inBins_2_NoCut_0: p_ij[11,2] + gamma[11] <= 98
 0_Supports_11_InSameBin1_3: - p_ij[0,3] + p_ij[11,3] - s[11,0] <= 0
 0_Supports_11_InSameBin2_3: p_ij[0,3] - p_ij[11,3] - s[11,0] <= 0
 Gamma0ForItems_11_inBins_3_NoCut_0: p_ij[11,3] + gamma[11] <= 98
 11_OnCutIn2_0: 10000 p_ij[11,2] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 11_OnCutIn3_0: 10000 p_ij[11,3] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 AbsZ111_1: - zi[11] + zi_hi[1] - k[11,1] <= 0
 AbsZ211_1: zi[11] - zi_hi[1] - k[11,1] <= 0
 1_Below_11: zi[11] - zi_hi[1] + k[11,1] + 310 m[11,1] <= 310
 1_NotBelow_11: - zi[11] + zi_hi[1] + k[11,1] - 310 m[11,1] <= 0
 HeightCompat_11_1: h[11,1] - k[11,1] <= 0
 NoHeightCompat_11_1: - 155 h[11,1] + k[11,1] <= 0
 Overlap_11_1: - x_p[1,11] - x_p[11,1] + o[11,1] <= 0
 NoOverlap_11_1: x_p[1,11] + x_p[11,1] - 2 o[11,1] <= 0
 3_Supports_11_1: - h[11,1] - o[11,1] - s[11,1] <= -1
 3_DoesntSupport_11_1: h[11,1] + o[11,1] + 2 s[11,1] <= 2
 beta1_11_1: - s[11,1] + beta1[11,1] <= 0
 beta2_11_1: - s[11,1] + beta2[11,1] <= 0
 sBeta1Flag_11_1: - s[11,1] + beta1[11,1] <= 0
 sBeta2Flag_11_1: - s[11,1] + beta2[11,1] <= 0
 Eta1Overlap_11_1: eta1[11,1] + beta1[11,1] <= 1
 Eta2Overlap_11_1: eta2[11,1] + beta2[11,1] <= 1
 Eta1Flag_11_1: xi_l[1] - xi_l[11] - 300 eta1[11,1] <= 0
 Eta2Flag_11_1: - xi_r[1] + xi_r[11] - 300 eta2[11,1] <= 0
 1_Supports_11_InSameBin1_0: - p_ij[1,0] + p_ij[11,0] - s[11,1] <= 0
 1_Supports_11_InSameBin2_0: p_ij[1,0] - p_ij[11,0] - s[11,1] <= 0
 Gamma0ForItems_11_inBins_0_NoCut_1: p_ij[11,0] + gamma[11] <= 1
 1_Supports_11_InSameBin1_1: - p_ij[1,1] + p_ij[11,1] - s[11,1] <= 0
 1_Supports_11_InSameBin2_1: p_ij[1,1] - p_ij[11,1] - s[11,1] <= 0
 Gamma0ForItems_11_inBins_1_NoCut_1: p_ij[11,1] + gamma[11] <= 1
 1_Supports_11_InSameBin1_2: - p_ij[1,2] + p_ij[11,2] - s[11,1] <= 0
 1_Supports_11_InSameBin2_2: p_ij[1,2] - p_ij[11,2] - s[11,1] <= 0
 Gamma0ForItems_11_inBins_2_NoCut_1: p_ij[11,2] + gamma[11] <= 98
 1_Supports_11_InSameBin1_3: - p_ij[1,3] + p_ij[11,3] - s[11,1] <= 0
 1_Supports_11_InSameBin2_3: p_ij[1,3] - p_ij[11,3] - s[11,1] <= 0
 Gamma0ForItems_11_inBins_3_NoCut_1: p_ij[11,3] + gamma[11] <= 98
 11_OnCutIn2_1: 10000 p_ij[11,2] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 11_OnCutIn3_1: 10000 p_ij[11,3] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 AbsZ111_2: - zi[11] + zi_hi[2] - k[11,2] <= 0
 AbsZ211_2: zi[11] - zi_hi[2] - k[11,2] <= 0
 2_Below_11: zi[11] - zi_hi[2] + k[11,2] + 310 m[11,2] <= 310
 2_NotBelow_11: - zi[11] + zi_hi[2] + k[11,2] - 310 m[11,2] <= 0
 HeightCompat_11_2: h[11,2] - k[11,2] <= 0
 NoHeightCompat_11_2: - 155 h[11,2] + k[11,2] <= 0
 Overlap_11_2: - x_p[2,11] - x_p[11,2] + o[11,2] <= 0
 NoOverlap_11_2: x_p[2,11] + x_p[11,2] - 2 o[11,2] <= 0
 3_Supports_11_2: - h[11,2] - o[11,2] - s[11,2] <= -1
 3_DoesntSupport_11_2: h[11,2] + o[11,2] + 2 s[11,2] <= 2
 beta1_11_2: - s[11,2] + beta1[11,2] <= 0
 beta2_11_2: - s[11,2] + beta2[11,2] <= 0
 sBeta1Flag_11_2: - s[11,2] + beta1[11,2] <= 0
 sBeta2Flag_11_2: - s[11,2] + beta2[11,2] <= 0
 Eta1Overlap_11_2: eta1[11,2] + beta1[11,2] <= 1
 Eta2Overlap_11_2: eta2[11,2] + beta2[11,2] <= 1
 Eta1Flag_11_2: xi_l[2] - xi_l[11] - 300 eta1[11,2] <= 0
 Eta2Flag_11_2: - xi_r[2] + xi_r[11] - 300 eta2[11,2] <= 0
 2_Supports_11_InSameBin1_0: - p_ij[2,0] + p_ij[11,0] - s[11,2] <= 0
 2_Supports_11_InSameBin2_0: p_ij[2,0] - p_ij[11,0] - s[11,2] <= 0
 Gamma0ForItems_11_inBins_0_NoCut_2: p_ij[11,0] + gamma[11] <= 1
 2_Supports_11_InSameBin1_1: - p_ij[2,1] + p_ij[11,1] - s[11,2] <= 0
 2_Supports_11_InSameBin2_1: p_ij[2,1] - p_ij[11,1] - s[11,2] <= 0
 Gamma0ForItems_11_inBins_1_NoCut_2: p_ij[11,1] + gamma[11] <= 1
 2_Supports_11_InSameBin1_2: - p_ij[2,2] + p_ij[11,2] - s[11,2] <= 0
 2_Supports_11_InSameBin2_2: p_ij[2,2] - p_ij[11,2] - s[11,2] <= 0
 Gamma0ForItems_11_inBins_2_NoCut_2: p_ij[11,2] + gamma[11] <= 98
 2_Supports_11_InSameBin1_3: - p_ij[2,3] + p_ij[11,3] - s[11,2] <= 0
 2_Supports_11_InSameBin2_3: p_ij[2,3] - p_ij[11,3] - s[11,2] <= 0
 Gamma0ForItems_11_inBins_3_NoCut_2: p_ij[11,3] + gamma[11] <= 98
 11_OnCutIn2_2: 10000 p_ij[11,2] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 11_OnCutIn3_2: 10000 p_ij[11,3] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 AbsZ111_3: - zi[11] + zi_hi[3] - k[11,3] <= 0
 AbsZ211_3: zi[11] - zi_hi[3] - k[11,3] <= 0
 3_Below_11: zi[11] - zi_hi[3] + k[11,3] + 310 m[11,3] <= 310
 3_NotBelow_11: - zi[11] + zi_hi[3] + k[11,3] - 310 m[11,3] <= 0
 HeightCompat_11_3: h[11,3] - k[11,3] <= 0
 NoHeightCompat_11_3: - 155 h[11,3] + k[11,3] <= 0
 Overlap_11_3: - x_p[3,11] - x_p[11,3] + o[11,3] <= 0
 NoOverlap_11_3: x_p[3,11] + x_p[11,3] - 2 o[11,3] <= 0
 3_Supports_11_3: - h[11,3] - o[11,3] - s[11,3] <= -1
 3_DoesntSupport_11_3: h[11,3] + o[11,3] + 2 s[11,3] <= 2
 beta1_11_3: - s[11,3] + beta1[11,3] <= 0
 beta2_11_3: - s[11,3] + beta2[11,3] <= 0
 sBeta1Flag_11_3: - s[11,3] + beta1[11,3] <= 0
 sBeta2Flag_11_3: - s[11,3] + beta2[11,3] <= 0
 Eta1Overlap_11_3: eta1[11,3] + beta1[11,3] <= 1
 Eta2Overlap_11_3: eta2[11,3] + beta2[11,3] <= 1
 Eta1Flag_11_3: xi_l[3] - xi_l[11] - 300 eta1[11,3] <= 0
 Eta2Flag_11_3: - xi_r[3] + xi_r[11] - 300 eta2[11,3] <= 0
 3_Supports_11_InSameBin1_0: - p_ij[3,0] + p_ij[11,0] - s[11,3] <= 0
 3_Supports_11_InSameBin2_0: p_ij[3,0] - p_ij[11,0] - s[11,3] <= 0
 Gamma0ForItems_11_inBins_0_NoCut_3: p_ij[11,0] + gamma[11] <= 1
 3_Supports_11_InSameBin1_1: - p_ij[3,1] + p_ij[11,1] - s[11,3] <= 0
 3_Supports_11_InSameBin2_1: p_ij[3,1] - p_ij[11,1] - s[11,3] <= 0
 Gamma0ForItems_11_inBins_1_NoCut_3: p_ij[11,1] + gamma[11] <= 1
 3_Supports_11_InSameBin1_2: - p_ij[3,2] + p_ij[11,2] - s[11,3] <= 0
 3_Supports_11_InSameBin2_2: p_ij[3,2] - p_ij[11,2] - s[11,3] <= 0
 Gamma0ForItems_11_inBins_2_NoCut_3: p_ij[11,2] + gamma[11] <= 98
 3_Supports_11_InSameBin1_3: - p_ij[3,3] + p_ij[11,3] - s[11,3] <= 0
 3_Supports_11_InSameBin2_3: p_ij[3,3] - p_ij[11,3] - s[11,3] <= 0
 Gamma0ForItems_11_inBins_3_NoCut_3: p_ij[11,3] + gamma[11] <= 98
 11_OnCutIn2_3: 10000 p_ij[11,2] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 11_OnCutIn3_3: 10000 p_ij[11,3] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 AbsZ111_4: - zi[11] + zi_hi[4] - k[11,4] <= 0
 AbsZ211_4: zi[11] - zi_hi[4] - k[11,4] <= 0
 4_Below_11: zi[11] - zi_hi[4] + k[11,4] + 310 m[11,4] <= 310
 4_NotBelow_11: - zi[11] + zi_hi[4] + k[11,4] - 310 m[11,4] <= 0
 HeightCompat_11_4: h[11,4] - k[11,4] <= 0
 NoHeightCompat_11_4: - 155 h[11,4] + k[11,4] <= 0
 Overlap_11_4: - x_p[4,11] - x_p[11,4] + o[11,4] <= 0
 NoOverlap_11_4: x_p[4,11] + x_p[11,4] - 2 o[11,4] <= 0
 3_Supports_11_4: - h[11,4] - o[11,4] - s[11,4] <= -1
 3_DoesntSupport_11_4: h[11,4] + o[11,4] + 2 s[11,4] <= 2
 beta1_11_4: - s[11,4] + beta1[11,4] <= 0
 beta2_11_4: - s[11,4] + beta2[11,4] <= 0
 sBeta1Flag_11_4: - s[11,4] + beta1[11,4] <= 0
 sBeta2Flag_11_4: - s[11,4] + beta2[11,4] <= 0
 Eta1Overlap_11_4: eta1[11,4] + beta1[11,4] <= 1
 Eta2Overlap_11_4: eta2[11,4] + beta2[11,4] <= 1
 Eta1Flag_11_4: xi_l[4] - xi_l[11] - 300 eta1[11,4] <= 0
 Eta2Flag_11_4: - xi_r[4] + xi_r[11] - 300 eta2[11,4] <= 0
 4_Supports_11_InSameBin1_0: - p_ij[4,0] + p_ij[11,0] - s[11,4] <= 0
 4_Supports_11_InSameBin2_0: p_ij[4,0] - p_ij[11,0] - s[11,4] <= 0
 Gamma0ForItems_11_inBins_0_NoCut_4: p_ij[11,0] + gamma[11] <= 1
 4_Supports_11_InSameBin1_1: - p_ij[4,1] + p_ij[11,1] - s[11,4] <= 0
 4_Supports_11_InSameBin2_1: p_ij[4,1] - p_ij[11,1] - s[11,4] <= 0
 Gamma0ForItems_11_inBins_1_NoCut_4: p_ij[11,1] + gamma[11] <= 1
 4_Supports_11_InSameBin1_2: - p_ij[4,2] + p_ij[11,2] - s[11,4] <= 0
 4_Supports_11_InSameBin2_2: p_ij[4,2] - p_ij[11,2] - s[11,4] <= 0
 Gamma0ForItems_11_inBins_2_NoCut_4: p_ij[11,2] + gamma[11] <= 98
 4_Supports_11_InSameBin1_3: - p_ij[4,3] + p_ij[11,3] - s[11,4] <= 0
 4_Supports_11_InSameBin2_3: p_ij[4,3] - p_ij[11,3] - s[11,4] <= 0
 Gamma0ForItems_11_inBins_3_NoCut_4: p_ij[11,3] + gamma[11] <= 98
 11_OnCutIn2_4: 10000 p_ij[11,2] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 11_OnCutIn3_4: 10000 p_ij[11,3] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 AbsZ111_5: - zi[11] + zi_hi[5] - k[11,5] <= 0
 AbsZ211_5: zi[11] - zi_hi[5] - k[11,5] <= 0
 5_Below_11: zi[11] - zi_hi[5] + k[11,5] + 310 m[11,5] <= 310
 5_NotBelow_11: - zi[11] + zi_hi[5] + k[11,5] - 310 m[11,5] <= 0
 HeightCompat_11_5: h[11,5] - k[11,5] <= 0
 NoHeightCompat_11_5: - 155 h[11,5] + k[11,5] <= 0
 Overlap_11_5: - x_p[5,11] - x_p[11,5] + o[11,5] <= 0
 NoOverlap_11_5: x_p[5,11] + x_p[11,5] - 2 o[11,5] <= 0
 3_Supports_11_5: - h[11,5] - o[11,5] - s[11,5] <= -1
 3_DoesntSupport_11_5: h[11,5] + o[11,5] + 2 s[11,5] <= 2
 beta1_11_5: - s[11,5] + beta1[11,5] <= 0
 beta2_11_5: - s[11,5] + beta2[11,5] <= 0
 sBeta1Flag_11_5: - s[11,5] + beta1[11,5] <= 0
 sBeta2Flag_11_5: - s[11,5] + beta2[11,5] <= 0
 Eta1Overlap_11_5: eta1[11,5] + beta1[11,5] <= 1
 Eta2Overlap_11_5: eta2[11,5] + beta2[11,5] <= 1
 Eta1Flag_11_5: xi_l[5] - xi_l[11] - 300 eta1[11,5] <= 0
 Eta2Flag_11_5: - xi_r[5] + xi_r[11] - 300 eta2[11,5] <= 0
 5_Supports_11_InSameBin1_0: - p_ij[5,0] + p_ij[11,0] - s[11,5] <= 0
 5_Supports_11_InSameBin2_0: p_ij[5,0] - p_ij[11,0] - s[11,5] <= 0
 Gamma0ForItems_11_inBins_0_NoCut_5: p_ij[11,0] + gamma[11] <= 1
 5_Supports_11_InSameBin1_1: - p_ij[5,1] + p_ij[11,1] - s[11,5] <= 0
 5_Supports_11_InSameBin2_1: p_ij[5,1] - p_ij[11,1] - s[11,5] <= 0
 Gamma0ForItems_11_inBins_1_NoCut_5: p_ij[11,1] + gamma[11] <= 1
 5_Supports_11_InSameBin1_2: - p_ij[5,2] + p_ij[11,2] - s[11,5] <= 0
 5_Supports_11_InSameBin2_2: p_ij[5,2] - p_ij[11,2] - s[11,5] <= 0
 Gamma0ForItems_11_inBins_2_NoCut_5: p_ij[11,2] + gamma[11] <= 98
 5_Supports_11_InSameBin1_3: - p_ij[5,3] + p_ij[11,3] - s[11,5] <= 0
 5_Supports_11_InSameBin2_3: p_ij[5,3] - p_ij[11,3] - s[11,5] <= 0
 Gamma0ForItems_11_inBins_3_NoCut_5: p_ij[11,3] + gamma[11] <= 98
 11_OnCutIn2_5: 10000 p_ij[11,2] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 11_OnCutIn3_5: 10000 p_ij[11,3] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 AbsZ111_6: - zi[11] + zi_hi[6] - k[11,6] <= 0
 AbsZ211_6: zi[11] - zi_hi[6] - k[11,6] <= 0
 6_Below_11: zi[11] - zi_hi[6] + k[11,6] + 310 m[11,6] <= 310
 6_NotBelow_11: - zi[11] + zi_hi[6] + k[11,6] - 310 m[11,6] <= 0
 HeightCompat_11_6: h[11,6] - k[11,6] <= 0
 NoHeightCompat_11_6: - 155 h[11,6] + k[11,6] <= 0
 Overlap_11_6: - x_p[6,11] - x_p[11,6] + o[11,6] <= 0
 NoOverlap_11_6: x_p[6,11] + x_p[11,6] - 2 o[11,6] <= 0
 3_Supports_11_6: - h[11,6] - o[11,6] - s[11,6] <= -1
 3_DoesntSupport_11_6: h[11,6] + o[11,6] + 2 s[11,6] <= 2
 beta1_11_6: - s[11,6] + beta1[11,6] <= 0
 beta2_11_6: - s[11,6] + beta2[11,6] <= 0
 sBeta1Flag_11_6: - s[11,6] + beta1[11,6] <= 0
 sBeta2Flag_11_6: - s[11,6] + beta2[11,6] <= 0
 Eta1Overlap_11_6: eta1[11,6] + beta1[11,6] <= 1
 Eta2Overlap_11_6: eta2[11,6] + beta2[11,6] <= 1
 Eta1Flag_11_6: xi_l[6] - xi_l[11] - 300 eta1[11,6] <= 0
 Eta2Flag_11_6: - xi_r[6] + xi_r[11] - 300 eta2[11,6] <= 0
 6_Supports_11_InSameBin1_0: - p_ij[6,0] + p_ij[11,0] - s[11,6] <= 0
 6_Supports_11_InSameBin2_0: p_ij[6,0] - p_ij[11,0] - s[11,6] <= 0
 Gamma0ForItems_11_inBins_0_NoCut_6: p_ij[11,0] + gamma[11] <= 1
 6_Supports_11_InSameBin1_1: - p_ij[6,1] + p_ij[11,1] - s[11,6] <= 0
 6_Supports_11_InSameBin2_1: p_ij[6,1] - p_ij[11,1] - s[11,6] <= 0
 Gamma0ForItems_11_inBins_1_NoCut_6: p_ij[11,1] + gamma[11] <= 1
 6_Supports_11_InSameBin1_2: - p_ij[6,2] + p_ij[11,2] - s[11,6] <= 0
 6_Supports_11_InSameBin2_2: p_ij[6,2] - p_ij[11,2] - s[11,6] <= 0
 Gamma0ForItems_11_inBins_2_NoCut_6: p_ij[11,2] + gamma[11] <= 98
 6_Supports_11_InSameBin1_3: - p_ij[6,3] + p_ij[11,3] - s[11,6] <= 0
 6_Supports_11_InSameBin2_3: p_ij[6,3] - p_ij[11,3] - s[11,6] <= 0
 Gamma0ForItems_11_inBins_3_NoCut_6: p_ij[11,3] + gamma[11] <= 98
 11_OnCutIn2_6: 10000 p_ij[11,2] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 11_OnCutIn3_6: 10000 p_ij[11,3] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 AbsZ111_7: - zi[11] + zi_hi[7] - k[11,7] <= 0
 AbsZ211_7: zi[11] - zi_hi[7] - k[11,7] <= 0
 7_Below_11: zi[11] - zi_hi[7] + k[11,7] + 310 m[11,7] <= 310
 7_NotBelow_11: - zi[11] + zi_hi[7] + k[11,7] - 310 m[11,7] <= 0
 HeightCompat_11_7: h[11,7] - k[11,7] <= 0
 NoHeightCompat_11_7: - 155 h[11,7] + k[11,7] <= 0
 Overlap_11_7: - x_p[7,11] - x_p[11,7] + o[11,7] <= 0
 NoOverlap_11_7: x_p[7,11] + x_p[11,7] - 2 o[11,7] <= 0
 3_Supports_11_7: - h[11,7] - o[11,7] - s[11,7] <= -1
 3_DoesntSupport_11_7: h[11,7] + o[11,7] + 2 s[11,7] <= 2
 beta1_11_7: - s[11,7] + beta1[11,7] <= 0
 beta2_11_7: - s[11,7] + beta2[11,7] <= 0
 sBeta1Flag_11_7: - s[11,7] + beta1[11,7] <= 0
 sBeta2Flag_11_7: - s[11,7] + beta2[11,7] <= 0
 Eta1Overlap_11_7: eta1[11,7] + beta1[11,7] <= 1
 Eta2Overlap_11_7: eta2[11,7] + beta2[11,7] <= 1
 Eta1Flag_11_7: xi_l[7] - xi_l[11] - 300 eta1[11,7] <= 0
 Eta2Flag_11_7: - xi_r[7] + xi_r[11] - 300 eta2[11,7] <= 0
 7_Supports_11_InSameBin1_0: - p_ij[7,0] + p_ij[11,0] - s[11,7] <= 0
 7_Supports_11_InSameBin2_0: p_ij[7,0] - p_ij[11,0] - s[11,7] <= 0
 Gamma0ForItems_11_inBins_0_NoCut_7: p_ij[11,0] + gamma[11] <= 1
 7_Supports_11_InSameBin1_1: - p_ij[7,1] + p_ij[11,1] - s[11,7] <= 0
 7_Supports_11_InSameBin2_1: p_ij[7,1] - p_ij[11,1] - s[11,7] <= 0
 Gamma0ForItems_11_inBins_1_NoCut_7: p_ij[11,1] + gamma[11] <= 1
 7_Supports_11_InSameBin1_2: - p_ij[7,2] + p_ij[11,2] - s[11,7] <= 0
 7_Supports_11_InSameBin2_2: p_ij[7,2] - p_ij[11,2] - s[11,7] <= 0
 Gamma0ForItems_11_inBins_2_NoCut_7: p_ij[11,2] + gamma[11] <= 98
 7_Supports_11_InSameBin1_3: - p_ij[7,3] + p_ij[11,3] - s[11,7] <= 0
 7_Supports_11_InSameBin2_3: p_ij[7,3] - p_ij[11,3] - s[11,7] <= 0
 Gamma0ForItems_11_inBins_3_NoCut_7: p_ij[11,3] + gamma[11] <= 98
 11_OnCutIn2_7: 10000 p_ij[11,2] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 11_OnCutIn3_7: 10000 p_ij[11,3] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 AbsZ111_8: - zi[11] + zi_hi[8] - k[11,8] <= 0
 AbsZ211_8: zi[11] - zi_hi[8] - k[11,8] <= 0
 8_Below_11: zi[11] - zi_hi[8] + k[11,8] + 310 m[11,8] <= 310
 8_NotBelow_11: - zi[11] + zi_hi[8] + k[11,8] - 310 m[11,8] <= 0
 HeightCompat_11_8: h[11,8] - k[11,8] <= 0
 NoHeightCompat_11_8: - 155 h[11,8] + k[11,8] <= 0
 Overlap_11_8: - x_p[8,11] - x_p[11,8] + o[11,8] <= 0
 NoOverlap_11_8: x_p[8,11] + x_p[11,8] - 2 o[11,8] <= 0
 3_Supports_11_8: - h[11,8] - o[11,8] - s[11,8] <= -1
 3_DoesntSupport_11_8: h[11,8] + o[11,8] + 2 s[11,8] <= 2
 beta1_11_8: - s[11,8] + beta1[11,8] <= 0
 beta2_11_8: - s[11,8] + beta2[11,8] <= 0
 sBeta1Flag_11_8: - s[11,8] + beta1[11,8] <= 0
 sBeta2Flag_11_8: - s[11,8] + beta2[11,8] <= 0
 Eta1Overlap_11_8: eta1[11,8] + beta1[11,8] <= 1
 Eta2Overlap_11_8: eta2[11,8] + beta2[11,8] <= 1
 Eta1Flag_11_8: xi_l[8] - xi_l[11] - 300 eta1[11,8] <= 0
 Eta2Flag_11_8: - xi_r[8] + xi_r[11] - 300 eta2[11,8] <= 0
 8_Supports_11_InSameBin1_0: - p_ij[8,0] + p_ij[11,0] - s[11,8] <= 0
 8_Supports_11_InSameBin2_0: p_ij[8,0] - p_ij[11,0] - s[11,8] <= 0
 Gamma0ForItems_11_inBins_0_NoCut_8: p_ij[11,0] + gamma[11] <= 1
 8_Supports_11_InSameBin1_1: - p_ij[8,1] + p_ij[11,1] - s[11,8] <= 0
 8_Supports_11_InSameBin2_1: p_ij[8,1] - p_ij[11,1] - s[11,8] <= 0
 Gamma0ForItems_11_inBins_1_NoCut_8: p_ij[11,1] + gamma[11] <= 1
 8_Supports_11_InSameBin1_2: - p_ij[8,2] + p_ij[11,2] - s[11,8] <= 0
 8_Supports_11_InSameBin2_2: p_ij[8,2] - p_ij[11,2] - s[11,8] <= 0
 Gamma0ForItems_11_inBins_2_NoCut_8: p_ij[11,2] + gamma[11] <= 98
 8_Supports_11_InSameBin1_3: - p_ij[8,3] + p_ij[11,3] - s[11,8] <= 0
 8_Supports_11_InSameBin2_3: p_ij[8,3] - p_ij[11,3] - s[11,8] <= 0
 Gamma0ForItems_11_inBins_3_NoCut_8: p_ij[11,3] + gamma[11] <= 98
 11_OnCutIn2_8: 10000 p_ij[11,2] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 11_OnCutIn3_8: 10000 p_ij[11,3] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 AbsZ111_9: - zi[11] + zi_hi[9] - k[11,9] <= 0
 AbsZ211_9: zi[11] - zi_hi[9] - k[11,9] <= 0
 9_Below_11: zi[11] - zi_hi[9] + k[11,9] + 310 m[11,9] <= 310
 9_NotBelow_11: - zi[11] + zi_hi[9] + k[11,9] - 310 m[11,9] <= 0
 HeightCompat_11_9: h[11,9] - k[11,9] <= 0
 NoHeightCompat_11_9: - 155 h[11,9] + k[11,9] <= 0
 Overlap_11_9: - x_p[9,11] - x_p[11,9] + o[11,9] <= 0
 NoOverlap_11_9: x_p[9,11] + x_p[11,9] - 2 o[11,9] <= 0
 3_Supports_11_9: - h[11,9] - o[11,9] - s[11,9] <= -1
 3_DoesntSupport_11_9: h[11,9] + o[11,9] + 2 s[11,9] <= 2
 beta1_11_9: - s[11,9] + beta1[11,9] <= 0
 beta2_11_9: - s[11,9] + beta2[11,9] <= 0
 sBeta1Flag_11_9: - s[11,9] + beta1[11,9] <= 0
 sBeta2Flag_11_9: - s[11,9] + beta2[11,9] <= 0
 Eta1Overlap_11_9: eta1[11,9] + beta1[11,9] <= 1
 Eta2Overlap_11_9: eta2[11,9] + beta2[11,9] <= 1
 Eta1Flag_11_9: xi_l[9] - xi_l[11] - 300 eta1[11,9] <= 0
 Eta2Flag_11_9: - xi_r[9] + xi_r[11] - 300 eta2[11,9] <= 0
 9_Supports_11_InSameBin1_0: - p_ij[9,0] + p_ij[11,0] - s[11,9] <= 0
 9_Supports_11_InSameBin2_0: p_ij[9,0] - p_ij[11,0] - s[11,9] <= 0
 Gamma0ForItems_11_inBins_0_NoCut_9: p_ij[11,0] + gamma[11] <= 1
 9_Supports_11_InSameBin1_1: - p_ij[9,1] + p_ij[11,1] - s[11,9] <= 0
 9_Supports_11_InSameBin2_1: p_ij[9,1] - p_ij[11,1] - s[11,9] <= 0
 Gamma0ForItems_11_inBins_1_NoCut_9: p_ij[11,1] + gamma[11] <= 1
 9_Supports_11_InSameBin1_2: - p_ij[9,2] + p_ij[11,2] - s[11,9] <= 0
 9_Supports_11_InSameBin2_2: p_ij[9,2] - p_ij[11,2] - s[11,9] <= 0
 Gamma0ForItems_11_inBins_2_NoCut_9: p_ij[11,2] + gamma[11] <= 98
 9_Supports_11_InSameBin1_3: - p_ij[9,3] + p_ij[11,3] - s[11,9] <= 0
 9_Supports_11_InSameBin2_3: p_ij[9,3] - p_ij[11,3] - s[11,9] <= 0
 Gamma0ForItems_11_inBins_3_NoCut_9: p_ij[11,3] + gamma[11] <= 98
 11_OnCutIn2_9: 10000 p_ij[11,2] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 11_OnCutIn3_9: 10000 p_ij[11,3] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 AbsZ111_10: - zi[11] + zi_hi[10] - k[11,10] <= 0
 AbsZ211_10: zi[11] - zi_hi[10] - k[11,10] <= 0
 10_Below_11: zi[11] - zi_hi[10] + k[11,10] + 310 m[11,10] <= 310
 10_NotBelow_11: - zi[11] + zi_hi[10] + k[11,10] - 310 m[11,10] <= 0
 HeightCompat_11_10: h[11,10] - k[11,10] <= 0
 NoHeightCompat_11_10: - 155 h[11,10] + k[11,10] <= 0
 Overlap_11_10: - x_p[10,11] - x_p[11,10] + o[11,10] <= 0
 NoOverlap_11_10: x_p[10,11] + x_p[11,10] - 2 o[11,10] <= 0
 3_Supports_11_10: - h[11,10] - o[11,10] - s[11,10] <= -1
 3_DoesntSupport_11_10: h[11,10] + o[11,10] + 2 s[11,10] <= 2
 beta1_11_10: - s[11,10] + beta1[11,10] <= 0
 beta2_11_10: - s[11,10] + beta2[11,10] <= 0
 sBeta1Flag_11_10: - s[11,10] + beta1[11,10] <= 0
 sBeta2Flag_11_10: - s[11,10] + beta2[11,10] <= 0
 Eta1Overlap_11_10: eta1[11,10] + beta1[11,10] <= 1
 Eta2Overlap_11_10: eta2[11,10] + beta2[11,10] <= 1
 Eta1Flag_11_10: xi_l[10] - xi_l[11] - 300 eta1[11,10] <= 0
 Eta2Flag_11_10: - xi_r[10] + xi_r[11] - 300 eta2[11,10] <= 0
 10_Supports_11_InSameBin1_0: - p_ij[10,0] + p_ij[11,0] - s[11,10] <= 0
 10_Supports_11_InSameBin2_0: p_ij[10,0] - p_ij[11,0] - s[11,10] <= 0
 Gamma0ForItems_11_inBins_0_NoCut_10: p_ij[11,0] + gamma[11] <= 1
 10_Supports_11_InSameBin1_1: - p_ij[10,1] + p_ij[11,1] - s[11,10] <= 0
 10_Supports_11_InSameBin2_1: p_ij[10,1] - p_ij[11,1] - s[11,10] <= 0
 Gamma0ForItems_11_inBins_1_NoCut_10: p_ij[11,1] + gamma[11] <= 1
 10_Supports_11_InSameBin1_2: - p_ij[10,2] + p_ij[11,2] - s[11,10] <= 0
 10_Supports_11_InSameBin2_2: p_ij[10,2] - p_ij[11,2] - s[11,10] <= 0
 Gamma0ForItems_11_inBins_2_NoCut_10: p_ij[11,2] + gamma[11] <= 98
 10_Supports_11_InSameBin1_3: - p_ij[10,3] + p_ij[11,3] - s[11,10] <= 0
 10_Supports_11_InSameBin2_3: p_ij[10,3] - p_ij[11,3] - s[11,10] <= 0
 Gamma0ForItems_11_inBins_3_NoCut_10: p_ij[11,3] + gamma[11] <= 98
 11_OnCutIn2_10: 10000 p_ij[11,2] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 11_OnCutIn3_10: 10000 p_ij[11,3] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 11_OnCutIn2_11: 10000 p_ij[11,2] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 11_OnCutIn3_11: 10000 p_ij[11,3] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 AbsZ111_12: - zi[11] + zi_hi[12] - k[11,12] <= 0
 AbsZ211_12: zi[11] - zi_hi[12] - k[11,12] <= 0
 12_Below_11: zi[11] - zi_hi[12] + k[11,12] + 310 m[11,12] <= 310
 12_NotBelow_11: - zi[11] + zi_hi[12] + k[11,12] - 310 m[11,12] <= 0
 HeightCompat_11_12: h[11,12] - k[11,12] <= 0
 NoHeightCompat_11_12: - 155 h[11,12] + k[11,12] <= 0
 Overlap_11_12: - x_p[11,12] - x_p[12,11] + o[11,12] <= 0
 NoOverlap_11_12: x_p[11,12] + x_p[12,11] - 2 o[11,12] <= 0
 3_Supports_11_12: - h[11,12] - o[11,12] - s[11,12] <= -1
 3_DoesntSupport_11_12: h[11,12] + o[11,12] + 2 s[11,12] <= 2
 beta1_11_12: - s[11,12] + beta1[11,12] <= 0
 beta2_11_12: - s[11,12] + beta2[11,12] <= 0
 sBeta1Flag_11_12: - s[11,12] + beta1[11,12] <= 0
 sBeta2Flag_11_12: - s[11,12] + beta2[11,12] <= 0
 Eta1Overlap_11_12: eta1[11,12] + beta1[11,12] <= 1
 Eta2Overlap_11_12: eta2[11,12] + beta2[11,12] <= 1
 Eta1Flag_11_12: - xi_l[11] + xi_l[12] - 300 eta1[11,12] <= 0
 Eta2Flag_11_12: xi_r[11] - xi_r[12] - 300 eta2[11,12] <= 0
 12_Supports_11_InSameBin1_0: p_ij[11,0] - p_ij[12,0] - s[11,12] <= 0
 12_Supports_11_InSameBin2_0: - p_ij[11,0] + p_ij[12,0] - s[11,12] <= 0
 Gamma0ForItems_11_inBins_0_NoCut_12: p_ij[11,0] + gamma[11] <= 1
 12_Supports_11_InSameBin1_1: p_ij[11,1] - p_ij[12,1] - s[11,12] <= 0
 12_Supports_11_InSameBin2_1: - p_ij[11,1] + p_ij[12,1] - s[11,12] <= 0
 Gamma0ForItems_11_inBins_1_NoCut_12: p_ij[11,1] + gamma[11] <= 1
 12_Supports_11_InSameBin1_2: p_ij[11,2] - p_ij[12,2] - s[11,12] <= 0
 12_Supports_11_InSameBin2_2: - p_ij[11,2] + p_ij[12,2] - s[11,12] <= 0
 Gamma0ForItems_11_inBins_2_NoCut_12: p_ij[11,2] + gamma[11] <= 98
 12_Supports_11_InSameBin1_3: p_ij[11,3] - p_ij[12,3] - s[11,12] <= 0
 12_Supports_11_InSameBin2_3: - p_ij[11,3] + p_ij[12,3] - s[11,12] <= 0
 Gamma0ForItems_11_inBins_3_NoCut_12: p_ij[11,3] + gamma[11] <= 98
 11_OnCutIn2_12: 10000 p_ij[11,2] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 11_OnCutIn3_12: 10000 p_ij[11,3] + 1.261904761904762 xi_l[11] + zi[11]
   + 10000 gamma[11] <= 20053
 GroundConstraint_12: zi[12] + 155 g[12] <= 155
 AbsZ112_0: - zi[12] + zi_hi[0] - k[12,0] <= 0
 AbsZ212_0: zi[12] - zi_hi[0] - k[12,0] <= 0
 0_Below_12: zi[12] - zi_hi[0] + k[12,0] + 310 m[12,0] <= 310
 0_NotBelow_12: - zi[12] + zi_hi[0] + k[12,0] - 310 m[12,0] <= 0
 HeightCompat_12_0: h[12,0] - k[12,0] <= 0
 NoHeightCompat_12_0: - 155 h[12,0] + k[12,0] <= 0
 Overlap_12_0: - x_p[0,12] - x_p[12,0] + o[12,0] <= 0
 NoOverlap_12_0: x_p[0,12] + x_p[12,0] - 2 o[12,0] <= 0
 3_Supports_12_0: - h[12,0] - o[12,0] - s[12,0] <= -1
 3_DoesntSupport_12_0: h[12,0] + o[12,0] + 2 s[12,0] <= 2
 beta1_12_0: - s[12,0] + beta1[12,0] <= 0
 beta2_12_0: - s[12,0] + beta2[12,0] <= 0
 sBeta1Flag_12_0: - s[12,0] + beta1[12,0] <= 0
 sBeta2Flag_12_0: - s[12,0] + beta2[12,0] <= 0
 Eta1Overlap_12_0: eta1[12,0] + beta1[12,0] <= 1
 Eta2Overlap_12_0: eta2[12,0] + beta2[12,0] <= 1
 Eta1Flag_12_0: xi_l[0] - xi_l[12] - 300 eta1[12,0] <= 0
 Eta2Flag_12_0: - xi_r[0] + xi_r[12] - 300 eta2[12,0] <= 0
 0_Supports_12_InSameBin1_0: - p_ij[0,0] + p_ij[12,0] - s[12,0] <= 0
 0_Supports_12_InSameBin2_0: p_ij[0,0] - p_ij[12,0] - s[12,0] <= 0
 Gamma0ForItems_12_inBins_0_NoCut_0: p_ij[12,0] + gamma[12] <= 1
 0_Supports_12_InSameBin1_1: - p_ij[0,1] + p_ij[12,1] - s[12,0] <= 0
 0_Supports_12_InSameBin2_1: p_ij[0,1] - p_ij[12,1] - s[12,0] <= 0
 Gamma0ForItems_12_inBins_1_NoCut_0: p_ij[12,1] + gamma[12] <= 1
 0_Supports_12_InSameBin1_2: - p_ij[0,2] + p_ij[12,2] - s[12,0] <= 0
 0_Supports_12_InSameBin2_2: p_ij[0,2] - p_ij[12,2] - s[12,0] <= 0
 Gamma0ForItems_12_inBins_2_NoCut_0: p_ij[12,2] + gamma[12] <= 98
 0_Supports_12_InSameBin1_3: - p_ij[0,3] + p_ij[12,3] - s[12,0] <= 0
 0_Supports_12_InSameBin2_3: p_ij[0,3] - p_ij[12,3] - s[12,0] <= 0
 Gamma0ForItems_12_inBins_3_NoCut_0: p_ij[12,3] + gamma[12] <= 98
 12_OnCutIn2_0: 10000 p_ij[12,2] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 12_OnCutIn3_0: 10000 p_ij[12,3] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 AbsZ112_1: - zi[12] + zi_hi[1] - k[12,1] <= 0
 AbsZ212_1: zi[12] - zi_hi[1] - k[12,1] <= 0
 1_Below_12: zi[12] - zi_hi[1] + k[12,1] + 310 m[12,1] <= 310
 1_NotBelow_12: - zi[12] + zi_hi[1] + k[12,1] - 310 m[12,1] <= 0
 HeightCompat_12_1: h[12,1] - k[12,1] <= 0
 NoHeightCompat_12_1: - 155 h[12,1] + k[12,1] <= 0
 Overlap_12_1: - x_p[1,12] - x_p[12,1] + o[12,1] <= 0
 NoOverlap_12_1: x_p[1,12] + x_p[12,1] - 2 o[12,1] <= 0
 3_Supports_12_1: - h[12,1] - o[12,1] - s[12,1] <= -1
 3_DoesntSupport_12_1: h[12,1] + o[12,1] + 2 s[12,1] <= 2
 beta1_12_1: - s[12,1] + beta1[12,1] <= 0
 beta2_12_1: - s[12,1] + beta2[12,1] <= 0
 sBeta1Flag_12_1: - s[12,1] + beta1[12,1] <= 0
 sBeta2Flag_12_1: - s[12,1] + beta2[12,1] <= 0
 Eta1Overlap_12_1: eta1[12,1] + beta1[12,1] <= 1
 Eta2Overlap_12_1: eta2[12,1] + beta2[12,1] <= 1
 Eta1Flag_12_1: xi_l[1] - xi_l[12] - 300 eta1[12,1] <= 0
 Eta2Flag_12_1: - xi_r[1] + xi_r[12] - 300 eta2[12,1] <= 0
 1_Supports_12_InSameBin1_0: - p_ij[1,0] + p_ij[12,0] - s[12,1] <= 0
 1_Supports_12_InSameBin2_0: p_ij[1,0] - p_ij[12,0] - s[12,1] <= 0
 Gamma0ForItems_12_inBins_0_NoCut_1: p_ij[12,0] + gamma[12] <= 1
 1_Supports_12_InSameBin1_1: - p_ij[1,1] + p_ij[12,1] - s[12,1] <= 0
 1_Supports_12_InSameBin2_1: p_ij[1,1] - p_ij[12,1] - s[12,1] <= 0
 Gamma0ForItems_12_inBins_1_NoCut_1: p_ij[12,1] + gamma[12] <= 1
 1_Supports_12_InSameBin1_2: - p_ij[1,2] + p_ij[12,2] - s[12,1] <= 0
 1_Supports_12_InSameBin2_2: p_ij[1,2] - p_ij[12,2] - s[12,1] <= 0
 Gamma0ForItems_12_inBins_2_NoCut_1: p_ij[12,2] + gamma[12] <= 98
 1_Supports_12_InSameBin1_3: - p_ij[1,3] + p_ij[12,3] - s[12,1] <= 0
 1_Supports_12_InSameBin2_3: p_ij[1,3] - p_ij[12,3] - s[12,1] <= 0
 Gamma0ForItems_12_inBins_3_NoCut_1: p_ij[12,3] + gamma[12] <= 98
 12_OnCutIn2_1: 10000 p_ij[12,2] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 12_OnCutIn3_1: 10000 p_ij[12,3] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 AbsZ112_2: - zi[12] + zi_hi[2] - k[12,2] <= 0
 AbsZ212_2: zi[12] - zi_hi[2] - k[12,2] <= 0
 2_Below_12: zi[12] - zi_hi[2] + k[12,2] + 310 m[12,2] <= 310
 2_NotBelow_12: - zi[12] + zi_hi[2] + k[12,2] - 310 m[12,2] <= 0
 HeightCompat_12_2: h[12,2] - k[12,2] <= 0
 NoHeightCompat_12_2: - 155 h[12,2] + k[12,2] <= 0
 Overlap_12_2: - x_p[2,12] - x_p[12,2] + o[12,2] <= 0
 NoOverlap_12_2: x_p[2,12] + x_p[12,2] - 2 o[12,2] <= 0
 3_Supports_12_2: - h[12,2] - o[12,2] - s[12,2] <= -1
 3_DoesntSupport_12_2: h[12,2] + o[12,2] + 2 s[12,2] <= 2
 beta1_12_2: - s[12,2] + beta1[12,2] <= 0
 beta2_12_2: - s[12,2] + beta2[12,2] <= 0
 sBeta1Flag_12_2: - s[12,2] + beta1[12,2] <= 0
 sBeta2Flag_12_2: - s[12,2] + beta2[12,2] <= 0
 Eta1Overlap_12_2: eta1[12,2] + beta1[12,2] <= 1
 Eta2Overlap_12_2: eta2[12,2] + beta2[12,2] <= 1
 Eta1Flag_12_2: xi_l[2] - xi_l[12] - 300 eta1[12,2] <= 0
 Eta2Flag_12_2: - xi_r[2] + xi_r[12] - 300 eta2[12,2] <= 0
 2_Supports_12_InSameBin1_0: - p_ij[2,0] + p_ij[12,0] - s[12,2] <= 0
 2_Supports_12_InSameBin2_0: p_ij[2,0] - p_ij[12,0] - s[12,2] <= 0
 Gamma0ForItems_12_inBins_0_NoCut_2: p_ij[12,0] + gamma[12] <= 1
 2_Supports_12_InSameBin1_1: - p_ij[2,1] + p_ij[12,1] - s[12,2] <= 0
 2_Supports_12_InSameBin2_1: p_ij[2,1] - p_ij[12,1] - s[12,2] <= 0
 Gamma0ForItems_12_inBins_1_NoCut_2: p_ij[12,1] + gamma[12] <= 1
 2_Supports_12_InSameBin1_2: - p_ij[2,2] + p_ij[12,2] - s[12,2] <= 0
 2_Supports_12_InSameBin2_2: p_ij[2,2] - p_ij[12,2] - s[12,2] <= 0
 Gamma0ForItems_12_inBins_2_NoCut_2: p_ij[12,2] + gamma[12] <= 98
 2_Supports_12_InSameBin1_3: - p_ij[2,3] + p_ij[12,3] - s[12,2] <= 0
 2_Supports_12_InSameBin2_3: p_ij[2,3] - p_ij[12,3] - s[12,2] <= 0
 Gamma0ForItems_12_inBins_3_NoCut_2: p_ij[12,3] + gamma[12] <= 98
 12_OnCutIn2_2: 10000 p_ij[12,2] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 12_OnCutIn3_2: 10000 p_ij[12,3] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 AbsZ112_3: - zi[12] + zi_hi[3] - k[12,3] <= 0
 AbsZ212_3: zi[12] - zi_hi[3] - k[12,3] <= 0
 3_Below_12: zi[12] - zi_hi[3] + k[12,3] + 310 m[12,3] <= 310
 3_NotBelow_12: - zi[12] + zi_hi[3] + k[12,3] - 310 m[12,3] <= 0
 HeightCompat_12_3: h[12,3] - k[12,3] <= 0
 NoHeightCompat_12_3: - 155 h[12,3] + k[12,3] <= 0
 Overlap_12_3: - x_p[3,12] - x_p[12,3] + o[12,3] <= 0
 NoOverlap_12_3: x_p[3,12] + x_p[12,3] - 2 o[12,3] <= 0
 3_Supports_12_3: - h[12,3] - o[12,3] - s[12,3] <= -1
 3_DoesntSupport_12_3: h[12,3] + o[12,3] + 2 s[12,3] <= 2
 beta1_12_3: - s[12,3] + beta1[12,3] <= 0
 beta2_12_3: - s[12,3] + beta2[12,3] <= 0
 sBeta1Flag_12_3: - s[12,3] + beta1[12,3] <= 0
 sBeta2Flag_12_3: - s[12,3] + beta2[12,3] <= 0
 Eta1Overlap_12_3: eta1[12,3] + beta1[12,3] <= 1
 Eta2Overlap_12_3: eta2[12,3] + beta2[12,3] <= 1
 Eta1Flag_12_3: xi_l[3] - xi_l[12] - 300 eta1[12,3] <= 0
 Eta2Flag_12_3: - xi_r[3] + xi_r[12] - 300 eta2[12,3] <= 0
 3_Supports_12_InSameBin1_0: - p_ij[3,0] + p_ij[12,0] - s[12,3] <= 0
 3_Supports_12_InSameBin2_0: p_ij[3,0] - p_ij[12,0] - s[12,3] <= 0
 Gamma0ForItems_12_inBins_0_NoCut_3: p_ij[12,0] + gamma[12] <= 1
 3_Supports_12_InSameBin1_1: - p_ij[3,1] + p_ij[12,1] - s[12,3] <= 0
 3_Supports_12_InSameBin2_1: p_ij[3,1] - p_ij[12,1] - s[12,3] <= 0
 Gamma0ForItems_12_inBins_1_NoCut_3: p_ij[12,1] + gamma[12] <= 1
 3_Supports_12_InSameBin1_2: - p_ij[3,2] + p_ij[12,2] - s[12,3] <= 0
 3_Supports_12_InSameBin2_2: p_ij[3,2] - p_ij[12,2] - s[12,3] <= 0
 Gamma0ForItems_12_inBins_2_NoCut_3: p_ij[12,2] + gamma[12] <= 98
 3_Supports_12_InSameBin1_3: - p_ij[3,3] + p_ij[12,3] - s[12,3] <= 0
 3_Supports_12_InSameBin2_3: p_ij[3,3] - p_ij[12,3] - s[12,3] <= 0
 Gamma0ForItems_12_inBins_3_NoCut_3: p_ij[12,3] + gamma[12] <= 98
 12_OnCutIn2_3: 10000 p_ij[12,2] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 12_OnCutIn3_3: 10000 p_ij[12,3] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 AbsZ112_4: - zi[12] + zi_hi[4] - k[12,4] <= 0
 AbsZ212_4: zi[12] - zi_hi[4] - k[12,4] <= 0
 4_Below_12: zi[12] - zi_hi[4] + k[12,4] + 310 m[12,4] <= 310
 4_NotBelow_12: - zi[12] + zi_hi[4] + k[12,4] - 310 m[12,4] <= 0
 HeightCompat_12_4: h[12,4] - k[12,4] <= 0
 NoHeightCompat_12_4: - 155 h[12,4] + k[12,4] <= 0
 Overlap_12_4: - x_p[4,12] - x_p[12,4] + o[12,4] <= 0
 NoOverlap_12_4: x_p[4,12] + x_p[12,4] - 2 o[12,4] <= 0
 3_Supports_12_4: - h[12,4] - o[12,4] - s[12,4] <= -1
 3_DoesntSupport_12_4: h[12,4] + o[12,4] + 2 s[12,4] <= 2
 beta1_12_4: - s[12,4] + beta1[12,4] <= 0
 beta2_12_4: - s[12,4] + beta2[12,4] <= 0
 sBeta1Flag_12_4: - s[12,4] + beta1[12,4] <= 0
 sBeta2Flag_12_4: - s[12,4] + beta2[12,4] <= 0
 Eta1Overlap_12_4: eta1[12,4] + beta1[12,4] <= 1
 Eta2Overlap_12_4: eta2[12,4] + beta2[12,4] <= 1
 Eta1Flag_12_4: xi_l[4] - xi_l[12] - 300 eta1[12,4] <= 0
 Eta2Flag_12_4: - xi_r[4] + xi_r[12] - 300 eta2[12,4] <= 0
 4_Supports_12_InSameBin1_0: - p_ij[4,0] + p_ij[12,0] - s[12,4] <= 0
 4_Supports_12_InSameBin2_0: p_ij[4,0] - p_ij[12,0] - s[12,4] <= 0
 Gamma0ForItems_12_inBins_0_NoCut_4: p_ij[12,0] + gamma[12] <= 1
 4_Supports_12_InSameBin1_1: - p_ij[4,1] + p_ij[12,1] - s[12,4] <= 0
 4_Supports_12_InSameBin2_1: p_ij[4,1] - p_ij[12,1] - s[12,4] <= 0
 Gamma0ForItems_12_inBins_1_NoCut_4: p_ij[12,1] + gamma[12] <= 1
 4_Supports_12_InSameBin1_2: - p_ij[4,2] + p_ij[12,2] - s[12,4] <= 0
 4_Supports_12_InSameBin2_2: p_ij[4,2] - p_ij[12,2] - s[12,4] <= 0
 Gamma0ForItems_12_inBins_2_NoCut_4: p_ij[12,2] + gamma[12] <= 98
 4_Supports_12_InSameBin1_3: - p_ij[4,3] + p_ij[12,3] - s[12,4] <= 0
 4_Supports_12_InSameBin2_3: p_ij[4,3] - p_ij[12,3] - s[12,4] <= 0
 Gamma0ForItems_12_inBins_3_NoCut_4: p_ij[12,3] + gamma[12] <= 98
 12_OnCutIn2_4: 10000 p_ij[12,2] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 12_OnCutIn3_4: 10000 p_ij[12,3] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 AbsZ112_5: - zi[12] + zi_hi[5] - k[12,5] <= 0
 AbsZ212_5: zi[12] - zi_hi[5] - k[12,5] <= 0
 5_Below_12: zi[12] - zi_hi[5] + k[12,5] + 310 m[12,5] <= 310
 5_NotBelow_12: - zi[12] + zi_hi[5] + k[12,5] - 310 m[12,5] <= 0
 HeightCompat_12_5: h[12,5] - k[12,5] <= 0
 NoHeightCompat_12_5: - 155 h[12,5] + k[12,5] <= 0
 Overlap_12_5: - x_p[5,12] - x_p[12,5] + o[12,5] <= 0
 NoOverlap_12_5: x_p[5,12] + x_p[12,5] - 2 o[12,5] <= 0
 3_Supports_12_5: - h[12,5] - o[12,5] - s[12,5] <= -1
 3_DoesntSupport_12_5: h[12,5] + o[12,5] + 2 s[12,5] <= 2
 beta1_12_5: - s[12,5] + beta1[12,5] <= 0
 beta2_12_5: - s[12,5] + beta2[12,5] <= 0
 sBeta1Flag_12_5: - s[12,5] + beta1[12,5] <= 0
 sBeta2Flag_12_5: - s[12,5] + beta2[12,5] <= 0
 Eta1Overlap_12_5: eta1[12,5] + beta1[12,5] <= 1
 Eta2Overlap_12_5: eta2[12,5] + beta2[12,5] <= 1
 Eta1Flag_12_5: xi_l[5] - xi_l[12] - 300 eta1[12,5] <= 0
 Eta2Flag_12_5: - xi_r[5] + xi_r[12] - 300 eta2[12,5] <= 0
 5_Supports_12_InSameBin1_0: - p_ij[5,0] + p_ij[12,0] - s[12,5] <= 0
 5_Supports_12_InSameBin2_0: p_ij[5,0] - p_ij[12,0] - s[12,5] <= 0
 Gamma0ForItems_12_inBins_0_NoCut_5: p_ij[12,0] + gamma[12] <= 1
 5_Supports_12_InSameBin1_1: - p_ij[5,1] + p_ij[12,1] - s[12,5] <= 0
 5_Supports_12_InSameBin2_1: p_ij[5,1] - p_ij[12,1] - s[12,5] <= 0
 Gamma0ForItems_12_inBins_1_NoCut_5: p_ij[12,1] + gamma[12] <= 1
 5_Supports_12_InSameBin1_2: - p_ij[5,2] + p_ij[12,2] - s[12,5] <= 0
 5_Supports_12_InSameBin2_2: p_ij[5,2] - p_ij[12,2] - s[12,5] <= 0
 Gamma0ForItems_12_inBins_2_NoCut_5: p_ij[12,2] + gamma[12] <= 98
 5_Supports_12_InSameBin1_3: - p_ij[5,3] + p_ij[12,3] - s[12,5] <= 0
 5_Supports_12_InSameBin2_3: p_ij[5,3] - p_ij[12,3] - s[12,5] <= 0
 Gamma0ForItems_12_inBins_3_NoCut_5: p_ij[12,3] + gamma[12] <= 98
 12_OnCutIn2_5: 10000 p_ij[12,2] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 12_OnCutIn3_5: 10000 p_ij[12,3] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 AbsZ112_6: - zi[12] + zi_hi[6] - k[12,6] <= 0
 AbsZ212_6: zi[12] - zi_hi[6] - k[12,6] <= 0
 6_Below_12: zi[12] - zi_hi[6] + k[12,6] + 310 m[12,6] <= 310
 6_NotBelow_12: - zi[12] + zi_hi[6] + k[12,6] - 310 m[12,6] <= 0
 HeightCompat_12_6: h[12,6] - k[12,6] <= 0
 NoHeightCompat_12_6: - 155 h[12,6] + k[12,6] <= 0
 Overlap_12_6: - x_p[6,12] - x_p[12,6] + o[12,6] <= 0
 NoOverlap_12_6: x_p[6,12] + x_p[12,6] - 2 o[12,6] <= 0
 3_Supports_12_6: - h[12,6] - o[12,6] - s[12,6] <= -1
 3_DoesntSupport_12_6: h[12,6] + o[12,6] + 2 s[12,6] <= 2
 beta1_12_6: - s[12,6] + beta1[12,6] <= 0
 beta2_12_6: - s[12,6] + beta2[12,6] <= 0
 sBeta1Flag_12_6: - s[12,6] + beta1[12,6] <= 0
 sBeta2Flag_12_6: - s[12,6] + beta2[12,6] <= 0
 Eta1Overlap_12_6: eta1[12,6] + beta1[12,6] <= 1
 Eta2Overlap_12_6: eta2[12,6] + beta2[12,6] <= 1
 Eta1Flag_12_6: xi_l[6] - xi_l[12] - 300 eta1[12,6] <= 0
 Eta2Flag_12_6: - xi_r[6] + xi_r[12] - 300 eta2[12,6] <= 0
 6_Supports_12_InSameBin1_0: - p_ij[6,0] + p_ij[12,0] - s[12,6] <= 0
 6_Supports_12_InSameBin2_0: p_ij[6,0] - p_ij[12,0] - s[12,6] <= 0
 Gamma0ForItems_12_inBins_0_NoCut_6: p_ij[12,0] + gamma[12] <= 1
 6_Supports_12_InSameBin1_1: - p_ij[6,1] + p_ij[12,1] - s[12,6] <= 0
 6_Supports_12_InSameBin2_1: p_ij[6,1] - p_ij[12,1] - s[12,6] <= 0
 Gamma0ForItems_12_inBins_1_NoCut_6: p_ij[12,1] + gamma[12] <= 1
 6_Supports_12_InSameBin1_2: - p_ij[6,2] + p_ij[12,2] - s[12,6] <= 0
 6_Supports_12_InSameBin2_2: p_ij[6,2] - p_ij[12,2] - s[12,6] <= 0
 Gamma0ForItems_12_inBins_2_NoCut_6: p_ij[12,2] + gamma[12] <= 98
 6_Supports_12_InSameBin1_3: - p_ij[6,3] + p_ij[12,3] - s[12,6] <= 0
 6_Supports_12_InSameBin2_3: p_ij[6,3] - p_ij[12,3] - s[12,6] <= 0
 Gamma0ForItems_12_inBins_3_NoCut_6: p_ij[12,3] + gamma[12] <= 98
 12_OnCutIn2_6: 10000 p_ij[12,2] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 12_OnCutIn3_6: 10000 p_ij[12,3] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 AbsZ112_7: - zi[12] + zi_hi[7] - k[12,7] <= 0
 AbsZ212_7: zi[12] - zi_hi[7] - k[12,7] <= 0
 7_Below_12: zi[12] - zi_hi[7] + k[12,7] + 310 m[12,7] <= 310
 7_NotBelow_12: - zi[12] + zi_hi[7] + k[12,7] - 310 m[12,7] <= 0
 HeightCompat_12_7: h[12,7] - k[12,7] <= 0
 NoHeightCompat_12_7: - 155 h[12,7] + k[12,7] <= 0
 Overlap_12_7: - x_p[7,12] - x_p[12,7] + o[12,7] <= 0
 NoOverlap_12_7: x_p[7,12] + x_p[12,7] - 2 o[12,7] <= 0
 3_Supports_12_7: - h[12,7] - o[12,7] - s[12,7] <= -1
 3_DoesntSupport_12_7: h[12,7] + o[12,7] + 2 s[12,7] <= 2
 beta1_12_7: - s[12,7] + beta1[12,7] <= 0
 beta2_12_7: - s[12,7] + beta2[12,7] <= 0
 sBeta1Flag_12_7: - s[12,7] + beta1[12,7] <= 0
 sBeta2Flag_12_7: - s[12,7] + beta2[12,7] <= 0
 Eta1Overlap_12_7: eta1[12,7] + beta1[12,7] <= 1
 Eta2Overlap_12_7: eta2[12,7] + beta2[12,7] <= 1
 Eta1Flag_12_7: xi_l[7] - xi_l[12] - 300 eta1[12,7] <= 0
 Eta2Flag_12_7: - xi_r[7] + xi_r[12] - 300 eta2[12,7] <= 0
 7_Supports_12_InSameBin1_0: - p_ij[7,0] + p_ij[12,0] - s[12,7] <= 0
 7_Supports_12_InSameBin2_0: p_ij[7,0] - p_ij[12,0] - s[12,7] <= 0
 Gamma0ForItems_12_inBins_0_NoCut_7: p_ij[12,0] + gamma[12] <= 1
 7_Supports_12_InSameBin1_1: - p_ij[7,1] + p_ij[12,1] - s[12,7] <= 0
 7_Supports_12_InSameBin2_1: p_ij[7,1] - p_ij[12,1] - s[12,7] <= 0
 Gamma0ForItems_12_inBins_1_NoCut_7: p_ij[12,1] + gamma[12] <= 1
 7_Supports_12_InSameBin1_2: - p_ij[7,2] + p_ij[12,2] - s[12,7] <= 0
 7_Supports_12_InSameBin2_2: p_ij[7,2] - p_ij[12,2] - s[12,7] <= 0
 Gamma0ForItems_12_inBins_2_NoCut_7: p_ij[12,2] + gamma[12] <= 98
 7_Supports_12_InSameBin1_3: - p_ij[7,3] + p_ij[12,3] - s[12,7] <= 0
 7_Supports_12_InSameBin2_3: p_ij[7,3] - p_ij[12,3] - s[12,7] <= 0
 Gamma0ForItems_12_inBins_3_NoCut_7: p_ij[12,3] + gamma[12] <= 98
 12_OnCutIn2_7: 10000 p_ij[12,2] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 12_OnCutIn3_7: 10000 p_ij[12,3] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 AbsZ112_8: - zi[12] + zi_hi[8] - k[12,8] <= 0
 AbsZ212_8: zi[12] - zi_hi[8] - k[12,8] <= 0
 8_Below_12: zi[12] - zi_hi[8] + k[12,8] + 310 m[12,8] <= 310
 8_NotBelow_12: - zi[12] + zi_hi[8] + k[12,8] - 310 m[12,8] <= 0
 HeightCompat_12_8: h[12,8] - k[12,8] <= 0
 NoHeightCompat_12_8: - 155 h[12,8] + k[12,8] <= 0
 Overlap_12_8: - x_p[8,12] - x_p[12,8] + o[12,8] <= 0
 NoOverlap_12_8: x_p[8,12] + x_p[12,8] - 2 o[12,8] <= 0
 3_Supports_12_8: - h[12,8] - o[12,8] - s[12,8] <= -1
 3_DoesntSupport_12_8: h[12,8] + o[12,8] + 2 s[12,8] <= 2
 beta1_12_8: - s[12,8] + beta1[12,8] <= 0
 beta2_12_8: - s[12,8] + beta2[12,8] <= 0
 sBeta1Flag_12_8: - s[12,8] + beta1[12,8] <= 0
 sBeta2Flag_12_8: - s[12,8] + beta2[12,8] <= 0
 Eta1Overlap_12_8: eta1[12,8] + beta1[12,8] <= 1
 Eta2Overlap_12_8: eta2[12,8] + beta2[12,8] <= 1
 Eta1Flag_12_8: xi_l[8] - xi_l[12] - 300 eta1[12,8] <= 0
 Eta2Flag_12_8: - xi_r[8] + xi_r[12] - 300 eta2[12,8] <= 0
 8_Supports_12_InSameBin1_0: - p_ij[8,0] + p_ij[12,0] - s[12,8] <= 0
 8_Supports_12_InSameBin2_0: p_ij[8,0] - p_ij[12,0] - s[12,8] <= 0
 Gamma0ForItems_12_inBins_0_NoCut_8: p_ij[12,0] + gamma[12] <= 1
 8_Supports_12_InSameBin1_1: - p_ij[8,1] + p_ij[12,1] - s[12,8] <= 0
 8_Supports_12_InSameBin2_1: p_ij[8,1] - p_ij[12,1] - s[12,8] <= 0
 Gamma0ForItems_12_inBins_1_NoCut_8: p_ij[12,1] + gamma[12] <= 1
 8_Supports_12_InSameBin1_2: - p_ij[8,2] + p_ij[12,2] - s[12,8] <= 0
 8_Supports_12_InSameBin2_2: p_ij[8,2] - p_ij[12,2] - s[12,8] <= 0
 Gamma0ForItems_12_inBins_2_NoCut_8: p_ij[12,2] + gamma[12] <= 98
 8_Supports_12_InSameBin1_3: - p_ij[8,3] + p_ij[12,3] - s[12,8] <= 0
 8_Supports_12_InSameBin2_3: p_ij[8,3] - p_ij[12,3] - s[12,8] <= 0
 Gamma0ForItems_12_inBins_3_NoCut_8: p_ij[12,3] + gamma[12] <= 98
 12_OnCutIn2_8: 10000 p_ij[12,2] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 12_OnCutIn3_8: 10000 p_ij[12,3] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 AbsZ112_9: - zi[12] + zi_hi[9] - k[12,9] <= 0
 AbsZ212_9: zi[12] - zi_hi[9] - k[12,9] <= 0
 9_Below_12: zi[12] - zi_hi[9] + k[12,9] + 310 m[12,9] <= 310
 9_NotBelow_12: - zi[12] + zi_hi[9] + k[12,9] - 310 m[12,9] <= 0
 HeightCompat_12_9: h[12,9] - k[12,9] <= 0
 NoHeightCompat_12_9: - 155 h[12,9] + k[12,9] <= 0
 Overlap_12_9: - x_p[9,12] - x_p[12,9] + o[12,9] <= 0
 NoOverlap_12_9: x_p[9,12] + x_p[12,9] - 2 o[12,9] <= 0
 3_Supports_12_9: - h[12,9] - o[12,9] - s[12,9] <= -1
 3_DoesntSupport_12_9: h[12,9] + o[12,9] + 2 s[12,9] <= 2
 beta1_12_9: - s[12,9] + beta1[12,9] <= 0
 beta2_12_9: - s[12,9] + beta2[12,9] <= 0
 sBeta1Flag_12_9: - s[12,9] + beta1[12,9] <= 0
 sBeta2Flag_12_9: - s[12,9] + beta2[12,9] <= 0
 Eta1Overlap_12_9: eta1[12,9] + beta1[12,9] <= 1
 Eta2Overlap_12_9: eta2[12,9] + beta2[12,9] <= 1
 Eta1Flag_12_9: xi_l[9] - xi_l[12] - 300 eta1[12,9] <= 0
 Eta2Flag_12_9: - xi_r[9] + xi_r[12] - 300 eta2[12,9] <= 0
 9_Supports_12_InSameBin1_0: - p_ij[9,0] + p_ij[12,0] - s[12,9] <= 0
 9_Supports_12_InSameBin2_0: p_ij[9,0] - p_ij[12,0] - s[12,9] <= 0
 Gamma0ForItems_12_inBins_0_NoCut_9: p_ij[12,0] + gamma[12] <= 1
 9_Supports_12_InSameBin1_1: - p_ij[9,1] + p_ij[12,1] - s[12,9] <= 0
 9_Supports_12_InSameBin2_1: p_ij[9,1] - p_ij[12,1] - s[12,9] <= 0
 Gamma0ForItems_12_inBins_1_NoCut_9: p_ij[12,1] + gamma[12] <= 1
 9_Supports_12_InSameBin1_2: - p_ij[9,2] + p_ij[12,2] - s[12,9] <= 0
 9_Supports_12_InSameBin2_2: p_ij[9,2] - p_ij[12,2] - s[12,9] <= 0
 Gamma0ForItems_12_inBins_2_NoCut_9: p_ij[12,2] + gamma[12] <= 98
 9_Supports_12_InSameBin1_3: - p_ij[9,3] + p_ij[12,3] - s[12,9] <= 0
 9_Supports_12_InSameBin2_3: p_ij[9,3] - p_ij[12,3] - s[12,9] <= 0
 Gamma0ForItems_12_inBins_3_NoCut_9: p_ij[12,3] + gamma[12] <= 98
 12_OnCutIn2_9: 10000 p_ij[12,2] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 12_OnCutIn3_9: 10000 p_ij[12,3] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 AbsZ112_10: - zi[12] + zi_hi[10] - k[12,10] <= 0
 AbsZ212_10: zi[12] - zi_hi[10] - k[12,10] <= 0
 10_Below_12: zi[12] - zi_hi[10] + k[12,10] + 310 m[12,10] <= 310
 10_NotBelow_12: - zi[12] + zi_hi[10] + k[12,10] - 310 m[12,10] <= 0
 HeightCompat_12_10: h[12,10] - k[12,10] <= 0
 NoHeightCompat_12_10: - 155 h[12,10] + k[12,10] <= 0
 Overlap_12_10: - x_p[10,12] - x_p[12,10] + o[12,10] <= 0
 NoOverlap_12_10: x_p[10,12] + x_p[12,10] - 2 o[12,10] <= 0
 3_Supports_12_10: - h[12,10] - o[12,10] - s[12,10] <= -1
 3_DoesntSupport_12_10: h[12,10] + o[12,10] + 2 s[12,10] <= 2
 beta1_12_10: - s[12,10] + beta1[12,10] <= 0
 beta2_12_10: - s[12,10] + beta2[12,10] <= 0
 sBeta1Flag_12_10: - s[12,10] + beta1[12,10] <= 0
 sBeta2Flag_12_10: - s[12,10] + beta2[12,10] <= 0
 Eta1Overlap_12_10: eta1[12,10] + beta1[12,10] <= 1
 Eta2Overlap_12_10: eta2[12,10] + beta2[12,10] <= 1
 Eta1Flag_12_10: xi_l[10] - xi_l[12] - 300 eta1[12,10] <= 0
 Eta2Flag_12_10: - xi_r[10] + xi_r[12] - 300 eta2[12,10] <= 0
 10_Supports_12_InSameBin1_0: - p_ij[10,0] + p_ij[12,0] - s[12,10] <= 0
 10_Supports_12_InSameBin2_0: p_ij[10,0] - p_ij[12,0] - s[12,10] <= 0
 Gamma0ForItems_12_inBins_0_NoCut_10: p_ij[12,0] + gamma[12] <= 1
 10_Supports_12_InSameBin1_1: - p_ij[10,1] + p_ij[12,1] - s[12,10] <= 0
 10_Supports_12_InSameBin2_1: p_ij[10,1] - p_ij[12,1] - s[12,10] <= 0
 Gamma0ForItems_12_inBins_1_NoCut_10: p_ij[12,1] + gamma[12] <= 1
 10_Supports_12_InSameBin1_2: - p_ij[10,2] + p_ij[12,2] - s[12,10] <= 0
 10_Supports_12_InSameBin2_2: p_ij[10,2] - p_ij[12,2] - s[12,10] <= 0
 Gamma0ForItems_12_inBins_2_NoCut_10: p_ij[12,2] + gamma[12] <= 98
 10_Supports_12_InSameBin1_3: - p_ij[10,3] + p_ij[12,3] - s[12,10] <= 0
 10_Supports_12_InSameBin2_3: p_ij[10,3] - p_ij[12,3] - s[12,10] <= 0
 Gamma0ForItems_12_inBins_3_NoCut_10: p_ij[12,3] + gamma[12] <= 98
 12_OnCutIn2_10: 10000 p_ij[12,2] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 12_OnCutIn3_10: 10000 p_ij[12,3] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 AbsZ112_11: - zi[12] + zi_hi[11] - k[12,11] <= 0
 AbsZ212_11: zi[12] - zi_hi[11] - k[12,11] <= 0
 11_Below_12: zi[12] - zi_hi[11] + k[12,11] + 310 m[12,11] <= 310
 11_NotBelow_12: - zi[12] + zi_hi[11] + k[12,11] - 310 m[12,11] <= 0
 HeightCompat_12_11: h[12,11] - k[12,11] <= 0
 NoHeightCompat_12_11: - 155 h[12,11] + k[12,11] <= 0
 Overlap_12_11: - x_p[11,12] - x_p[12,11] + o[12,11] <= 0
 NoOverlap_12_11: x_p[11,12] + x_p[12,11] - 2 o[12,11] <= 0
 3_Supports_12_11: - h[12,11] - o[12,11] - s[12,11] <= -1
 3_DoesntSupport_12_11: h[12,11] + o[12,11] + 2 s[12,11] <= 2
 beta1_12_11: - s[12,11] + beta1[12,11] <= 0
 beta2_12_11: - s[12,11] + beta2[12,11] <= 0
 sBeta1Flag_12_11: - s[12,11] + beta1[12,11] <= 0
 sBeta2Flag_12_11: - s[12,11] + beta2[12,11] <= 0
 Eta1Overlap_12_11: eta1[12,11] + beta1[12,11] <= 1
 Eta2Overlap_12_11: eta2[12,11] + beta2[12,11] <= 1
 Eta1Flag_12_11: xi_l[11] - xi_l[12] - 300 eta1[12,11] <= 0
 Eta2Flag_12_11: - xi_r[11] + xi_r[12] - 300 eta2[12,11] <= 0
 11_Supports_12_InSameBin1_0: - p_ij[11,0] + p_ij[12,0] - s[12,11] <= 0
 11_Supports_12_InSameBin2_0: p_ij[11,0] - p_ij[12,0] - s[12,11] <= 0
 Gamma0ForItems_12_inBins_0_NoCut_11: p_ij[12,0] + gamma[12] <= 1
 11_Supports_12_InSameBin1_1: - p_ij[11,1] + p_ij[12,1] - s[12,11] <= 0
 11_Supports_12_InSameBin2_1: p_ij[11,1] - p_ij[12,1] - s[12,11] <= 0
 Gamma0ForItems_12_inBins_1_NoCut_11: p_ij[12,1] + gamma[12] <= 1
 11_Supports_12_InSameBin1_2: - p_ij[11,2] + p_ij[12,2] - s[12,11] <= 0
 11_Supports_12_InSameBin2_2: p_ij[11,2] - p_ij[12,2] - s[12,11] <= 0
 Gamma0ForItems_12_inBins_2_NoCut_11: p_ij[12,2] + gamma[12] <= 98
 11_Supports_12_InSameBin1_3: - p_ij[11,3] + p_ij[12,3] - s[12,11] <= 0
 11_Supports_12_InSameBin2_3: p_ij[11,3] - p_ij[12,3] - s[12,11] <= 0
 Gamma0ForItems_12_inBins_3_NoCut_11: p_ij[12,3] + gamma[12] <= 98
 12_OnCutIn2_11: 10000 p_ij[12,2] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 12_OnCutIn3_11: 10000 p_ij[12,3] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 12_OnCutIn2_12: 10000 p_ij[12,2] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 12_OnCutIn3_12: 10000 p_ij[12,3] + 1.261904761904762 xi_l[12] + zi[12]
   + 10000 gamma[12] <= 20053
 FlaggingConstraint_0_0: - p_ij[0,0] + u_j[0] >= 0
 FlaggingConstraint_0_1: - p_ij[0,1] + u_j[1] >= 0
 FlaggingConstraint_0_2: - p_ij[0,2] + u_j[2] >= 0
 FlaggingConstraint_0_3: - p_ij[0,3] + u_j[3] >= 0
 FlaggingConstraint_1_0: - p_ij[1,0] + u_j[0] >= 0
 FlaggingConstraint_1_1: - p_ij[1,1] + u_j[1] >= 0
 FlaggingConstraint_1_2: - p_ij[1,2] + u_j[2] >= 0
 FlaggingConstraint_1_3: - p_ij[1,3] + u_j[3] >= 0
 FlaggingConstraint_2_0: - p_ij[2,0] + u_j[0] >= 0
 FlaggingConstraint_2_1: - p_ij[2,1] + u_j[1] >= 0
 FlaggingConstraint_2_2: - p_ij[2,2] + u_j[2] >= 0
 FlaggingConstraint_2_3: - p_ij[2,3] + u_j[3] >= 0
 FlaggingConstraint_3_0: - p_ij[3,0] + u_j[0] >= 0
 FlaggingConstraint_3_1: - p_ij[3,1] + u_j[1] >= 0
 FlaggingConstraint_3_2: - p_ij[3,2] + u_j[2] >= 0
 FlaggingConstraint_3_3: - p_ij[3,3] + u_j[3] >= 0
 FlaggingConstraint_4_0: - p_ij[4,0] + u_j[0] >= 0
 FlaggingConstraint_4_1: - p_ij[4,1] + u_j[1] >= 0
 FlaggingConstraint_4_2: - p_ij[4,2] + u_j[2] >= 0
 FlaggingConstraint_4_3: - p_ij[4,3] + u_j[3] >= 0
 FlaggingConstraint_5_0: - p_ij[5,0] + u_j[0] >= 0
 FlaggingConstraint_5_1: - p_ij[5,1] + u_j[1] >= 0
 FlaggingConstraint_5_2: - p_ij[5,2] + u_j[2] >= 0
 FlaggingConstraint_5_3: - p_ij[5,3] + u_j[3] >= 0
 FlaggingConstraint_6_0: - p_ij[6,0] + u_j[0] >= 0
 FlaggingConstraint_6_1: - p_ij[6,1] + u_j[1] >= 0
 FlaggingConstraint_6_2: - p_ij[6,2] + u_j[2] >= 0
 FlaggingConstraint_6_3: - p_ij[6,3] + u_j[3] >= 0
 FlaggingConstraint_7_0: - p_ij[7,0] + u_j[0] >= 0
 FlaggingConstraint_7_1: - p_ij[7,1] + u_j[1] >= 0
 FlaggingConstraint_7_2: - p_ij[7,2] + u_j[2] >= 0
 FlaggingConstraint_7_3: - p_ij[7,3] + u_j[3] >= 0
 FlaggingConstraint_8_0: - p_ij[8,0] + u_j[0] >= 0
 FlaggingConstraint_8_1: - p_ij[8,1] + u_j[1] >= 0
 FlaggingConstraint_8_2: - p_ij[8,2] + u_j[2] >= 0
 FlaggingConstraint_8_3: - p_ij[8,3] + u_j[3] >= 0
 FlaggingConstraint_9_0: - p_ij[9,0] + u_j[0] >= 0
 FlaggingConstraint_9_1: - p_ij[9,1] + u_j[1] >= 0
 FlaggingConstraint_9_2: - p_ij[9,2] + u_j[2] >= 0
 FlaggingConstraint_9_3: - p_ij[9,3] + u_j[3] >= 0
 FlaggingConstraint_10_0: - p_ij[10,0] + u_j[0] >= 0
 FlaggingConstraint_10_1: - p_ij[10,1] + u_j[1] >= 0
 FlaggingConstraint_10_2: - p_ij[10,2] + u_j[2] >= 0
 FlaggingConstraint_10_3: - p_ij[10,3] + u_j[3] >= 0
 FlaggingConstraint_11_0: - p_ij[11,0] + u_j[0] >= 0
 FlaggingConstraint_11_1: - p_ij[11,1] + u_j[1] >= 0
 FlaggingConstraint_11_2: - p_ij[11,2] + u_j[2] >= 0
 FlaggingConstraint_11_3: - p_ij[11,3] + u_j[3] >= 0
 FlaggingConstraint_12_0: - p_ij[12,0] + u_j[0] >= 0
 FlaggingConstraint_12_1: - p_ij[12,1] + u_j[1] >= 0
 FlaggingConstraint_12_2: - p_ij[12,2] + u_j[2] >= 0
 FlaggingConstraint_12_3: - p_ij[12,3] + u_j[3] >= 0
 Fragile_0_0: s[0,0] <= 13
 Fragile_0_1: s[0,1] <= 13
 Fragile_0_2: s[0,2] <= 13
 Fragile_0_3: s[0,3] <= 0
 Fragile_0_4: s[0,4] <= 13
 Fragile_0_5: s[0,5] <= 13
 Fragile_0_6: s[0,6] <= 13
 Fragile_0_7: s[0,7] <= 13
 Fragile_0_8: s[0,8] <= 13
 Fragile_0_9: s[0,9] <= 0
 Fragile_0_10: s[0,10] <= 13
 Fragile_0_11: s[0,11] <= 13
 Fragile_0_12: s[0,12] <= 0
 Fragile_1_0: s[1,0] <= 13
 Fragile_1_1: s[1,1] <= 13
 Fragile_1_2: s[1,2] <= 13
 Fragile_1_3: s[1,3] <= 0
 Fragile_1_4: s[1,4] <= 13
 Fragile_1_5: s[1,5] <= 13
 Fragile_1_6: s[1,6] <= 13
 Fragile_1_7: s[1,7] <= 13
 Fragile_1_8: s[1,8] <= 13
 Fragile_1_9: s[1,9] <= 0
 Fragile_1_10: s[1,10] <= 13
 Fragile_1_11: s[1,11] <= 13
 Fragile_1_12: s[1,12] <= 0
 Fragile_2_0: s[2,0] <= 13
 Fragile_2_1: s[2,1] <= 13
 Fragile_2_2: s[2,2] <= 13
 Fragile_2_3: s[2,3] <= 0
 Fragile_2_4: s[2,4] <= 13
 Fragile_2_5: s[2,5] <= 13
 Fragile_2_6: s[2,6] <= 13
 Fragile_2_7: s[2,7] <= 13
 Fragile_2_8: s[2,8] <= 13
 Fragile_2_9: s[2,9] <= 0
 Fragile_2_10: s[2,10] <= 13
 Fragile_2_11: s[2,11] <= 13
 Fragile_2_12: s[2,12] <= 0
 Fragile_3_0: s[3,0] <= 13
 Fragile_3_1: s[3,1] <= 13
 Fragile_3_2: s[3,2] <= 13
 Fragile_3_3: s[3,3] <= 0
 Fragile_3_4: s[3,4] <= 13
 Fragile_3_5: s[3,5] <= 13
 Fragile_3_6: s[3,6] <= 13
 Fragile_3_7: s[3,7] <= 13
 Fragile_3_8: s[3,8] <= 13
 Fragile_3_9: s[3,9] <= 0
 Fragile_3_10: s[3,10] <= 13
 Fragile_3_11: s[3,11] <= 13
 Fragile_3_12: s[3,12] <= 0
 Fragile_4_0: s[4,0] <= 13
 Fragile_4_1: s[4,1] <= 13
 Fragile_4_2: s[4,2] <= 13
 Fragile_4_3: s[4,3] <= 0
 Fragile_4_4: s[4,4] <= 13
 Fragile_4_5: s[4,5] <= 13
 Fragile_4_6: s[4,6] <= 13
 Fragile_4_7: s[4,7] <= 13
 Fragile_4_8: s[4,8] <= 13
 Fragile_4_9: s[4,9] <= 0
 Fragile_4_10: s[4,10] <= 13
 Fragile_4_11: s[4,11] <= 13
 Fragile_4_12: s[4,12] <= 0
 Fragile_5_0: s[5,0] <= 13
 Fragile_5_1: s[5,1] <= 13
 Fragile_5_2: s[5,2] <= 13
 Fragile_5_3: s[5,3] <= 0
 Fragile_5_4: s[5,4] <= 13
 Fragile_5_5: s[5,5] <= 13
 Fragile_5_6: s[5,6] <= 13
 Fragile_5_7: s[5,7] <= 13
 Fragile_5_8: s[5,8] <= 13
 Fragile_5_9: s[5,9] <= 0
 Fragile_5_10: s[5,10] <= 13
 Fragile_5_11: s[5,11] <= 13
 Fragile_5_12: s[5,12] <= 0
 Fragile_6_0: s[6,0] <= 13
 Fragile_6_1: s[6,1] <= 13
 Fragile_6_2: s[6,2] <= 13
 Fragile_6_3: s[6,3] <= 0
 Fragile_6_4: s[6,4] <= 13
 Fragile_6_5: s[6,5] <= 13
 Fragile_6_6: s[6,6] <= 13
 Fragile_6_7: s[6,7] <= 13
 Fragile_6_8: s[6,8] <= 13
 Fragile_6_9: s[6,9] <= 0
 Fragile_6_10: s[6,10] <= 13
 Fragile_6_11: s[6,11] <= 13
 Fragile_6_12: s[6,12] <= 0
 Fragile_7_0: s[7,0] <= 13
 Fragile_7_1: s[7,1] <= 13
 Fragile_7_2: s[7,2] <= 13
 Fragile_7_3: s[7,3] <= 0
 Fragile_7_4: s[7,4] <= 13
 Fragile_7_5: s[7,5] <= 13
 Fragile_7_6: s[7,6] <= 13
 Fragile_7_7: s[7,7] <= 13
 Fragile_7_8: s[7,8] <= 13
 Fragile_7_9: s[7,9] <= 0
 Fragile_7_10: s[7,10] <= 13
 Fragile_7_11: s[7,11] <= 13
 Fragile_7_12: s[7,12] <= 0
 Fragile_8_0: s[8,0] <= 13
 Fragile_8_1: s[8,1] <= 13
 Fragile_8_2: s[8,2] <= 13
 Fragile_8_3: s[8,3] <= 0
 Fragile_8_4: s[8,4] <= 13
 Fragile_8_5: s[8,5] <= 13
 Fragile_8_6: s[8,6] <= 13
 Fragile_8_7: s[8,7] <= 13
 Fragile_8_8: s[8,8] <= 13
 Fragile_8_9: s[8,9] <= 0
 Fragile_8_10: s[8,10] <= 13
 Fragile_8_11: s[8,11] <= 13
 Fragile_8_12: s[8,12] <= 0
 Fragile_9_0: s[9,0] <= 13
 Fragile_9_1: s[9,1] <= 13
 Fragile_9_2: s[9,2] <= 13
 Fragile_9_3: s[9,3] <= 0
 Fragile_9_4: s[9,4] <= 13
 Fragile_9_5: s[9,5] <= 13
 Fragile_9_6: s[9,6] <= 13
 Fragile_9_7: s[9,7] <= 13
 Fragile_9_8: s[9,8] <= 13
 Fragile_9_9: s[9,9] <= 0
 Fragile_9_10: s[9,10] <= 13
 Fragile_9_11: s[9,11] <= 13
 Fragile_9_12: s[9,12] <= 0
 Fragile_10_0: s[10,0] <= 13
 Fragile_10_1: s[10,1] <= 13
 Fragile_10_2: s[10,2] <= 13
 Fragile_10_3: s[10,3] <= 0
 Fragile_10_4: s[10,4] <= 13
 Fragile_10_5: s[10,5] <= 13
 Fragile_10_6: s[10,6] <= 13
 Fragile_10_7: s[10,7] <= 13
 Fragile_10_8: s[10,8] <= 13
 Fragile_10_9: s[10,9] <= 0
 Fragile_10_10: s[10,10] <= 13
 Fragile_10_11: s[10,11] <= 13
 Fragile_10_12: s[10,12] <= 0
 Fragile_11_0: s[11,0] <= 13
 Fragile_11_1: s[11,1] <= 13
 Fragile_11_2: s[11,2] <= 13
 Fragile_11_3: s[11,3] <= 0
 Fragile_11_4: s[11,4] <= 13
 Fragile_11_5: s[11,5] <= 13
 Fragile_11_6: s[11,6] <= 13
 Fragile_11_7: s[11,7] <= 13
 Fragile_11_8: s[11,8] <= 13
 Fragile_11_9: s[11,9] <= 0
 Fragile_11_10: s[11,10] <= 13
 Fragile_11_11: s[11,11] <= 13
 Fragile_11_12: s[11,12] <= 0
 Fragile_12_0: s[12,0] <= 13
 Fragile_12_1: s[12,1] <= 13
 Fragile_12_2: s[12,2] <= 13
 Fragile_12_3: s[12,3] <= 0
 Fragile_12_4: s[12,4] <= 13
 Fragile_12_5: s[12,5] <= 13
 Fragile_12_6: s[12,6] <= 13
 Fragile_12_7: s[12,7] <= 13
 Fragile_12_8: s[12,8] <= 13
 Fragile_12_9: s[12,9] <= 0
 Fragile_12_10: s[12,10] <= 13
 Fragile_12_11: s[12,11] <= 13
 Fragile_12_12: s[12,12] <= 0
 Perishable_radioactive_0: p_ij[10,0] + p_ij[12,0] <= 1
 Perishable_radioactive_1: p_ij[10,1] + p_ij[12,1] <= 1
 Perishable_radioactive_2: p_ij[10,2] + p_ij[12,2] <= 1
 Perishable_radioactive_3: p_ij[10,3] + p_ij[12,3] <= 1
Bounds
Binaries
 p_ij[0,0] p_ij[0,1] p_ij[0,2] p_ij[0,3] p_ij[1,0] p_ij[1,1] p_ij[1,2]
 p_ij[1,3] p_ij[2,0] p_ij[2,1] p_ij[2,2] p_ij[2,3] p_ij[3,0] p_ij[3,1]
 p_ij[3,2] p_ij[3,3] p_ij[4,0] p_ij[4,1] p_ij[4,2] p_ij[4,3] p_ij[5,0]
 p_ij[5,1] p_ij[5,2] p_ij[5,3] p_ij[6,0] p_ij[6,1] p_ij[6,2] p_ij[6,3]
 p_ij[7,0] p_ij[7,1] p_ij[7,2] p_ij[7,3] p_ij[8,0] p_ij[8,1] p_ij[8,2]
 p_ij[8,3] p_ij[9,0] p_ij[9,1] p_ij[9,2] p_ij[9,3] p_ij[10,0] p_ij[10,1]
 p_ij[10,2] p_ij[10,3] p_ij[11,0] p_ij[11,1] p_ij[11,2] p_ij[11,3]
 p_ij[12,0] p_ij[12,1] p_ij[12,2] p_ij[12,3] u_j[0] u_j[1] u_j[2] u_j[3]
 x_p[0,0] x_p[0,1] x_p[0,2] x_p[0,3] x_p[0,4] x_p[0,5] x_p[0,6] x_p[0,7]
 x_p[0,8] x_p[0,9] x_p[0,10] x_p[0,11] x_p[0,12] x_p[1,0] x_p[1,1] x_p[1,2]
 x_p[1,3] x_p[1,4] x_p[1,5] x_p[1,6] x_p[1,7] x_p[1,8] x_p[1,9] x_p[1,10]
 x_p[1,11] x_p[1,12] x_p[2,0] x_p[2,1] x_p[2,2] x_p[2,3] x_p[2,4] x_p[2,5]
 x_p[2,6] x_p[2,7] x_p[2,8] x_p[2,9] x_p[2,10] x_p[2,11] x_p[2,12] x_p[3,0]
 x_p[3,1] x_p[3,2] x_p[3,3] x_p[3,4] x_p[3,5] x_p[3,6] x_p[3,7] x_p[3,8]
 x_p[3,9] x_p[3,10] x_p[3,11] x_p[3,12] x_p[4,0] x_p[4,1] x_p[4,2] x_p[4,3]
 x_p[4,4] x_p[4,5] x_p[4,6] x_p[4,7] x_p[4,8] x_p[4,9] x_p[4,10] x_p[4,11]
 x_p[4,12] x_p[5,0] x_p[5,1] x_p[5,2] x_p[5,3] x_p[5,4] x_p[5,5] x_p[5,6]
 x_p[5,7] x_p[5,8] x_p[5,9] x_p[5,10] x_p[5,11] x_p[5,12] x_p[6,0] x_p[6,1]
 x_p[6,2] x_p[6,3] x_p[6,4] x_p[6,5] x_p[6,6] x_p[6,7] x_p[6,8] x_p[6,9]
 x_p[6,10] x_p[6,11] x_p[6,12] x_p[7,0] x_p[7,1] x_p[7,2] x_p[7,3] x_p[7,4]
 x_p[7,5] x_p[7,6] x_p[7,7] x_p[7,8] x_p[7,9] x_p[7,10] x_p[7,11] x_p[7,12]
 x_p[8,0] x_p[8,1] x_p[8,2] x_p[8,3] x_p[8,4] x_p[8,5] x_p[8,6] x_p[8,7]
 x_p[8,8] x_p[8,9] x_p[8,10] x_p[8,11] x_p[8,12] x_p[9,0] x_p[9,1] x_p[9,2]
 x_p[9,3] x_p[9,4] x_p[9,5] x_p[9,6] x_p[9,7] x_p[9,8] x_p[9,9] x_p[9,10]
 x_p[9,11] x_p[9,12] x_p[10,0] x_p[10,1] x_p[10,2] x_p[10,3] x_p[10,4]
 x_p[10,5] x_p[10,6] x_p[10,7] x_p[10,8] x_p[10,9] x_p[10,10] x_p[10,11]
 x_p[10,12] x_p[11,0] x_p[11,1] x_p[11,2] x_p[11,3] x_p[11,4] x_p[11,5]
 x_p[11,6] x_p[11,7] x_p[11,8] x_p[11,9] x_p[11,10] x_p[11,11] x_p[11,12]
 x_p[12,0] x_p[12,1] x_p[12,2] x_p[12,3] x_p[12,4] x_p[12,5] x_p[12,6]
 x_p[12,7] x_p[12,8] x_p[12,9] x_p[12,10] x_p[12,11] x_p[12,12] z_p[0,0]
 z_p[0,1] z_p[0,2] z_p[0,3] z_p[0,4] z_p[0,5] z_p[0,6] z_p[0,7] z_p[0,8]
 z_p[0,9] z_p[0,10] z_p[0,11] z_p[0,12] z_p[1,0] z_p[1,1] z_p[1,2] z_p[1,3]
 z_p[1,4] z_p[1,5] z_p[1,6] z_p[1,7] z_p[1,8] z_p[1,9] z_p[1,10] z_p[1,11]
 z_p[1,12] z_p[2,0] z_p[2,1] z_p[2,2] z_p[2,3] z_p[2,4] z_p[2,5] z_p[2,6]
 z_p[2,7] z_p[2,8] z_p[2,9] z_p[2,10] z_p[2,11] z_p[2,12] z_p[3,0] z_p[3,1]
 z_p[3,2] z_p[3,3] z_p[3,4] z_p[3,5] z_p[3,6] z_p[3,7] z_p[3,8] z_p[3,9]
 z_p[3,10] z_p[3,11] z_p[3,12] z_p[4,0] z_p[4,1] z_p[4,2] z_p[4,3] z_p[4,4]
 z_p[4,5] z_p[4,6] z_p[4,7] z_p[4,8] z_p[4,9] z_p[4,10] z_p[4,11] z_p[4,12]
 z_p[5,0] z_p[5,1] z_p[5,2] z_p[5,3] z_p[5,4] z_p[5,5] z_p[5,6] z_p[5,7]
 z_p[5,8] z_p[5,9] z_p[5,10] z_p[5,11] z_p[5,12] z_p[6,0] z_p[6,1] z_p[6,2]
 z_p[6,3] z_p[6,4] z_p[6,5] z_p[6,6] z_p[6,7] z_p[6,8] z_p[6,9] z_p[6,10]
 z_p[6,11] z_p[6,12] z_p[7,0] z_p[7,1] z_p[7,2] z_p[7,3] z_p[7,4] z_p[7,5]
 z_p[7,6] z_p[7,7] z_p[7,8] z_p[7,9] z_p[7,10] z_p[7,11] z_p[7,12] z_p[8,0]
 z_p[8,1] z_p[8,2] z_p[8,3] z_p[8,4] z_p[8,5] z_p[8,6] z_p[8,7] z_p[8,8]
 z_p[8,9] z_p[8,10] z_p[8,11] z_p[8,12] z_p[9,0] z_p[9,1] z_p[9,2] z_p[9,3]
 z_p[9,4] z_p[9,5] z_p[9,6] z_p[9,7] z_p[9,8] z_p[9,9] z_p[9,10] z_p[9,11]
 z_p[9,12] z_p[10,0] z_p[10,1] z_p[10,2] z_p[10,3] z_p[10,4] z_p[10,5]
 z_p[10,6] z_p[10,7] z_p[10,8] z_p[10,9] z_p[10,10] z_p[10,11] z_p[10,12]
 z_p[11,0] z_p[11,1] z_p[11,2] z_p[11,3] z_p[11,4] z_p[11,5] z_p[11,6]
 z_p[11,7] z_p[11,8] z_p[11,9] z_p[11,10] z_p[11,11] z_p[11,12] z_p[12,0]
 z_p[12,1] z_p[12,2] z_p[12,3] z_p[12,4] z_p[12,5] z_p[12,6] z_p[12,7]
 z_p[12,8] z_p[12,9] z_p[12,10] z_p[12,11] z_p[12,12] r[0,0,0] r[0,0,1]
 r[0,1,0] r[0,1,1] r[1,0,0] r[1,0,1] r[1,1,0] r[1,1,1] r[2,0,0] r[2,0,1]
 r[2,1,0] r[2,1,1] r[3,0,0] r[3,0,1] r[3,1,0] r[3,1,1] r[4,0,0] r[4,0,1]
 r[4,1,0] r[4,1,1] r[5,0,0] r[5,0,1] r[5,1,0] r[5,1,1] r[6,0,0] r[6,0,1]
 r[6,1,0] r[6,1,1] r[7,0,0] r[7,0,1] r[7,1,0] r[7,1,1] r[8,0,0] r[8,0,1]
 r[8,1,0] r[8,1,1] r[9,0,0] r[9,0,1] r[9,1,0] r[9,1,1] r[10,0,0] r[10,0,1]
 r[10,1,0] r[10,1,1] r[11,0,0] r[11,0,1] r[11,1,0] r[11,1,1] r[12,0,0]
 r[12,0,1] r[12,1,0] r[12,1,1] g[0] g[1] g[2] g[3] g[4] g[5] g[6] g[7] g[8]
 g[9] g[10] g[11] g[12] h[0,0] h[0,1] h[0,2] h[0,3] h[0,4] h[0,5] h[0,6]
 h[0,7] h[0,8] h[0,9] h[0,10] h[0,11] h[0,12] h[1,0] h[1,1] h[1,2] h[1,3]
 h[1,4] h[1,5] h[1,6] h[1,7] h[1,8] h[1,9] h[1,10] h[1,11] h[1,12] h[2,0]
 h[2,1] h[2,2] h[2,3] h[2,4] h[2,5] h[2,6] h[2,7] h[2,8] h[2,9] h[2,10]
 h[2,11] h[2,12] h[3,0] h[3,1] h[3,2] h[3,3] h[3,4] h[3,5] h[3,6] h[3,7]
 h[3,8] h[3,9] h[3,10] h[3,11] h[3,12] h[4,0] h[4,1] h[4,2] h[4,3] h[4,4]
 h[4,5] h[4,6] h[4,7] h[4,8] h[4,9] h[4,10] h[4,11] h[4,12] h[5,0] h[5,1]
 h[5,2] h[5,3] h[5,4] h[5,5] h[5,6] h[5,7] h[5,8] h[5,9] h[5,10] h[5,11]
 h[5,12] h[6,0] h[6,1] h[6,2] h[6,3] h[6,4] h[6,5] h[6,6] h[6,7] h[6,8]
 h[6,9] h[6,10] h[6,11] h[6,12] h[7,0] h[7,1] h[7,2] h[7,3] h[7,4] h[7,5]
 h[7,6] h[7,7] h[7,8] h[7,9] h[7,10] h[7,11] h[7,12] h[8,0] h[8,1] h[8,2]
 h[8,3] h[8,4] h[8,5] h[8,6] h[8,7] h[8,8] h[8,9] h[8,10] h[8,11] h[8,12]
 h[9,0] h[9,1] h[9,2] h[9,3] h[9,4] h[9,5] h[9,6] h[9,7] h[9,8] h[9,9]
 h[9,10] h[9,11] h[9,12] h[10,0] h[10,1] h[10,2] h[10,3] h[10,4] h[10,5]
 h[10,6] h[10,7] h[10,8] h[10,9] h[10,10] h[10,11] h[10,12] h[11,0] h[11,1]
 h[11,2] h[11,3] h[11,4] h[11,5] h[11,6] h[11,7] h[11,8] h[11,9] h[11,10]
 h[11,11] h[11,12] h[12,0] h[12,1] h[12,2] h[12,3] h[12,4] h[12,5] h[12,6]
 h[12,7] h[12,8] h[12,9] h[12,10] h[12,11] h[12,12] o[0,0] o[0,1] o[0,2]
 o[0,3] o[0,4] o[0,5] o[0,6] o[0,7] o[0,8] o[0,9] o[0,10] o[0,11] o[0,12]
 o[1,0] o[1,1] o[1,2] o[1,3] o[1,4] o[1,5] o[1,6] o[1,7] o[1,8] o[1,9]
 o[1,10] o[1,11] o[1,12] o[2,0] o[2,1] o[2,2] o[2,3] o[2,4] o[2,5] o[2,6]
 o[2,7] o[2,8] o[2,9] o[2,10] o[2,11] o[2,12] o[3,0] o[3,1] o[3,2] o[3,3]
 o[3,4] o[3,5] o[3,6] o[3,7] o[3,8] o[3,9] o[3,10] o[3,11] o[3,12] o[4,0]
 o[4,1] o[4,2] o[4,3] o[4,4] o[4,5] o[4,6] o[4,7] o[4,8] o[4,9] o[4,10]
 o[4,11] o[4,12] o[5,0] o[5,1] o[5,2] o[5,3] o[5,4] o[5,5] o[5,6] o[5,7]
 o[5,8] o[5,9] o[5,10] o[5,11] o[5,12] o[6,0] o[6,1] o[6,2] o[6,3] o[6,4]
 o[6,5] o[6,6] o[6,7] o[6,8] o[6,9] o[6,10] o[6,11] o[6,12] o[7,0] o[7,1]
 o[7,2] o[7,3] o[7,4] o[7,5] o[7,6] o[7,7] o[7,8] o[7,9] o[7,10] o[7,11]
 o[7,12] o[8,0] o[8,1] o[8,2] o[8,3] o[8,4] o[8,5] o[8,6] o[8,7] o[8,8]
 o[8,9] o[8,10] o[8,11] o[8,12] o[9,0] o[9,1] o[9,2] o[9,3] o[9,4] o[9,5]
 o[9,6] o[9,7] o[9,8] o[9,9] o[9,10] o[9,11] o[9,12] o[10,0] o[10,1]
 o[10,2] o[10,3] o[10,4] o[10,5] o[10,6] o[10,7] o[10,8] o[10,9] o[10,10]
 o[10,11] o[10,12] o[11,0] o[11,1] o[11,2] o[11,3] o[11,4] o[11,5] o[11,6]
 o[11,7] o[11,8] o[11,9] o[11,10] o[11,11] o[11,12] o[12,0] o[12,1] o[12,2]
 o[12,3] o[12,4] o[12,5] o[12,6] o[12,7] o[12,8] o[12,9] o[12,10] o[12,11]
 o[12,12] s[0,0] s[0,1] s[0,2] s[0,3] s[0,4] s[0,5] s[0,6] s[0,7] s[0,8]
 s[0,9] s[0,10] s[0,11] s[0,12] s[1,0] s[1,1] s[1,2] s[1,3] s[1,4] s[1,5]
 s[1,6] s[1,7] s[1,8] s[1,9] s[1,10] s[1,11] s[1,12] s[2,0] s[2,1] s[2,2]
 s[2,3] s[2,4] s[2,5] s[2,6] s[2,7] s[2,8] s[2,9] s[2,10] s[2,11] s[2,12]
 s[3,0] s[3,1] s[3,2] s[3,3] s[3,4] s[3,5] s[3,6] s[3,7] s[3,8] s[3,9]
 s[3,10] s[3,11] s[3,12] s[4,0] s[4,1] s[4,2] s[4,3] s[4,4] s[4,5] s[4,6]
 s[4,7] s[4,8] s[4,9] s[4,10] s[4,11] s[4,12] s[5,0] s[5,1] s[5,2] s[5,3]
 s[5,4] s[5,5] s[5,6] s[5,7] s[5,8] s[5,9] s[5,10] s[5,11] s[5,12] s[6,0]
 s[6,1] s[6,2] s[6,3] s[6,4] s[6,5] s[6,6] s[6,7] s[6,8] s[6,9] s[6,10]
 s[6,11] s[6,12] s[7,0] s[7,1] s[7,2] s[7,3] s[7,4] s[7,5] s[7,6] s[7,7]
 s[7,8] s[7,9] s[7,10] s[7,11] s[7,12] s[8,0] s[8,1] s[8,2] s[8,3] s[8,4]
 s[8,5] s[8,6] s[8,7] s[8,8] s[8,9] s[8,10] s[8,11] s[8,12] s[9,0] s[9,1]
 s[9,2] s[9,3] s[9,4] s[9,5] s[9,6] s[9,7] s[9,8] s[9,9] s[9,10] s[9,11]
 s[9,12] s[10,0] s[10,1] s[10,2] s[10,3] s[10,4] s[10,5] s[10,6] s[10,7]
 s[10,8] s[10,9] s[10,10] s[10,11] s[10,12] s[11,0] s[11,1] s[11,2] s[11,3]
 s[11,4] s[11,5] s[11,6] s[11,7] s[11,8] s[11,9] s[11,10] s[11,11] s[11,12]
 s[12,0] s[12,1] s[12,2] s[12,3] s[12,4] s[12,5] s[12,6] s[12,7] s[12,8]
 s[12,9] s[12,10] s[12,11] s[12,12] eta1[0,0] eta1[0,1] eta1[0,2] eta1[0,3]
 eta1[0,4] eta1[0,5] eta1[0,6] eta1[0,7] eta1[0,8] eta1[0,9] eta1[0,10]
 eta1[0,11] eta1[0,12] eta1[1,0] eta1[1,1] eta1[1,2] eta1[1,3] eta1[1,4]
 eta1[1,5] eta1[1,6] eta1[1,7] eta1[1,8] eta1[1,9] eta1[1,10] eta1[1,11]
 eta1[1,12] eta1[2,0] eta1[2,1] eta1[2,2] eta1[2,3] eta1[2,4] eta1[2,5]
 eta1[2,6] eta1[2,7] eta1[2,8] eta1[2,9] eta1[2,10] eta1[2,11] eta1[2,12]
 eta1[3,0] eta1[3,1] eta1[3,2] eta1[3,3] eta1[3,4] eta1[3,5] eta1[3,6]
 eta1[3,7] eta1[3,8] eta1[3,9] eta1[3,10] eta1[3,11] eta1[3,12] eta1[4,0]
 eta1[4,1] eta1[4,2] eta1[4,3] eta1[4,4] eta1[4,5] eta1[4,6] eta1[4,7]
 eta1[4,8] eta1[4,9] eta1[4,10] eta1[4,11] eta1[4,12] eta1[5,0] eta1[5,1]
 eta1[5,2] eta1[5,3] eta1[5,4] eta1[5,5] eta1[5,6] eta1[5,7] eta1[5,8]
 eta1[5,9] eta1[5,10] eta1[5,11] eta1[5,12] eta1[6,0] eta1[6,1] eta1[6,2]
 eta1[6,3] eta1[6,4] eta1[6,5] eta1[6,6] eta1[6,7] eta1[6,8] eta1[6,9]
 eta1[6,10] eta1[6,11] eta1[6,12] eta1[7,0] eta1[7,1] eta1[7,2] eta1[7,3]
 eta1[7,4] eta1[7,5] eta1[7,6] eta1[7,7] eta1[7,8] eta1[7,9] eta1[7,10]
 eta1[7,11] eta1[7,12] eta1[8,0] eta1[8,1] eta1[8,2] eta1[8,3] eta1[8,4]
 eta1[8,5] eta1[8,6] eta1[8,7] eta1[8,8] eta1[8,9] eta1[8,10] eta1[8,11]
 eta1[8,12] eta1[9,0] eta1[9,1] eta1[9,2] eta1[9,3] eta1[9,4] eta1[9,5]
 eta1[9,6] eta1[9,7] eta1[9,8] eta1[9,9] eta1[9,10] eta1[9,11] eta1[9,12]
 eta1[10,0] eta1[10,1] eta1[10,2] eta1[10,3] eta1[10,4] eta1[10,5]
 eta1[10,6] eta1[10,7] eta1[10,8] eta1[10,9] eta1[10,10] eta1[10,11]
 eta1[10,12] eta1[11,0] eta1[11,1] eta1[11,2] eta1[11,3] eta1[11,4]
 eta1[11,5] eta1[11,6] eta1[11,7] eta1[11,8] eta1[11,9] eta1[11,10]
 eta1[11,11] eta1[11,12] eta1[12,0] eta1[12,1] eta1[12,2] eta1[12,3]
 eta1[12,4] eta1[12,5] eta1[12,6] eta1[12,7] eta1[12,8] eta1[12,9]
 eta1[12,10] eta1[12,11] eta1[12,12] eta2[0,0] eta2[0,1] eta2[0,2]
 eta2[0,3] eta2[0,4] eta2[0,5] eta2[0,6] eta2[0,7] eta2[0,8] eta2[0,9]
 eta2[0,10] eta2[0,11] eta2[0,12] eta2[1,0] eta2[1,1] eta2[1,2] eta2[1,3]
 eta2[1,4] eta2[1,5] eta2[1,6] eta2[1,7] eta2[1,8] eta2[1,9] eta2[1,10]
 eta2[1,11] eta2[1,12] eta2[2,0] eta2[2,1] eta2[2,2] eta2[2,3] eta2[2,4]
 eta2[2,5] eta2[2,6] eta2[2,7] eta2[2,8] eta2[2,9] eta2[2,10] eta2[2,11]
 eta2[2,12] eta2[3,0] eta2[3,1] eta2[3,2] eta2[3,3] eta2[3,4] eta2[3,5]
 eta2[3,6] eta2[3,7] eta2[3,8] eta2[3,9] eta2[3,10] eta2[3,11] eta2[3,12]
 eta2[4,0] eta2[4,1] eta2[4,2] eta2[4,3] eta2[4,4] eta2[4,5] eta2[4,6]
 eta2[4,7] eta2[4,8] eta2[4,9] eta2[4,10] eta2[4,11] eta2[4,12] eta2[5,0]
 eta2[5,1] eta2[5,2] eta2[5,3] eta2[5,4] eta2[5,5] eta2[5,6] eta2[5,7]
 eta2[5,8] eta2[5,9] eta2[5,10] eta2[5,11] eta2[5,12] eta2[6,0] eta2[6,1]
 eta2[6,2] eta2[6,3] eta2[6,4] eta2[6,5] eta2[6,6] eta2[6,7] eta2[6,8]
 eta2[6,9] eta2[6,10] eta2[6,11] eta2[6,12] eta2[7,0] eta2[7,1] eta2[7,2]
 eta2[7,3] eta2[7,4] eta2[7,5] eta2[7,6] eta2[7,7] eta2[7,8] eta2[7,9]
 eta2[7,10] eta2[7,11] eta2[7,12] eta2[8,0] eta2[8,1] eta2[8,2] eta2[8,3]
 eta2[8,4] eta2[8,5] eta2[8,6] eta2[8,7] eta2[8,8] eta2[8,9] eta2[8,10]
 eta2[8,11] eta2[8,12] eta2[9,0] eta2[9,1] eta2[9,2] eta2[9,3] eta2[9,4]
 eta2[9,5] eta2[9,6] eta2[9,7] eta2[9,8] eta2[9,9] eta2[9,10] eta2[9,11]
 eta2[9,12] eta2[10,0] eta2[10,1] eta2[10,2] eta2[10,3] eta2[10,4]
 eta2[10,5] eta2[10,6] eta2[10,7] eta2[10,8] eta2[10,9] eta2[10,10]
 eta2[10,11] eta2[10,12] eta2[11,0] eta2[11,1] eta2[11,2] eta2[11,3]
 eta2[11,4] eta2[11,5] eta2[11,6] eta2[11,7] eta2[11,8] eta2[11,9]
 eta2[11,10] eta2[11,11] eta2[11,12] eta2[12,0] eta2[12,1] eta2[12,2]
 eta2[12,3] eta2[12,4] eta2[12,5] eta2[12,6] eta2[12,7] eta2[12,8]
 eta2[12,9] eta2[12,10] eta2[12,11] eta2[12,12] beta1[0,0] beta1[0,1]
 beta1[0,2] beta1[0,3] beta1[0,4] beta1[0,5] beta1[0,6] beta1[0,7]
 beta1[0,8] beta1[0,9] beta1[0,10] beta1[0,11] beta1[0,12] beta1[1,0]
 beta1[1,1] beta1[1,2] beta1[1,3] beta1[1,4] beta1[1,5] beta1[1,6]
 beta1[1,7] beta1[1,8] beta1[1,9] beta1[1,10] beta1[1,11] beta1[1,12]
 beta1[2,0] beta1[2,1] beta1[2,2] beta1[2,3] beta1[2,4] beta1[2,5]
 beta1[2,6] beta1[2,7] beta1[2,8] beta1[2,9] beta1[2,10] beta1[2,11]
 beta1[2,12] beta1[3,0] beta1[3,1] beta1[3,2] beta1[3,3] beta1[3,4]
 beta1[3,5] beta1[3,6] beta1[3,7] beta1[3,8] beta1[3,9] beta1[3,10]
 beta1[3,11] beta1[3,12] beta1[4,0] beta1[4,1] beta1[4,2] beta1[4,3]
 beta1[4,4] beta1[4,5] beta1[4,6] beta1[4,7] beta1[4,8] beta1[4,9]
 beta1[4,10] beta1[4,11] beta1[4,12] beta1[5,0] beta1[5,1] beta1[5,2]
 beta1[5,3] beta1[5,4] beta1[5,5] beta1[5,6] beta1[5,7] beta1[5,8]
 beta1[5,9] beta1[5,10] beta1[5,11] beta1[5,12] beta1[6,0] beta1[6,1]
 beta1[6,2] beta1[6,3] beta1[6,4] beta1[6,5] beta1[6,6] beta1[6,7]
 beta1[6,8] beta1[6,9] beta1[6,10] beta1[6,11] beta1[6,12] beta1[7,0]
 beta1[7,1] beta1[7,2] beta1[7,3] beta1[7,4] beta1[7,5] beta1[7,6]
 beta1[7,7] beta1[7,8] beta1[7,9] beta1[7,10] beta1[7,11] beta1[7,12]
 beta1[8,0] beta1[8,1] beta1[8,2] beta1[8,3] beta1[8,4] beta1[8,5]
 beta1[8,6] beta1[8,7] beta1[8,8] beta1[8,9] beta1[8,10] beta1[8,11]
 beta1[8,12] beta1[9,0] beta1[9,1] beta1[9,2] beta1[9,3] beta1[9,4]
 beta1[9,5] beta1[9,6] beta1[9,7] beta1[9,8] beta1[9,9] beta1[9,10]
 beta1[9,11] beta1[9,12] beta1[10,0] beta1[10,1] beta1[10,2] beta1[10,3]
 beta1[10,4] beta1[10,5] beta1[10,6] beta1[10,7] beta1[10,8] beta1[10,9]
 beta1[10,10] beta1[10,11] beta1[10,12] beta1[11,0] beta1[11,1] beta1[11,2]
 beta1[11,3] beta1[11,4] beta1[11,5] beta1[11,6] beta1[11,7] beta1[11,8]
 beta1[11,9] beta1[11,10] beta1[11,11] beta1[11,12] beta1[12,0] beta1[12,1]
 beta1[12,2] beta1[12,3] beta1[12,4] beta1[12,5] beta1[12,6] beta1[12,7]
 beta1[12,8] beta1[12,9] beta1[12,10] beta1[12,11] beta1[12,12] beta2[0,0]
 beta2[0,1] beta2[0,2] beta2[0,3] beta2[0,4] beta2[0,5] beta2[0,6]
 beta2[0,7] beta2[0,8] beta2[0,9] beta2[0,10] beta2[0,11] beta2[0,12]
 beta2[1,0] beta2[1,1] beta2[1,2] beta2[1,3] beta2[1,4] beta2[1,5]
 beta2[1,6] beta2[1,7] beta2[1,8] beta2[1,9] beta2[1,10] beta2[1,11]
 beta2[1,12] beta2[2,0] beta2[2,1] beta2[2,2] beta2[2,3] beta2[2,4]
 beta2[2,5] beta2[2,6] beta2[2,7] beta2[2,8] beta2[2,9] beta2[2,10]
 beta2[2,11] beta2[2,12] beta2[3,0] beta2[3,1] beta2[3,2] beta2[3,3]
 beta2[3,4] beta2[3,5] beta2[3,6] beta2[3,7] beta2[3,8] beta2[3,9]
 beta2[3,10] beta2[3,11] beta2[3,12] beta2[4,0] beta2[4,1] beta2[4,2]
 beta2[4,3] beta2[4,4] beta2[4,5] beta2[4,6] beta2[4,7] beta2[4,8]
 beta2[4,9] beta2[4,10] beta2[4,11] beta2[4,12] beta2[5,0] beta2[5,1]
 beta2[5,2] beta2[5,3] beta2[5,4] beta2[5,5] beta2[5,6] beta2[5,7]
 beta2[5,8] beta2[5,9] beta2[5,10] beta2[5,11] beta2[5,12] beta2[6,0]
 beta2[6,1] beta2[6,2] beta2[6,3] beta2[6,4] beta2[6,5] beta2[6,6]
 beta2[6,7] beta2[6,8] beta2[6,9] beta2[6,10] beta2[6,11] beta2[6,12]
 beta2[7,0] beta2[7,1] beta2[7,2] beta2[7,3] beta2[7,4] beta2[7,5]
 beta2[7,6] beta2[7,7] beta2[7,8] beta2[7,9] beta2[7,10] beta2[7,11]
 beta2[7,12] beta2[8,0] beta2[8,1] beta2[8,2] beta2[8,3] beta2[8,4]
 beta2[8,5] beta2[8,6] beta2[8,7] beta2[8,8] beta2[8,9] beta2[8,10]
 beta2[8,11] beta2[8,12] beta2[9,0] beta2[9,1] beta2[9,2] beta2[9,3]
 beta2[9,4] beta2[9,5] beta2[9,6] beta2[9,7] beta2[9,8] beta2[9,9]
 beta2[9,10] beta2[9,11] beta2[9,12] beta2[10,0] beta2[10,1] beta2[10,2]
 beta2[10,3] beta2[10,4] beta2[10,5] beta2[10,6] beta2[10,7] beta2[10,8]
 beta2[10,9] beta2[10,10] beta2[10,11] beta2[10,12] beta2[11,0] beta2[11,1]
 beta2[11,2] beta2[11,3] beta2[11,4] beta2[11,5] beta2[11,6] beta2[11,7]
 beta2[11,8] beta2[11,9] beta2[11,10] beta2[11,11] beta2[11,12] beta2[12,0]
 beta2[12,1] beta2[12,2] beta2[12,3] beta2[12,4] beta2[12,5] beta2[12,6]
 beta2[12,7] beta2[12,8] beta2[12,9] beta2[12,10] beta2[12,11] beta2[12,12]
 gamma[0] gamma[1] gamma[2] gamma[3] gamma[4] gamma[5] gamma[6] gamma[7]
 gamma[8] gamma[9] gamma[10] gamma[11] gamma[12] m[0,0] m[0,1] m[0,2]
 m[0,3] m[0,4] m[0,5] m[0,6] m[0,7] m[0,8] m[0,9] m[0,10] m[0,11] m[0,12]
 m[1,0] m[1,1] m[1,2] m[1,3] m[1,4] m[1,5] m[1,6] m[1,7] m[1,8] m[1,9]
 m[1,10] m[1,11] m[1,12] m[2,0] m[2,1] m[2,2] m[2,3] m[2,4] m[2,5] m[2,6]
 m[2,7] m[2,8] m[2,9] m[2,10] m[2,11] m[2,12] m[3,0] m[3,1] m[3,2] m[3,3]
 m[3,4] m[3,5] m[3,6] m[3,7] m[3,8] m[3,9] m[3,10] m[3,11] m[3,12] m[4,0]
 m[4,1] m[4,2] m[4,3] m[4,4] m[4,5] m[4,6] m[4,7] m[4,8] m[4,9] m[4,10]
 m[4,11] m[4,12] m[5,0] m[5,1] m[5,2] m[5,3] m[5,4] m[5,5] m[5,6] m[5,7]
 m[5,8] m[5,9] m[5,10] m[5,11] m[5,12] m[6,0] m[6,1] m[6,2] m[6,3] m[6,4]
 m[6,5] m[6,6] m[6,7] m[6,8] m[6,9] m[6,10] m[6,11] m[6,12] m[7,0] m[7,1]
 m[7,2] m[7,3] m[7,4] m[7,5] m[7,6] m[7,7] m[7,8] m[7,9] m[7,10] m[7,11]
 m[7,12] m[8,0] m[8,1] m[8,2] m[8,3] m[8,4] m[8,5] m[8,6] m[8,7] m[8,8]
 m[8,9] m[8,10] m[8,11] m[8,12] m[9,0] m[9,1] m[9,2] m[9,3] m[9,4] m[9,5]
 m[9,6] m[9,7] m[9,8] m[9,9] m[9,10] m[9,11] m[9,12] m[10,0] m[10,1]
 m[10,2] m[10,3] m[10,4] m[10,5] m[10,6] m[10,7] m[10,8] m[10,9] m[10,10]
 m[10,11] m[10,12] m[11,0] m[11,1] m[11,2] m[11,3] m[11,4] m[11,5] m[11,6]
 m[11,7] m[11,8] m[11,9] m[11,10] m[11,11] m[11,12] m[12,0] m[12,1] m[12,2]
 m[12,3] m[12,4] m[12,5] m[12,6] m[12,7] m[12,8] m[12,9] m[12,10] m[12,11]
 m[12,12]
End
