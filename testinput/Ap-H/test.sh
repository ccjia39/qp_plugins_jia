rm -r H.ezfio
source /home/cjia/quantum_package/qp2/quantum_package.rc

qp create_ezfio -b "cc-pvtz" -m 2 H.xyz -o H.ezfio
qp run scf #| tee log_scf
qp set cippres finput_cippres h.xml
#qp set electrons elec_alpha_num 1
#qp set electrons elec_beta_num 0
#qp set cippres ici1 1
#qp set cippres ici2 1

qp run cippres_gencsf
qp set cippres ifcsf 1
qp set cippres ici1 1
qp set cippres ici2 1
qp set cippres n_sta_cistate_analysis 5
qp run cippres_runci
