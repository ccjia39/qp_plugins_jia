rm -r He.ezfio
source /home/cjia/quantum_package/qp2/quantum_package.rc


qp create_ezfio -b "cc-pvdz" He.xyz -o He.ezfio
qp run scf #| tee log_scf
qp set cippres finput_cippres he.xml
#qp set electrons elec_alpha_num 1
#qp set electrons elec_beta_num 0

qp run cippres_gencsf
qp set cippres ifcsf 1
qp set cippres ici1 1
qp set cippres ici2 1
qp set cippres n_sta_cistate_analysis 5
qp run cippres_runci

#qp set cippres finput_coll coll_input.xml
#qp set cippres ici1 1
#qp set cippres ici2 1
#qp run cippres_setup_collision
