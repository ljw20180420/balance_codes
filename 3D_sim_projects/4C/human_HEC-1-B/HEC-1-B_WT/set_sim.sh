cd /home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/human_HEC-1-B/HEC-1-B_WT
/home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/set_range_and_CTCF_ver3/build/set_range_and_ctcf_ver3 -start 140160700 -end 140920000 \
-C4s_rep1 4C/4C_HEC-1-B_a12_WT_rep1.bedGraph,4C/4C_HEC-1-B_HS5-1_WT_rep1.bedGraph -C4s_rep2 4C/4C_HEC-1-B_a12_WT_rep2.bedGraph,4C/4C_HEC-1-B_HS5-1_WT_rep2.bedGraph \
-CTCF_peaks chip/HEC-1-B_WT_CTCF.bam_peaks.txt -nipbls ~/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/human_HEC-1-B/nipbl/Hec1bwt_nipbl_treat_afterfiting_chr5.bdg \
-mut_inform mut_inform -chromosome chr5 -stiffnesses 2 -densities 0.2 -processivities 400 -separations 400 -D3_per_D1s 1250 -capture_ratios 2.000000 -record_1D 50000
cd /home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/human_HEC-1-B/HEC-1-B_WT/run1
/home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/loop_extrusion_3d_sim_ver5/build/loop_extrusion_3d_sim_ver5 \
-pos_file /home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/human_HEC-1-B/HEC-1-B_WT/pos_file -THR_MAX 1 &
cd /home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/human_HEC-1-B/HEC-1-B_WT/run2
/home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/loop_extrusion_3d_sim_ver5/build/loop_extrusion_3d_sim_ver5 \
-pos_file /home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/human_HEC-1-B/HEC-1-B_WT/pos_file -THR_MAX 1 &
cd /home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/human_HEC-1-B/HEC-1-B_WT/run3
/home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/loop_extrusion_3d_sim_ver5/build/loop_extrusion_3d_sim_ver5 \
-pos_file /home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/human_HEC-1-B/HEC-1-B_WT/pos_file -THR_MAX 1 &
cd /home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/human_HEC-1-B/HEC-1-B_WT/run4
/home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/loop_extrusion_3d_sim_ver5/build/loop_extrusion_3d_sim_ver5 \
-pos_file /home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/human_HEC-1-B/HEC-1-B_WT/pos_file -THR_MAX 1 &
cd /home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/human_HEC-1-B/HEC-1-B_WT/run5
/home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/loop_extrusion_3d_sim_ver5/build/loop_extrusion_3d_sim_ver5 \
-pos_file /home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/human_HEC-1-B/HEC-1-B_WT/pos_file -THR_MAX 1 &
cd /home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/human_HEC-1-B/HEC-1-B_WT/run6
/home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/loop_extrusion_3d_sim_ver5/build/loop_extrusion_3d_sim_ver5 \
-pos_file /home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/human_HEC-1-B/HEC-1-B_WT/pos_file -THR_MAX 1 &
cd /home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/human_HEC-1-B/HEC-1-B_WT/run7
/home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/loop_extrusion_3d_sim_ver5/build/loop_extrusion_3d_sim_ver5 \
-pos_file /home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/human_HEC-1-B/HEC-1-B_WT/pos_file -THR_MAX 1 &
wait
