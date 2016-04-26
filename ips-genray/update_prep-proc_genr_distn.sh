#BH: 090927
#BH: 130912 update
#BH: 131021 update
#my_ips_trunk/components/rf/genray/src/
#Execute in my_ips_trunk/components/rf/genray/src/
#Copies necessary files for prep/proc genray from ~/genray_cswim_svn/trunk

PWD_RETURN=$PWD
HM=$HOME
GENR_DISTN=$HM'/genray_cswim_svn/trunk'
cd $GENR_DISTN
cp -p read_write_genray_input_prep.f bcast.f $PWD_RETURN
cp -p \
adj_nml.i           lsc_approach_nml.i                     name_tokamak.i  \
cone_nml.i          name_adj.i                             name_uniform_mesh_profiles.i  \
dinit_nml.i         name_eccone.i                          one_nml.i  \
edge_prof.i         name_edge_prof_nml.i                   onetwo_nml.i  \
edge_prof_nml.i     name_genr.i                            output_nml.i  \
edge_prof_no_nml.i  name_grill.i                           param.i  \
emissa_nml.i        name.i                                 rkutta.i  \
grill_nml.i         name_lsc_approach_nml.i                scatnper_nml.i  \
ions_nml.i          name_non_uniform_mesh_profiles_line.i  six_nml.i  \
$PWD_RETURN/genray_includes
cd $PWD_RETURN

