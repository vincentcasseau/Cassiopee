import KCore.Dist as Dist
from KCore.config import *
#==============================================================================
# Fichiers c++
#==============================================================================
cpp_srcs = []
#==============================================================================
cpp_srcs1 = ['XCore/scotch/library_mesh_io_habo.c', 'XCore/scotch/vgraph_separate_gp.c', 'XCore/scotch/vgraph_separate_es.c', 'XCore/scotch/arch_vcmplt.c', 'XCore/scotch/library_mesh_io_scot_f.c', 'XCore/scotch/library_version.c', 'XCore/scotch/parser_yy.c', 'XCore/scotch/mesh_graph.c', 'XCore/scotch/library_dgraph_order_gather.c', 'XCore/scotch/common_memory.c', 'XCore/scotch/kgraph_map_ex.c', 'XCore/scotch/dgraph_fold_comm.c', 'XCore/scotch/library_dgraph_order.c', 'XCore/scotch/library_dgraph.c', 'XCore/scotch/library_dgraph_order_io_block_f.c', 'XCore/scotch/hmesh_induce.c', 'XCore/scotch/dgraph_induce.c', 'XCore/scotch/dgraph_gather.c', 'XCore/scotch/kgraph_map_ml.c', 'XCore/scotch/hdgraph_order_si.c', 'XCore/scotch/arch_vhcub.c', 'XCore/scotch/graph_clone.c', 'XCore/scotch/hmesh_order_gr.c', 'XCore/scotch/dgraph_band.c', 'XCore/scotch/graph_coarsen_edge.c', 'XCore/scotch/kdgraph_map_rb.c', 'XCore/scotch/library_dgraph_redist_f.c', 'XCore/scotch/kgraph_check.c', 'XCore/scotch/mesh_io_scot.c', 'XCore/scotch/hmesh_mesh.c', 'XCore/scotch/vgraph_separate_df.c', 'XCore/scotch/vdgraph_separate_zr.c', 'XCore/scotch/library_version_f.c', 'XCore/scotch/kgraph_map_st.c', 'XCore/scotch/dgraph_coarsen.c', 'XCore/scotch/library_dgraph_induce.c', 'XCore/scotch/mesh_io.c', 'XCore/scotch/library_dgraph_gather.c', 'XCore/scotch/library_dgraph_grow.c', 'XCore/scotch/arch_mesh.c', 'XCore/scotch/graph_diam.c', 'XCore/scotch/library_dgraph_f.c', 'XCore/scotch/library_graph_io_habo_f.c', 'XCore/scotch/graph_coarsen.c', 'XCore/scotch/common_file_decompress.c', 'XCore/scotch/bgraph_bipart_zr.c', 'XCore/scotch/library_dgraph_map_view.c', 'XCore/scotch/library_graph_map.c', 'XCore/scotch/dgraph_band_grow.c', 'XCore/scotch/mesh_check.c', 'XCore/scotch/library_dgraph_map_f.c', 'XCore/scotch/mapping.c', 'XCore/scotch/kdgraph_map_rb_map.c', 'XCore/scotch/vmesh_separate_gr.c', 'XCore/scotch/common_sort.c', 'XCore/scotch/library_memory_f.c', 'XCore/scotch/wgraph_store.c', 'XCore/scotch/vdgraph_store.c', 'XCore/scotch/dgraph_io_save.c', 'XCore/scotch/arch_cmplt.c', 'XCore/scotch/arch.c', 'XCore/scotch/fibo.c', 'XCore/scotch/kdgraph.c', 'XCore/scotch/vmesh.c', 'XCore/scotch/library_dgraph_order_io_block.c', 'XCore/scotch/dgraph_match_sync_coll.c', 'XCore/scotch/graph_induce.c', 'XCore/scotch/dgraph_scatter.c', 'XCore/scotch/bdgraph_check.c', 'XCore/scotch/hgraph.c', 'XCore/scotch/hmesh_order_bl.c', 'XCore/scotch/library_graph_base.c', 'XCore/scotch/graph_io_habo.c', 'XCore/scotch/library_graph_order_f.c', 'XCore/scotch/bgraph_bipart_gp.c', 'XCore/scotch/dgraph_ghst.c', 'XCore/scotch/kgraph_store.c', 'XCore/scotch/library_graph_coarsen.c', 'XCore/scotch/library_mesh_graph.c', 'XCore/scotch/dgraph_check.c', 'XCore/scotch/library_arch_f.c', 'XCore/scotch/dgraph_fold_dup.c', 'XCore/scotch/library_dgraph_band_f.c', 'XCore/scotch/mapping_io.c', 'XCore/scotch/library_dgraph_map.c', 'XCore/scotch/vdgraph_separate_df.c', 'XCore/scotch/hdgraph_order_nd.c', 'XCore/scotch/wgraph_part_fm.c', 'XCore/scotch/common_error.c', 'XCore/scotch/hmesh_check.c', 'XCore/scotch/dgraph_fold.c', 'XCore/scotch/hmesh_order_hd.c', 'XCore/scotch/bdgraph_bipart_bd.c', 'XCore/scotch/dgraph_redist.c', 'XCore/scotch/hgraph_order_bl.c', 'XCore/scotch/dgraph_allreduce.c', 'XCore/scotch/vgraph.c', 'XCore/scotch/dgraph_gather_all.c', 'XCore/scotch/library_graph_f.c', 'XCore/scotch/library_dgraph_io_save_f.c', 'XCore/scotch/library_graph_io_habo.c', 'XCore/scotch/kdgraph_gather.c', 'XCore/scotch/bgraph_bipart_df.c', 'XCore/scotch/graph_match_scan.c', 'XCore/scotch/library_dgraph_scatter_f.c', 'XCore/scotch/kdgraph_map_st.c', 'XCore/scotch/vgraph_separate_ml.c', 'XCore/scotch/dgraph.c', 'XCore/scotch/hmesh_order_st.c', 'XCore/scotch/library_graph_diam.c', 'XCore/scotch/hgraph_order_hd.c', 'XCore/scotch/library_graph_check.c', 'XCore/scotch/library_dmapping.c', 'XCore/scotch/hall_order_hd.c', 'XCore/scotch/kgraph_map_cp.c', 'XCore/scotch/library_dgraph_halo.c', 'XCore/scotch/vmesh_separate_ml.c', 'XCore/scotch/graph.c', 'XCore/scotch/library_dgraph_redist.c', 'XCore/scotch/vdgraph_separate_sq.c', 'XCore/scotch/comm.c', 'XCore/scotch/hgraph_order_kp.c', 'XCore/scotch/dmapping.c', 'XCore/scotch/library_graph_map_io.c', 'XCore/scotch/hgraph_order_st.c', 'XCore/scotch/library_geom.c', 'XCore/scotch/hmesh_order_hf.c', 'XCore/scotch/dgraph_match_sync_ptop.c', 'XCore/scotch/library_parser.c', 'XCore/scotch/library_mesh.c', 'XCore/scotch/vgraph_separate_st.c', 'XCore/scotch/bdgraph_store.c', 'XCore/scotch/bgraph_check.c', 'XCore/scotch/order_io.c', 'XCore/scotch/vmesh_separate_st.c', 'XCore/scotch/library_graph_map_view.c', 'XCore/scotch/hdgraph.c', 'XCore/scotch/library_dgraph_order_gather_f.c', 'XCore/scotch/bdgraph_bipart_zr.c', 'XCore/scotch/library_graph.c', 'XCore/scotch/vgraph_separate_gg.c', 'XCore/scotch/library_dgraph_order_tree_dist_f.c', 'XCore/scotch/kgraph_map_fm.c', 'XCore/scotch/graph_io_chac.c', 'XCore/scotch/library_dgraph_order_io.c', 'XCore/scotch/vdgraph_gather_all.c', 'XCore/scotch/hgraph_order_hf.c', 'XCore/scotch/vmesh_separate_gg.c', 'XCore/scotch/graph_band.c', 'XCore/scotch/hall_order_hf.c', 'XCore/scotch/kgraph_map_rb_map.c', 'XCore/scotch/wgraph_part_zr.c', 'XCore/scotch/mesh_coarsen.c', 'XCore/scotch/bgraph_bipart_df_loop.c', 'XCore/scotch/vdgraph_separate_ml.c', 'XCore/scotch/library_graph_map_io_f.c', 'XCore/scotch/kgraph_map_bd.c', 'XCore/scotch/library_mapping.c', 'XCore/scotch/library_dgraph_io_load.c', 'XCore/scotch/library_strat.c', 'XCore/scotch/wgraph.c', 'XCore/scotch/graph_check.c', 'XCore/scotch/library_dorder.c', 'XCore/scotch/dgraph_halo.c', 'XCore/scotch/vmesh_check.c', 'XCore/scotch/library_dgraph_order_tree_dist.c', 'XCore/scotch/library_dgraph_stat_f.c', 'XCore/scotch/library_graph_io_chac.c', 'XCore/scotch/graph_list.c', 'XCore/scotch/bgraph_bipart_ml.c', 'XCore/scotch/hmesh_order_hx.c', 'XCore/scotch/library_mesh_io_scot.c', 'XCore/scotch/bgraph_bipart_ex.c', 'XCore/scotch/wgraph_part_gp.c', 'XCore/scotch/hmesh.c', 'XCore/scotch/library_mesh_io_habo_f.c', 'XCore/scotch/library_dgraph_coarsen_f.c', 'XCore/scotch/library_graph_io_mmkt_f.c', 'XCore/scotch/hdgraph_order_sq.c', 'XCore/scotch/mesh_io_habo.c', 'XCore/scotch/dorder_io_block.c', 'XCore/scotch/hmesh_order_cp.c', 'XCore/scotch/vgraph_check.c', 'XCore/scotch/vdgraph_separate_st.c', 'XCore/scotch/bdgraph_gather_all.c', 'XCore/scotch/arch_deco2.c', 'XCore/scotch/dgraph_view.c', 'XCore/scotch/common_file_compress.c', 'XCore/scotch/graph_ielo.c', 'XCore/scotch/bdgraph_bipart_df.c', 'XCore/scotch/library_dgraph_check_f.c', 'XCore/scotch/hmesh_hgraph.c', 'XCore/scotch/hgraph_order_hx.c', 'XCore/scotch/library_graph_io_scot_f.c', 'XCore/scotch/hall_order_hx.c', 'XCore/scotch/bgraph_bipart_st.c', 'XCore/scotch/bgraph_store.c', 'XCore/scotch/library_geom_f.c', 'XCore/scotch/hdgraph_fold.c', 'XCore/scotch/library_common_f.c', 'XCore/scotch/hgraph_order_cp.c', 'XCore/scotch/kdgraph_map_rb_part.c', 'XCore/scotch/common.c', 'XCore/scotch/library_graph_color_f.c', 'XCore/scotch/hmesh_order_si.c', 'XCore/scotch/library_arch.c', 'XCore/scotch/library_dgraph_halo_f.c', 'XCore/scotch/arch_cmpltw.c', 'XCore/scotch/library_error_exit.c', 'XCore/scotch/arch_deco.c', 'XCore/scotch/library_graph_induce.c', 'XCore/scotch/bgraph_bipart_gg.c', 'XCore/scotch/dgraph_match_scan.c', 'XCore/scotch/library_dgraph_induce_f.c', 'XCore/scotch/arch_sub.c', 'XCore/scotch/common_string.c', 'XCore/scotch/dgraph_halo_fill.c', 'XCore/scotch/graph_io_mmkt.c', 'XCore/scotch/library_dgraph_order_io_f.c', 'XCore/scotch/library_graph_order.c', 'XCore/scotch/library_graph_base_f.c', 'XCore/scotch/library_dgraph_map_view_f.c', 'XCore/scotch/library_memory.c', 'XCore/scotch/library_dgraph_check.c', 'XCore/scotch/parser.c', 'XCore/scotch/hgraph_order_si.c', 'XCore/scotch/library_graph_io_chac_f.c', 'XCore/scotch/kgraph_band.c', 'XCore/scotch/wgraph_part_rb.c', 'XCore/scotch/library_graph_coarsen_f.c', 'XCore/scotch/dorder_gather.c', 'XCore/scotch/bdgraph_bipart_sq.c', 'XCore/scotch/library_mesh_order.c', 'XCore/scotch/common_thread.c', 'XCore/scotch/dgraph_match.c', 'XCore/scotch/library_graph_map_view_f.c', 'XCore/scotch/library_graph_part_ovl.c', 'XCore/scotch/vgraph_separate_fm.c', 'XCore/scotch/vmesh_store.c', 'XCore/scotch/common_integer.c', 'XCore/scotch/vgraph_separate_th.c', 'XCore/scotch/arch_torus.c', 'XCore/scotch/library_random_f.c', 'XCore/scotch/library_dgraph_order_perm_f.c', 'XCore/scotch/library_error.c', 'XCore/scotch/hdgraph_order_st.c', 'XCore/scotch/graph_io_scot.c', 'XCore/scotch/dorder_tree_dist.c', 'XCore/scotch/library_dgraph_stat.c', 'XCore/scotch/vmesh_separate_fm.c', 'XCore/scotch/dorder_perm.c', 'XCore/scotch/vgraph_store.c', 'XCore/scotch/library_graph_io_mmkt.c', 'XCore/scotch/dorder_io.c', 'XCore/scotch/library_dgraph_order_perm.c', 'XCore/scotch/kgraph_map_df_loop.c', 'XCore/scotch/vgraph_separate_bd.c', 'XCore/scotch/dgraph_io_load.c', 'XCore/scotch/mesh_induce_sepa.c', 'XCore/scotch/hgraph_induce.c', 'XCore/scotch/hgraph_order_cc.c', 'XCore/scotch/library_errcom.c', 'XCore/scotch/dorder.c', 'XCore/scotch/kgraph_map_df.c', 'XCore/scotch/kgraph_map_rb_part.c', 'XCore/scotch/library_graph_diam_f.c', 'XCore/scotch/library_graph_map_f.c', 'XCore/scotch/library_dgraph_order_f.c', 'XCore/scotch/bdgraph_bipart_ml.c', 'XCore/scotch/library_dgraph_coarsen.c', 'XCore/scotch/bdgraph_bipart_ex.c', 'XCore/scotch/library_graph_io_scot.c', 'XCore/scotch/library_mesh_graph_f.c', 'XCore/scotch/hdgraph_check.c', 'XCore/scotch/order.c', 'XCore/scotch/hmesh_order_nd.c', 'XCore/scotch/dorder_io_tree.c', 'XCore/scotch/dmapping_io.c', 'XCore/scotch/common_file.c', 'XCore/scotch/library_dgraph_gather_f.c', 'XCore/scotch/wgraph_part_ml.c', 'XCore/scotch/hgraph_check.c', 'XCore/scotch/order_check.c', 'XCore/scotch/dummysizes.c', 'XCore/scotch/kgraph.c', 'XCore/scotch/dgraph_match_check.c', 'XCore/scotch/vgraph_separate_vw.c', 'XCore/scotch/kgraph_map_rb.c', 'XCore/scotch/library_dgraph_io_save.c', 'XCore/scotch/bgraph.c', 'XCore/scotch/hdgraph_induce.c', 'XCore/scotch/hdgraph_gather.c', 'XCore/scotch/library_dgraph_band.c', 'XCore/scotch/bdgraph_bipart_st.c', 'XCore/scotch/arch_hcub.c', 'XCore/scotch/library_dgraph_scatter.c', 'XCore/scotch/library_graph_color.c', 'XCore/scotch/hgraph_order_nd.c', 'XCore/scotch/arch_tleaf.c', 'XCore/scotch/library_order.c', 'XCore/scotch/bdgraph.c', 'XCore/scotch/gain.c', 'XCore/scotch/library_mesh_order_f.c', 'XCore/scotch/library_graph_induce_f.c', 'XCore/scotch/wgraph_part_st.c', 'XCore/scotch/vdgraph_separate_bd.c', 'XCore/scotch/hmesh_order_gp.c', 'XCore/scotch/common_stub.c', 'XCore/scotch/library_mesh_f.c', 'XCore/scotch/library_parser_f.c', 'XCore/scotch/bgraph_bipart_fm.c', 'XCore/scotch/vgraph_separate_zr.c', 'XCore/scotch/geom.c', 'XCore/scotch/vdgraph.c', 'XCore/scotch/mesh.c', 'XCore/scotch/library_graph_part_ovl_f.c', 'XCore/scotch/library_random.c', 'XCore/scotch/vmesh_separate_zr.c', 'XCore/scotch/parser_ll.c', 'XCore/scotch/library_dgraph_io_load_f.c', 'XCore/scotch/wgraph_check.c', 'XCore/scotch/vdgraph_check.c', 'XCore/scotch/wgraph_part_gg.c', 'XCore/scotch/arch_dist.c', 'XCore/scotch/bgraph_bipart_bd.c', 'XCore/scotch/graph_match.c', 'XCore/scotch/hgraph_induce_edge.c', 'XCore/scotch/library_graph_check_f.c', 'XCore/scotch/graph_io.c', 'XCore/scotch/hgraph_order_gp.c', 'XCore/scotch/graph_base.c']
#==============================================================================
cpp_srcs2 = ['XCore/paradigma/pdm.c', 'XCore/paradigma/ppart/pdm_hilbert.c', 'XCore/paradigma/ppart/pdm_part.c', 'XCore/paradigma/ppart/pdm_part_coarse_mesh.c', 'XCore/paradigma/ppart/pdm_part_geom.c', 'XCore/paradigma/ppart/pdm_part_renum.c', 'XCore/paradigma/ext_wrapper/pdm_ext_wrapper.c', 'XCore/paradigma/io/pdm_error.c', 'XCore/paradigma/io/pdm_printf.c', 'XCore/paradigma/mesh/pdm_cellface_orient.c', 'XCore/paradigma/mesh/pdm_dmesh_nodal.c', 'XCore/paradigma/mesh/pdm_elt_parent_find.c', 'XCore/paradigma/mesh/pdm_geom_elem.c', 'XCore/paradigma/mesh/pdm_graph_bound.c', 'XCore/paradigma/mesh/pdm_line.c', 'XCore/paradigma/mesh/pdm_mesh_check.c', 'XCore/paradigma/mesh/pdm_mesh_nodal.c', 'XCore/paradigma/mesh/pdm_part_bound.c', 'XCore/paradigma/mesh/pdm_plane.c', 'XCore/paradigma/mesh/pdm_points_merge.c', 'XCore/paradigma/mesh/pdm_polygon.c', 'XCore/paradigma/mesh/pdm_surf_mesh.c', 'XCore/paradigma/mesh/pdm_surf_part.c', 'XCore/paradigma/mesh/pdm_triangle.c', 'XCore/paradigma/meshgen/pdm_dcube_gen.c', 'XCore/paradigma/meshgen/pdm_poly_surf_gen.c', 'XCore/paradigma/mpi_wrapper/mpi/pdm_mpi.c', 'XCore/paradigma/mpi_wrapper/mpi/pdm_mpi_ext_dependencies.c', 'XCore/paradigma/mpi_wrapper/no_mpi/pdm_no_mpi.c', 'XCore/paradigma/mpi_wrapper/no_mpi/pdm_no_mpi_ext_dependencies.c', 'XCore/paradigma/pario/pdm_file_par.c', 'XCore/paradigma/pario/pdm_file_seq.c', 'XCore/paradigma/pario/pdm_io.c', 'XCore/paradigma/pario/pdm_io_tab.c', 'XCore/paradigma/pario/pdm_writer.c', 'XCore/paradigma/pario/pdm_writer_ensight.c', 'XCore/paradigma/pario/pdm_writer_ensight_case.c', 'XCore/paradigma/struct/pdm_binary_search.c', 'XCore/paradigma/struct/pdm_block_to_block.c', 'XCore/paradigma/struct/pdm_block_to_part.c', 'XCore/paradigma/struct/pdm_box.c', 'XCore/paradigma/struct/pdm_box_tree.c', 'XCore/paradigma/struct/pdm_cuthill.c', 'XCore/paradigma/struct/pdm_dbbtree.c', 'XCore/paradigma/struct/pdm_distrib.c', 'XCore/paradigma/struct/pdm_global_mean.c', 'XCore/paradigma/struct/pdm_gnum.c', 'XCore/paradigma/struct/pdm_handles.c', 'XCore/paradigma/struct/pdm_hash_tab.c', 'XCore/paradigma/struct/pdm_morton.c', 'XCore/paradigma/struct/pdm_octree.c', 'XCore/paradigma/struct/pdm_octree_seq.c', 'XCore/paradigma/struct/pdm_order.c', 'XCore/paradigma/struct/pdm_part_graph.c', 'XCore/paradigma/struct/pdm_part_to_block.c', 'XCore/paradigma/struct/pdm_sort.c', 'XCore/paradigma/util/pdm_fortran_to_c_string.c', 'XCore/paradigma/util/pdm_mpi_node_first_rank.c', 'XCore/paradigma/util/pdm_remove_blank.c', 'XCore/paradigma/util/pdm_timer.c']
#==============================================================================