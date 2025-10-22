cmake_minimum_required(VERSION 3.10)
project(mdmc LANGUAGES Fortran)

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Define source files
set(SOURCES
    algorithms/mc_control.f90
    algorithms/md_control.f90
    algorithms/md_control_time_corr.f90
    algorithms/md_gridsearch_control.f90
    algorithms/mdmc_control.f90
    algorithms/mdmc_control_time_corr.f90
    functions/g_d_rt_fom.f90
    functions/lennard_jones.f90
    functions/rdf_fom.f90
    functions/s_q_omega_fom.f90
    functions/s_q_time_fom.f90
    readers/alloc_q_omega_arrays.f90
    readers/fnc_constraint_reader.f90
    readers/g_d_readers.f90
    readers/s_qo_reader.f90
    readers/s_qt_reader.f90
    readers/structure_reader.f90
    xmlf90-1.2/sax/flib_sax.f90
    xmlf90-1.2/sax/m_buffer.f90
    xmlf90-1.2/sax/m_charset.f90
    xmlf90-1.2/sax/m_converters.f90
    xmlf90-1.2/sax/m_debug.f90
    xmlf90-1.2/sax/m_dictionary.f90
    xmlf90-1.2/sax/m_elstack.f90
    xmlf90-1.2/sax/m_entities.f90
    xmlf90-1.2/sax/m_fsm.f90
    xmlf90-1.2/sax/m_io.f90
    xmlf90-1.2/sax/m_reader.f90
    xmlf90-1.2/sax/m_xml_error.f90
    xmlf90-1.2/sax/m_xml_parser.f90
    xmlf90-1.2/wxml/flib_wxml.f90
    xmlf90-1.2/wxml/m_wxml_buffer.f90
    xmlf90-1.2/wxml/m_wxml_core.f90
    xmlf90-1.2/wxml/m_wxml_dictionary.f90
    xmlf90-1.2/wxml/m_wxml_elstack.f90
    xmlf90-1.2/wxml/m_wxml_text.f90
    xmlf90-1.2/xpath/flib_xpath.f90
    xmlf90-1.2/xpath/m_path.f90
    common_block.f90
    configuration.f90
    control_containers.f90
    converters.f90
    fnc_constraint.f90
    fom_readers.f90
    func_param.f90
    func_params_wrapper.f90
    function.f90
    handler.f90
    handler_structure_part.f90
    histogram.f90
    histogram_cutdown.f90
    md_properties.f90
    mdmc.f90
    phasespace.f90
    rdf.f90
    s_q_omega.f90
    s_q_time.f90
    structure.f90
    structure_nn_methods.f90
    structure_type_definition.f90
    tic_toc.f90
    time_corr_algorithm.f90
    time_corr_hist_container.f90
    various_constants.f90
)

# Create executable
add_executable(mdmc ${SOURCES})

# Set build configurations
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")