Data in `Van_Well_thesis_Ag_data.xml` are as published in van Well thesis' or in van Well et al., 1985
(van Well, A. A., Verkerk, P., de Graaf, L. A., Suck, J.-B., & Copley, J. R. D. (1985). Density fluctuations
in liquid argon: Coherent dynamic structure factor along the 120-K isotherm obtained by neutron scattering.
Physical Review A, 31(5), 3391–3414. https://doi.org/10.1103/PhysRevA.31.3391 )

The script used to convert the table in thoese files into .xml file is `converted_scanned_images_into_xml_format.m`.
This script also creates `Well_s_q_omega_Ag_data_unsymmetrised.xml`, i.e. an unsymmetrised version of the T=120K, ρ=0.0176Å-3 argon S(Q, ω) data from (van Well et al., 1985)
multiplied by exp(h_bar*ω/(2*k_b*T)). This multiplication is performed to ‘un-symmetrise’ the symmetrised S(Q, ω) data, see Eq. (5) in (van Well et al., 1985).
