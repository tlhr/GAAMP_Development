#!/bin/bash
INSTALLDIR='../bin'
mkdir -p $INSTALLDIR
cp SIZES $INSTALLDIR

CC='g++'
FC='ifort'
CFLAGS='-g  -Wno-unused-variable -Wno-write-strings'
NLOPTLIB='/home/huanglei/tools/mynlopt/nlopt-2.2.4/.libs/libnlopt.a'

$CC $CFLAGS -o $INSTALLDIR/1d-org 1D-fitting_org.cpp ff.cpp $NLOPTLIB
$CC $CFLAGS -o $INSTALLDIR/1d-rotamer-fitting 1D-rotamer-fitting.cpp ff.cpp $NLOPTLIB 
$CC $CFLAGS -o $INSTALLDIR/org-rotamer-E 1D-rotamer-fitting_org.cpp ff.cpp $NLOPTLIB
$CC $CFLAGS -o $INSTALLDIR/acceptor mol_H_acceptor.cpp ff.cpp $NLOPTLIB -static
$CC $CFLAGS -o $INSTALLDIR/acceptor_para mol_H_acceptor_para.cpp ff.cpp $NLOPTLIB -static
$CC $CFLAGS -o $INSTALLDIR/add-tip3 add_tip3.cpp -static
$FC -o $INSTALLDIR/cgrid ConnollyGrid.f Surface.f -static
$CC $CFLAGS -o $INSTALLDIR/check-b0-theta0 check_b0_theta0.cpp ff.cpp -static
$CC $CFLAGS -o $INSTALLDIR/check_lj modify_lj.cpp -static
$CC $CFLAGS -o $INSTALLDIR/clustering-phi clustering.cpp -static
$CC $CFLAGS -o $INSTALLDIR/cut_LJ14 cut-LJ14.cpp -static
$CC $CFLAGS -o $INSTALLDIR/gen_soft_list gen_soft_list.cpp ff.cpp $NLOPTLIB -static
$CC $CFLAGS -o $INSTALLDIR/donor mol_H_donor.cpp ff.cpp $NLOPTLIB -static
$CC $CFLAGS -o $INSTALLDIR/donor_para mol_H_donor_para.cpp ff.cpp $NLOPTLIB -static
$CC $CFLAGS -o $INSTALLDIR/equiv_atom equiv_atom.cpp -static
$CC $CFLAGS -o $INSTALLDIR/exclude_H2 exclude-H2.cpp -static
$CC $CFLAGS -o $INSTALLDIR/fitcharge fitcharge.cpp ff.cpp $NLOPTLIB -static
$CC $CFLAGS -o $INSTALLDIR/fitcharge-again fitcharge_again.cpp ff.cpp $NLOPTLIB -static
$CC $CFLAGS -o $INSTALLDIR/gen-esp gen-esp.cpp -static
$CC $CFLAGS -o $INSTALLDIR/gen_xpsf gen_xpsf.cpp -static
$CC $CFLAGS -o $INSTALLDIR/mm_pes mm_pes.cpp ff.cpp $NLOPTLIB -static
$CC $CFLAGS -o $INSTALLDIR/mm_pes_large mm_pes_large.cpp ff.cpp $NLOPTLIB -static
$CC $CFLAGS -o $INSTALLDIR/pdb_to_crd pdb_to_crd.cpp ff.cpp -static
$CC $CFLAGS -o $INSTALLDIR/prep_cgenff extract_param.cpp ff.cpp -static
$CC $CFLAGS -o $INSTALLDIR/qm-1d-scan QM_1D_scan.cpp ff.cpp $NLOPTLIB -static
$CC $CFLAGS -o $INSTALLDIR/qm-1d-scan_large QM_1D_scan_large.cpp ff.cpp $NLOPTLIB -static
$CC $CFLAGS -o $INSTALLDIR/qm-1d-scan_large_para QM_1D_scan_large_para.cpp ff.cpp $NLOPTLIB -static
$CC $CFLAGS -o $INSTALLDIR/qm-1d-scan-para QM_1D_scan_para.cpp ff.cpp $NLOPTLIB -static
$CC $CFLAGS -o $INSTALLDIR/qm-1d-scan-single QM_1D_scan-single.cpp ff.cpp $NLOPTLIB -static
$CC $CFLAGS -o $INSTALLDIR/qm-1d-scan-single_para QM_1D_scan-single_para.cpp ff.cpp $NLOPTLIB -static
$CC $CFLAGS -o $INSTALLDIR/qm-rotamer-scan QM_rotamer_scan.cpp ff.cpp $NLOPTLIB -static
$CC $CFLAGS -openmp -o $INSTALLDIR/qm-rotamer-scan-large QM_rotamer_scan_large.cpp ff.cpp $NLOPTLIB -static
$CC $CFLAGS -o $INSTALLDIR/update-tor-para update-torsion-para.cpp ff.cpp $NLOPTLIB -static
$CC $CFLAGS -o $INSTALLDIR/update-xpsf update-xpsf.cpp -static

$FC -o $INSTALLDIR/cgrid drude-ConnollyGrid.f drude-Surface.f -static
$CC $CFLAGS -o $INSTALLDIR/drude-1d-fitting drude-1D-fitting.cpp ff.cpp $NLOPTLIB
$CC $CFLAGS -o $INSTALLDIR/drude-1d-rotamer-fitting drude-1D-rotamer-fitting.cpp ff.cpp $NLOPTLIB
$CC $CFLAGS -o $INSTALLDIR/drude-fitcharge_again drude-fitcharge_again.cpp ff.cpp $NLOPTLIB -static
$CC $CFLAGS -o $INSTALLDIR/drude-fitcharge drude-fitcharge.cpp ff.cpp $NLOPTLIB -static
$CC $CFLAGS -o $INSTALLDIR/drude-gen-esp drude-gen-esp.cpp -static
$CC $CFLAGS -o $INSTALLDIR/drude-gen_xpsf drude-gen_xpsf.cpp -static
$CC $CFLAGS -o $INSTALLDIR/drude-modify_lj drude-modify_lj.cpp -static
$CC $CFLAGS -o $INSTALLDIR/drude-ff drude-modify_rtf.cpp ff.cpp $NLOPTLIB -static
$CC $CFLAGS -o $INSTALLDIR/drude-scale-alpha drude-scale-alpha.cpp -static
$CC $CFLAGS -o $INSTALLDIR/drude-scale-thole drude-scale-thole.cpp -static
$CC $CFLAGS -o $INSTALLDIR/drude-scale-vdw-r6 drude-scale-vdw-r6.cpp -static
$CC $CFLAGS -o $INSTALLDIR/drude-update-torsion-para drude-update-torsion-para.cpp ff.cpp $NLOPTLIB -static

