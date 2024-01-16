#!bash/bin

#pdb4amber stest.pdb > stest.amber.pdb

#apply the forcefiled
#tleap -f leap.in
#parmchk2 -i etest.mol2 -f mol2 -o etest.frcmod

#resize the box size
cpptraj -p etest.prmtop -y etest.ncrst -i resize_box.cpptraj
#cd $AMBERHOME/dat/leap/parm/


#pmemd -O -i ./em1_Prot-Lip.in -p ./etest.prmtop -c ./resize_etest.ncrst -o etest_em_1.out -r etest_em_1.ncrst &

#pmemd -O -i ./em2_Prot-Lip.in -p ./etest.prmtop -c ./resize_etest.ncrst -o etest_em_1.out -r etest_em_1.ncrst &i


pmemd -O -i ./eq_Prot-Lip.in -p ./etest.prmtop -c resize_etest.ncrst  -o etest.out -r etest_em.ncrst -x etest_eq_1.nc
