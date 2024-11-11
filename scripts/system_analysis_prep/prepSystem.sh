#!/bin/sh

SCRIPTS=.
alias vmd="/Applications/VMD_1.9.4a57-arm64-Rev12.app/Contents/Resources/VMD.app/Contents/MacOS/VMD"

function StripWaters () 
{
  local sys="${1}"
  local dst_dir=./RSMC_stripped

  mkdir -p ${dst_dir}

  vmd -dispdev text -e ${SCRIPTS}/strip_waters_dsv.tcl -args          \
    ../systems/${sys}.psf             ../systems/${sys}.pdb                            \
    ./target.stk                ${dst_dir}/stripped.psf               \
    ${dst_dir}/first_prot.pdb         ${sys}                                \
    ${dst_dir}
}


function AlignProt ()
{
  local src_dir=./RSMC_stripped
  local dst_dir=./RSMC_stripped_and_aligned

  mkdir -p ${dst_dir}

  vmd -dispdev text -e ${SCRIPTS}/orient.tcl -args                        \
    ${src_dir}/first_prot.pdb     ${src_dir}/first_prot.oriented.pdb

  for dcd in ${src_dir}/*.dcd; do
    dcd=$(basename ${dcd})
    vmd -dispdev text -e ${SCRIPTS}/align_prot.tcl -args                   \
      ${src_dir}/stripped.psf                 ${src_dir}/${dcd}            \
      ${src_dir}/first_prot.oriented.pdb      ${dst_dir}/${dcd}
    echo "Finished with ${dcd}!"
  done 
}


function getPBC ()
{
  local src_dir=./RSMC_stripped
  local dst_dir=.

  touch ${dst_dir}/"RSMC_PBC.txt"

  for dcd in ${src_dir}/*.dcd; do
    dcd=$(basename ${dcd})
    vmd -dispdev text -e ${SCRIPTS}/getPBCSize.tcl -args                   \
    ${src_dir}/stripped.psf                 ${src_dir}/${dcd}              \
    ${dst_dir}/"RSMC_PBC.txt"

    echo "Finished with ${dcd}!"
  done
}


function CreateDatabase ()
{
  local src_dir=./RSMC_stripped
  local dst_dir=.

  touch ${dst_dir}/"RSCM_contacts.txt"

  for dcd in ${src_dir}/*.dcd; do
    dcd=$(basename ${dcd})
    vmd -dispdev text -e getLipidProtInter-txt.tcl -args                   \
    ${src_dir}/stripped.psf                 ${src_dir}/${dcd}              \
    ${dst_dir}/"RSCM_contacts.txt"

    echo "Finished with ${dcd}!"
  done

  python createDatabase.py
}

basename='smallChargedMembrane_2-ions-wats'

StripWaters $basename
AlignProt
GetPBC
CreateDatabase

