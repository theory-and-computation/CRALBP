#!/bin/sh

SCRIPTS=.
alias vmd="/Applications/VMD_1.9.4a57-arm64-Rev12.app/Contents/Resources/VMD.app/Contents/MacOS/VMD"

function CalcDensityDX ()
{

  local src_dir=/Users/danielsantos/Downloads/CRALBP_analysis/SMC_stripped_and_aligned_2_4-10
  local dst_dir=SMC_volmap_RTA

  mkdir -p ${dst_dir}

  for dcd in ${src_dir}/*.dcd; do
    dcd=${dcd##*/}
    dcd=${dcd%.*}
    echo ${dcd}

    vmd -dispdev text -e ${SCRIPTS}/calcRTAVolmap.tcl -args                 \
      ${src_dir}/stripped.psf         ${src_dir}/${dcd}.dcd                  \
      ${dst_dir}/tmp.dx             ${dst_dir}/${dcd}.dx
  done
}

function CombineDX ()
{

  local src_dir=SMC_volmap_RTA
  local dst_dir=SMC_volmap_RTA_sum
  local sys_name=${1}
  local start_traj=${2}

  mkdir -p ${dst_dir}

  prev_dx='0'

  for dx in ${src_dir}/*.dx; do
    dx_name=${dx##*/}
    dx_name=${dx_name%.*}

    #if is gonna depend on the start of the stripped file number
    if [[ "$dx_name" == "${sys_name}-${start_traj}-stripped" ]]; then
      cp ${src_dir}/${dx_name}.dx  ${dst_dir}/sum-${dx_name}.dx
    else
      vmd -dispdev text -e ${SCRIPTS}/combine_dx.tcl -args                  \
        ${src_dir}/${dx_name}.dx        ${dst_dir}/sum-${prev_dx}.dx        \
        ${dst_dir}/sum-${dx_name}.dx
    fi

    prev_dx=${dx_name}
  done
}

basename="smallChargedMembrane_2-ions-wats"
trajectory_start_number="00004"

CalcDensityDX $basename
CombineDX $basename $trajectory_start_number
