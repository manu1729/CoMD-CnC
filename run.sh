echo $1
echo $2
export OCR_CONFIG=./ocr-config/x86-pthread-x86/mach-hc-$1w.cfg
#export OCR_CONFIG=/home/mshanth/software/cnc-ocr/examples/CoMD-CnC/ocr-config/x86-pthread-fsim/mach-1block.cfg
./comd.exe -N $2
