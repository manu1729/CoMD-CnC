echo $1
echo $2
export OCR_CONFIG=./ocr-config/x86-pthread-x86/mach-hc-$1w.cfg
./comd.exe -N $2
