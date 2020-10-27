# Download example data
DATADIR=tmpdata
IWAVEATLASURL="https://cloudstor.aarnet.edu.au/plus/s/DnOYpT3ZPnqjNsS/download"
IWAVEATLASFILE=NWS_2km_GLORYS_hex_2013_2014_SSHBC_Harmonics.nc

# IWAVE atlas from a cloudstor public repo
cd $DATADIR
echo "Downloading internal wave atlas netcdf data ~160 GB"
echo $IWAVEATLASURL $IWAVEATLASFILE
curl $IWAVEATLASURL -o $IWAVEATLASFILE
