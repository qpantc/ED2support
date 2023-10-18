library(rhdf5)

cruncep_hdf5 <- H5Fopen('./data/1950JAN.h5')
cruncep_hdf5

flux_hdf5 <- H5Fopen('./sites/FLX_GF-Guy2004_2020/FLX_GF-Guy2004_2020_2004JAN.h5')
flux_hdf5

paracou_hdf5 <- H5Fopen("./data/Paracou_2004JAN.h5")
paracou_hdf5


# 0  dlwrf H5I_DATASET  FLOAT 124 x 1 x 1
summary(paracou_hdf5$dlwrf)
summary(flux_hdf5$dlwrf)
# 1  hgt   H5I_DATASET  FLOAT 124 x 1 x 1
summary(paracou_hdf5$hgt)
summary(flux_hdf5$hgt)
# 2  nbdsf H5I_DATASET  FLOAT 124 x 1 x 1
summary(paracou_hdf5$nbdsf)
summary(flux_hdf5$nbdsf)
# 3  nddsf H5I_DATASET  FLOAT 124 x 1 x 1
summary(paracou_hdf5$nddsf)
summary(flux_hdf5$nddsf)
# 4  prate H5I_DATASET  FLOAT 124 x 1 x 1
summary(paracou_hdf5$prate)
summary(flux_hdf5$prate)
# 5  pres  H5I_DATASET  FLOAT 124 x 1 x 1
summary(paracou_hdf5$pres)
summary(flux_hdf5$pres)
# 6  sh    H5I_DATASET  FLOAT 124 x 1 x 1
summary(paracou_hdf5$sh)
summary(flux_hdf5$sh)
# 7  tmp   H5I_DATASET  FLOAT 124 x 1 x 1
summary(paracou_hdf5$tmp)
summary(flux_hdf5$tmp)
# 8  ugrd  H5I_DATASET  FLOAT 124 x 1 x 1
summary(paracou_hdf5$ugrd)
summary(flux_hdf5$ugrd)
# 9  vbdsf H5I_DATASET  FLOAT 124 x 1 x 1
summary(paracou_hdf5$vbdsf)
summary(flux_hdf5$vbdsf)
# 10 vddsf H5I_DATASET  FLOAT 124 x 1 x 1
summary(paracou_hdf5$vddsf)
summary(flux_hdf5$vddsf)
# 11 vgrd  H5I_DATASET  FLOAT 124 x 1 x 1
summary(paracou_hdf5$vgrd)
summary(flux_hdf5$vgrd)

