#pack.nec<-c("sfsmisc","shapefiles","leaps","raster","maptools","rgdal","sp","gstat","geoR","MASS","plyr","RPostgreSQL","RSAGA")
#install.packages(pack.nec)

library(shapefiles)
library(leaps)
library(raster)
#library(maptools)
#library(rgdal)
library(sp)
library(gstat)
library(geoR)
library(MASS)
library(plyr)
library(RPostgreSQL)
library(RSAGA)
library(sfsmisc)
library(sf) # nouveau maptools
library(RSAGA)
library(ggplot2)


chemin='data/CODE_KRIGE_NEW'
setwd(chemin)

env <- rsaga.env(path = "/saga-9.7.2_x64")
print(env)

x.coord.grille=seq(640000,839990, by=500)
y.coord.grille=seq(6199990, 6360000, by=500)
grille<-as.data.frame(xy.grid(x.coord.grille, y.coord.grille))
id<-as.data.frame(seq(1, nrow(grille), by=1))
grille<-cbind(id, grille)
names(grille)<-c('id','x','y')


#CHARGEMENT DU MNT

mnt15<-raster("data/MNT_T15.tif")
crs(mnt15) <- CRS("+init=epsg:2154")  # Attribution de la projection
names(mnt15) <- "mnt"

# Sauvegarde du MNT dans un format SAGA (fichier .sgrd)
writeRaster(x = mnt15, filename = "mnt15.sgrd", format = "SAGA", overwrite = TRUE)

# Exécution du module MRVBF (Module 8 du tool ta_morphometry)
#rsaga.geoprocessor("ta_morphometry", module = 8,
#                   param = list(DEM       = "mnt15.sgrd",
#                                MRVBF     = "MRVBF.sgrd",
#                                MRRTF     = "MRRTF.sgrd",
#                                CLASSIFY  = FALSE,
#                                T_SLOPE   = 100,
#                                T_PCTL_V  = 0.4,
#                                T_PCTL_R  = 0.35,
#                                P_SLOPE   = 4,
#                                P_PCTL    = 3,
#                                MAX_RES   = 100),
#                   env = env)
#
## Chargement du résultat MRVBF
MRVBF <- raster("MRVBF.sdat")
names(MRVBF) <- "mrvbf"

#writeRaster(x = MRVBF, filename = "MRVBF.tif", format = "GTiff", overwrite = TRUE)

pente<-terrain(mnt15,opt="slope",unit="degrees")
aspect<-terrain(mnt15,opt="aspect",unit="degrees")
# Charger le MNT original
mnt <- raster("mnt15.sgrd")

mnt_df <- as.data.frame(mnt, xy = TRUE)

# Visualiser avec ggplot
ggplot(data = mnt_df, aes(x = x, y = y, fill = mnt)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma", name = "Altitude (m)", begin = 0.2, end = 1) +
  coord_fixed() +
  labs(title = "",
       x = "Longitude", y = "Latitude") +
  theme_minimal()# Réduire la résolution par agrégation (par exemple, par un facteur de 4)




library(ggplot2)

# Convertir le raster en data frame
pente_df <- as.data.frame(pente, xy = TRUE)

# Visualiser avec ggplot
ggplot(data = pente_df, aes(x = x, y = y, fill = slope)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma",name = "Pente (en degrés)", begin = 0.2, end = 1) +
  coord_fixed() +
  labs(title = "",
       x = "Longitude", y = "Latitude") +
  theme_minimal()# Réduire la résolution par agrégation (par exemple, par un facteur de 4)

aspect_df <- as.data.frame(aspect, xy = TRUE)

# Visualiser avec ggplot
ggplot(data = aspect_df, aes(x = x, y = y, fill = aspect)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma",name = "Orientations (en degrés)", begin = 0.2, end = 1) +
  coord_fixed() +
  labs(title = "",
       x = "Longitude", y = "Latitude") +
  theme_minimal()# Réduire la résolution par agrégation (par exemple, par un facteur de 4)





# Ce facteur dépend de tes données et du compromis entre vitesse et précision.
#mnt_lowres <- aggregate(mnt, fact = 4, fun = mean)
#
## Définir la mer comme altitude <= 0 sur le MNT à basse résolution
#mer_mask <- mnt_lowres <= 0
#
## Transformer le masque : conserver 1 pour la mer, mettre NA pour la terre
#mer_mask2 <- mer_mask
#mer_mask2[mer_mask == 0] <- NA
#
## Calculer le raster de distance
#distmer <- distance(mer_mask2)

# Sauvegarder le résultat
#writeRaster(distmer, "distmer_lowres.tif", format = "GTiff", overwrite = TRUE)
distmer<-raster("distmer_lowres.tif")
# Visualiser
plot(distmer)

distmer_df <- as.data.frame(distmer, xy = TRUE)

ggplot(data = distmer_df, aes(x = x, y = y, fill = distmer_lowres)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Distance à la mer (m)") +
  coord_fixed() +
  labs(title = "Carte des distances à la mer",
       x = "Longitude", y = "Latitude") +
  theme_minimal()# Réduire la résolution par agrégation (par exemple, par un facteur de 4)


distmer.rs<-resample(distmer, mnt15, method='ngb')
MRVBF.rs<-resample(MRVBF, mnt15, method='ngb')

mrvbf_df <- as.data.frame(MRVBF, xy = TRUE)

ggplot(data = mrvbf_df, aes(x = x, y = y, fill = mrvbf)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma",name = "Indice de vallée", begin = 0.2, end = 1) +
  coord_fixed() +
  labs(title = "",
       x = "Longitude", y = "Latitude") +
  theme_minimal()# Réduire la résolution par agrégation (par exemple, par un facteur de 4)



covar<-stack(mnt15,pente,aspect,distmer.rs,MRVBF.rs)

#writeRaster(covar, filename = "covar_stack.grd", overwrite = TRUE)


covar.mnt15.pente.aspect<-stack(mnt15,pente,aspect)
covar.distmer<-stack(distmer)
covar.mrvbf<-stack(MRVBF)


grille.covar<-extract(covar, grille[,2:3], method='simple')

grille.mnt15.pente.aspect<-extract(covar.mnt15.pente.aspect, grille[,2:3], method='simple')
grille.distmer<-extract(covar.distmer, grille[,2:3], method='simple')
grille.mrvbf<-extract(covar.mrvbf, grille[,2:3], method='simple')


grille.tot.test<-cbind(as.data.frame(grille), grille.covar)
grille.tot<-cbind(as.data.frame(grille), grille.mnt15.pente.aspect, grille.distmer, grille.mrvbf)


write.table(grille.tot.test, 'grille-covar-test.txt', row.names=F)
write.table(grille.tot, 'grille-covar.txt', row.names=F)


stations<-read.table('data/resultats_2010-2020.csv', header=T, sep=',')
stations.covar<-extract(covar, stations[,c("X","Y")], method='simple')
stations.mnt15.pente.aspect<-extract(covar.mnt15.pente.aspect, stations[,c("X","Y")], method='simple')
stations.distmer<-extract(covar.distmer, stations[,c("X","Y")], method='simple')
stations.mrvbf<-extract(covar.mrvbf, stations[,c("X","Y")], method='simple')

stations.tot<-cbind(as.data.frame(stations), stations.mnt15.pente.aspect, stations.distmer, stations.mrvbf)
stations.tot.test<-cbind(as.data.frame(stations), stations.covar)

write.table(stations.tot, 'stations-covar.txt', row.names=F)
write.table(stations.tot.test, 'stations-covar-test.txt', row.names=F)
