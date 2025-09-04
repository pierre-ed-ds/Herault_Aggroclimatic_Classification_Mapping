# Carte interactive des unités agroclimatiques de l’Hérault

Ce projet propose une visualisation interactive des unités agroclimatiques du département de l’Hérault, réalisée avec **R** et le package **leaflet**.  
L’objectif est de représenter la variabilité agroclimatique locale en tenant compte des données météorologiques et des caractéristiques géomorphologiques du territoire.

---

## Données utilisées
- Données météorologiques de **Météo-France**  
- Données départementales de l’Hérault  
- Variables géomorphologiques :
  - Modèle Numérique de Terrain (MNT)  
  - Multi-Resolution Valley Bottom Flatness (MRVBF)  
  - Pente  
  - Orientation  
  - Distance à la mer  

---

## Méthodologie
1. **Classification des stations météo**  
   Réalisée avec la méthode de clustering **ClustGeo**, intégrant à la fois la proximité géographique et la similitude climatique.  

2. **Interpolation spatiale**  
   Utilisation de **forêts aléatoires (Random Forest)** pour interpoler les unités agroclimatiques sur le territoire.  
   Les variables géomorphologiques listées ci-dessus ont été intégrées comme prédicteurs.  

3. **Période étudiée**  
   Analyse réalisée sur trois décennies : **1990 à 2020**.  

---

## Résultat
Une **carte interactive leaflet** permettant d’explorer :  
- Les unités agroclimatiques de l’Hérault  
- La distribution spatiale issue de l’interpolation  
- Les informations associées aux stations et variables explicatives  

---

## Technologies
- **Langage** : R  
- **Packages principaux** :  
  - `leaflet` (cartographie interactive)  
  - `clustGeo` (classification hiérarchique contrainte)  
  - `randomForest` (modèle d’interpolation)  
  - `raster`, `sf` (traitement spatial)  

---

## Utilisation
1. Cloner le dépôt  
   ```bash
   git clone https://github.com/ton-profil/units-agroclimatiques-herault.git
   cd units-agroclimatiques-herault
   ```
## Aperçu

![Carte interactive des unités agroclimatiques](images/carte.png)
