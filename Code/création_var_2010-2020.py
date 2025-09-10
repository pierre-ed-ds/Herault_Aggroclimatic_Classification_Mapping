
#%%
import pandas as pd

def saison(mois):
    if mois in [12, 1, 2]:
        return 'Hiver'
    elif mois in [6, 7, 8]:
        return 'Été'
    else:
        return 'Autre'


#%%

def calcul_indicateurs(df, station):
    # Supprimer les stations avec des valeurs manquantes
    stations_incompletes = df[df.isnull().any(axis=1)]['CODE_STATION'].unique()
    df = df[~df['CODE_STATION'].isin(stations_incompletes)]

    # 2. Supprimer les stations qui n'ont pas exactement 120 lignes (10 ans * 12 mois)
    station_counts = df['CODE_STATION'].value_counts()
    stations_avec_120_lignes = station_counts[station_counts > 100].index
    df = df[df['CODE_STATION'].isin(stations_avec_120_lignes)]

    print(f"{len(stations_incompletes)} stations supprimées pour valeurs manquantes")
    print(f"{sum(station_counts <= 100)} stations supprimées pour séries incomplètes (≠ 120 mois)")
            
    df['SAISON'] = df['MOIS'].apply(saison)

    # Température moyenne annuelle par station
    temp_moy_annuelle = df.groupby(['CODE_STATION', 'ANNEE'])[['TEMP_MAX_MOY', 'TEMP_MIN_MOY']]\
                          .mean().mean(axis=1).groupby('CODE_STATION').mean()

    # Cumul pluviométrique annuel moyen
    cumul_pluv_annuel = df.groupby(['CODE_STATION', 'ANNEE'])['CUMUL_PRECIP']\
                           .sum().groupby('CODE_STATION').mean()

    # Précipitations saisonnières
    precip_hiver = df[df['SAISON'] == 'Hiver'].groupby(['CODE_STATION', 'ANNEE'])['CUMUL_PRECIP']\
                     .sum().groupby('CODE_STATION').mean()
    precip_ete = df[df['SAISON'] == 'Été'].groupby(['CODE_STATION', 'ANNEE'])['CUMUL_PRECIP']\
                   .sum().groupby('CODE_STATION').mean()

    # Températures saisonnières
    temp_min_hiver = df[df['SAISON'] == 'Hiver'].groupby('CODE_STATION')['TEMP_MIN_MOY'].mean()
    temp_max_hiver = df[df['SAISON'] == 'Hiver'].groupby('CODE_STATION')['TEMP_MAX_MOY'].mean()
    temp_min_ete = df[df['SAISON'] == 'Été'].groupby('CODE_STATION')['TEMP_MIN_MOY'].mean()
    temp_max_ete = df[df['SAISON'] == 'Été'].groupby('CODE_STATION')['TEMP_MAX_MOY'].mean()

    # Mois secs
    df['MOIS_SEC'] = df.apply(lambda row: row['CUMUL_PRECIP'] < 2 * ((row['TEMP_MAX_MOY'] + row['TEMP_MIN_MOY']) / 2), axis=1)
    frequence_mois_secs = df.groupby('CODE_STATION')['MOIS_SEC'].mean()

    # Mois très secs
    df['MOIS_TRES_SEC'] = df.apply(lambda row: 1.5 * row['CUMUL_PRECIP'] < 2 * ((row['TEMP_MAX_MOY'] + row['TEMP_MIN_MOY']) / 2), axis=1)
    frequence_mois_tres_secs = df.groupby('CODE_STATION')['MOIS_TRES_SEC'].mean()

    # Assemblage des résultats
    result = pd.concat([
        temp_moy_annuelle.rename('temp_moy_annuelle'),
        cumul_pluv_annuel.rename('cumul_pluv_annuel'),
        precip_hiver.rename('precip_hiver'),
        precip_ete.rename('precip_ete'),
        temp_min_hiver.rename('temp_min_hiver'),
        temp_max_hiver.rename('temp_max_hiver'),
        temp_min_ete.rename('temp_min_ete'),
        temp_max_ete.rename('temp_max_ete'),
        frequence_mois_secs.rename('frequence_mois_secs'),
        frequence_mois_tres_secs.rename('frequence_mois_tres_secs')
    ], axis=1)

    # Ajout des coordonnées des stations
    result = result.merge(station[['CODE_STATION', 'X', 'Y']], on='CODE_STATION', how='left')

    return result
#%%

df = pd.read_csv('big_2010-2020.csv')

station = pd.read_csv('sig_clim_stations.csv')
resultats = calcul_indicateurs(df, station)


# %%
resultats
#resultats.to_csv('resultats_2010-2020.csv', index=False)
# %%

df2 = pd.read_csv('big_1990-2000.csv')
resultats2 = calcul_indicateurs(df2, station)
# %%
resultats2
#enregistrer les résultats dans un fichier CSV
#resultats2.to_csv('resultats_1990-2000.csv', index=False)
# %%
df3 = pd.read_csv('big_2000-2010.csv')
resultats3 = calcul_indicateurs(df3, station)

# %%
resultats3
#resultats3.to_csv('resultats_2000-2010.csv', index=False)
# %%
####################################################
#################### VAR SAFRAN ####################
####################################################
def transforme_donnees_safran(fichier_csv, annee_debut=1990, annee_fin=2000):
    
    # Lecture du fichier
    df = pd.read_csv(fichier_csv, sep=";")
    
    # Conversion des dates
    df["DATE_MESURE"] = pd.to_datetime(df["DATE_MESURE"], format="%d/%m/%Y")
    df["ANNEE"] = df["DATE_MESURE"].dt.year
    df["MOIS"] = df["DATE_MESURE"].dt.month
    
    # Filtrage des années
    df_filtre = df[(df["ANNEE"] >= annee_debut) & (df["ANNEE"] < annee_fin)]
    
    # Agrégation mensuelle
    result = df_filtre.groupby(["ID_POINT", "X_L93", "Y_L93", "ANNEE", "MOIS"]).agg(
        temp_max_moy = ("TSUP_H_Q", "mean"),
        temp_min_moy = ("TINF_H_Q", "mean"),
        cumul_precip = ("PRELIQ_Q", "sum"),
        temp_max_max = ("TSUP_H_Q", "max"),
        temp_min_min = ("TINF_H_Q", "min")
    ).reset_index()
    
    # Tri
    result = result.sort_values(by=["ID_POINT", "ANNEE", "MOIS"])
    
    return result
#%%

saf = pd.read_csv('Safran.csv', sep=";")

#%%

saf90_2000 = transforme_donnees_safran("Safran.csv", 1990, 2000)

#%%
# Afficher un aperçu
print(saf90_2000.head())
# %%
def calcul_indicateurs_saf(df):
    # Supprimer les stations avec des valeurs manquantes
    stations_incompletes = df[df.isnull().any(axis=1)]['ID_POINT'].unique()
    df = df[~df['ID_POINT'].isin(stations_incompletes)]

    # Supprimer les stations avec < 100 mois
    station_counts = df['ID_POINT'].value_counts()
    stations_avec_serie_complete = station_counts[station_counts > 100].index
    df = df[df['ID_POINT'].isin(stations_avec_serie_complete)]

    print(f"{len(stations_incompletes)} stations supprimées pour valeurs manquantes")
    print(f"{sum(station_counts <= 100)} stations supprimées pour séries incomplètes (≠ 120 mois)")

    # Définir la saison
    df['SAISON'] = df['MOIS'].apply(saison)

    # Moyenne annuelle des températures
    temp_moy_annuelle = df.groupby(['ID_POINT', 'ANNEE'])[['temp_max_moy', 'temp_min_moy']].mean() \
                          .mean(axis=1).groupby('ID_POINT').mean()

    # Cumul annuel de précipitations
    cumul_pluv_annuel = df.groupby(['ID_POINT', 'ANNEE'])['cumul_precip'] \
                           .sum().groupby('ID_POINT').mean()

    # Précipitations saisonnières moyennes
    precip_hiver = df[df['SAISON'] == 'Hiver'].groupby(['ID_POINT', 'ANNEE'])['cumul_precip'] \
                     .sum().groupby('ID_POINT').mean()
    precip_ete = df[df['SAISON'] == 'Été'].groupby(['ID_POINT', 'ANNEE'])['cumul_precip'] \
                   .sum().groupby('ID_POINT').mean()

    # Températures saisonnières
    temp_min_hiver = df[df['SAISON'] == 'Hiver'].groupby('ID_POINT')['temp_min_moy'].mean()
    temp_max_hiver = df[df['SAISON'] == 'Hiver'].groupby('ID_POINT')['temp_max_moy'].mean()
    temp_min_ete = df[df['SAISON'] == 'Été'].groupby('ID_POINT')['temp_min_moy'].mean()
    temp_max_ete = df[df['SAISON'] == 'Été'].groupby('ID_POINT')['temp_max_moy'].mean()

    # Mois secs (précipitations < 2 * Tmoy)
    df['MOIS_SEC'] = df.apply(lambda row: row['cumul_precip'] < 2 * ((row['temp_max_moy'] + row['temp_min_moy']) / 2), axis=1)
    frequence_mois_secs = df.groupby('ID_POINT')['MOIS_SEC'].mean()

    # Mois très secs (précipitations < 1.5 * Tmoy)
    df['MOIS_TRES_SEC'] = df.apply(lambda row: row['cumul_precip'] < 1.5 * ((row['temp_max_moy'] + row['temp_min_moy']) / 2), axis=1)
    frequence_mois_tres_secs = df.groupby('ID_POINT')['MOIS_TRES_SEC'].mean()

    # Coordonnées moyennes par station
    coords = df.groupby('ID_POINT')[['X_L93', 'Y_L93']].mean()

    # Assemblage final
    result = pd.concat([
        temp_moy_annuelle.rename('temp_moy_annuelle'),
        cumul_pluv_annuel.rename('cumul_pluv_annuel'),
        precip_hiver.rename('precip_hiver'),
        precip_ete.rename('precip_ete'),
        temp_min_hiver.rename('temp_min_hiver'),
        temp_max_hiver.rename('temp_max_hiver'),
        temp_min_ete.rename('temp_min_ete'),
        temp_max_ete.rename('temp_max_ete'),
        frequence_mois_secs.rename('frequence_mois_secs'),
        frequence_mois_tres_secs.rename('frequence_mois_tres_secs'),
        coords
    ], axis=1)

    return result.reset_index()
# %%

resultats_saf_90_2000 = calcul_indicateurs_saf(saf90_2000)
# %%

resultats_saf_90_2000

resultats_saf_90_2000 = resultats_saf_90_2000.rename(columns={
    "ID_POINT": "CODE_STATION",
    "X_L93": "X",
    "Y_L93": "Y"
})
# %%

# Enregistrer les résultats dans un fichier CSV
#resultats_saf_90_2000.to_csv('resultats_saf_90_2000.csv', index=False)
# %%

saf00_10 = transforme_donnees_safran("Safran.csv", 2000, 2010)

resultats_saf_00_10 = calcul_indicateurs_saf(saf00_10)
#%%
resultats_saf_00_10.head()
resultats_saf_00_10 = resultats_saf_00_10.rename(columns={
    "ID_POINT": "CODE_STATION",
    "X_L93": "X",
    "Y_L93": "Y"
})
# %%
#resultats_saf_00_10.to_csv('resultats_saf_00_10.csv', index=False)
# %%

saf10_20 = transforme_donnees_safran("Safran.csv", 2010, 2020)
resultats_saf_10_20 = calcul_indicateurs_saf(saf10_20)

#%%
resultats_saf_10_20.head()
resultats_saf_10_20 = resultats_saf_10_20.rename(columns={
    "ID_POINT": "CODE_STATION",
    "X_L93": "X",
    "Y_L93": "Y"
})

#resultats_saf_10_20.to_csv('resultats_saf_10_20.csv', index=False)
# %%
